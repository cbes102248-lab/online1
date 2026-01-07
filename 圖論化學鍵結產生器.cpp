#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#define MAX_ATOMS 50
typedef struct {
    char symbol[3];
    char name[15];
    float weight, mp, bp, electroneg;//電負度
    int valence_e;//原價電子數
    bool is_metal;//金屬?
    int octet_req;;//八隅
} Element;
Element element_dict[] = {
    {"H", "Hydrogen", 1.008, -259.1, -252.9, 2.20, 1, false, 2},
    {"He", "Helium", 4.003, -272.2, -268.9, 0.00, 2, false, 2},
    {"Li", "Lithium", 6.941, 180.5, 1342.0, 0.98, 1, true, 8},
    {"Be", "Beryllium", 9.012, 1287.0, 2471.0, 1.57, 2, true, 8},
    {"B", "Boron", 10.811, 2075.0, 4000.0, 2.04, 3, false, 8},
    {"C", "Carbon", 12.011, 3550.0, 4827.0, 2.55, 4, false, 8},
    {"N", "Nitrogen", 14.007, -210.0, -195.8, 3.04, 5, false, 8},
    {"O", "Oxygen", 15.999, -218.8, -183.0, 3.44, 6, false, 8},
    {"F", "Fluorine", 18.998, -219.6, -188.1, 3.98, 7, false, 8},
    {"Ne", "Neon", 20.180, -248.6, -246.1, 0.00, 8, false, 8},
    {"Na", "Sodium", 22.990, 97.7, 883.0, 0.93, 1, true, 8},
    {"Mg", "Magnesium", 24.305, 650.0, 1090.0, 1.31, 2, true, 8},
    {"Al", "Aluminum", 26.982, 660.3, 2519.0, 1.61, 3, true, 8},
    {"Si", "Silicon", 28.085, 1414.0, 3265.0, 1.90, 4, false, 8},
    {"P", "Phosphorus", 30.974, 44.1, 280.5, 2.19, 5, false, 8},
    {"S", "Sulfur", 32.065, 115.2, 444.6, 2.58, 6, false, 8},
    {"Cl", "Chlorine", 35.453, -101.5, -34.0, 3.16, 7, false, 8},
    {"Ar", "Argon", 39.948, -189.3, -185.8, 0.00, 8, false, 8},
    {"K", "Potassium", 39.098, 63.4, 759.0, 0.82, 1, true, 8}
};
typedef struct { int type_idx, rem_e; } AtomNode;//分子字典定義
typedef struct { AtomNode nodes[MAX_ATOMS]; int adj[MAX_ATOMS][MAX_ATOMS], total; } MolGraph;//分子中每個原子附值
int find_el(char* sym) {//查找原子字典
    for (int i = 0; i < 19; i++) if (strcmp(element_dict[i].symbol, sym) == 0) return i;
    return -1;
}
void auto_bond(MolGraph *mg) {//計算全單鍵
    for (int i = 0; i < mg->total; i++) mg->nodes[i].rem_e = element_dict[mg->nodes[i].type_idx].valence_e;
    for (int i = 0; i < mg->total; i++) {
        for (int j = i + 1; j < mg->total; j++) {
            if (mg->adj[i][j] >= 1) { mg->nodes[i].rem_e--; mg->nodes[j].rem_e--; }
        }
    }
    for (int i = 0; i < mg->total; i++) {//嘗試多鍵
        while (mg->nodes[i].rem_e > 0) {
            bool added = false;
            for (int j = 0; j < mg->total; j++) {
                if (mg->adj[i][j] >= 1 && mg->nodes[j].rem_e > 0) {
                    mg->adj[i][j]++; mg->adj[j][i]++; mg->nodes[i].rem_e--; mg->nodes[j].rem_e--;
                    added = true; break;
                }
            }
            if (!added) break;
        }
    }
}
void analyze(MolGraph *mg) {
    int chg;
    printf("\n輸入分子總電荷 (如 SO4 2- 輸入 -2, NH4+ 輸入 1): ");
    scanf("%d", &chg);
    mg->nodes[0].rem_e -= chg; // 將電荷電子池計入中心原子
    auto_bond(mg);
    bool has_HB = false; float max_en = -1, min_en = 10;
    printf("\n%-5s %-5s %-8s %-10s %-8s %-15s\n", "編號", "原子", "總鍵數", "鍵級分配", "LP對", "空間形狀");
    for (int i = 0; i < mg->total; i++) {
        Element e = element_dict[mg->nodes[i].type_idx];
        if (e.electroneg > max_en) max_en = e.electroneg;
        if (e.electroneg < min_en && e.electroneg > 0) min_en = e.electroneg;
        int b_sum = 0, nbr = 0; char desc[32] = "";
        for (int j = 0; j < mg->total; j++) {
            if (mg->adj[i][j] > 0) {
                b_sum += mg->adj[i][j]; nbr++;
                Element ne = element_dict[mg->nodes[j].type_idx];
                if (strcmp(e.symbol, "H") == 0 && (strchr("FON", ne.symbol[0]))) has_HB = true;
                char tmp[8]; sprintf(tmp, "%d-", mg->adj[i][j]); strcat(desc, tmp);
            }
        }
        int lp = mg->nodes[i].rem_e / 2, sn = nbr + lp;//判斷LP
        char *sh;
        sh = (sn == 4) ? (lp == 0 ? "四面體" : (lp == 1 ? "三角錐" : "角形")) : //判斷形狀
           (sn == 3) ? (lp == 0 ? "平面三角" : "角形") : 
           (sn == 2) ? "線性" : "末端";
        printf("%-7d %-6s %-9d %-11s %-8d %-15s\n", i, e.symbol, b_sum, desc, lp, sh);
    }
    float de = max_en - min_en;
    printf("\n[作用力] 極性強度: %.2f | 作用力: %s%s倫敦分散力\n", de, has_HB?"氫鍵, ":"", de>0.5?"偶極-偶極, ":"");//作用力
}
int main() {
    MolGraph mg; memset(mg.adj, 0, sizeof(mg.adj));
    char f[100]; int n = 0;//定義正列
    printf("輸入化學式 (如 (C)(1)(O)(2)): "); scanf("%s", f);
    char *p = f;
    while (*p) {
        if (*p == '(') {//查找原子字典與數量
            p++; 
            char s[10] = {0};
            int i = 0;
            while (*p && *p != ')') s[i++] = *p++; 
            if (*p == ')') p++; 
            int idx = find_el(s);
            if (idx != -1) {
                int c = 1;
                if (*p == '(') {
                    p++; 
                    char nb[10] = {0};
                    int j = 0;
                    while (*p && *p != ')') nb[j++] = *p++;
                    if (*p == ')') p++; 
                    c = atoi(nb);
                }
                for (int k = 0; k < c; k++) {
                    if (n < MAX_ATOMS) mg.nodes[n++].type_idx = idx;
                }
            }
        } else {
            p++; 
        }
    }
    mg.total = n;
    printf("\n原子編號清單:\n");
    for(int i=0; i<n; i++) printf("[%d]: %s  ", i, element_dict[mg.nodes[i].type_idx].symbol);
    printf("\n\n輸入連接關係 (A B), 輸入 -1 -1 結束:\n");
    int a, b;
    while (scanf("%d %d", &a, &b) && a != -1) if (a<n && b<n) { mg.adj[a][b] = 1; mg.adj[b][a] = 1; }
    analyze(&mg);
    return 0;
}