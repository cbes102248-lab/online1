// script.js

document.addEventListener('DOMContentLoaded', () => {
    console.log("化學視界網站已加載！");
});
// script.js 完整元素資料
const elements = [
    // 第一週期
    { num: 1, symbol: "H", name: "氫", row: 1, col: 1 },
    { num: 2, symbol: "He", name: "氦", row: 1, col: 18 },
    // 第二週期
    { num: 3, symbol: "Li", name: "鋰", row: 2, col: 1 },
    { num: 4, symbol: "Be", name: "鈹", row: 2, col: 2 },
    { num: 5, symbol: "B", name: "硼", row: 2, col: 13 },
    { num: 6, symbol: "C", name: "碳", row: 2, col: 14 },
    { num: 7, symbol: "N", name: "氮", row: 2, col: 15 },
    { num: 8, symbol: "O", name: "氧", row: 2, col: 16 },
    { num: 9, symbol: "F", name: "氟", row: 2, col: 17 },
    { num: 10, symbol: "Ne", name: "氖", row: 2, col: 18 },
    // 第三週期
    { num: 11, symbol: "Na", name: "鈉", row: 3, col: 1 },
    { num: 12, symbol: "Mg", name: "鎂", row: 3, col: 2 },
    { num: 13, symbol: "Al", name: "鋁", row: 3, col: 13 },
    { num: 14, symbol: "Si", name: "矽", row: 3, col: 14 },
    { num: 15, symbol: "P", name: "磷", row: 3, col: 15 },
    { num: 16, symbol: "S", name: "硫", row: 3, col: 16 },
    { num: 17, symbol: "Cl", name: "氯", row: 3, col: 17 },
    { num: 18, symbol: "Ar", name: "氬", row: 3, col: 18 },
    // 第四週期
    { num: 19, symbol: "K", name: "鉀", row: 4, col: 1 },
    { num: 20, symbol: "Ca", name: "鈣", row: 4, col: 2 },
    { num: 21, symbol: "Sc", name: "鈧", row: 4, col: 3 },
    { num: 22, symbol: "Ti", name: "鈦", row: 4, col: 4 },
    { num: 23, symbol: "V", name: "釩", row: 4, col: 5 },
    { num: 24, symbol: "Cr", name: "鉻", row: 4, col: 6 },
    { num: 25, symbol: "Mn", name: "錳", row: 4, col: 7 },
    { num: 26, symbol: "Fe", name: "鐵", row: 4, col: 8 },
    { 
    num: 27, 
    symbol: "Co", 
    name: "鈷", 
    row: 4, 
    col: 9, 
    desc: "曾被視為礦坑中的妖怪（Kobold），是維生素 B-12 的核心。",
    link: "periodic-table/article27-cobalt.html", 
    history: "名稱來自德國地下怪物kobold，因常伴隨砷出產而含劇毒，。",
    usage: "用於製造藍色或具有磁性的合金，曾經鋰電池的正極",
    funFact: "鈷也在隕石中被發現過"
},
    { num: 28, symbol: "Ni", name: "鎳", row: 4, col: 10 },
    { num: 29, symbol: "Cu", name: "銅", row: 4, col: 11 },
    { num: 30, symbol: "Zn", name: "鋅", row: 4, col: 12 },
    { num: 31, symbol: "Ga", name: "鎵", row: 4, col: 13 },
    { num: 32, symbol: "Ge", name: "鍺", row: 4, col: 14 },
    { num: 33, symbol: "As", name: "砷", row: 4, col: 15 },
    { num: 34, symbol: "Se", name: "硒", row: 4, col: 16 },
    { num: 35, symbol: "Br", name: "溴", row: 4, col: 17 },
    { num: 36, symbol: "Kr", name: "氪", row: 4, col: 18 },
    // 第五週期
    { num: 37, symbol: "Rb", name: "銣", row: 5, col: 1 },
    { num: 38, symbol: "Sr", name: "鍶", row: 5, col: 2 },
    { num: 39, symbol: "Y", name: "釔", row: 5, col: 3 },
    { num: 40, symbol: "Zr", name: "鋯", row: 5, col: 4 },
    { num: 41, symbol: "Nb", name: "鈮", row: 5, col: 5 },
    { num: 42, symbol: "Mo", name: "鉬", row: 5, col: 6 },
    { num: 43, symbol: "Tc", name: "鎝", row: 5, col: 7 },
    { num: 44, symbol: "Ru", name: "釕", row: 5, col: 8 },
    { num: 45, symbol: "Rh", name: "銠", row: 5, col: 9 },
    { num: 46, symbol: "Pd", name: "鈀", row: 5, col: 10 },
    { num: 47, symbol: "Ag", name: "銀", row: 5, col: 11 },
    { num: 48, symbol: "Cd", name: "鎘", row: 5, col: 12 },
    { num: 49, symbol: "In", name: "銦", row: 5, col: 13 },
    { num: 50, symbol: "Sn", name: "錫", row: 5, col: 14 },
    {num: 51, 
        symbol: "Sb", 
        name: "銻", 
        row: 5, 
        col: 15, 
        desc: "傳說中剋死僧侶的『反僧侶』金屬。", 
        link: "periodic-table/article51-antimony.html",
        history: "名稱來自 'Anti-monk'，因早期鍊金僧侶頻繁中毒而得名。",
        usage: "主要用於阻燃劑、電池合金，古代曾作眼影。",
        funFact: "17世紀流行的『恆久藥丸』排泄後可洗淨重複使用。"
    },
    { num: 52, symbol: "Te", name: "碲", row: 5, col: 16 },
    { num: 53, symbol: "I", name: "碘", row: 5, col: 17 },
    { num: 54, symbol: "Xe", name: "氙", row: 5, col: 18 },
    // 第六週期
    { num: 55, symbol: "Cs", name: "銫", row: 6, col: 1 },
    { num: 56, symbol: "Ba", name: "鋇", row: 6, col: 2 },
    { num: 57, symbol: "La", name: "鑭", row: 9, col: 4 }, // 鑭系開始
    { num: 58, symbol: "Ce", name: "鈰", row: 9, col: 5 },
    { num: 59, symbol: "Pr", name: "鐠", row: 9, col: 6 },
    { num: 60, symbol: "Nd", name: "釹", row: 9, col: 7 },
    { num: 61, symbol: "Pm", name: "鉕", row: 9, col: 8 },
    { num: 62, symbol: "Sm", name: "釤", row: 9, col: 9 },
    { num: 63, symbol: "Eu", name: "銪", row: 9, col: 10 },
    { num: 64, symbol: "Gd", name: "釓", row: 9, col: 11 },
    { num: 65, symbol: "Tb", name: "鋱", row: 9, col: 12 },
    { num: 66, symbol: "Dy", name: "鏑", row: 9, col: 13 },
    { num: 67, symbol: "Ho", name: "鈥", row: 9, col: 14 },
    { num: 68, symbol: "Er", name: "鉺", row: 9, col: 15 },
    { num: 69, symbol: "Tm", name: "銩", row: 9, col: 16 },
    { num: 70, symbol: "Yb", name: "鐿", row: 9, col: 17 },
    { num: 71, symbol: "Lu", name: "鑥", row: 9, col: 18 },
    { num: 72, symbol: "Hf", name: "鉿", row: 6, col: 4 },
    { num: 73, symbol: "Ta", name: "鉭", row: 6, col: 5 },
    { num: 74, symbol: "W", name: "鎢", row: 6, col: 6 },
    { num: 75, symbol: "Re", name: "錸", row: 6, col: 7 },
    { num: 76, symbol: "Os", name: "鋨", row: 6, col: 8 },
    { num: 77, symbol: "Ir", name: "銥", row: 6, col: 9 },
    { num: 78, symbol: "Pt", name: "鉑", row: 6, col: 10 },
    { num: 79, symbol: "Au", name: "金", row: 6, col: 11 },
    { num: 80, symbol: "Hg", name: "汞", row: 6, col: 12 },
    { num: 81, symbol: "Tl", name: "鉈", row: 6, col: 13 },
    { num: 82, symbol: "Pb", name: "鉛", row: 6, col: 14 },
    { num: 83, symbol: "Bi", name: "鉍", row: 6, col: 15 },
    { num: 84, symbol: "Po", name: "釙", row: 6, col: 16 },
    { num: 85, symbol: "At", name: "砈", row: 6, col: 17 },
    { num: 86, symbol: "Rn", name: "氡", row: 6, col: 18 },
    // 第七週期
    { num: 87, symbol: "Fr", name: "鍅", row: 7, col: 1 },
    { num: 88, symbol: "Ra", name: "鐳", row: 7, col: 2 },
    { num: 89, symbol: "Ac", name: "錒", row: 10, col: 4 }, // 錒系開始
    { num: 90, symbol: "Th", name: "釷", row: 10, col: 5 },
    { num: 91, symbol: "Pa", name: "鏷", row: 10, col: 6 },
    { num: 92, symbol: "U", name: "鈾", row: 10, col: 7 },
    { num: 93, symbol: "Np", name: "錼", row: 10, col: 8 },
    { num: 94, symbol: "Pu", name: "鈽", row: 10, col: 9 },
    { num: 95, symbol: "Am", name: "鋂", row: 10, col: 10 },
    { num: 96, symbol: "Cm", name: "鋦", row: 10, col: 11 },
    { num: 97, symbol: "Bk", name: "鉳", row: 10, col: 12 },
    { num: 98, symbol: "Cf", name: "鉲", row: 10, col: 13 },
    { num: 99, symbol: "Es", name: "鑀", row: 10, col: 14 },
    { num: 100, symbol: "Fm", name: "鐨", row: 10, col: 15 },
    { num: 101, symbol: "Md", name: "鍆", row: 10, col: 16 },
    { num: 102, symbol: "No", name: "鍩", row: 10, col: 17 },
    { num: 103, symbol: "Lr", name: "鐒", row: 10, col: 18 },
    { num: 104, symbol: "Rf", name: "鑪", row: 7, col: 4 },
    { num: 105, symbol: "Db", name: "𨧀", row: 7, col: 5 },
    { num: 106, symbol: "Sg", name: "𨭎", row: 7, col: 6 },
    { num: 107, symbol: "Bh", name: "𨨏", row: 7, col: 7 },
    { num: 108, symbol: "Hs", name: "𨭆", row: 7, col: 8 },
    { num: 109, symbol: "Mt", name: "䥑", row: 7, col: 9 },
    { num: 110, symbol: "Ds", name: "鐽", row: 7, col: 10 },
    { num: 111, symbol: "Rg", name: "錀", row: 7, col: 11 },
    { num: 112, symbol: "Cn", name: "鎶", row: 7, col: 12 },
    { num: 113, symbol: "Nh", name: "鉨", row: 7, col: 13 },
    { num: 114, symbol: "Fl", name: "鈇", row: 7, col: 14 },
    { num: 115, symbol: "Mc", name: "鏌", row: 7, col: 15 },
    { num: 116, symbol: "Lv", name: "鉝", row: 7, col: 16 },
    { num: 117, symbol: "Ts", name: "鿬", row: 7, col: 17 },
    { num: 118, symbol: "Og", name: "鿫", row: 7, col: 18 }
];
document.addEventListener('DOMContentLoaded', () => {
    // --- 週期表生成邏輯 ---
    const table = document.getElementById('periodicTable');
    const extraRows = document.getElementById('extraRows');
    const modal = document.getElementById('elementModal');
    if (table) {// 如果目前頁面有週期表容器才執行 (避免在首頁報錯)
        elements.forEach(el => {
            const box = document.createElement('div');
            box.className = 'element-box';
            box.style.gridColumn = el.col;
            box.style.gridRow = el.row;
            box.innerHTML = `
                <div class="element-number">${el.num}</div>
                <div class="element-symbol">${el.symbol}</div>
                <div class="element-name">${el.name}</div>
            `;
            box.onclick = () => showModal(el);// 點擊方塊開啟彈窗
            if (el.row > 7) {
                box.style.gridRow = el.row - 8; 
                box.style.gridColumn = el.col;
                if (extraRows) extraRows.appendChild(box);
            } else {
                box.style.gridRow = el.row;
                box.style.gridColumn = el.col;
                table.appendChild(box);
            }
        });
    }
    function showModal(el) {//彈窗顯示
        if (!modal) return;
        document.getElementById('modalName').innerText = `${el.name} (${el.symbol})`;
        document.getElementById('modalNumber').innerText = el.num;
        // 優先顯示詳細 desc，沒有的話顯示基本訊息
        document.getElementById('modalDesc').innerHTML = `
            ${el.desc || "更多資訊編寫中..."}<br><br>
            ${el.history ? `<strong>歷史：</strong>${el.history}<br>` : ''}
            ${el.usage ? `<strong>用途：</strong>${el.usage}<br>` : ''}
            ${el.funFact ? `<i style="color: #7f8c8d;">冷知識：${el.funFact}</i>` : ''}
        `;
        
        const linkContainer = document.getElementById('modalLinkContainer');
        linkContainer.innerHTML = ''; 
        if (el.link) {
            const btn = document.createElement('a');
            btn.href = el.link;
            btn.className = 'btn';
            btn.style.marginTop = '15px';
            btn.style.display = 'inline-block';
            btn.innerText = '閱讀深度專題 →';
            linkContainer.appendChild(btn);
        }
        modal.style.display = 'flex';
    }
    // --- 關閉彈窗邏輯 ---
    const closeBtn = document.querySelector('.close-button');
    if (closeBtn) {
        closeBtn.onclick = () => modal.style.display = 'none';
    }
    window.onclick = (e) => { if (e.target == modal) modal.style.display = 'none'; };
    // --- 3. 特殊互動區：首頁卡片提示 (放在最下面) ---
    const antimonyCardLink = document.querySelector('a[href="article-antimony.html"]');
    if (antimonyCardLink) {
        const antimonyCard = antimonyCardLink.parentElement;
        antimonyCard.addEventListener('mouseenter', () => {
            console.log("提示：你知道銻丸曾被重複使用了幾代人嗎？");
        });
    }
});
const elementData = {
    "Sb": {
        name: "銻 (Antimony)",
        number: 51,
        history: "名稱傳說來自 'Anti-monk' (反僧侶)，因早期鍊金僧侶頻繁接觸導致中毒死亡而得名。",
        usage: "主要用於塑膠阻燃劑、鉛酸電池強化合金，古代曾被當作眼影使用。",
        funFact: "雖然它有毒，但在17世紀曾流行過『銻丸』，這種藥丸排泄後可以洗淨回收重複使用，被稱為恆久藥丸。"
    }
};
document.addEventListener('DOMContentLoaded', () => {
    const antimonyCard = document.querySelector('a[href="article-antimony.html"]').parentElement;
    antimonyCard.addEventListener('mouseenter', () => {
        console.log("提示：你知道銻丸曾被重複使用了幾代人嗎？");
    });
});