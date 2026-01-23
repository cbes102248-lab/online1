// script.js
document.addEventListener('DOMContentLoaded', () => {
    const musicTrigger = document.getElementById('backmusic');
    const audio = document.getElementById('bgMusic');
    // 檢查瀏覽器是否已啟動音樂
    if (sessionStorage.getItem('hasStartedMusic') === 'true') {
        if (musicTrigger) musicTrigger.style.display = 'none'; // 直接隱藏，不要有動畫
        audio.play().catch(() => { // 嘗試播放
            // 如果失敗，通常是因為新頁面也需要一次互動，這時可以讓遮罩再次出現
            if (musicTrigger) musicTrigger.style.display = 'flex';
        });
    }
    if (musicTrigger && audio) { // 兩個函式同時為true觸發
        musicTrigger.addEventListener('click', () => { // 確認點擊
            audio.play(); // 音樂播放
            musicTrigger.classList.add('hidden'); // 隱藏此函式
            sessionStorage.setItem('hasStartedMusic', 'true'); // 標記使用者點過了
        });
    }
});
document.addEventListener('click', (e) => { // 確認點擊函式
    const link = e.target.closest('a');
    // 確保同網域、非錨點、非新視窗才攔截
    if (link && link.hostname === window.location.hostname && !link.hash && link.target !== "_blank") {
        e.preventDefault(); // 阻止瀏覽器跳轉
        const url = link.href;
        loadPage(url); // 執行自定義載入函式
    }
});
function loadPage(url) { // 自定義網頁載入函式
    fetch(url)
        .then(response => response.text())
        .then(html => {
            const parser = new DOMParser();
            const doc = parser.parseFromString(html, 'text/html');
            const newMain = doc.getElementById('main-content');
            const currentMain = document.getElementById('main-content');
            const heroHeader = document.querySelector('.hero'); // 抓取首頁大標頭
            if (newMain && currentMain) {
                // 處理 Header首頁顯示分頁隱藏
                // 增加對 /online1/index.html 的判斷
                if (url.includes('index.html') || url.endsWith('/') || url.endsWith('online1/')) {
                    if (heroHeader) heroHeader.style.display = 'block';
                } else {
                    if (heroHeader) heroHeader.style.display = 'none';
                }
                currentMain.className = newMain.className; // 替換main內容
                currentMain.innerHTML = newMain.innerHTML;
                history.pushState({ path: url }, '', url); // 更新URL與標題
                document.title = doc.title;
                setTimeout(reInitPageScripts, 50); // 延遲執行腳本初始化(預留載入時間50ms)
                window.scrollTo(0, 0);
            }
        })
        .catch(err => {
            console.error("載入失敗:", err);
            window.location.href = url; // 失敗時的保險機制
        });
}
function initPeriodicTable() { // 週期表生成函式
    const table = document.getElementById('periodicTable');
    const extraRows = document.getElementById('extraRows');
    const modal = document.getElementById('elementModal');
    // 獲取目前的顯示模式 (由 HTML 的 select 決定)
    const modeSelect = document.getElementById('viewMode');
    const mode = modeSelect ? modeSelect.value : 'standard';
    if (!table) return; // 如果這頁沒週期表，就不跑
    console.log(`偵測到週期表容器，模式：${mode}，開始渲染...`);
    table.innerHTML = ''; // 清空舊的，避免重複生成
    if (extraRows) extraRows.innerHTML = '';
    elements.forEach(el => {
        const box = document.createElement('div');
        box.className = 'element-box';
        // --- 根據顯示模式決定樣式與顯示內容 ---
        let displayText = el.name; // 預設顯示中文名稱
        let bgColor = '';
        let bgImg = '';
        if (mode === 'standard') {
            box.classList.add(el.category); // 套用元素性質分類顏色 (如 nonmetal)
            // 標準模式下顯示原子圖片
            if (el.image) {
                bgImg = `linear-gradient(rgba(255,255,255,0.7), rgba(255,255,255,0.7)), url('${el.image}')`;
            }
        } 
        else if (mode === 'electronegativity') {
            displayText = el.electronegativity || 'N/A';
            // 電負度越高顏色越深 (以橘紅色系為基礎)
            const alpha = el.electronegativity ? (el.electronegativity / 4.5) : 0.1;
            bgColor = `rgba(255, 87, 34, ${alpha})`;
        } 
        else if (mode === 'flame') {
            // 如果資料中有焰色，則強化邊框與陰影
            if (el.flame && el.flame !== 'none') {
                box.style.boxShadow = `inset 0 0 15px ${el.flame}`;
                box.style.borderColor = el.flame;
            }
        } 
        else if (mode === 'year') {
            displayText = el.year <= 0 ? '古代' : el.year;
            // 越近代發現的元素顏色越淺，古代最深
            const brightness = el.year <= 0 ? 80 : Math.max(30, 100 - (2026 - el.year) / 10);
            bgColor = `hsl(200, 70%, ${brightness}%)`;
        }
        // 套用動態計算的樣式
        if (bgColor) box.style.backgroundColor = bgColor;
        if (bgImg) {
            box.style.backgroundImage = bgImg;
            box.style.backgroundSize = 'cover';
        }
        // 設定網格位置
        if (el.row > 7) {
            box.style.gridRow = el.row - 8; 
            box.style.gridColumn = el.col;
            if (extraRows) extraRows.appendChild(box);
        } else {
            box.style.gridRow = el.row;
            box.style.gridColumn = el.col;
            if (table) table.appendChild(box);
        }
        box.innerHTML = `
            <div class="element-number">${el.num}</div>
            <div class="element-symbol">${el.symbol}</div>
            <div class="element-name" style="font-size: 0.7rem;">${displayText}</div>
        `;
        box.onclick = () => showModal(el); // 當點擊呼叫框框
    });
    const closeBtn = document.querySelector('.close-button'); // 重新綁定關閉按鈕事件
    if (closeBtn && modal) {
        closeBtn.onclick = () => modal.style.display = 'none';
        window.onclick = (e) => { if (e.target == modal) modal.style.display = 'none'; };
    }
}
function reInitPageScripts() {
    console.log("執行換頁後的腳本初始化...");
    initPeriodicTable(); // 重新生成週期表
    // --- PyScript 支援區塊 ---
    const hasPyScript = document.querySelector('script[type="py"]'); // 檢查新內容是否有 Python 腳本
    if (hasPyScript && window.pyscript) {
        console.log("偵測到 Python 分析器，PyScript 正在待命...");
    }
    // 重新綁定首頁卡片提示
    const antimonyLink = document.querySelector('a[href="article-antimony.html"]');
    if (antimonyLink) {
        const antimonyCard = antimonyLink.parentElement;
        antimonyCard.addEventListener('mouseenter', () => {
            console.log("提示：你知道銻丸曾被重複使用了幾代人嗎？");
        });
    }
}
document.addEventListener('DOMContentLoaded', () => {
    reInitPageScripts(); // 初次進入首頁或週期表頁時執行一次腳本初始化
    try {
        if (typeof initMusicControl === 'function') {
            initMusicControl();
        }
    } catch (e) {
        console.log("音樂控制初始化提示:", e);
    }
});
const elements = [//週期表表格位置與彈窗訊息
    // 第一週期
    { 
        num: 1, symbol: "H", name: "氫", row: 1, col: 1,
        category: "nonmetal" , electronegativity: 2.20, year: 1766, flame: "none",image: "images/atom/H.jpg",
     },
    {
        num: 2, symbol: "He", name: "氦", row: 1, col: 18 ,
        category: "noble-gas", electronegativity: 0, year: 1868, flame: "none",image: "images/atom/He.jpg",
    },
    // 第二週期
    { 
        num: 3, symbol: "Li", name: "鋰", row: 2, col: 1 ,
        category:  "alkali-metal", electronegativity: 0.98, year: 1817, flame:(220,20,60) ,image: "images/atom/Li.jpg",
    },
    { 
        num: 4, symbol: "Be", name: "鈹", row: 2, col: 2 ,
        category: "alkaline-earth-metal", electronegativity: 1.57, year: 1798, flame: "none",image: "images/atom/Be.jpg",

    },
    { 
        num: 5, symbol: "B", name: "硼", row: 2, col: 13 ,
        category: "metalloid", electronegativity: 2.04, year: 1808, flame:"none" ,image: "images/atom/B.jpg",

    },
    { 
        num: 6, symbol: "C", name: "碳", row: 2, col: 14 ,
        category: "nonmetal", electronegativity: 2.55, year: 0, flame:"none" ,image: "images/atom/C.jpg",

    },
    { 
        num: 7, symbol: "N", name: "氮", row: 2, col: 15 ,
        category: "nonmetal", electronegativity: 3.04, year: 1772, flame:"none" ,image: "images/atom/N.jpg",

    },
    { 
        num: 8, symbol: "O", name: "氧", row: 2, col: 16 ,
        category: "nonmetal", electronegativity: 3.44, year: 1774, flame:"none" ,image: "images/atom/O.jpg",

    },
    { 
        num: 9, symbol: "F", name: "氟", row: 2, col: 17 ,
        category: "nonmetal", electronegativity: 3.98, year: 1886, flame:"none" ,image: "images/atom/F.jpg",

    },
    { 
        num: 10, symbol: "Ne", name: "氖", row: 2, col: 18 ,
        category: "noble-gas", electronegativity: 0, year: 1898, flame:"none" ,image: "images/atom/Ne.jpg",

    },
    // 第三週期
    { 
        num: 11, symbol: "Na", name: "鈉", row: 3, col: 1 ,
        category: "alkali-metal", electronegativity: 0.93, year: 1807, flame: (255,215,0),image: "images/atom/Na.jpg",

    },
    { 
        num: 12, symbol: "Mg", name: "鎂", row: 3, col: 2 ,
        category: "alkaline-earth-metal", electronegativity: 1.31, year: 1808, flame:"none" ,image: "images/atom/Mg.jpg",

    },
    { 
        num: 13, symbol: "Al", name: "鋁", row: 3, col: 13 ,
        category: "post-transition-metal", electronegativity: 1.61, year: 1825, flame:"none" ,image: "images/atom/Al.jpg",

    },
    { 
        num: 14, symbol: "Si", name: "矽", row: 3, col: 14 ,
        category: "metalloid", electronegativity: 1.90, year: 1824, flame:"none" ,image: "images/atom/Si.jpg",

    },
    { 
        num: 15, symbol: "P", name: "磷", row: 3, col: 15 ,
        category: "nonmetal", electronegativity: 2.19, year: 1669, flame:"none" ,image: "images/atom/P.jpg",

    },
    { 
        num: 16, symbol: "S", name: "硫", row: 3, col: 16 ,
        category: "nonmetal", electronegativity: 2.58, year: 0, flame:"none" ,image: "images/atom/S.jpg",

    },
    { 
        num: 17, symbol: "Cl", name: "氯", row: 3, col: 17 ,
        category: "nonmetal", electronegativity: 3.16, year: 1774, flame:"none" ,image: "images/atom/Cl.jpg",

    },
    { 
        num: 18, symbol: "Ar", name: "氬", row: 3, col: 18 ,
        category: "noble-gas", electronegativity: 0, year: 1894, flame:"none" ,image: "images/atom/Ar.jpg",

    },
    // 第四週期
    { 
        num: 19, symbol: "K", name: "鉀", row: 4, col: 1 ,
        category: "alkali-metal", electronegativity: 0.82, year: 1807, flame:(216,191,216) ,image: "images/atom/K.jpg",

    },
    { 
        num: 20, symbol: "Ca", name: "鈣", row: 4, col: 2 ,
        category: "alkaline-earth-metal", electronegativity: 1.00, year: 1808, flame:(255,69,0) ,image: "images/atom/Ca.jpg",

    },
    { 
        num: 21, symbol: "Sc", name: "鈧", row: 4, col: 3 ,
        category: "transition-metal", electronegativity: 1.36, year: 1879, flame:"none" ,image: "images/atom/Sc.jpg",

    },
    { 
        num: 22, symbol: "Ti", name: "鈦", row: 4, col: 4 ,
        category: "transition-metal", electronegativity: 1.54, year: 1791, flame:"none" ,image: "images/atom/Ti.jpg",

    },
    { 
        num: 23, symbol: "V", name: "釩", row: 4, col: 5 ,
        category: "transition-metal", electronegativity: 1.63, year: 1801, flame:"none" ,image: "images/atom/V.jpg",

    },
    { 
        num: 24, symbol: "Cr", name: "鉻", row: 4, col: 6 ,
        category: "transition-metal", electronegativity: 1.66, year: 1797, flame:"none" ,image: "images/atom/Cr.jpg",

    },
    { 
        num: 25, symbol: "Mn", name: "錳", row: 4, col: 7 ,
        category: "transition-metal", electronegativity: 1.55, year: 1774, flame:"none" ,image: "images/atom/Mn.jpg",

    },
    { 
        num: 26, symbol: "Fe", name: "鐵", row: 4, col: 8 ,
        category: "transition-metal", electronegativity: 1.83, year: 0, flame: (255,165,0),image: "images/atom/Fe.jpg",

    },
    { 
        num: 27, symbol: "Co", name: "鈷", row: 4, col: 9 ,
        category: "transition-metal", electronegativity: 1.88, year: 1735, flame:"none" ,image: "images/atom/Co.jpg",
        desc: "曾被視為礦坑中的妖怪（Kobold），是維生素 B-12 的核心。",
        link: "periodic-table/article27-cobalt.html", 
        history: "名稱來自德國地下怪物kobold，因常伴隨砷出產而含劇毒，。",
        usage: "用於製造藍色或具有磁性的合金，曾經鋰電池的正極",
        funFact: "鈷也在隕石中被發現過"
    },
    { 
        num: 28, symbol: "Ni", name: "鎳", row: 4, col: 10 ,
        category: "transition-metal", electronegativity: 1.91, year: 1751, flame:"none" ,image: "images/atom/Ni.jpg",

    },
    {
        num: 29, symbol: "Cu", name: "銅", row: 4, col: 11 ,
        category: "transition-metal", electronegativity: 1.90, year: 0, flame:(0,255,255) ,image: "images/atom/Cu.jpg",

    },
    { 
        num: 30, symbol: "Zn", name: "鋅", row: 4, col: 12 ,
        category: "transition-metal", electronegativity: 1.65, year: 1746, flame:"none" ,image: "images/atom/Zn.jpg",

    },
    { 
        num: 31, symbol: "Ga", name: "鎵", row: 4, col: 13 ,
        category: "post-transition-metal", electronegativity: 1.81, year: 1875, flame:"none" ,image: "images/atom/Ga.jpg",

    },
    { 
        num: 32, symbol: "Ge", name: "鍺", row: 4, col: 14 ,
        category: "metalloid", electronegativity: 2.01, year: 1886, flame:"none" ,image: "images/atom/Ge.jpg",

    },
    { 
        num: 33, symbol: "As", name: "砷", row: 4, col: 15 ,
        category: "metalloid", electronegativity: 2.18, year: 0, flame:"none" ,image: "images/atom/As.jpg",

    },
    { 
        num: 34, symbol: "Se", name: "硒", row: 4, col: 16 ,
        category: "nonmetal", electronegativity: 2.55, year: 1817, flame:"none" ,image: "images/atom/Se.jpg",

    },
    { 
        num: 35, symbol: "Br", name: "溴", row: 4, col: 17 ,
        category: "nonmetal", electronegativity: 2.96, year: 1826, flame:"none" ,image: "images/atom/Br.jpg",

    },
    { 
        num: 36, symbol: "Kr", name: "氪", row: 4, col: 18 ,
        category: "noble-gas", electronegativity: 3.00, year: 1898, flame:"none" ,image: "images/atom/Kr.jpg",

    },
    // 第五週期
    { 
        num: 37, symbol: "Rb", name: "銣", row: 5, col: 1 ,
        category: "alkali-metal", electronegativity: 0.82, year: 1861, flame: (148,0,211) ,image: "images/atom/Rb.jpg",

    },
    { 
        num: 38, symbol: "Sr", name: "鍶", row: 5, col: 2 ,
        category: "alkaline-earth-metal", electronegativity: 0.95, year: 1790, flame:(255,0,0) ,image: "images/atom/Sr.jpg",

    },
    { 
        num: 39, symbol: "Y", name: "釔", row: 5, col: 3 ,
        category: "transition-metal", electronegativity: 1.22, year: 1794, flame:"none" ,image: "images/atom/Y.jpg",

    },
    { 
        num: 40, symbol: "Zr", name: "鋯", row: 5, col: 4 ,
        category: "transition-metal", electronegativity: 1.33, year: 1789, flame:"none" ,image: "images/atom/Zr.jpg",

    },
    { 
        num: 41, symbol: "Nb", name: "鈮", row: 5, col: 5 ,
        category: "transition-metal", electronegativity: 1.60, year: 1801, flame:"none" ,image: "images/atom/Nb.jpg",

    },
    { 
        num: 42, symbol: "Mo", name: "鉬", row: 5, col: 6 ,
        category: "transition-metal", electronegativity: 2.16, year: 1778, flame:"none" ,image: "images/atom/Mo.jpg",

    },
    { 
        num: 43, symbol: "Tc", name: "鎝", row: 5, col: 7 ,
        category: "transition-metal", electronegativity: 1.90, year: 1937, flame:"none" ,image: "images/atom/Tc.jpg",

    },
    { 
        num: 44, symbol: "Ru", name: "釕", row: 5, col: 8 ,
        category: "transition-metal", electronegativity: 2.20, year: 1844, flame:"none" ,image: "images/atom/Ru.jpg",

    },
    { 
        num: 45, symbol: "Rh", name: "銠", row: 5, col: 9 ,
        category: "transition-metal", electronegativity: 2.28, year: 1803, flame:"none" ,image: "images/atom/Rh.jpg",

    },
    { 
        num: 46, symbol: "Pd", name: "鈀", row: 5, col: 10 ,
        category: "transition-metal", electronegativity: 2.20, year: 1803, flame:"none" ,image: "images/atom/Pd.jpg",

    },
    { 
        num: 47, symbol: "Ag", name: "銀", row: 5, col: 11 ,
        category: "transition-metal", electronegativity: 1.93, year: 0, flame:"none" ,image: "images/atom/Ag.jpg",

    },
    { 
        num: 48, symbol: "Cd", name: "鎘", row: 5, col: 12 ,
        category: "transition-metal", electronegativity: 1.69, year: 1817, flame:"none" ,image: "images/atom/Cd.jpg",

    },
    { 
        num: 49, symbol: "In", name: "銦", row: 5, col: 13 ,
        category: "post-transition-metal", electronegativity: 1.78, year: 1863, flame:(75,0,130) ,image: "images/atom/In.jpg",

    },
    { 
        num: 50, symbol: "Sn", name: "錫", row: 5, col: 14 ,
        category: "post-transition-metal", electronegativity: 1.96, year: 0, flame:"none" ,image: "images/atom/Sn.jpg",

    },
    {
        num: 51, symbol: "Sb", name: "銻", row: 5, col: 15 ,
        category: "metalloid", electronegativity: 2.05, year: 0, flame:"none" ,image: "images/atom/Sb.jpg",
        desc: "傳說中剋死僧侶的『反僧侶』金屬。", 
        link: "periodic-table/article51-antimony.html",
        history: "名稱來自 'Anti-monk'，因早期鍊金僧侶頻繁中毒而得名。",
        usage: "主要用於阻燃劑、電池合金，古代曾作眼影。",
        funFact: "17世紀流行的『恆久藥丸』排泄後可洗淨重複使用。"
    },
    { 
        num: 52, symbol: "Te", name: "碲", row: 5, col: 16 ,
        category: "metalloid", electronegativity: 2.10, year: 1782, flame:"none" ,image: "images/atom/Te.jpg",

    },
    { 
        num: 53, symbol: "I", name: "碘", row: 5, col: 17 ,
        category: "nonmetal", electronegativity: 2.66, year: 1811, flame:"none" ,image: "images/atom/I.jpg",

    },
    { 
        num: 54, symbol: "Xe", name: "氙", row: 5, col: 18 ,
        category: "noble-gas", electronegativity: 2.60, year: 1898, flame:"none" ,image: "images/atom/Xe.jpg",

    },
    // 第六週期
    { 
        num: 55, symbol: "Cs", name: "銫", row: 6, col: 1 ,
        category: "alkali-metal", electronegativity: 0.79, year: 1860, flame: (0,0,255) ,image: "images/atom/Cs.jpg",

    },
    { 
        num: 56, symbol: "Ba", name: "鋇", row: 6, col: 2 ,
        category: "alkaline-earth-metal", electronegativity: 0.89, year: 1808, flame: (144,238,144) ,image: "images/atom/Ba.jpg",

    },
    { 
        num: 57, symbol: "La", name: "鑭", row: 9, col: 4 ,
        category: "lanthanide", electronegativity: 1.10, year: 1839, flame:"none" ,image: "images/atom/La.jpg",

    }, // 鑭系開始
    { 
        num: 58, symbol: "Ce", name: "鈰", row: 9, col: 5 ,
        category: "lanthanide", electronegativity: 1.12, year: 1803, flame:"none" ,image: "images/atom/Ce.jpg",

    },
    { 
        num: 59, symbol: "Pr", name: "鐠", row: 9, col: 6 ,
        category: "lanthanide", electronegativity: 1.13, year: 1885, flame:"none" ,image: "images/atom/Pr.jpg",

    },
    { 
        num: 60, symbol: "Nd", name: "釹", row: 9, col: 7 ,
        category: "lanthanide", electronegativity: 1.14, year: 1885, flame:"none" ,image: "images/atom/Nd.jpg",

    },
    { 
        num: 61, symbol: "Pm", name: "鉕", row: 9, col: 8 ,
        category: "lanthanide", electronegativity: 1.13, year: 1945, flame:"none" ,image: "images/atom/Pm.jpg",

    },
    {
        num: 62, symbol: "Sm", name: "釤", row: 9, col: 9 ,
        category: "lanthanide", electronegativity: 1.17, year: 1879, flame:"none" ,image: "images/atom/Sm.jpg",

    },
    { 
        num: 63, symbol: "Eu", name: "銪", row: 9, col: 10 ,
        category: "lanthanide", electronegativity: 1.20, year: 1901, flame:"none" ,image: "images/atom/Eu.jpg",

    },
    { 
        num: 64, symbol: "Gd", name: "釓", row: 9, col: 11 ,
        category: "lanthanide", electronegativity: 1.20, year: 1880, flame:"none" ,image: "images/atom/Gd.jpg",

    },
    { 
        num: 65, symbol: "Tb", name: "鋱", row: 9, col: 12 ,
        category: "lanthanide", electronegativity: 1.10, year: 1843, flame:"none" ,image: "images/atom/Tb.jpg",

    },
    { 
        num: 66, symbol: "Dy", name: "鏑", row: 9, col: 13 ,
        category: "lanthanide", electronegativity: 1.22, year: 1886, flame:"none" ,image: "images/atom/Dy.jpg",

    },
    { 
        num: 67, symbol: "Ho", name: "鈥", row: 9, col: 14 ,
        category: "lanthanide", electronegativity: 1.23, year: 1878, flame:"none" ,image: "images/atom/Ho.jpg",

    },
    { 
        num: 68, symbol: "Er", name: "鉺", row: 9, col: 15 ,
        category: "lanthanide", electronegativity: 1.24, year: 1842, flame:"none" ,image: "images/atom/Er.jpg",

    },
    { 
        num: 69, symbol: "Tm", name: "銩", row: 9, col: 16 ,
        category: "lanthanide", electronegativity: 1.25, year: 1879, flame:"none" ,image: "images/atom/Tm.jpg",

    },
    { 
        num: 70, symbol: "Yb", name: "鐿", row: 9, col: 17 ,
        category: "lanthanide", electronegativity: 1.10, year: 1878, flame:"none" ,image: "images/atom/Yb.jpg",

    },
    { 
        num: 71, symbol: "Lu", name: "鑥", row: 9, col: 18 ,
        category: "lanthanide", electronegativity: 1.27, year: 1907, flame:"none" ,image: "images/atom/Lu.jpg",

    },
    { 
        num: 72, symbol: "Hf", name: "鉿", row: 6, col: 4 ,
        category: "transition-metal", electronegativity: 1.30, year: 1923, flame:"none" ,image: "images/atom/Hf.jpg",

    },
    { 
        num: 73, symbol: "Ta", name: "鉭", row: 6, col: 5 ,
        category: "transition-metal", electronegativity: 1.50, year:1802 , flame:"none" ,image: "images/atom/Ta.jpg",

    },
    { 
        num: 74, symbol: "W", name: "鎢", row: 6, col: 6 ,
        category: "transition-metal", electronegativity: 2.36, year: 1783, flame:"none" ,image: "images/atom/W.jpg",

    },
    { 
        num: 75, symbol: "Re", name: "錸", row: 6, col: 7 ,
        category: "transition-metal", electronegativity: 1.90, year: 1925, flame:"none" ,image: "images/atom/Re.jpg",

    },
    { 
        num: 76, symbol: "Os", name: "鋨", row: 6, col: 8 ,
        category: "transition-metal", electronegativity: 2.20, year: 1803, flame:"none" ,image: "images/atom/Os.jpg",

    },
    { 
        num: 77, symbol: "Ir", name: "銥", row: 6, col: 9 ,
        category: "transition-metal", electronegativity: 2.20, year: 1803, flame:"none" ,image: "images/atom/Ir.jpg",

    },
    { 
        num: 78, symbol: "Pt", name: "鉑", row: 6, col: 10 ,
        category: "transition-metal", electronegativity: 2.28, year: 1735, flame:"none" ,image: "images/atom/Pt.jpg",

    },
    { 
        num: 79, symbol: "Au", name: "金", row: 6, col: 11 ,
        category: "transition-metal", electronegativity: 2.54, year: 0, flame:"none" ,image: "images/atom/Au.jpg",

    },
    { 
        num: 80, symbol: "Hg", name: "汞", row: 6, col: 12 ,
        category: "transition-metal", electronegativity: 2.00, year: 0, flame:"none" ,image: "images/atom/Hg.jpg",

    },
    { 
        num: 81, symbol: "Tl", name: "鉈", row: 6, col: 13 ,
        category: "post-transition-metal", electronegativity: 1.62, year: 1861, flame:(0,255,0) ,image: "images/atom/Tl.jpg",

    },
    { 
        num: 82, symbol: "Pb", name: "鉛", row: 6, col: 14 ,
        category: "post-transition-metal", electronegativity: 2.33, year: 0, flame:"none" ,image: "images/atom/Pb.jpg",

    },
    { 
        num: 83, symbol: "Bi", name: "鉍", row: 6, col: 15 ,
        category: "post-transition-metal", electronegativity: 2.02, year: 0, flame:"none" ,image: "images/atom/Bi.jpg",

    },
    { 
        num: 84, symbol: "Po", name: "釙", row: 6, col: 16 ,
        category: "metalloid", electronegativity: 2.00, year: 1898, flame:"none" ,image: "images/atom/Po.jpg",

    },
    { 
        num: 85, symbol: "At", name: "砈", row: 6, col: 17 ,
        category: "nonmetal", electronegativity: 2.20, year: 1940, flame:"none" ,image: "images/atom/At.jpg",

    },
    { 
        num: 86, symbol: "Rn", name: "氡", row: 6, col: 18 ,
        category: "noble-gas", electronegativity: 0, year: 1900, flame:"none" ,image: "images/atom/Rn.jpg",

    },
    // 第七週期
    { 
        num: 87, symbol: "Fr", name: "鍅", row: 7, col: 1 ,
        category: "alkali-metal", electronegativity: 0.70, year: 1939, flame:"none" ,image: "images/atom/Fr.jpg",

    },
    { 
        num: 88, symbol: "Ra", name: "鐳", row: 7, col: 2 ,
        category: "alkaline-earth-metal", electronegativity: 0.90, year: 1898, flame:(255,0,0) ,image: "images/atom/Ra.jpg",

    },
    { 
        num: 89, symbol: "Ac", name: "錒", row: 10, col: 4 ,
        category: "actinide", electronegativity: 1.10, year: 1899, flame:"none" ,image: "images/atom/Ac.jpg",

    }, // 錒系開始
    { 
        num: 90, symbol: "Th", name: "釷", row: 10, col: 5 ,
        category: "actinide", electronegativity: 1.30, year: 1828, flame:"none" ,image: "images/atom/Th.jpg",

    },
    { 
        num: 91, symbol: "Pa", name: "鏷", row: 10, col: 6 ,
        category: "actinide", electronegativity: 1.50, year: 1913, flame:"none" ,image: "images/atom/Pa.jpg",

    },
    { 
        num: 92, symbol: "U", name: "鈾", row: 10, col: 7 ,
        category: "actinide", electronegativity: 1.38, year: 1789, flame:"none" ,image: "images/atom/U.jpg",

    },
    { 
        num: 93, symbol: "Np", name: "錼", row: 10, col: 8 ,
        category: "actinide", electronegativity: 1.36, year: 1940, flame:"none" ,image: "images/atom/Np.jpg",

    },
    { 
        num: 94, symbol: "Pu", name: "鈽", row: 10, col: 9 ,
        category: "actinide", electronegativity: 1.28, year: 1940, flame:"none" ,image: "images/atom/Pu.jpg",

    },
    { 
        num: 95, symbol: "Am", name: "鋂", row: 10, col: 10 ,
        category: "actinide", electronegativity: 1.13, year: 1944, flame:"none" ,image: "images/atom/Am.jpg",

    },
    { 
        num: 96, symbol: "Cm", name: "鋦", row: 10, col: 11 ,
        category: "actinide", electronegativity: 1.28, year: 1944, flame:"none" ,image: "images/atom/Cm.jpg",

    },
    { 
        num: 97, symbol: "Bk", name: "鉳", row: 10, col: 12 ,
        category: "actinide", electronegativity: 1.30, year: 1949, flame:"none" ,image: "images/atom/Bk.jpg",

    },
    { 
        num: 98, symbol: "Cf", name: "鉲", row: 10, col: 13 ,
        category: "actinide", electronegativity: 1.30, year: 1950, flame:"none" ,image:"images/atom/Cf.jpg" ,

    },
    { 
        num: 99, symbol: "Es", name: "鑀", row: 10, col: 14 ,
        category: "actinide", electronegativity: 1.30, year: 1952, flame:"none" ,image: "images/atom/Es.jpg",

    },
    { 
        num: 100, symbol: "Fm", name: "鐨", row: 10, col: 15 ,
        category: "actinide", electronegativity: 1.30, year: 1953, flame:"none" ,image: "images/atom/Fm.jpg",

    },
    { 
        num: 101, symbol: "Md", name: "鍆", row: 10, col: 16 ,
        category: "actinide", electronegativity: 1.30, year: 1955, flame:"none" ,image: "images/atom/Md.jpg",

    },
    { 
        num: 102, symbol: "No", name: "鍩", row: 10, col: 17 ,
        category: "actinide", electronegativity: 1.30, year: 1958, flame:"none" ,image: "images/atom/No.jpg",

    },
    { 
        num: 103, symbol: "Lr", name: "鐒", row: 10, col: 18 ,
        category: "actinide", electronegativity: 1.30, year: 1961, flame:"none" ,image: "images/atom/Lr.jpg",

    },
    { 
        num: 104, symbol: "Rf", name: "鑪", row: 7, col: 4 ,
        category: "transition-metal", electronegativity: 0, year: 1964, flame:"none" ,image: "images/atom/Rf.jpg",

    },
    { 
        num: 105, symbol: "Db", name: "𨧀", row: 7, col: 5 ,
        category: "transition-metal", electronegativity: 0, year: 1967, flame:"none" ,image: "images/atom/Db.jpg",

    },
    { 
        num: 106, symbol: "Sg", name: "𨭎", row: 7, col: 6 ,
        category: "transition-metal", electronegativity: 0, year: 1974, flame:"none" ,image: "images/atom/Sg.jpg",

    },
    { 
        num: 107, symbol: "Bh", name: "𨨏", row: 7, col: 7 ,
        category: "transition-metal", electronegativity: 0, year: 1981, flame:"none" ,image: "images/atom/Bh.jpg",

    },
    { 
        num: 108, symbol: "Hs", name: "𨭆", row: 7, col: 8 ,
        category: "transition-metal", electronegativity: 0, year: 1984, flame:"none" ,image: "images/atom/Hs.jpg",

    },
    { 
        num: 109, symbol: "Mt", name: "䥑", row: 7, col: 9 ,
        category: "transition-metal", electronegativity: 0, year: 1982, flame:"none" ,image: "images/atom/Mt.jpg",

    },
    { 
        num: 110, symbol: "Ds", name: "鐽", row: 7, col: 10 ,
        category: "transition-metal", electronegativity: 0, year: 1994, flame:"none" ,image: "images/atom/Ds.jpg",

    },
    { 
        num: 111, symbol: "Rg", name: "錀", row: 7, col: 11 ,
        category: "transition-metal", electronegativity: 0, year: 1994, flame:"none" ,image: "images/atom/Rg.jpg",

    },
    { 
        num: 112, symbol: "Cn", name: "鎶", row: 7, col: 12 ,
        category: "transition-metal", electronegativity: 0, year: 1996, flame:"none" ,image: "images/atom/Cn.jpg",

    },
    { 
        num: 113, symbol: "Nh", name: "鉨", row: 7, col: 13 ,
        category: "post-transition-metal", electronegativity: 0, year: 2004, flame:"none" ,image: "images/atom/Nh.jpg",

    },
    { 
        num: 114, symbol: "Fl", name: "鈇", row: 7, col: 14 ,
        category: "post-transition-metal", electronegativity: 0, year: 1999, flame:"none" ,image: "images/atom/Fl.jpg",

    },
    { 
        num: 115, symbol: "Mc", name: "鏌", row: 7, col: 15 ,
        category: "post-transition-metal", electronegativity: 0, year: 2003, flame:"none" ,image: "images/atom/Mc.jpg",

    },
    { 
        num: 116, symbol: "Lv", name: "鉝", row: 7, col: 16 ,
        category: "post-transition-metal", electronegativity: 0, year: 2000, flame:"none" ,image: "images/atom/Lv.jpg",

    },
    { 
        num: 117, symbol: "Ts", name: "鿬", row: 7, col: 17 ,
        category: "nonmetal", electronegativity: 2.0, year: 2010, flame:"none" ,image: "images/atom/Ts.jpg",

    },
    { 
        num: 118, symbol: "Og", name: "鿫", row: 7, col: 18 ,
        category: "noble-gas", electronegativity: 0, year: 2006, flame:"none" ,image: "images/atom/Og.jpg",

    }
];
// 彈窗函式
function showModal(el) {
    const modal = document.getElementById('elementModal');
    if (!modal) return;
    
    document.getElementById('modalName').innerText = `${el.name} (${el.symbol})`;
    document.getElementById('modalNumber').innerText = el.num;
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
        btn.innerText = '閱讀深度專題 →';
        linkContainer.appendChild(btn);
    }
    modal.style.display = 'flex';
}