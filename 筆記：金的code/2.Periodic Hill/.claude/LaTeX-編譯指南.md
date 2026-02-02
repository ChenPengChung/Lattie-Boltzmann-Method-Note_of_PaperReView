# LaTeX 編譯指南

## 快速修正常見問題

### 1. 字體錯誤修正

**問題：** `fontspec package requires XeTeX or LuaTeX`

**解決：** 使用 `xelatex` 而非 `pdflatex`
```bash
xelatex your_file.tex
```

### 2. 中文字體缺失（DFKai-SB not found）

**問題：** macOS 系統上缺少 Windows 字體 DFKai-SB

**解決方案 1（推薦）：** 修改 .tex 檔案，將字體改為 macOS 內建字體

在你的 .tex 檔案中找到類似這樣的行：
```latex
\setCJKmonofont{DFKai-SB}
```

改為：
```latex
\setCJKmonofont{Kaiti TC}
```

**解決方案 2：** 使用其他 macOS 字體
- `Kaiti TC`（楷體）
- `Songti TC`（宋體）
- `Heiti TC`（黑體）

### 3. 等寬字體無法顯示中文

**問題：** 程式碼區塊中的中文變成空白或方框

**原因：** `lmmono10-regular` (Latin Modern Mono) 不支援中文字符

**解決：** 確保正確設定 CJK 等寬字體
```latex
\setCJKmonofont{Kaiti TC}  % macOS
% 或
\setCJKmonofont{DFKai-SB}  % Windows
```

### 4. 交叉引用未定義（Undefined references）

**問題：** `Warning: There were undefined references`

**解決：** 執行兩次編譯
```bash
xelatex your_file.tex
xelatex your_file.tex  # 第二次更新交叉引用
```

## 完整編譯流程

### 基本編譯
```bash
cd "/Volumes/Seagate/14.Periodic Hill"
xelatex "4.Periodic Hill_evolution.tex"
```

### 完整編譯（含書籤、目錄）
```bash
xelatex "4.Periodic Hill_evolution.tex"
xelatex "4.Periodic Hill_evolution.tex"
```

### 含參考文獻的編譯
```bash
xelatex your_file.tex
bibtex your_file
xelatex your_file.tex
xelatex your_file.tex
```

## 檢查編譯錯誤

### 查看錯誤詳細資訊
```bash
# 查看完整 log
less "4.Periodic Hill_evolution.log"

# 搜尋錯誤
grep -i "error" "4.Periodic Hill_evolution.log"

# 搜尋警告
grep -i "warning" "4.Periodic Hill_evolution.log"

# 搜尋缺失字符
grep -i "missing character" "4.Periodic Hill_evolution.log"
```

## 推薦的 .tex 檔案字體設定（macOS）

```latex
\documentclass[12pt]{article}
\usepackage[a4paper,margin=2cm]{geometry}
\usepackage{fontspec}
\usepackage{xeCJK}

% 英文字體
\setmainfont{Times New Roman}

% 中文字體（macOS）
\setCJKmainfont[
    BoldFont={Kaiti TC Bold},
    ItalicFont={Kaiti TC},
    BoldItalicFont={Kaiti TC Bold}
]{Kaiti TC}

% CJK 等寬字體（重要！）
\setCJKmonofont{Kaiti TC}

% 其餘套件...
\usepackage{amsmath}
\usepackage{listings}
% ...
```

## listings 套件中文設定

若程式碼區塊包含中文註解：

```latex
\usepackage{listings}
\usepackage{xcolor}

\lstset{
    basicstyle=\ttfamily\small,
    commentstyle=\color{green!50!black},
    keywordstyle=\color{blue},
    stringstyle=\color{red},
    numbers=left,
    numberstyle=\tiny\color{gray},
    stepnumber=1,
    numbersep=5pt,
    backgroundcolor=\color{white},
    showspaces=false,
    showstringspaces=false,
    showtabs=false,
    frame=single,
    tabsize=4,
    captionpos=b,
    breaklines=true,
    breakatwhitespace=false,
    escapeinside={(*@}{@*)},  % 允許 LaTeX 指令
    columns=flexible,
    keepspaces=true
}
```

## 常見套件依賴

確保已安裝以下 LaTeX 套件：

**基礎：**
- `xeCJK` - 中文支援
- `fontspec` - 字體管理
- `geometry` - 頁面設定

**數學：**
- `amsmath` - 數學方程式
- `physics` - 物理符號
- `amsthm` - 定理環境

**圖表：**
- `graphicx` - 圖片插入
- `tikz` - 繪圖
- `tcolorbox` - 彩色文字框
- `float` - 浮動體控制

**其他：**
- `listings` - 程式碼區塊
- `hyperref` - 超連結與書籤
- `booktabs` - 專業表格

## 效能優化

### 草稿模式（快速預覽）
```latex
\documentclass[12pt,draft]{article}
```

### 跳過圖片載入
```latex
\usepackage[draft]{graphicx}
```

### 使用 -interaction 選項
```bash
# 遇到錯誤不停止，繼續編譯
xelatex -interaction=nonstopmode your_file.tex

# 遇到錯誤立即停止
xelatex -interaction=errorstopmode your_file.tex
```

## 疑難排解

### 問題：編譯卡住不動

**可能原因：** 等待使用者輸入

**解決：**
1. 按 Ctrl+C 中斷
2. 加上 `-interaction=nonstopmode` 選項
3. 檢查 .log 檔案找出錯誤行

### 問題：PDF 書籤亂碼

**解決：** 確保使用 `hyperref` 套件並加上：
```latex
\usepackage{hyperref}
\hypersetup{
    pdfencoding=unicode,
    colorlinks=true,
    linkcolor=blue
}
```

### 問題：目錄頁碼錯誤

**解決：** 執行兩次編譯更新頁碼

---

**提示：** 使用 Claude Code 的 LBM LaTeX Helper skill 可獲得更多編輯建議！
