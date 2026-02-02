# LaTeX 拼字自動修正設定指南

## 功能說明
1. **VS Code 外掛自動修正**: 儲存 `.tex` 檔案時自動修正拼字
2. **系統排程自動修正**: 每 5 分鐘自動掃描並修正所有 `.tex` 檔案的拼字錯誤

---

## 1. VS Code 設定 (已完成)

### 已安裝的設定
- `settings.json` 已設定 Code Spell Checker
- 儲存時自動修正拼字錯誤
- Run On Save 外掛設定 (可選)

### 建議安裝的 VS Code 外掛
1. **Code Spell Checker** (`streetsidesoftware.code-spell-checker`)
   - 即時拼字檢查
   - 儲存時自動修正

2. **Run On Save** (`emeraldwalk.RunOnSave`) - 可選
   - 儲存 `.tex` 時執行 spellfix_tex.sh

### 安裝外掛指令
```bash
code --install-extension streetsidesoftware.code-spell-checker
code --install-extension emeraldwalk.RunOnSave
```

---

## 2. macOS launchd 排程設定

### ⚠️ 重要：授予磁碟存取權限

由於 macOS 安全性限制，launchd 執行外接硬碟上的腳本需要授予權限：

#### 方法 A：授予 bash 完整磁碟存取權限
1. 開啟「系統設定」→「隱私與安全性」→「完整磁碟存取權限」
2. 點擊「+」按鈕
3. 按 `Cmd+Shift+G`，輸入 `/bin/bash`
4. 選擇 `bash` 並加入清單
5. 重新載入排程

#### 方法 B：將腳本複製到本機 (推薦)
```bash
# 建立本機腳本目錄
mkdir -p ~/scripts

# 複製腳本
cp /Volumes/Seagate/Lattie-Boltzmann-Method-Note_of_PaperReView/scripts/spellfix_tex.sh ~/scripts/

# 修改腳本中的路徑 (如果需要)
nano ~/scripts/spellfix_tex.sh
```

### 安裝排程
```bash
cd /Volumes/Seagate/Lattie-Boltzmann-Method-Note_of_PaperReView/scripts
./install_launchd.sh
```

### 卸載排程
```bash
./uninstall_launchd.sh
```

### 檢查排程狀態
```bash
launchctl list | grep spellfix
```

### 查看日誌
```bash
tail -f /tmp/spellfix-tex.log
tail -f /tmp/spellfix-tex.err
```

### 手動觸發執行
```bash
launchctl start com.user.spellfix-tex
```

---

## 3. 手動執行拼字修正

如果排程有問題，可以直接執行：
```bash
/Volumes/Seagate/Lattie-Boltzmann-Method-Note_of_PaperReView/scripts/spellfix_tex.sh
```

---

## 檔案清單

```
scripts/
├── spellfix_tex.sh              # 拼字修正腳本
├── com.user.spellfix-tex.plist  # launchd 排程設定
├── install_launchd.sh           # 安裝排程腳本
├── uninstall_launchd.sh         # 卸載排程腳本
└── SPELLFIX_README.md           # 本說明文件
```
