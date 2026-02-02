#!/usr/bin/env bash
# 安裝 launchd 排程任務（每 5 分鐘自動修正 .tex 拼字）
# 使用方式: ./install_launchd.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PLIST_NAME="com.user.spellfix-tex.plist"
PLIST_SRC="$SCRIPT_DIR/$PLIST_NAME"
PLIST_DST="$HOME/Library/LaunchAgents/$PLIST_NAME"

echo "=== 安裝 LaTeX 拼字自動修正排程 ==="

# 確保腳本有執行權限
chmod +x "$SCRIPT_DIR/spellfix_tex.sh"
echo "✓ 已設定 spellfix_tex.sh 執行權限"

# 複製 plist 到 LaunchAgents
mkdir -p "$HOME/Library/LaunchAgents"
cp "$PLIST_SRC" "$PLIST_DST"
echo "✓ 已複製 plist 到 ~/Library/LaunchAgents/"

# 卸載舊的（如果存在）
launchctl unload "$PLIST_DST" 2>/dev/null || true

# 載入新的排程
launchctl load "$PLIST_DST"
echo "✓ 已載入 launchd 排程"

echo ""
echo "=== 安裝完成！ ==="
echo "排程將每 5 分鐘自動執行拼字修正"
echo ""
echo "查看狀態: launchctl list | grep spellfix"
echo "查看日誌: tail -f /tmp/spellfix-tex.log"
echo "卸載排程: ./uninstall_launchd.sh"
