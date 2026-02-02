#!/usr/bin/env bash
# 卸載 launchd 排程任務
# 使用方式: ./uninstall_launchd.sh

set -euo pipefail

PLIST_NAME="com.user.spellfix-tex.plist"
PLIST_DST="$HOME/Library/LaunchAgents/$PLIST_NAME"

echo "=== 卸載 LaTeX 拼字自動修正排程 ==="

if [ -f "$PLIST_DST" ]; then
    launchctl unload "$PLIST_DST" 2>/dev/null || true
    rm -f "$PLIST_DST"
    echo "✓ 已卸載並移除排程"
else
    echo "⚠ 排程不存在，無需卸載"
fi

echo ""
echo "=== 卸載完成！ ==="
