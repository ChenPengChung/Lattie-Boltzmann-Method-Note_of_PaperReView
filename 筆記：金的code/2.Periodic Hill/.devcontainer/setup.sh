#!/usr/bin/env bash
set -e

echo "==> Installing system packages (C++ toolchain + full LaTeX)..."
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  build-essential cmake ninja-build gdb \
  clang clang-format clang-tidy cppcheck \
  pkg-config git curl ca-certificates \
  texlive-full latexmk

echo "==> Installing AI CLIs (Claude Code + Codex)..."
sudo npm install -g @anthropic-ai/claude-code @openai/codex

echo "==> Done."
