#!/usr/bin/env bash
# Riker Engine — DELIBERATE publish.
#
# The post-commit auto-push hook was DISABLED 2026-06-29 (renamed to
# .git/hooks/post-commit.disabled). Committing no longer publishes anything.
# Pushing `main` to the public mirror (web) + GitHub (origin) is now an EXPLICIT
# act: run this script only when you mean to publish, and review what it prints
# before confirming.
set -euo pipefail
cd "$(dirname "$0")/.."

branch="$(git rev-parse --abbrev-ref HEAD)"
if [ "$branch" != "main" ]; then
  echo "Refusing: you are on '$branch', not 'main'. Checkout main first."
  exit 1
fi

echo "=== About to publish  main -> web (Vultr mirror) + origin (GitHub) ==="
echo "Local main : $(git rev-parse --short main)"
echo "Mirror has : $(curl -fsS --max-time 15 https://rikerengine.quickaffordablesites.com/STATE.json 2>/dev/null | grep -oE '"commit_short": *"[^"]*"' || echo '?')"
echo
echo "Commits that will become PUBLIC (web/main..main):"
git log --oneline web/main..main 2>/dev/null || echo "  (cannot compute offline; review 'git log' manually)"
echo
read -rp "Type 'publish' to push these to the public mirror + GitHub: " ans
[ "$ans" = "publish" ] || { echo "Aborted — nothing pushed."; exit 1; }

git push --follow-tags web main
git push --follow-tags origin main
echo "Published main to web + origin."
