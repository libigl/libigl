#!/bin/sh
#
# To enable this hook, issue:
#
#     ln -s ../../scripts/pre-commit.sh .git/hooks/pre-commit
#
RED='\033[0;31m'
NC='\033[0m'

## Throw error if any files contain "std::__1"

# list of changed files
CHANGED=$(git diff --cached --name-only --diff-filter=ACM)
# Ignore this file!
CHANGED=$(echo "$CHANGED" | grep -v "scripts/pre-commit.sh")
# Changed files containing the namespace "std::__1"
STD1=$(grep -Il "std::__1" $CHANGED)
if  [ $? -eq 0 ]; then
  STD1=$(echo "$STD1" | sed -e "s/^/  /")
  >&2 echo "[pre-commit hook] Error: Commit prohibited.

The following files contain the offensive \"std::__1\" namespace: 

${RED}$STD1${NC}

Consider issueing:

    sed -i '' -e \"s/std::__1/std/g\"" $STD1

  exit 1
else
  exit 0
fi
