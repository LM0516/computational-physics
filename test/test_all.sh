#!/bin/bash

# Parse flags
NOCONFIRM=false
for arg in "$@"; do
    case $arg in
    --noconfirm) NOCONFIRM=true ;;
    esac
done

# Navigate to the repository root relative to this script
cd "$(dirname "$0")/.." || exit 1

echo "Looking for .jl files in exercises/01 to exercises/06..."
echo ""

# Find all .jl files in folders 01 to 06, sort them, and loop through them
# mapfile -d '' files < <(find 0[1-6]* -type f -name "*.jl" -print0 2>/dev/null | sort -z)
while IFS= read -r -d '' f; do
    files+=("$f")
done < <(find exercises/0[1-6]* -type f -name "*.jl" -print0 2>/dev/null | sort -z)

for file in "${files[@]}"; do
    # Visual Separator
    echo "============================================================"
    echo "Running: $file"
    echo "============================================================"

    # Run the julia file
    # If --noconfirm, pipe endless newlines into julia to auto-confirm any readline() calls
    if [ "$NOCONFIRM" = true ]; then
        yes "" | julia --project=. "$file"
    else
        julia --project=. "$file"
    fi

    # Pause for user input (mimics readline())
    echo ""
    if [ "$NOCONFIRM" = false ]; then
        read -p "Press Enter to run the next file (or Ctrl+C to exit)..." </dev/tty
    fi
done

echo "Done! No more Julia files to run."
