#!/bin/bash

# 1. Check if a directory was provided
if [ -z "$1" ]; then
    echo "Error: Please specify a directory."
    echo "Usage: ./run_julia_folder.sh \"path/to/folder\""
    exit 1
fi

SEARCH_DIR="$1"

# 2. Check if the directory exists
if [ ! -d "$SEARCH_DIR" ]; then
    echo "Error: Directory '$SEARCH_DIR' not found."
    exit 1
fi

echo "Looking for .jl files in: $SEARCH_DIR"
echo ""

# 3. Find all .jl files, sort them, and loop through them
# We use 'find' to handle recursion and 'read' to handle spaces in filenames safely
find "$SEARCH_DIR" -type f -name "*.jl" | sort | while read -r file; do

    # Visual Separator
    echo "============================================================"
    echo "Running: $file"
    echo "============================================================"

    # 4. Run the julia file
    # We use quotes "$file" to handle spaces in the folder names
    julia "$file"

    # 5. Pause for user input (mimics readline())
    echo "" # Empty line for spacing
    read -p "Press Enter to run the next file (or Ctrl+C to exit)..." </dev/tty

done

echo "Done! No more Julia files in this folder."
