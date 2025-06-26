#!/usr/bin/env bash

# -------------------------------------------------------------------------------
# Script Name : count_file_extensions.sh
# Description : Recursively count and list all file extensions under a given directory.
# Usage       : 
#   bash count_file_extensions.sh
#   chmod +x count_file_extensions.sh && ./count_file_extensions.sh
# -------------------------------------------------------------------------------

# --- Configuration -------------------------------------------------------------

# Set this to the directory you want to scan
TARGET_DIR="/gpfs/chencao/ysbioinfor/project/proteohubProject/pqtl_share/"

# --- Main Logic ----------------------------------------------------------------

# Step 1: Verify that the target directory exists
if [ ! -d "$TARGET_DIR" ]; then
  echo "Error: Directory '$TARGET_DIR' does not exist or is not a valid directory."
  exit 1
fi

echo "Scanning directory: $TARGET_DIR"
echo "----------------------------------------------------"
printf "%-8s %s\n" "Count" "Extension"
echo "----------------------------------------------------"

# Step 2: Execute the core pipeline to extract and count file extensions
#
# Pipeline breakdown:
#   1. find "$TARGET_DIR" -type f -name '*.*'
#      - Locate all files (-type f) whose names contain a dot.
#   2. sed 's/.*\.//'
#      - Strip everything up to the last dot, leaving only the extension.
#   3. tr '[:upper:]' '[:lower:]'
#      - Normalize extensions to lowercase (e.g. .JPG → .jpg).
#   4. sort
#      - Sort the list of extensions alphabetically.
#   5. uniq -c
#      - Count occurrences of each extension.
#   6. sort -nr
#      - Sort the results numerically in descending order.
#
find "$TARGET_DIR" -type f -name '*.*' \
  | sed 's/.*\.//' \
  | tr '[:upper:]' '[:lower:]' \
  | sort \
  | uniq -c \
  | sort -nr \
  | awk '{ printf "%-8s %s\n", $1, $2 }'

echo "----------------------------------------------------"
echo "Scan complete."