#!/usr/bin/env bash
#===============================================================================
# sync_pqtl.sh
#
# Description:
#   1. Mirror the entire `web_out/pqtl/` directory tree to `pqtl_share/`,
#      excluding all `output/` subdirectories.
#   2. Copy only those files in each `output/` subdirectory whose names contain "2".
#   3. In the target `output/` directories, remove the character "2" from each filename.
#   4. Decompress all `*.txt.gz` files in the target tree and delete the original archives.
#===============================================================================

set -Eeuo pipefail

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
SOURCE_DIR="/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl"
TARGET_DIR="/gpfs/chencao/ysbioinfor/project/proteohubProject/pqtl_share"

#-------------------------------------------------------------------------------
# Step 1: Mirror everything except the "output/" subdirectories
#-------------------------------------------------------------------------------
echo "=== Step 1: Sync all files except 'output/' directories ==="
rsync -av \
    --prune-empty-dirs \        # avoid creating empty output/ directories
    --exclude='*/output/*' \    # exclude all contents of any output/ subdirectory
    "${SOURCE_DIR}/" \
    "${TARGET_DIR}/"

#-------------------------------------------------------------------------------
# Step 2: Copy only files whose names contain "2" from each output/ subtree
#-------------------------------------------------------------------------------
echo "=== Step 2: Copy files containing '2' in any output/ subdirectory ==="
cd "${SOURCE_DIR}"
find . \
    -type f \
    -path '*/output/*' \
    -name '*2*' \
    -print0 | \
  rsync -av \
    --files-from=- \
    --from0 \
    "${SOURCE_DIR}/" \
    "${TARGET_DIR}/"

#-------------------------------------------------------------------------------
# Step 3: Remove the character "2" from filenames in all target output/ dirs
#-------------------------------------------------------------------------------
echo "=== Step 3: Strip '2' from filenames in target output/ directories ==="
cd "${TARGET_DIR}"
# Requires Perl's rename utility (Debian/Ubuntu package: rename; RHEL/CentOS: perl-File-Rename)
find . \
    -type f \
    -path '*/output/*' \
    -name '*2*' \
    -exec rename 's/2//g' {} +

# Alternative if 'rename' is unavailable:
# find . -type f -path '*/output/*' -name '*2*' -print0 | \
# while IFS= read -r -d '' file; do
#   dir=$(dirname "$file")
#   base=$(basename "$file")
#   mv "$file" "$dir/${base//2/}"
# done

#-------------------------------------------------------------------------------
# Step 4: Decompress all .txt.gz files (removing the .gz automatically)
#-------------------------------------------------------------------------------
echo "=== Step 4: Decompress all .txt.gz files in target directory ==="
find "${TARGET_DIR}" \
    -type f \
    -name '*.txt.gz' \
    -exec gunzip {} +

echo ">>> All tasks completed successfully!"