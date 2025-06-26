#!/usr/bin/env bash
# -----------------------------------------------------------------------
# Batch tarball generator for pQTL protein sub-projects
# -----------------------------------------------------------------------
# This script compresses all first-level directories whose names start
# with the letters a–z into individual *.tar.gz archives.  Each archive
# is stamped with the current date.
#
# Author : Hydraulik
# Updated: 2025-06-22
# -----------------------------------------------------------------------

# Move to the project root; abort if the directory does not exist
cd /gpfs/chencao/ysbioinfor/project/proteohubProject/pqtl_share || exit 1

# Shell options
shopt -s nullglob   # Allow empty globs (unmatched patterns expand to nothing)
shopt -s nocaseglob # Make globbing case-insensitive (matches A/, B/, …)

# Date stamp used in the archive names (YYYYMMDD)
today=$(date +%Y%m%d)

# Iterate over the alphabet
for letter in {a..z}; do
    # Collect first-level directories that start with the current letter
    dirs=( "${letter}"*/ )

    # Skip this letter if no matching directories were found
    (( ${#dirs[@]} == 0 )) && continue

    # Target tarball name: e.g. ProteoNexus_pQTL_protein_a_20250101.tar
    tarball="ProteoNexus_pQTL_protein_${letter}_${today}.tar"

    # Report progress
    echo "Packing ${dirs[*]} -> ${tarball}"

    tar -cf "${tarball}" -C . "${dirs[@]}"
done