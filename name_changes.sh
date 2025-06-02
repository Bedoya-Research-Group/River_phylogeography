#!/bin/bash

# Change to the directory where your files are
#cd /path/to/your/files

# Read each line of the mapping file
while IFS=$'\t' read -r old_prefix new_prefix; do
  for file in ${old_prefix}*; do
    # Construct new filename
    new_file=${file/$old_prefix/$new_prefix}
    echo "Renaming: $file -> $new_file"
    mv "$file" "$new_file"
  done
done < name_map.tsv