#!/bin/bash

# Check the number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 directory"
    exit 1
fi

directory="$1"
output_file="list.txt"

# List all files in the directory and save to the output file
alien_find "$directory" "AO2D.root" > "$output_file"
root cp.C

echo "File names listed and saved to $output_file"
