#!/bin/bash

# Collect ECNODE experiment IDs from FACTORBOOK MEME ChIP-seq Catalog and download its corrsponding Json files

# URL, Input/Output files
url="https://storage.googleapis.com/gcp.wenglab.org/factorbook_chipseq_meme_motifs.tsv"
metadata_file=$(basename $url)
output_file="uniq_experimentID.txt"
output_folder="Json/experiments"

# Create the output folder if it doesn't exist
mkdir -p $output_folder

# Check if the metadata file already exists
if [ ! -f $metadata_file ]; then
    echo "Downloading TSV file..."
    curl -o $metadata_file $url
else
    echo "TSV file already exists. Skipping download."
fi

# Extract ENCODE experiment IDs
awk -F'\t' 'NR > 1 {print $4}' $metadata_file | sort | uniq > $output_file

echo "Unique IDs have been written to $output_file"

# Download corresponding JSON files using the unique IDs
cat $output_file | while read -r id; do
    json_file="$output_folder/$id.json"
    if [ ! -f $json_file ]; then
        echo "Downloading JSON for $id..."
        curl -L -H "Accept: application/json" "https://www.encodeproject.org/files/$id/" > "$json_file"
    else
        echo "JSON for $id already exists. Skipping download."
    fi
done

echo "JSON files have been checked and downloaded to $output_folder where necessary."
