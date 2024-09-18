#!/usr/bin/env python

import requests
import gzip
import shutil
import os

# Function to download file from URL
def download_file(url, save_path):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(save_path, 'wb') as file:
            file.write(response.content)

# Function to decompress the gzipped file
def decompress_gz(gz_path, decompressed_path):
    with gzip.open(gz_path, 'rb') as f_in:
        with open(decompressed_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# Function to process MEME file into individual motifs
def process_meme_file(file_path):
    motifs = []
    current_motif = []
    header = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('MOTIF'):
                if current_motif:
                    motifs.append(current_motif)
                    current_motif = []
                current_motif.append(line.strip())
            elif line.startswith('MEME') or line.startswith('ALPHABET'):
                header.append(line.strip())
            elif current_motif:  # Append lines only if we are processing a motif
                current_motif.append(line.strip())
        
        # Append last motif
        if current_motif:
            motifs.append(current_motif)
    
    return motifs, header

# Function to write motifs to individual .meme format files without appending bigBed accession
def write_motifs_to_files(motifs, header, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    for motif in motifs:
        first_line = motif[0].split()
        motif_id = first_line[1]  # Extract motif ID from the first line
        file_name = os.path.join(output_folder, f"{motif_id}.meme")
        with open(file_name, 'w') as f:
            # Write the header
            f.write('\n'.join(header) + '\n\n')
            # Write the motif content
            f.write('\n'.join(motif))

# URL to download the file
url = 'https://downloads.wenglab.org/factorbook-download/complete-factorbook-catalog.meme.gz'
gz_save_path = 'complete-factorbook-catalog.meme.gz'
decompressed_save_path = 'complete-factorbook-catalog.meme'

# Download the file
download_file(url, gz_save_path)

# Decompress the file
decompress_gz(gz_save_path, decompressed_save_path)

# Process the decompressed MEME file
motifs, header = process_meme_file(decompressed_save_path)

# Output folder name
output_folder = 'pwm'

# Write motifs to individual .meme format files without appending bigBed accession
write_motifs_to_files(motifs, header, output_folder)

print(f"Motifs written to .meme files in '{output_folder}' folder.")
