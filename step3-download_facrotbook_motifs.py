#!/usr/local/bin/python

'''
Download Factorbook motifs.

Usage: python3 step3-download_facrotbook_motifs.py bedtobigBed_mapping.txt fb_motifs
'''

import os
import sys
import pandas as pd
import requests
import gzip

def download_files(file_to_download, download_folder):
    # Load the data from the specified file
    try:
        df = pd.read_csv(file_to_download, sep='\t')
    except FileNotFoundError:
        print(f"Error: File '{file_to_download}' not found.")
        return
    except pd.errors.EmptyDataError:
        print(f"Error: File '{file_to_download}' is empty or cannot be read as CSV.")
        return

    # Create the download folder if it doesn't exist
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)

    base_url = "https://screen-beta-api.wenglab.org/factorbook_downloads/hq-occurrences/"
    deleted_files = []

    # Iterate through the BED_Consensus column to construct URLs and download files
    for bed_consensus in df['BED_Consensus']:
        local_filename = os.path.join(download_folder, f"{bed_consensus}.gz")
        if os.path.exists(local_filename):
            print(f"File already exists: {local_filename}, skipping download.")
            continue
        
        file_url = f"{base_url}{bed_consensus}.gz"
        try:
            with requests.get(file_url, stream=True) as r:
                r.raise_for_status()
                with open(local_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            print(f"Downloaded: {local_filename}")
            
            # Verify if the downloaded file is a valid gz file
            try:
                with gzip.open(local_filename, 'rb') as gz_file:
                    gz_file.read(1)  # Try reading one byte from the gz file to ensure it's valid
            except (OSError, gzip.BadGzipFile):
                print(f"Invalid gz file, deleting: {local_filename}")
                os.remove(local_filename)
                deleted_files.append(local_filename)
        
        except requests.exceptions.HTTPError as err:
            print(f"HTTP Error for {file_url}: {err}")
        except Exception as err:
            print(f"Error downloading {file_url}: {err}")

    # Log deleted files
    if deleted_files:
        with open('deleted_gz_files.txt', 'w') as log_file:
            for filename in deleted_files:
                log_file.write(f"{filename}\n")
        print("Deleted invalid gz files have been logged in 'deleted_gz_files.txt'.")
    else:
        print("No invalid gz files were found.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script_name.py <file_to_download> <download_directory>")
        sys.exit(1)

    file_to_download = sys.argv[1]
    download_directory = sys.argv[2]

    download_files(file_to_download, download_directory)
