#!/usr/local/bin/python

'''
Process ENCODE metadata to map BED files to bigBed files for all Factorbook motfs.

'''

import pandas as pd
import sys
import os
import subprocess
import tempfile
import json
import shutil

def process_biological_replicates(metadata):
    # Extract all unique biological replicates from the filtered metadata
    unique_replicates = metadata.iloc[:, 5].unique()
    
    # Flatten and collect all unique replicates
    all_replicates = set()
    for value in unique_replicates:
        replicates = str(value).split(',')
        for rep in replicates:
            all_replicates.add(rep.strip())
    
    # Sort and join to create the filter string
    filter_list = sorted(all_replicates, key=int)  # Sort numerically
    filter_str = ', '.join(filter_list)
    
    # Return filtered metadata
    return metadata[metadata.iloc[:, 5] == filter_str]

def verify_file_origin(bed_files, bigbed_files, output_data, experiment_id):
    matched = False
    for _, bigbed_file in bigbed_files.iterrows():
        derived_from_paths = [path.strip() for path in bigbed_file[6].split(',')]
        for path in derived_from_paths:
            derived_file_accession = path.split('/')[2]
            matching_bed_file = bed_files[bed_files.iloc[:, 0] == derived_file_accession]
            if not matching_bed_file.empty:
                output_data.append([
                    experiment_id,
                    bigbed_file[0],  # bigBed accession
                    matching_bed_file.iloc[0, 0],  # BED accession
                    bigbed_file[6]  # bigBed derived from column
                ])
                matched = True
    return matched

def metadata_filtering(metadata, experiment_id, output_data, empty_experiment_ids):
    # Filter metadata for the experiment ID (5th column, index 4)
    filtered_metadata = metadata[metadata.iloc[:, 4] == experiment_id]
    #print(f"Filtered by Experiment ID ({experiment_id}):")
    #print(filtered_metadata)
    #print("\n")


    # Further filter for 'IDR threshold peaks' output type (4th column, index 3)
    filtered_metadata = filtered_metadata[filtered_metadata.iloc[:, 3] == 'IDR thresholded peaks']
    #print("Filtered by 'IDR thresholded peaks' Output Type:")
    #print(filtered_metadata)
    #print("\n")

    # Further filter for biological replicates
    filtered_metadata = process_biological_replicates(filtered_metadata)
    

    # Check if the filtered metadata is empty
    if filtered_metadata.empty:
        #print(f"Experiment ID {experiment_id} not found after filtering.\n")
        empty_experiment_ids.append(experiment_id)
        return

    # Separate BED and bigBed files (2nd column, index 1)
    bed_files = filtered_metadata[filtered_metadata.iloc[:, 1] == 'bed']
    bigbed_files = filtered_metadata[filtered_metadata.iloc[:, 1] == 'bigBed']

    # Verify derived from relationship and store matching pairs in output_data
    if not verify_file_origin(bed_files, bigbed_files, output_data, experiment_id):
        print(f"BED and bigBed do not match for Experiment ID {experiment_id}\n")

def derived_from_QC(file):
    # Read the file into a DataFrame
    df = pd.read_csv(file, delimiter='\t')
    
    # Function to extract the relevant part from the 4th column
    def extract_id(path):
        return path.split('/')[-2]
    
    # Apply the extraction function to the 4th column
    df['Parsed Column 4'] = df['bigBed Derived From'].apply(extract_id)
    
    # Compare the 3rd and the parsed 4th columns
    mismatches = df[df['BED Accession'] != df['Parsed Column 4']]
    
    # Print non-matching rows if there are any mismatches
    if not mismatches.empty:
        print("Non-matching rows:")
        print(mismatches)
    else:
        print("bigBed to bed mapping is verified.")

def download_json_files(empty_ids_file):
    with open(empty_ids_file, 'r') as f:
        experiment_ids = [line.strip() for line in f]

    # Create a temporary directory for downloading JSON files in the current directory
    tmp_dir = tempfile.mkdtemp(prefix="experiment_json_", dir=".")

    # Download JSON files using parallel curl commands
    for experiment_id in experiment_ids:
        output_path = os.path.join(tmp_dir, f"{experiment_id}.json")
        curl_command = f'curl -L -H "Accept: application/json" https://www.encodeproject.org/files/{experiment_id}/ > {output_path}'
        subprocess.run(curl_command, shell=True, check=True)

    print(f"Downloaded JSON files to {tmp_dir}")
    return tmp_dir

def update_empty_ids_with_json_info(empty_ids_file, tmp_dir):
    with open(empty_ids_file, 'r') as f:
        experiment_ids = [line.strip() for line in f]

    with open(empty_ids_file, 'w') as f:
        for experiment_id in experiment_ids:
            json_file_path = os.path.join(tmp_dir, f"{experiment_id}.json")
            if os.path.exists(json_file_path):
                with open(json_file_path, 'r') as json_file:
                    data = json.load(json_file)
                    biosample_summary = data.get('biosample_summary', '')
                    status = data.get('status', '')
                    
                    if 'Mus musculus' in biosample_summary:
                        f.write(f"{experiment_id}\t{biosample_summary}\n")
                    else:
                        f.write(f"{experiment_id}\t{status}\n")
            else:
                f.write(f"{experiment_id}\tJSON file not found\n")

    print(f"Updated {empty_ids_file} with biosample_summary and status information.")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <metadata_file> <experiment_ids_file> <output_file>")
        sys.exit(1)

    metadata_file = sys.argv[1]
    experiment_ids_file = sys.argv[2]
    output_file = sys.argv[3]

    # Load the metadata
    metadata = pd.read_csv(metadata_file, sep='\t')

    # Load the experiment IDs
    with open(experiment_ids_file, 'r') as f:
        experiment_ids = [line.strip() for line in f]
    
    # List to store output data
    output_data = []
    # List to store empty experiment IDs
    empty_experiment_ids = []

    #experiment_id = 'ENCSR000EAD'

    # For each experiment ID, call the function and process the lines
    for experiment_id in experiment_ids:
        metadata_filtering(metadata, experiment_id, output_data, empty_experiment_ids)

    # Print the number of empty experiment IDs
    print(f"Total number of empty experiment IDs: {len(empty_experiment_ids)}")

    # Convert the output data to a DataFrame and save it to a file
    output_df = pd.DataFrame(output_data, columns=['Experiment ID', 'bigBed Accession', 'BED Accession', 'bigBed Derived From'])
    output_df.to_csv(output_file, sep='\t', index=False)

    # Save the empty experiment IDs to a file
    empty_ids_file = "empty_ids.txt"
    with open(empty_ids_file, 'w') as f:
        for eid in empty_experiment_ids:
            f.write(f"{eid}\n")

    # Perform QC check on the output file
    derived_from_QC(output_file)

    # Download JSON files for the empty experiment IDs
    tmp_dir = download_json_files(empty_ids_file)

    # Update empty_ids.txt with information from the downloaded JSON files
    update_empty_ids_with_json_info(empty_ids_file, tmp_dir)

    # Delete the temporary directory
    shutil.rmtree(tmp_dir)
    print(f"Deleted temporary directory {tmp_dir}")

