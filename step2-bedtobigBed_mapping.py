#!/usr/local/bin/python

'''
Process ENCODE metadata to map BED files to bigBed files for all ENCODE experiments used in calaloging Factorbook motifs. 
Usage: python3 step2-bedtobigBed_mapping.py filtered_metadata.txt ENCODE_experimentIDs.txt factorbook_chipseq_meme_motifs.tsv bedtobigBed_mapping.txt 
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

    # Further filter for 'IDR threshold peaks' output type (4th column, index 3)
    filtered_metadata = filtered_metadata[filtered_metadata.iloc[:, 3] == 'IDR thresholded peaks']

    # Further filter for biological replicates
    filtered_metadata = process_biological_replicates(filtered_metadata)
    
    # Check if the filtered metadata is empty
    if filtered_metadata.empty:
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
    
    def extract_id(path):
        return path.split('/')[-2]
    
    # Process derived from column
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

    # Download JSON files
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

def merge_files(bedtobigBed_df, factorbook_file, output_file):
    # Load the factorbook data
    factorbook_df = pd.read_csv(factorbook_file, sep='\t')

    # Filter out rows where 'Animal' is 'Mus musculus'
    factorbook_df = factorbook_df[factorbook_df.iloc[:, 0] != 'Mus musculus']

    # Merge the data on 'Experiment ID' using column indices
    merged_df = pd.merge(factorbook_df, bedtobigBed_df.iloc[:, [0, 1, 2]], left_on=factorbook_df.columns[3], right_on=bedtobigBed_df.columns[0], how='left')

    # Combine BED Accession and Consensus into a new column
    merged_df['BED_Consensus'] = merged_df.iloc[:, len(merged_df.columns) - 1] + '_' + merged_df.iloc[:, 4]

    # Arrange final set of columns in order
    final_df = merged_df.loc[:, [factorbook_df.columns[0], factorbook_df.columns[1], factorbook_df.columns[2], factorbook_df.columns[3], bedtobigBed_df.columns[0], bedtobigBed_df.columns[1], bedtobigBed_df.columns[2], factorbook_df.columns[4], factorbook_df.columns[5], 'BED_Consensus']]

    # Rename columns
    final_df.columns = ['Animal', 'Biosample', 'Target', 'Experiment_ID', 'Experiment_ID', 'bigBed_Accession', 'BED_Accession', 'Consensus', 'Sample_Type', 'BED_Consensus']

    # Write output
    final_df.to_csv(output_file, sep='\t', index=False)

    print(f"Output file '{output_file}' has been created successfully.")
    return final_df

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python script.py <metadata_file> <experiment_ids_file> <factorbook_file> <output_file>")
        sys.exit(1)

    metadata_file = sys.argv[1]
    experiment_ids_file = sys.argv[2]
    factorbook_file = sys.argv[3]
    output_file = sys.argv[4]

    # Load the metadata
    metadata = pd.read_csv(metadata_file, sep='\t')

    # Load the experiment IDs
    with open(experiment_ids_file, 'r') as f:
        experiment_ids = [line.strip() for line in f]
    
    # List to store output data
    output_data = []

    # List to store empty experiment IDs
    empty_experiment_ids = []

    # For each experiment ID, call metadata_filtering function
    for experiment_id in experiment_ids:
        metadata_filtering(metadata, experiment_id, output_data, empty_experiment_ids)

    # Print the number of empty experiment IDs
    print(f"Total number of empty experiment IDs: {len(empty_experiment_ids)}")

    # Convert the output data to a DataFrame
    output_df = pd.DataFrame(output_data, columns=['Experiment ID', 'bigBed Accession', 'BED Accession', 'bigBed Derived From'])

    # Save the empty experiment IDs to a file
    empty_ids_file = "empty_ids.txt"
    with open(empty_ids_file, 'w') as f:
        for eid in empty_experiment_ids:
            f.write(f"{eid}\n")

    # Perform derived from QC check on the output DataFrame
    output_df.to_csv('temp_bedtobigBed_mapping.txt', sep='\t', index=False)
    derived_from_QC('temp_bedtobigBed_mapping.txt')

    # Download JSON files for the empty experiment IDs
    tmp_dir = download_json_files(empty_ids_file)

    # Update empty_ids.txt with information from the downloaded JSON files
    update_empty_ids_with_json_info(empty_ids_file, tmp_dir)

    # Delete the temporary directory
    shutil.rmtree(tmp_dir)
    print(f"Deleted temporary directory {tmp_dir}")

    # Merge files and download data based on BED_Consensus column
    merged_df = merge_files(output_df, factorbook_file, output_file)
