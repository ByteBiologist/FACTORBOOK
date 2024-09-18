#!/usr/local/bin/python

import subprocess
import os
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import requests
import csv

def run_tomtom(command: str, output_dir: str):
    """Run a Tomtom command using subprocess if output directory doesn't exist."""
    if not os.path.exists(output_dir):
        subprocess.run(command, shell=True)
    else:
        print(f"Skipping Tomtom run, output directory '{output_dir}' already exists.")

def process_tomtom_file(file_path: str) -> dict:
    # Load the data into a DataFrame
    df = pd.read_csv(file_path, delimiter='\t')

    # Convert 'q-value' to numeric, including scientific notation, coercing errors to NaN
    df['q-value'] = pd.to_numeric(df['q-value'], errors='coerce')

    # Drop rows where 'q-value' is NaN
    df = df.dropna(subset=['q-value'])

    # Initialize an empty dictionary to store the results
    query_target_dict = {}

    # Group by 'Query_ID'
    grouped = df.groupby('Query_ID')

    for query_id, group in grouped:
        # Find the row with the smallest q-value for this Query_ID
        min_q_value_row = group.loc[group['q-value'].idxmin()]

        # Check if Query_consensus is part of Target_consensus or vice versa
        if min_q_value_row['Query_consensus'] in min_q_value_row['Target_consensus'] or min_q_value_row['Target_consensus'] in min_q_value_row['Query_consensus']:
            # Extract relevant information
            target_id = min_q_value_row['Target_ID']
            q_value = min_q_value_row['q-value']
            target_consensus = min_q_value_row['Target_consensus']
            query_consensus = min_q_value_row['Query_consensus']
            
            # Store the results in the dictionary
            query_target_dict[query_id] = {
                'Target_ID': target_id,
                'q-value': q_value,
                'Target_consensus': target_consensus,
                'Query_consensus': query_consensus
            }
    print(f"Processed {len(query_target_dict)} entries.")
    return query_target_dict


def fetch_jaspar_info_matrix(matrix_id: str) -> str:
    """Fetch the transcription factor name and species from JASPAR using the Matrix ID."""
    url = f"https://jaspar.genereg.net/api/v1/matrix/{matrix_id}/"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        tf_name = data.get("name", "NA")
        species_info = data.get("species", [])
        species = species_info[0].get("name", "NA") if species_info else "NA"
        return tf_name, species
    return "NA", "NA"

def verify_hocomoco_tf_name(hocomoco_target: str, hocomoco_tf_name: str, metadata_file: str) -> str:
    """Verify the HOCOMOCO TF name using the H12CORE_motifs.tsv file."""
    df = pd.read_csv(metadata_file, delimiter='\t')
    
    #  Find the row where the full HOCOMOCO target matches
    matching_row = df[df['Motif'] == hocomoco_target]
    #matching_row = df[df['Motif'].str.contains(hocomoco_tf_name, regex=False)]
    
    if not matching_row.empty:
        # Extract the UniProt ID (human) and process it
        uniprot_id_human = matching_row.iloc[0, 8]
        tf_human = uniprot_id_human.replace('_HUMAN', '')
        
        # Compare gene_human with hocomoco_tf_name and return uniprot human ID
        if tf_human == hocomoco_tf_name:
            return uniprot_id_human
    return "NA"

#def merge_tomtom_results(hocomoco_file: str, jaspar_file: str) -> dict:
def merge_tomtom_results(hocomoco_file: str, jaspar_file: str, metadata_file: str) -> dict:

    """Merge Tomtom results from HOCOMOCO and JASPAR into a dictionary."""
    hocomoco_dict = process_tomtom_file(hocomoco_file)
    jaspar_dict = process_tomtom_file(jaspar_file)

    merged_dict = {}
    all_query_ids = set(hocomoco_dict.keys()).union(set(jaspar_dict.keys()))

    for query_id in all_query_ids:
        hocomoco_info = hocomoco_dict.get(query_id, {})
        jaspar_info = jaspar_dict.get(query_id, {})

        hocomoco_target = hocomoco_info.get('Target_ID', 'no_match')
        hocomoco_q_value = hocomoco_info.get('q-value', 'NA')
        hocomoco_consensus = hocomoco_info.get('Target_consensus', 'NA')
        hocomoco_query = hocomoco_info.get('Query_consensus', 'NA')

        jaspar_target = jaspar_info.get('Target_ID', 'no_match')
        jaspar_q_value = jaspar_info.get('q-value', 'NA')
        jaspar_consensus = jaspar_info.get('Target_consensus', 'NA')
        jaspar_query = jaspar_info.get ('Query_consensus', 'NA')

        hocomoco_tf_name = hocomoco_target.split('.')[0] if hocomoco_target != 'no_match' else 'NA'
        uniprot_id_human = verify_hocomoco_tf_name(hocomoco_target, hocomoco_tf_name, metadata_file)
        
        #jaspar_tf_name = fetch_jaspar_info_matrix(jaspar_target) if jaspar_target != 'no_match' else 'NA'
        jaspar_tf_name, jaspar_species = fetch_jaspar_info_matrix(jaspar_target) if jaspar_target != 'no_match' else ('NA', 'NA')

        merged_dict[query_id] = {
            "HOCOMOCO_ID": hocomoco_target,
            "HOCOMOCO_q-value": hocomoco_q_value,
            "HOCOMOCO_query_consensus": hocomoco_query,
            "HOCOMOCO_Target_consensus": hocomoco_consensus,
            "HOCOMOCO_TF_name": hocomoco_tf_name,
            "HOCOMOCO_species": uniprot_id_human,
            "JASPAR_ID": jaspar_target,
            "JASPAR_q-value": jaspar_q_value,
            "JASPAR_query_consensus": jaspar_query,
            "JASPAR_Target_consensus": jaspar_consensus,
            "JASPAR_TF_name": jaspar_tf_name,
            "JASPAR_species": jaspar_species
        }

    return merged_dict

def write_dict_to_file(merged_dict: dict, output_file: str, non_human_file: str):
    """Write the merged dictionary to a tab-separated file with a header."""
    header = [
        "query_ID", "HOCOMOCO_ID", "HOCOMOCO_q-value", "HOCOMOCO_query_consensus", "HOCOMOCO_Target_consensus", "HOCOMOCO_TF_name", "HOCOMOCO_species", "JASPAR_ID", "JASPAR_q-value", "JASPAR_query_consensus", "JASPAR_Target_consensus", "JASPAR_TF_name", "JASPAR_species"
    ]
    
    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=header, delimiter='\t')
        writer.writeheader()
        
        for query_id, values in merged_dict.items():
            row = {**{"query_ID": query_id}, **values}
            writer.writerow(row)

def main():
    hocomoco_command = "./tomtom -no-ssc -oc tomtom_HOCOMOCO -verbosity 1 -min-overlap 10 -dist pearson -thresh 0.05 complete-factorbook-catalog.meme db/H12CORE_meme_format.meme"
    jaspar_command = "./tomtom -no-ssc -oc tomtom_JASPAR -verbosity 1 -min-overlap 10 -dist pearson -thresh 0.05 complete-factorbook-catalog.meme db/JASPAR2022_CORE_non-redundant_v2.meme"

    # Run Tomtom commands in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=3) as executor:
        executor.submit(run_tomtom, hocomoco_command, 'tomtom_HOCOMOCO')
        executor.submit(run_tomtom, jaspar_command, 'tomtom_JASPAR')

    # Wait for the Tomtom runs to complete before proceeding
    executor.shutdown(wait=True)

    # Process the results
    hocomoco_file = 'tomtom_HOCOMOCO/tomtom.tsv'
    jaspar_file = 'tomtom_JASPAR/tomtom.tsv'
    hocomoco_metadata = 'H12CORE_motifs.tsv'
    result_dict = merge_tomtom_results(hocomoco_file, jaspar_file, hocomoco_metadata)

    # Define output file names
    output_file = 'tomtom_mapping.txt'
    non_human_file = 'jaspar_non_human.txt'

    # Write the results to files
    write_dict_to_file(result_dict, output_file, non_human_file)

if __name__ == "__main__":
    main()
