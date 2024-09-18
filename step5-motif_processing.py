#!/usr/bin/env python

import pandas as pd
import gzip
import subprocess
import os

def get_reference(sequenceID, fasta_file):
    output = subprocess.check_output(["samtools", "faidx", fasta_file, sequenceID]).decode("utf-8")
    lines = output.split('\n')
    ref_sequence = ''.join(lines[1:])
    ref_sequence = ref_sequence.upper()  # Convert to uppercase
    return ref_sequence

def get_reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    reverse_complement = ''.join(complement.get(base, base) for base in reversed(sequence))
    return reverse_complement

def process_bedtobigBed_mapping(mapping_file_path):
    """Process bedtobigBed_mapping.txt and create lookup dictionaries."""
    input_data = pd.read_csv(mapping_file_path, sep='\t')
    
    bed_consensus_to_target = dict(zip(input_data['BED_Consensus'], input_data['Target']))
    bed_consensus_to_experiment_id = dict(zip(input_data['BED_Consensus'], input_data['Experiment_ID']))
    bed_consensus_to_bigbed_accession = dict(zip(input_data['BED_Consensus'], input_data['bigBed_Accession']))
    bed_consensus_to_bed_accession = dict(zip(input_data['BED_Consensus'], input_data['BED_Accession']))
    
    return bed_consensus_to_target, bed_consensus_to_experiment_id, bed_consensus_to_bigbed_accession, bed_consensus_to_bed_accession

def process_tomtom_mapping(tomtom_file_path, query_ID):
    """Process tomtom_mapping.txt and return relevant data for a given query_ID."""
    tomtom_data = pd.read_csv(tomtom_file_path, sep='\t')
    
    # Find the row that matches the query_ID
    tomtom_row = tomtom_data[tomtom_data['query_ID'] == query_ID]

    if tomtom_row.empty:
        # Return NA for all fields if no match is found
        return {
            'HOCOMOCO_pwm_ID': 'NA',
            'HOCOMOCO_q-value': 'NA',
            'HOCOMOCO_Target_consensus': 'NA',
            'HOCOMOCO_TF_name': 'NA',
            'JASPAR_pwm_ID': 'NA',
            'JASPAR_q-value': 'NA',
            'JASPAR_Target_consensus': 'NA',
            'JASPAR_TF_name': 'NA'
        }
    
    # Extract HOCOMOCO information regardless of JASPAR species
    hocomoco_pwm_id = tomtom_row['HOCOMOCO_ID'].values[0] if pd.notna(tomtom_row['HOCOMOCO_ID'].values[0]) else 'NA'
    hocomoco_q_value = tomtom_row['HOCOMOCO_q-value'].values[0] if pd.notna(tomtom_row['HOCOMOCO_q-value'].values[0]) else 'NA'
    hocomoco_target_consensus = tomtom_row['HOCOMOCO_Target_consensus'].values[0] if pd.notna(tomtom_row['HOCOMOCO_Target_consensus'].values[0]) else 'NA'
    hocomoco_tf_name = tomtom_row['HOCOMOCO_TF_name'].values[0] if pd.notna(tomtom_row['HOCOMOCO_TF_name'].values[0]) else 'NA'
    
    # Only extract JASPAR information if the species is Homo sapiens
    if tomtom_row['JASPAR_species'].values[0] == 'Homo sapiens':
        jaspar_pwm_id = tomtom_row['JASPAR_ID'].values[0] if pd.notna(tomtom_row['JASPAR_ID'].values[0]) else 'NA'
        jaspar_q_value = tomtom_row['JASPAR_q-value'].values[0] if pd.notna(tomtom_row['JASPAR_q-value'].values[0]) else 'NA'
        jaspar_target_consensus = tomtom_row['JASPAR_Target_consensus'].values[0] if pd.notna(tomtom_row['JASPAR_Target_consensus'].values[0]) else 'NA'
        jaspar_tf_name = tomtom_row['JASPAR_TF_name'].values[0] if pd.notna(tomtom_row['JASPAR_TF_name'].values[0]) else 'NA'
    else:
        # If the species is not Homo sapiens, set JASPAR fields to 'NA'
        jaspar_pwm_id = 'NA'
        jaspar_q_value = 'NA'
        jaspar_target_consensus = 'NA'
        jaspar_tf_name = 'NA'
        print('Non human', tomtom_row['JASPAR_species'].values[0])  
    return {
        'HOCOMOCO_pwm_ID': hocomoco_pwm_id,
        'HOCOMOCO_q-value': hocomoco_q_value,
        'HOCOMOCO_Target_consensus': hocomoco_target_consensus,
        'HOCOMOCO_TF_name': hocomoco_tf_name,
        'JASPAR_pwm_ID': jaspar_pwm_id,
        'JASPAR_q-value': jaspar_q_value,
        'JASPAR_Target_consensus': jaspar_target_consensus,
        'JASPAR_TF_name': jaspar_tf_name
    }

def prepare_bed_metadata(bed_file_path, bed_mapping_dicts, tomtom_file_path):
    """Prepare metadata for the BED file using lookup data from mappings."""
    bed_file_name = os.path.basename(bed_file_path)
    bed_data = pd.read_csv(bed_file_path, sep='\t', header=None, names=['#chrom', 'chromStart', 'chromEnd', 'strand', 'FB_score'])

    bed_consensus_to_target, bed_consensus_to_experiment_id, bed_consensus_to_bigbed_accession, bed_consensus_to_bed_accession = bed_mapping_dicts

    bed_file_key = bed_file_name.split('.')[0]
    consensus_sequence = bed_file_name.split('_')[1].split('.')[0]

    experiment_ID = bed_consensus_to_experiment_id.get(bed_file_key, 'Unknown')
    bigBed_Accession = bed_consensus_to_bigbed_accession.get(bed_file_key, 'Unknown')
    BED_Accession = bed_consensus_to_bed_accession.get(bed_file_key, 'Unknown')

    # Combine experiment_ID, bigBed_Accession, and BED_Accession into ENCODE_exp_bigBed_bed_IDs
    bed_data['ENCODE_exp_bigBed_bed_IDs'] = f"{experiment_ID}-{bigBed_Accession}-{BED_Accession}"

    bed_data['FB_ChIP-Seq_target'] = bed_consensus_to_target.get(bed_file_key, 'Unknown')
    bed_data['FB_consensus_sequence'] = consensus_sequence
    bed_data['FB_pwm_ID'] = bed_file_key

    # Get tomtom info
    tomtom_info = process_tomtom_mapping(tomtom_file_path, f"{experiment_ID}_{consensus_sequence}")
    
    bed_data['HOCOMOCO_pwm_ID'] = tomtom_info['HOCOMOCO_pwm_ID']
    bed_data['HOCOMOCO_q-value'] = tomtom_info['HOCOMOCO_q-value']
    bed_data['HOCOMOCO_Target_consensus'] = tomtom_info['HOCOMOCO_Target_consensus']
    hocomoco_tf_name = tomtom_info['HOCOMOCO_TF_name']
    
    bed_data['JASPAR_pwm_ID'] = tomtom_info['JASPAR_pwm_ID']
    bed_data['JASPAR_q-value'] = tomtom_info['JASPAR_q-value']
    bed_data['JASPAR_Target_consensus'] = tomtom_info['JASPAR_Target_consensus']
    jaspar_tf_name = tomtom_info['JASPAR_TF_name']
    
    # Combine HOCOMOCO_TF_name and JASPAR_TF_name into TF_name
    tf_names = []
    if hocomoco_tf_name != 'NA':
        tf_names.append(hocomoco_tf_name)
    if jaspar_tf_name != 'NA':
        tf_names.append(jaspar_tf_name)
    bed_data['hocomoco_jaspar_TF_name'] = ','.join(tf_names) if tf_names else 'NA'
    
    return bed_data

def process_and_write_bed(bed_data, fasta_file, unfiltered_output_directory, filtered_output_directory):
    """Process sequences, apply filtering, and write both filtered and unfiltered BED files to separate directories."""
    
    # Ensure the DataFrame has correct column names
    bed_data.rename(columns={
        'consensus_sequence': 'FB_consensus_sequence',
        'ENCODE_lookup': 'ENCODE_exp_bigBed_bed_IDs',
        'ChIP-Seq_target': 'FB_ChIP-Seq_target'
    }, inplace=True)

    # Add the binding_sequence column
    bed_data['binding_sequence'] = ''

    # Rearrange the columns in the correct new order
    bed_data = bed_data[['#chrom', 'chromStart', 'chromEnd', 'hocomoco_jaspar_TF_name', 'FB_score', 'strand', 'binding_sequence',
                         'FB_pwm_ID', 'FB_consensus_sequence', 'FB_ChIP-Seq_target',
                         'HOCOMOCO_pwm_ID', 'HOCOMOCO_Target_consensus', 'HOCOMOCO_q-value',
                         'JASPAR_pwm_ID', 'JASPAR_Target_consensus', 'JASPAR_q-value',
                         'ENCODE_exp_bigBed_bed_IDs']]

    bed_data.sort_values(by=['#chrom', 'chromStart', 'chromEnd'], inplace=True)

    # Create directories if they don't exist
    os.makedirs(unfiltered_output_directory, exist_ok=True)
    os.makedirs(filtered_output_directory, exist_ok=True)

    # Write the unfiltered data (all rows) to the unfiltered_output_directory
    bed_file_key = bed_data['FB_pwm_ID'].iloc[0]
    unfiltered_output_file_name = f"{unfiltered_output_directory}/{bed_file_key}_sorted.bed.gz"

    with gzip.open(unfiltered_output_file_name, "wt", encoding="utf-8") as output_file:
        output_file.write('\t'.join(bed_data.columns) + '\n')
        for index, row in bed_data.iterrows():
            start = int(row['chromStart']) + 1  # Convert from 0-based to 1-based
            end = int(row['chromEnd'])
            strand = row['strand']
            sequenceID = f"{row['#chrom']}:{start}-{end}"

            if strand == '+':
                binding_sequence = get_reference(sequenceID, fasta_file)
            elif strand == '-':
                ref_sequence = get_reference(sequenceID, fasta_file)
                binding_sequence = get_reverse_complement(ref_sequence)

            bed_data.at[index, 'binding_sequence'] = binding_sequence

        for index, row in bed_data.iterrows():
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f'Processed and saved unfiltered BED file to {unfiltered_output_file_name}')
    
    # Apply filtering: keep only rows where 'FB_score' column is less than 1e-4 and make an explicit copy
    bed_data_filtered = bed_data[bed_data['FB_score'] < 1e-4].copy()

    # Check if the filtered DataFrame is empty
    if bed_data_filtered.empty:
        # Log the empty file information
        with open('empty_motif_files.log', 'a') as log_file:
            log_file.write(f"{bed_file_key}\n")
        print(f"Motif: {bed_file_key} is empty after filtering.")
        return  # Skip further processing for this file

    # Write the filtered data to the filtered_output_directory
    filtered_output_file_name = f"{filtered_output_directory}/{bed_file_key}_filtered.bed.gz"
    
    with gzip.open(filtered_output_file_name, "wt", encoding="utf-8") as output_file:
        output_file.write('\t'.join(bed_data_filtered.columns) + '\n')
        for index, row in bed_data_filtered.iterrows():
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f'Processed and saved filtered BED file to {filtered_output_file_name}')


def main(input_directory, bed_mapping_file, tomtom_file_path, fasta_file, unfiltered_output_directory, filtered_output_directory):
    bed_mapping_dicts = process_bedtobigBed_mapping(bed_mapping_file)

    for bed_file in os.listdir(input_directory):
        if bed_file.endswith('.gz'):
            bed_file_path = os.path.join(input_directory, bed_file)
            bed_data = prepare_bed_metadata(bed_file_path, bed_mapping_dicts, tomtom_file_path)
            process_and_write_bed(bed_data, fasta_file, unfiltered_output_directory, filtered_output_directory)

if __name__ == "__main__":
    input_directory = "fb_motifs"
    bed_mapping_file = "bedtobigBed_mapping.txt"
    tomtom_file_path = "tomtom_mapping.txt"
    fasta_file = "hg38.fa"
    unfiltered_output_directory = "fb_motifs_processed"
    filtered_output_directory = "fb_motifs_filtered"

    main(input_directory, bed_mapping_file, tomtom_file_path, fasta_file, unfiltered_output_directory, filtered_output_directory)

