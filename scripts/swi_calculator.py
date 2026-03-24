import pandas as pd
import numpy as np
import re
import os
import argparse

# Define the weights and constants
weights = {'A': 0.8356471476582918,
           'C': 0.5208088354857734,
           'E': 0.9876987431418378,
           'D': 0.9079044671339564,
           'G': 0.7997168496420723,
           'F': 0.5849790194237692,
           'I': 0.6784124413866582,
           'H': 0.8947913996466419,
           'K': 0.9267104557513497,
           'M': 0.6296623675420369,
           'L': 0.6554221515081433,
           'N': 0.8597433107431216,
           'Q': 0.789434648348208,
           'P': 0.8235328714705341,
           'S': 0.7440908318492778,
           'R': 0.7712466317693457,
           'T': 0.8096922697856334,
           'W': 0.6374678690957594,
           'V': 0.7357837119163659,
           'Y': 0.6112801822947587}

A = 81.0581
B = -62.7775

## Function to read FASTA file
#def fasta_reader(file):
#    '''Converts .fasta to a pandas dataframe with accession as index and sequence in a column 'sequence' '''
#    valid = re.compile('^[ACEDGFIHKMLNQPSRTWVY]+$')
#    with open(file, 'r') as f:
#        lines = f.read().splitlines()
#    
#    data = {'Accession': [], 'Sequence': []}
#    for line in lines:
#        if line.startswith('>'):
#            accession = line[1:]
#            sequence = ''
#        else:
#            sequence = line.strip()
#            if valid.match(sequence):
#                data['Accession'].append(accession)
#                data['Sequence'].append(sequence)
#    
#    fasta_df = pd.DataFrame(data)
#    return fasta_df

def fasta_reader(file):
    """
    Reads a multi-line FASTA file and returns a DataFrame with:
    - 'Accession': header line (without '>')
    - 'Sequence': full concatenated sequence
    """
    valid = re.compile('^[ACEDGFIHKMLNQPSRTWVY]+$', re.IGNORECASE)
    data = {'Accession': [], 'Sequence': []}
    
    accession = None
    sequence_lines = []

    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            if line.startswith('>'):
                # Save the previous sequence
                if accession and sequence_lines:
                    full_seq = ''.join(sequence_lines).upper()
                    if valid.match(full_seq):
                        data['Accession'].append(accession)
                        data['Sequence'].append(full_seq)
                # Start new record
                accession = line[1:].strip()
                sequence_lines = []
            else:
                sequence_lines.append(line)

        # Don't forget the last one!
        if accession and sequence_lines:
            full_seq = ''.join(sequence_lines).upper()
            if valid.match(full_seq):
                data['Accession'].append(accession)
                data['Sequence'].append(full_seq)

    return pd.DataFrame(data)

# Function to compute SWI
def compute_swi(df):
    '''Computes the Solubility-Weighted Index and Probability of Solubility for each sequence in the dataframe.'''
    df['SWI'] = df['Sequence'].apply(lambda x: np.mean([weights[i] for i in x]))
    df['Prob. of Solubility'] = 1 / (1 + np.exp(-(A * df['SWI'] + B)))
    return df

def main():
    parser = argparse.ArgumentParser(description='Calculate Solubility-Weighted Index (SWI) from FASTA file.')
    parser.add_argument('fasta_file', type=str, help='Path to the input FASTA file.')
    parser.add_argument('-o', '--output', type=str, default='swi_results_sorted.csv', help='Path to the output CSV file.')
    
    args = parser.parse_args()
    
    # Read the FASTA file
    df = fasta_reader(args.fasta_file)
    
    # Compute SWI and Probability of Solubility
    df = compute_swi(df)
    
    # Sort the DataFrame by the SWI column
    df = df.sort_values(by='SWI', ascending=False)
    
    # Save the sorted results to a CSV file
    df.to_csv(args.output, index=False)
    print(f'Sorted results saved to {args.output}')

if __name__ == '__main__':
    main()
