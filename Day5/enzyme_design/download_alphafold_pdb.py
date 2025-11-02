import requests
import os
import pandas as pd
# Function to download a PDB file for a given UniProt ID
def download_pdb(uniprot_id, database_version='v6', output_dir='pdb_files_new'):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Construct AlphaFold DB URL
    alphafold_id = f'AF-{uniprot_id}-F1'
    model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_id}-model_{database_version}.pdb'
    print(model_url)
    output_path = os.path.join(output_dir, f'{alphafold_id}.pdb')
    
    try:
        # Send request to download the PDB file
        response = requests.get(model_url)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            print(f'Successfully downloaded {alphafold_id}.pdb to {output_path}')
            return True
        else:
            print(f'Failed to download {alphafold_id}.pdb (Status code: {response.status_code})')
            return False
    except Exception as e:
        print(f'Error downloading {alphafold_id}.pdb: {str(e)}')
        return False

# Main function to process UniProt IDs
def main():
    
    df = pd.read_csv("pfam14428_f100_id.csv")
    uniprot_ids = list(set(df['id']))
    
    # Download PDB files for each UniProt ID
    database_version = 'v6'  # Adjust if a newer version is available
    successful = 0
    failed = 0
    
    for uniprot_id in uniprot_ids:
        if download_pdb(uniprot_id, database_version):
            successful += 1
        else:
            failed += 1
    
    # Print summary
    print(f'\nDownload Summary: {successful} successful, {failed} failed')

if __name__ == '__main__':
    main()