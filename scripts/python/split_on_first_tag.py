import os
import re
from pdb import set_trace
import click
from tqdm import tqdm

@click.command()
@click.option('--complete_clusters', required=True, type=click.Path(exists=True), help='Path to complete clusters file(s)')
@click.option('--output_dir', required=True, type=click.Path(), help='Directory to store the output cluster files')
def main(complete_clusters, output_dir):
    line_count = 0
    conditions = set()
    
    print("Collecting ROUND1 barcodes")
    with open(complete_clusters, 'r') as clusters:
        for line in clusters:
            # Get a line count for a progress bar + strip barcodes from cluster lines
            line_count += 1
            cluster_barcode = line.strip('\n').split('\t', 1)[0]
            barcodes = cluster_barcode.split('.')[:-1]

            # Get only the ROUND1 barcodes which could be in different indices and multiple indices
            round_1_barcodes = [barcode for barcode in barcodes if barcode.startswith("ROUND1_")]

            # Strip alphanumeric suffixes at the end of control/experimental tags (e.g. ROUND1_CNTRL_A10 --> ROUND1_CNTRL)
            round1_tags = ['_'.join(b.split('_')[:-1]) for b in round_1_barcodes]

            for tag in round1_tags:
                conditions.add(tag.split('_')[1])
    
    # For each condition, store a list of clusters with each condition (e.g. CNTRL, TORIN, etc)
    output_dict = {condition: [] for condition in conditions}

    # Splitting clusters based on the first tag
    with open(complete_clusters, 'r') as clusters:
        for line in tqdm(clusters, total=line_count):
            cluster_barcode = line.strip('\n').split('\t', 1)[0]
            barcodes = cluster_barcode.split('.')[:-1]

            # Last index corresponds to the ROUND1 tags
            first_barcode = barcodes[-1]

            # Define a regular expression pattern to extract the middle part
            pattern = r'_([^_]+)_'
            search_results = re.search(pattern, first_barcode)
            condition = search_results.group(1)
            output_dict[condition].append(line)

    print("Writing each cluster condition to disk")

    os.makedirs(output_dir, exist_ok=True)

    for condition, lines in output_dict.items():
        print(f"Writing {condition} to disk")
        filename = os.path.basename(complete_clusters)
        filename = filename.replace('.complete.clusters', f'.{condition}.clusters')
        filepath = os.path.join(output_dir, filename)

        with open(filepath, 'w') as f:
            f.write(''.join(lines))
            
if __name__ == '__main__':
    main()
