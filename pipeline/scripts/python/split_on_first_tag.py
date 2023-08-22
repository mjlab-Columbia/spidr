from pdb import set_trace
import click
from tqdm import tqdm

def split_by_first_tag(complete_clusters, control_output, experimental_output, groups):
    line_count = 0
    barcode_set = set()
    
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
            cntrl_vs_exp = ['_'.join(b.split('_')[:-1]) for b in round_1_barcodes]

            for tag in cntrl_vs_exp:
                barcode_set.add(tag)
    
    # TODO: Add helpful print statement to this assertion
    assert len(barcode_set) == groups
    
    # Splitting clusters based on the first tag
    with open(complete_clusters, 'r') as clusters, \
    open(control_output, 'wt') as control_out, \
    open(experimental_output, 'wt') as experimental_out:
        for line in tqdm(clusters, total=line_count):
            cluster_barcode = line.strip('\n').split('\t', 1)[0]
            barcodes = cluster_barcode.split('.')[:-1]

            # Last index corresponds to the ROUND1 tags
            first_barcode = barcodes[-1]
            
            if first_barcode.startswith('ROUND1_CNTRL'):
                control_out.write(line)
            elif first_barcode.startswith('ROUND1_TORIN'):
                experimental_out.write(line)

@click.command()
@click.option('--complete_clusters', required=True, type=click.Path(exists=True), help='Path to complete clusters file(s)')
@click.option('--control_output', required=True, type=click.Path(), help='Path to control output file')
@click.option('--experimental_output', required=True, type=click.Path(), help='Path to experimental output file')
@click.option('--groups', required=False, type=int, default=2, help='Number of expected different first barcodes (default assumes 1 control and 1 experiments)')
def main(complete_clusters, control_output, experimental_output, groups):
    split_by_first_tag(complete_clusters, control_output, experimental_output, groups)

if __name__ == '__main__':
    main()
