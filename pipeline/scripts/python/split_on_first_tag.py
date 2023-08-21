import click
from tqdm import tqdm

def split_by_first_tag(complete_clusters, control_output, experimental_output, groups):
    line_count = 0
    barcode_set = {}
    
    with open(complete_clusters, 'r') as clusters:
        for line in clusters:
            # Get a line count for a progress bar
            line_count += 1
            
            # Keep track of all possible first barcodes
            cluster_barcode = line.strip('\n').split('\t', 1)[0]
            barcodes = cluster_barcode.split('.')[:-1]
            barcode_set.add(barcodes[0])
    
    assert len(barcode_set) == groups, "--groups does not match the number of barcodes found\n" + print(barcode_set)
    
    # Your code for splitting clusters based on the first tag goes here
    with open(complete_clusters, 'r') as clusters, \
    open(control_output, 'wt') as control_out, \
    open(experimental_output, 'wt') as experimental_out:
        for line in tqdm(line, total=line_count):
            cluster_barcode = line.strip('\n').split('\t', 1)[0]
            barcodes = cluster_barcode.split('.')[:-1]
            first_barcode = barcodes[0]
            
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
