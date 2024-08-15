from pdb import set_trace
import os
import click
import pandas as pd
from datetime import timedelta
from tqdm import tqdm
from collections import defaultdict

# Function to convert h:m:s to seconds
def hms_to_seconds(hms_str):
    h, m, s = map(int, hms_str.split(':'))
    return h * 3600 + m * 60 + s

@click.command()
@click.option('--log-dir', default='.', help='Path to the directory containing log files')
def consolidate_benchmarks(log_dir):
    # Initialize a DataFrame to store consolidated data
    columns = ['rule_name', 'sample_id', 'max_rss', 'max_vms', 'max_uss', 'max_pss', 'io_in', 'io_out', 'mean_load', 'cpu_time', 'time_seconds']
    consolidated_data = pd.DataFrame(columns=columns)

    filedict = defaultdict(lambda: defaultdict(list))

    # Iterate through log files
    logs = os.listdir(log_dir)
    for filename in tqdm(logs, total=len(logs), desc="Parsing filenames"):
        if filename.endswith('.tsv'):
            parts = filename.split('.')
            if len(parts) == 3:  # Format: <sample_id>.<rule name>.tsv
                sample_id, rule_name = parts[0], parts[1]
                split_id = None
            elif len(parts) == 4:  # Format: <sample_id>.<split id>.<rule name>.tsv
                sample_id, split_id, rule_name = parts[0], parts[1], parts[2]

            if len(filedict[sample_id][rule_name]) == 0:
                filedict[sample_id][rule_name] = [log for log in logs if log.startswith(sample_id) and rule_name in log]

    consolidated_data = []
    for sample_id, rule_dict in filedict.items():
        for rule_name, files in rule_dict.items():
            file_paths = [os.path.join(log_dir, f) for f in files]
            data = [pd.read_csv(filepath, sep='\t') for filepath in file_paths]
            df = pd.concat(data, axis=0, join='outer')
            df.drop('h:m:s', axis=1, inplace=True)
            df['rule_name'] = [rule_name for _ in range(df.shape[0])]
            df['sample_id'] = [sample_id for _ in range(df.shape[0])]
            consolidated_data.append(df)

    consolidated_data = pd.concat(consolidated_data, axis=0, join='outer')

    # Write consolidated data to a file
    consolidated_data_path = os.path.join("workup", "consolidated_benchmarks.tsv")
    consolidated_data.to_csv(consolidated_data_path, sep='\t', index=False)
    print("Consolidation complete. Consolidated data saved to 'consolidated_benchmark.tsv'.")

if __name__ == '__main__':
    consolidate_benchmarks()

