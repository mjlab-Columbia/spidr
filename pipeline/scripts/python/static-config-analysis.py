import sys
from pdb import set_trace
import pandas as pd
import yaml
import click
import json
from colorama import Fore, Back, Style, init

# Initialize colorama
init(autoreset=True)

def check_barcodes(barcode_df, config_conditions):
    round1_tags = barcode_df[barcode_df['tag-name'].str.contains('ROUND1')]['tag-name']
    conditions = round1_tags.str.split('_').apply(lambda x: x[1])
    barcode_conditions = set(conditions)

    try:
        assert barcode_conditions == set(config_conditions)
    except AssertionError:
        print(f"{Fore.RED}The conditions in the barcode configuration don't match the ones listed in the pipeline configuration")
        print(f"{Fore.RED}The conditions you listed in the pipeline configuration are:")
        print(json.dumps(config_conditions, indent=4), "\n\n")

        print(f"{Fore.RED}However, only the following conditions were found in ROUND1 barcodes:")
        print(json.dumps(list(barcode_conditions), indent=4), '\n\n')
        return

    print(f"{Fore.GREEN}Barcodes look good :)")

def check_rounds_format(rounds_df, config_conditions):
    round1_tags = rounds_df[rounds_df['tag-name'].str.contains('ROUND1')]['tag-name']
    conditions = round1_tags.str.split('_').apply(lambda x: x[1])
    rounds_conditions = set(conditions)

    try:
        assert rounds_conditions == set(config_conditions)
    except AssertionError:
        print(f"{Fore.RED}The conditions in the split-pool rounds configuration don't match the ones listed in the pipeline configuration")
        print(f"{Fore.RED}The conditions you listed in the pipeline configuration are:")
        print(json.dumps(config_conditions, indent=4), "\n\n")

        print(f"{Fore.RED}However, the ROUND1 barcodes look like this:")
        print(json.dumps(list(rounds_conditions), indent=4), '\n\n')
        return

    print(f"{Fore.GREEN}Round formats look good :)")

@click.command()
@click.option('--snakemake_config', required=True, type=click.Path(exists=True))
@click.option('--fix', required=False, type=bool, default=False)   
def main(snakemake_config, fix):
    with open(snakemake_config, 'r') as f:
        smk_config = yaml.load(
                f.read(), 
                Loader=yaml.CLoader
        )
    
    # Names based on SPRITE GitHub Wiki 
    # https://github.com/GuttmanLab/sprite-pipeline/wiki/1.-Barcode-Identification#configuration-file
    try:
        barcode_ids = pd.read_csv(
                smk_config['bID'], 
                sep='\t', 
                skiprows=3, 
                skipfooter=1,
                engine='python',
                names=['tag-category', 'tag-name', 'tag-sequence', 'tag-error-tolerance']
        ) 
    except FileNotFoundError:
        print(f"{Fore.RED}Looks like the value for the `bID` doesn't point to a correct path.")
        print(f"{Fore.RED}Please check the `bID` key in {snakemake_config}")
        sys.exit()

    # Only `tag-position` and `tag-name` are used in the pipeline
    try:
        rounds_format = pd.read_csv(
                smk_config['rounds_format'], 
                sep='\t',
                names=['tag-position', 'tag-name', 'tag-sequence', 'tag-error-tolerance']
        ) 
    except FileNotFoundError:
        print(f"{Fore.RED}Looks like the value for the `rounds_format` doesn't point to a correct path.")
        print(f"{Fore.RED}Please check the `rounds_format` key in {snakemake_config}")
        sys.exit()

    check_barcodes(barcode_ids, smk_config['conditions'])
    check_rounds_format(rounds_format, smk_config['conditions'])


if __name__ == "__main__":
    main()
