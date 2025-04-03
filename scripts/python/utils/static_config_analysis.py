import sys
import os
import pandas as pd
import yaml
import click
import json

# Initialize colorama
init(autoreset=True)


def prepend_to_file(filename, lines):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)

        if isinstance(lines, list):
            lines = '\n'.join(lines)
            f.write(lines.rstrip('\r\n') + '\n' + content)


def fix_barcodes(barcode_df):
    df = barcode_df.copy()

    # Issue: ROUND1 tags don't contain a condition in the middle
    # Solution: add a single condition called BLANK
    def add_blank_tag(tag):
        split_tags = tag.split('_')
        split_tags.insert(1, 'BLANK')
        return '_'.join(split_tags)

    round1_rows = df[df['tag-name'].str.contains('ROUND1')].copy()
    round1_rows['tag-name'] = round1_rows['tag-name'].apply(lambda x: add_blank_tag(x))

    mask = df['tag-name'].str.contains('ROUND1')
    df.loc[mask, 'tag-name'] = round1_rows['tag-name'].values

    # TODO: Add more automatic fixes

    return df


def fix_rounds_format(rounds_df):
    df = rounds_df.copy()

    # Issue: ROUND1 tags don't contain a condition in the middle
    # Solution: add a single condition called BLANK
    def add_blank_tag(tag):
        split_tags = tag.split('_')
        split_tags.insert(1, 'BLANK')
        return '_'.join(split_tags)

    round1_rows = df[df['tag-name'].str.contains('ROUND1')].copy()
    round1_rows['tag-name'] = round1_rows['tag-name'].apply(lambda x: add_blank_tag(x))

    mask = df['tag-name'].str.contains('ROUND1')
    df.loc[mask, 'tag-name'] = round1_rows['tag-name'].values

    # TODO: Add more automatic fixes

    return df


def check_barcodes(barcode_df, config, fix, barcodes_path):
    round1_tags = barcode_df[barcode_df['tag-name'].str.contains('ROUND1')]['tag-name']
    conditions = round1_tags.str.split('_').apply(lambda x: x[1])
    barcode_conditions = set(conditions)

    config_conditions = config['conditions']

    try:
        assert barcode_conditions == set(config_conditions)
    except AssertionError:
        if not fix:
            print(f"{Fore.RED}The conditions in the barcode configuration don't match the ones listed in the pipeline configuration")
            print(f"{Fore.RED}The conditions you listed in the pipeline configuration are:")
            print(json.dumps(config_conditions, indent=4), "\n\n")

            print(f"{Fore.RED}However, only the following conditions were found in ROUND1 barcodes:")
            print(json.dumps(list(barcode_conditions), indent=4), '\n\n')
            print(f"{Fore.YELLOW}If you would like to try to automtically fix these issues, rerun the command with --fix True")
            return None
        else:
            fixed_barcode_df = fix_barcodes(barcode_df)

            barcodes_basename = os.path.basename(barcodes_path)
            new_barcodes_basename = ''.join(barcodes_basename.split('.')[:-1]) + "_FIXED" + ".txt"
            new_barcodes_path = os.path.join(os.path.dirname(barcodes_path), new_barcodes_basename)

            fixed_barcode_df.to_csv(new_barcodes_path, sep='\t', index=False, header=False)

            with open(config['bID'], 'r') as f:
                preamble = f.read().split('\n')[:3]

            with open(new_barcodes_path, 'r') as f:
                original_content = f.read()

            with open(new_barcodes_path, 'w') as f:
                for line in preamble:
                    f.write(line + '\n')

                f.write(original_content)

            return new_barcodes_path

    print(f"{Fore.GREEN}Barcodes look good :)")


def check_rounds_format(rounds_df, config, fix, rounds_path):
    round1_tags = rounds_df[rounds_df['tag-name'].str.contains('ROUND1')]['tag-name']
    conditions = round1_tags.str.split('_').apply(lambda x: x[1])
    rounds_conditions = set(conditions)

    config_conditions = config['conditions']

    try:
        assert rounds_conditions == set(config_conditions)
    except AssertionError:
        if not fix:
            print(f"{Fore.RED}The conditions in the barcode configuration don't match the ones listed in the pipeline configuration")
            print(f"{Fore.RED}The conditions you listed in the pipeline configuration are:")
            print(json.dumps(config_conditions, indent=4), "\n\n")

            print(f"{Fore.RED}However, only the following conditions were found in ROUND1 barcodes:")
            print(json.dumps(list(rounds_conditions), indent=4), '\n\n')
            print(f"{Fore.YELLOW}If you would like to try to automtically fix these issues, rerun the command with --fix True")
            return None
        else:
            fixed_rounds_df = fix_rounds_format(rounds_df)

            rounds_basename = os.path.basename(rounds_path)
            new_rounds_basename = ''.join(rounds_basename.split('.')[:-1]) + "_FIXED" + ".txt"
            new_rounds_path = os.path.join(os.path.dirname(rounds_path), new_rounds_basename)

            fixed_rounds_df.to_csv(new_rounds_path, sep='\t', index=False, header=False)
            return new_rounds_path

    print(f"{Fore.GREEN}Round formats look good :)")


def precheck_config(smk_config):
    barcode_config_path = smk_config['bID']
    rounds_format_path = smk_config['rounds_format']
    barcode = os.path.basename(barcode_config_path).split('.')[0]
    rounds = os.path.basename(rounds_format_path).split('.')[0]

    if barcode.endswith('_FIXED') or rounds.endswith('_FIXED'):
        print(f'{Fore.YELLOW}Ignoring --fix because barcode config and/or rounds format have been fixed')
        print(f"{Fore.YELLOW}barcodes fixed: {barcode.endswith('_FIXED')}")
        print(f"{Fore.YELLOW}rounds format fixed: {rounds.endswith('_FIXED')}")

        return False

    return True


@click.command()
@click.option('--snakemake_config', required=True, type=click.Path(exists=True), show_default=True)
@click.option('--fix', required=False, type=bool, default=False, show_default=True)
def main(snakemake_config, fix):

    # Load the YAML
    with open(snakemake_config, 'r') as f:
        smk_config = yaml.load(
            f.read(),
            Loader=yaml.CLoader
        )

    # Avoid fixing a file that's already been fixed
    if fix:
        fix = precheck_config(smk_config)

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

    barcode_results = check_barcodes(barcode_ids, smk_config, fix, smk_config['bID'])
    rounds_results = check_rounds_format(rounds_format, smk_config, fix, smk_config['rounds_format'])

    if barcode_results is not None:
        print(f"{Fore.GREEN}Fixed barcode config file located at: {barcode_results}")
        smk_config['bID'] = barcode_results

    if rounds_results is not None:
        print(f"{Fore.GREEN}Fixed rounds format file located at: {rounds_results}")
        smk_config['rounds_format'] = rounds_results

    if (barcode_results is not None) or (rounds_results is not None):
        smk_config['conditions'] = ['BLANK']

        print(f"{Fore.GREEN}Updated {snakemake_config} with new values\n")
        with open(snakemake_config, 'w') as f:
            yaml.dump(smk_config, f)


if __name__ == "__main__":
    main()
