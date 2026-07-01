import sys
import re
import csv
import click

def stream_errors(file_stream, include_header, sort_by_rule):
    fieldnames = ['name of the rule', 'jobid', 'log']
    
    # Configure csv.DictWriter to stream directly to stdout
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter='\t')
    
    if include_header:
        writer.writeheader()

    # Regex patterns to capture relevant fields within an error block
    rule_pattern = re.compile(r"Error in rule\s+(\w+):")
    jobid_pattern = re.compile(r"jobid:\s*(\d+)")
    log_pattern = re.compile(r"log:\s*([^\s\n\(\)]+)")

    current_rule = None
    current_jobid = None
    current_log = None
    inside_error_block = False

    # If sorting is enabled, cache results in a list; otherwise, write instantly
    collected_errors = []

    def handle_row(rule, jobid, log):
        row = {
            'name of the rule': rule or '',
            'jobid': jobid or '',
            'log': log or ''
        }
        if sort_by_rule:
            collected_errors.append(row)
        else:
            writer.writerow(row)
            sys.stdout.flush()

    for line in file_stream:
        # Check for the start of a new error block
        rule_match = rule_pattern.search(line)
        if rule_match:
            # If we were already tracking a block, emit/collect it before starting the next one
            if inside_error_block and (current_rule or current_jobid or current_log):
                handle_row(current_rule, current_jobid, current_log)
            
            # Reset and track the new block
            current_rule = rule_match.group(1)
            current_jobid = None
            current_log = None
            inside_error_block = True
            continue

        if inside_error_block:
            # Look for jobid
            jobid_match = jobid_pattern.search(line)
            if jobid_match:
                current_jobid = jobid_match.group(1)
                continue
            
            # Look for log path
            log_match = log_pattern.search(line)
            if log_match:
                current_log = log_match.group(1)
                continue

            # End of a Snakemake block typically indicated by empty line or specific progress text
            if line.strip() == "" or "Finished job" in line:
                handle_row(current_rule, current_jobid, current_log)
                current_rule = current_jobid = current_log = None
                inside_error_block = False

    # Catch any trailing error block at the absolute end of the stream
    if inside_error_block and (current_rule or current_jobid or current_log):
        handle_row(current_rule, current_jobid, current_log)

    # If sorting is requested, sort by 'name of the rule' and output all rows now
    if sort_by_rule:
        collected_errors.sort(key=lambda x: x['name of the rule'])
        for row in collected_errors:
            writer.writerow(row)
        sys.stdout.flush()

@click.command()
@click.argument('log_file', type=click.File('r', encoding='utf-8'), default='spidr.log')
@click.option('--header', '-h', is_flag=True, default=False, help='Include TSV column headers.')
@click.option('--sort', '-s', is_flag=True, default=False, help='Sort rows alphabetically by the name of the rule.')
def main(log_file, header, sort):
    """Extracts Snakemake rule errors from LOG_FILE (default: spidr.log) into a streamable TSV."""
    try:
        stream_errors(log_file, header, sort)
    except BrokenPipeError:
        # Gracefully handle downstream termination (like piping into `head`)
        sys.stderr.close()

if __name__ == '__main__':
    main()
