from os import path, listdir, system
from argparse import PARSER, ArgumentParser

# Argument parser to decide where to look for inputs, outputs, GTF database, and logs
parser = ArgumentParser()
parser.add_argument("--input", help='Directory containing bed files')
parser.add_argument("--output", help='Directory to write outputs of annotator')
parser.add_argument("--logs", help='Directory to write log files')
parser.add_argument("--gtfdb", help='Location of gtfdb')
args = parser.parse_args()

if __name__ == "__main__":
    # Turn command line arguments into absolute paths
    input = path.abspath(args.input)
    output = path.abspath(args.output)
    gtfdb = path.abspath(args.gtfdb)
    logs = path.abspath(args.logs)
    
    # Create lists of absolute paths for inputs and outputs
    bed_files = [path.join(input, f) for f in listdir(input)]
    output_files = [path.join(output, path.basename(bed).split(".")[:-1][0]) + ".out" for bed in bed_files]

    # Iterate through each combination of bed file and output file and submit a job
    for bed, out in zip(bed_files, output_files):
        # Take the bed file name, strip the .bed extension and replace it with .log
        logfile = path.basename(bed).split('.')[:-1][0] + ".log"
        
        # Turn logfile into an absolute path
        logpath = path.join(logs, logfile)
        
        # Build command in components
        command = "sbatch "
        command += f"-A mjlab "
        command += f"-o {logpath} "
        command += f"--export=INPUT={bed},OUTPUT={out},GTFDB={gtfdb} "
        command += f"template.sh"
        system(command)
