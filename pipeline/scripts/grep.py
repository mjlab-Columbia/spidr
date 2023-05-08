import yaml
import os
import subprocess

with open("../envs/sprite.yaml", "r") as f:
    packages = yaml.safe_load(f)

for package in packages['dependencies']:
    package = package.split("=")[0]
    if package not in ["python", "pigz", "samtools", "seqkit"]:
        output = subprocess.check_output(f"grep -r 'import {package}' ./")
        print(output.decode().strip())
