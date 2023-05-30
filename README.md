# Analysis and Code for SPIDR Paper

Repository to host all the analysis, code, and some data for the SPIDR paper.

## Setup
### Conda Environment 
Run 

```
conda env create -f environment.yml
```

then activate the environment with the following command

```
conda activate spidr
```

### Directory/File Structure
By default only the compressed versions of the annotation files are committed to GitHub because there is a strict 100mb file size limit imposed by GitHub. Here is what the directory tree looks like.

```
spidr/
├── README.md
├── analysis
│   ├── calculation.ipynb
│   ├── scripts
│   └── visualization.ipynb
├── annotated
│   ├── Summary_info_encode_Suppl_Data_4.xlsx
│   ├── compressed
│   └── uncompressed
├── environment.yml
├── figures
│   ├── heatmap_spidr_vs_encode_annotator_heatmap.png
│   └── heatmap_spidr_vs_encode_supplementary_heatmap.png
├── hpc
│   ├── submit_jobs.py
│   └── template.sh
└── output
    ├── encode_percent_by_region.csv
    ├── encode_supp_percent_by_region.csv
    └── spidr_percent_by_region.csv
```

To use the annotation files, decompress them with whatever utility you would like and move them into the uncompressed folder. On MacOS decompressing can be done by simply opening the zip files in Finder. The `spidr/annotated/uncompressed/` directory may not exist at first, but you can simply create at the same level as compressed (i.e. where it currently sits in the directory tree). You can re-run any analyses by running `calculation.ipynb` and `visualization.ipynb`. `calculation.ipynb` creates a matrix of percentages of each annotation type for the given analysis across RBPs. Each column represents an RBP and each row is an annotation type (e.g. 5utr). Script versions of each portion of code can be found in `spidr/analysis/scripts/`. The scripts are intended for "headless" usage (i.e. on a server or other environment where you don't have a screen).
