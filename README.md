# kmerSV
kmerSV: a visualization and annotation tool for structural variants (SVs). Users can input both a reference and a target sequence, which will generate a SV plot complete with genomic annotations. Additionally, kmerSV offers support for visualizing SVs against a pangenomic reference. Specifcally, users can utilize the pangenome tool such as odgi tool (https://github.com/pangenome/odgi) to extract all paths within a specified region of the pangenome graph. Once these samples are extracted and saved in a fasta file, they can be seamlessly integrated into the kmerSV for visualization.

## Installation
Download repository
```
git clone https://github.com/sgtc-stanford/kmerSV.git
```

kmerSV is optimized for:
- **Python** version 3.9.0 or higher
- **numpy** version 1.20.1 or higher
- **matplotlib** version 3.7.0 or higher
- **pandas** version 1.2.0 or higher

To directly install these dependencies, run the following command:

```
pip install -r requirements.txt
```
## General Usage of kmerSV

### 1. **Paired Sequences**

To generate an SV plot for a pair of sequences, execute the following command with `kmerSV_plot.py`:

```bash
python kmerSV_plot.py [-h] -i INPUT -o OUTPUT -k KMER [-a ANNOTATION] [-c CHR] [-s START] [-e END]
```
#### Optional Arguments:

- `-h, --help` : Show help message and exit.
- `-r REFERENCE, --reference REFERENCE` : Define the reference sequence.
- `-i INPUT, --input INPUT` : Define the target sequence.
- `-o OUTPUT, --output OUTPUT` : Specify the output SV plot.
- `-k KMER, --kmer KMER` : Specify the kmer length.
- `-a ANNOTATION, --annotation ANNOTATION` : Provide the Genomic Annotation File.
- `-c CHR, --chr CHR` : Define the chromosome region.
- `-s START, --start START` : Set the start position.
- `-e END, --end END` : Set the end position.
- `-f TEXT, --text TEXT` : Output a text file with detailed information.

### 2. **Pangenome Reference**

To generate an SV plot using a pangenome reference, execute the following command with `kmerSV_pangenome.py`:

```bash
python kmerSV_pangenome.py [-h] -i INPUT -o OUTPUT -k KMER [-a ANNOTATION] [-c CHR] [-s START] [-e END]
```
#### Optional Arguments:

- `-h, --help`                                : Show help message and exit.
- `-i INPUT, --input INPUT`                   : Provide the input FASTA file with all sample paths.
- `-o OUTPUT, --output OUTPUT`                : Specify the output SV plot using pangenome reference.
- `-k KMER, --kmer KMER`                      : Specify the kmer length.
- `-a ANNOTATION, --annotation ANNOTATION`    : Provide the Genomic Annotation File.
- `-c CHR, --chr CHR`                         : Define the chromosome region.
- `-s START, --start START`                   : Set the start position.
- `-e END, --end END`                         : Set the end position.

## Example Usage of kmerSV

This section presents the example usages of `kmerSV`. We will delve into its application for both paired sequences and a pangenome reference.

### Datasets

All datasets used in these examples are housed in the `data/test` directory:

1. **chr5.fasta**:
   - Description: Represents the region `GRCH38.chr5:141169588-141184761`.

2. **hg02080.fasta**:
   - Description: Contains data for the region `HG02080#2#JAHEOV010000100.1:14419686-14450465`.

3. **chr19_pan.fasta**:
   - Description: This file was created using the `odgi` tool (version 0.6.2). It has been used to extract paths within the region `GRCh38.chr19:4512541 - 4513161` from the [pangenome graph](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr19.full.og).
   - Details: The file aggregates a total of 90 haplotypes. For our purposes, we use the first entry in the fasta file (`GRCh38.chr19:4512541 - 4513161`) as the primary pangenome reference.

### Annotations

The associated annotation data, found in `data/annotation/annotation.bed`, originates from the Matched Annotation from NCBI and EMBL-EBI (MANE).


### Visualization and Annotation with Paired Sequences

To visualize and annotate the SVs using a pair of sequences, use the following command:

```bash
python kmerSV_plot.py -r data/test/chr5.fasta -i data/test/hg02080.fasta -o SV_plot.png -a data/annotation/annotation.bed -c chr5 -s 141169588 -e 141184761 -k 31 -f SV_info.txt
```

**Output:**
- **SV Plot:** SV_plot.png
- **SV Details:** SV_info.txt


To visuliaze and annotate the SVs with the pangenome reference, use the following command::

```bash
python3 kmerSV_pangenome.py -i data/test/chr19_pan.fasta -o SV_pan_plot.png -c chr19 -s 4512541 -e 4513161 -k 31
```
**Output:**
- **SV Plot:** SV_pan_plot.png
