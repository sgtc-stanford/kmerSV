# kmerSV
kmerSV: a visualization and annotation tool for structural variants (SVs). Users can input both a reference and a target sequence, which will generate a SV plot complete with genomic annotations. Additionally, kmerSV offers support for visualizing SVs against a pangenomic reference. Specifcally, users can utilize the pangenome tool such as odgi tool (https://github.com/pangenome/odgi) to extract all paths within a specified region of the pangenome graph. Once these samples are extracted and saved in a fasta file, they can be seamlessly integrated into the kmerSV for visualization.

## Installation Requirements

kmerSV is optimized for:
- **Python** version 3.9.0 or higher
- **numpy** version 1.20.1 or higher
- **matplotlib** version 3.7.0 or higher
- **pandas** version 1.2.0 or higher

To seamlessly install these dependencies, run the following command:

```
pip install -r requirements.txt
```
## Usage

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

