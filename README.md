# SENSV

## Installation

### Step 1. Install required packages
```
# config for conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install minimap2 samtools pigz grabix pypy survivor

# for python
pip install pandas
pip install scipy
pip install pysam

```

### Step 2. Clone the repository

```
git clone https://github.com/HKU-BAL/SENSV.git
```

### Step 3. Fill in the paths for the required file

In config.ini, change the path for samtools and minimap2 if they are not available in PATH.
A GRCh37 reference file is also needed. If you do not have it in advance, you can download it with the following commands
```
curl ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > hs37d5.fa.gz
gzip -d hs37d5.fa.gz
samtools faidx hs37d5.fa
```
```
[common]
samtools = <path_of_samtools>
minimap2 = <path_of_minimap2>
ref_37 = <path_of_GRCh37_ref>
```
Please also make sure that the reference index is also available in <path_to_GRCh37_ref>.fai

In config.sh, change the path for bgzip, pigz and grabix if they are not available in PATH.

```
PIGZ = <path_of_pigz>
GRABIX = <path_of_grabix>
BGZIP = <path_of_bgzip>
```

In merge_sv.sh, change the path for survivor if it is not available in PATH.

```
survivor = <path_of_survivor>
```

## Usage

You will need a fastq file of the sample to run SENSV.

python SENSV.py [options]

```
Required Arguments:
-sample_name - Name of the sample
-fastq - The path to the reads, either gziped or raw
-output_prefix - Output prefix for all intermediate files and final output, preferably inside a folder.

Optional Arguments:
-min_sv_size - Minimum SV size to be called
-max_sv_size - Maximum SV size to be called
-target_sv_type - [DUP/DEL/DUP,DEL], default is DUP,DEL
```

After the process is done, the result bed2 file will be available in <output_prefix>_final.result.
