# SENSV

## Installation

### Option 1. Bioconda
```
# make sure channels are added in conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment named "sensv-env"
conda create -n sensv-env -c bioconda sensv python=3.7
conda activate sensv-env

# post setup
cd $(dirname `which sensv`) && make post_link && cd -

# run sensv like this afterwards
sensv --help
```

### Option 2. Build and anaconda virtual environment step by step
```
# add conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environemnt named "sensv-env"
conda create -n sensv-env python=3.7 minimap2=2.17 samtools=1.7 pigz=2.3.4 grabix=0.1.8 pypy3.6=7.3.0 survivor=1.0.6 pandas=1.0.1 scipy=1.4.1 pysam=0.15.3 htslib=1.10.2 intervaltree=3.0.2 vcflib=1.0.0 pyfaidx=0.5.8
conda activate sensv-env

# clone repo
git clone https://github.com/HKU-BAL/SENSV.git

# setup sensv
cd SENSV
make
export PATH=`pwd`":$PATH"

# download data files and uncompress it
curl http://www.bio8.cs.hku.hk/sensv/data_v1.0.1.tar.gz --output data.tar.gz
tar -xf data.tar.gz

# run sensv like this afterwards
sensv --help
```

## After installation

### Download GRCh37 reference file
```
curl ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > hs37d5.fa.gz

gzip -d hs37d5.fa.gz

# make sure that the reference index is also available in <path_to_GRCh37_ref>.fai
samtools faidx hs37d5.fa
```


## Usage

You will need a fastq file of the sample to run SENSV.

sensv [options]

```
Required Arguments:
  -sample_name - Name of the sample
  -fastq - The path to the reads, either gziped or raw
  -ref - Reference fasta file
  -output_prefix - Output prefix for all intermediate files and final output, preferably inside a folder.

Optional Arguments:
  -min_sv_size - Minimum SV size to be called
  -max_sv_size - Maximum SV size to be called
```

After the process is done, the result bed2 file will be available in <output_prefix>_final.vcf.
