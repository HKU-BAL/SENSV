## Installation

### Required Utilities
```
minimap2
samtools
pigz
bgzip
grabix
```

### Required Python packages
```
pysam
```

### Required R packages
```
ggplot2
magrittr
reshape2
optparse
dplyr
devtools
cttobin/ggthemr
```

Can install by running the following R scripts by super-users:
```
install.packages('ggplot2')
install.packages('magrittr')
install.packages('reshape2')
install.packages('optparse')
install.packages('dplyr')

install.packages('devtools')
devtools::install_github('cttobin/ggthemr')
```

## Required data files (for depth analysis)
```
1. depth_edward/mask_table.csv
2. depth_edward/ref_mask_22samples_absdepth_norm_sexchr_nomask.csv
```

## Configuration file
```
1. config.sh (config file for shell scripts)
2. config.ini (config file for python scripts)
```

1. In config.sh,
```
SCRIPT_PATH=
PIGZ=
GRABIX=
BGZIP=
```

2. In config.ini,
```
[common]
minimap2 = [path]/minimap2
smatools = [path]/samtools
ref = [path]/hs37d5.fa

[default_value]
minSvSize = 10000
maxSvSize = 0
svType = DEL

[realignSv]
aligner = minimap2 
```

## Execution
```
python realignSv.py -sample_name <sample_name> -fastq <fastq> -output_prefix <output_prefix> [-min_svsize <min_svsize>]
```

Please note that the script will detect the existence of resulting files of a particular step. It starts at the step where the reuslting files of that step do not exist. Hence, it is necessary to remove the resulting files manually to allow re-creation.

## Generation of deletion region/breakpoints
Step 1: 
Input: bam file 
Command: ./depth.hjyu.sh -b "bamfile path" 
Output: depth_df_sample.csv 

Step 2:
Input: output from step 1 
Command: Rscript normalize_sample.R -i depth_df_sample path 
Output: depth_sample_mask_norm.csv 

Step 3: 
Input: output from step 2 
Command: Rscript find_abnormal_region.R -i depth_sample_mask_norm.csv -r /home/hjyu/sv/ref_mask_20samples_1.csv -t 1.25 -w 24 -d 40 
Output: report.csv 
