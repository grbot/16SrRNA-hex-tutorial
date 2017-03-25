# 16S analysis part of Analysing the Microbiome Workshop 2017

The dataset we will be working are the practice dataset from the [H3ABioNet 16S rDNA diversity analysis SOP](http://www.h3abionet.org/tools-and-resources/sops/16s-rrna-diversity-analysis). The source data can be accessed [here](http://h3data.cbio.uct.ac.za/assessments/16SrRNADiversityAnalysis/practice) but for our purposes it is already on the cluster and stored here:`/global/mb/amw/dog_stool_samples`

The table below contains the metadata associated with the dog stool samples. There are three dogs which are treated with increased percentage of a compound in their diet: 5 different treatments (0-4, representing an increased percentage of a compound in their diet).

Sample | Dog | Treatment
------ | --- | ---------
Dog1 | B | 2
Dog2 | G | 3
Dog3 | K | 3
Dog8 | B | 4
Dog9 | G | 0
Dog10 | K | 4
Dog15 | B | 1
Dog16 | G | 4
Dog17 | K | 0
Dog22 | B | 3
Dog23 | G | 1
Dog24 | K | 2
Dog29 | B | 0
Dog30 | G | 2
Dog31 | K | 1

## Do some local setup
### Maybe pull this code here

### Setup some PATHS


### Add pipeline here

Raw reads --> FastQC 
    |-------> Merge
    
## Lets do some QC on the raw data

To be able to run all tools in the tutorial be sure the bionf module is loaded
```bash
module load bioinf
```
### Run FastQC
```bash
fastq_dir=/global/mb/amw/dog_stool_samples
fastqc_dir=/global/mb/amw/run/process/fastqc

mkdir $fastqc_dir
fastqc --extract -f fastq -o $fastqc_dir -t 12 fastq_dir/*
```
### Combine FastQC reports
```bash
fastqc_dir=/global/mb/amw/run/process/fastqc
/global/mb/amw/run/code/fastqc_combine/fastqc_combine.pl -v --out $fastqc_dir --skip --files "$fastqc_dir/*_fastqc"
```
Now lets view the reports.

## Run the UPARSE pipeline
### First rename the read headers so that they are compatible with the UPARSE pipline
```bash
renamed_dir="/global/mb/amw/run/process/usearch/renamed"
mkdir renamed
while read sid_fastq_pair; do sid=`echo $sid_fastq_pair | awk -F ' ' '{print $1}'`; fastq_r1=`echo $sid_fastq_pair | awk -F ' ' '{print $2}'`; fastq_r2=`echo $sid_fastq_pair | awk -F ' ' '{print $3}'`; fastq_r1_renamed=$renamed_dir"/"$(basename $fastq_r1); fastq_r2_renamed=$renamed_dir"/"$(basename $fastq_r2); /global/mb/amw/run/code/rename_fastq_headers.sh $sid $fastq_r1 $fastq_r2 $fastq_r1_renamed $fastq_r2_renamed;done < /global/mb/amw/run/code/sid.fastq_pair.list
```
This will take about 10 minutes to run. Lets have a look at the headers once done.

### Merge the paired reads
```bash
fastq_maxdiffs=3
merged_dir="/global/mb/amw/run/process/usearch/merged"
mkdir /global/mb/amw/run/process/usearch/merged

while read sid_fastq_pair; do sid=`echo $sid_fastq_pair | awk -F ' ' '{print $1}'`; fastq_r1=`echo $sid_fastq_pair | awk -F ' ' '{print $2}'`; fastq_r2=`echo $sid_fastq_pair | awk -F ' ' '{print $3}'`; fastq_r1_renamed=$renamed_dir"/"$(basename $fastq_r1); fastq_r2_renamed=$renamed_dir"/"$(basename $fastq_r2); usearch -fastq_mergepairs $fastq_r1_renamed -reverse $fastq_r2_renamed -fastq_maxdiffs $fastq_maxdiffs -fastqout $merged_dir"/"$sid".merged.fastq";done < /global/mb/amw/run/code/sid.fastq_pair.list
```







