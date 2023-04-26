# HgvsGo

Program For Analyzing Hgvs for SNV and Small INDEL

## How to use it? 

### Step1: Download The repo
```
git clone https://github.com/SoloEdward/HgvsGo.git
cd ./HgvsGo/
```

### Step2: Downlaod and Prepare the Human Genome
```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
gunzip  GRCh37_latest_genomic.fna.gz 
python parse_genome.py
```

### Step3: Download RNA Sequences for all transcripts 
```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_rna.fna.gz
gunzip  GRCh37_latest_rna.fna.gz
```

### Step4: Run The Program
```
chmod 755 ./HgvsGo
./HgvsGo ./GRCh37_latest_rna.fna.gz ./human.genome.fa ./refseq.select.hg19.parsed.txt demo.input.txt demo.output.txt
```

## Input Format 
Input is a tab-delimited txt file. It must contain at least four columns with chrom,pos,ref,alt as their header. See demo.input.txt as an example.

## Output Format
Output is also a tab-delimited txt file. It will contian all the columns in the input, and will write the hgvs result after each line. See demo.output.txt as an example.

