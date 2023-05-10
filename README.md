# HgvsGo

Program For Analyzing "c." and "p." hgvs for SNV and Small INDEL after variant calling.

## Why HgvsGo?

Firstly, HgvsGo is invented for clinical usage.
Secondly, HgvsGo is accurate. It applies 3' rule to both "c." and "p." hgvs.
Thirdly, HgvsGo is fast. It takes only 20 seconds on my mac to annotate >1,300,000 variants (variants from clinvar database which is download from ncbi).

## How to use it? 

### Step1: Download The repo
```
git clone https://github.com/SoloEdward/HgvsGo.git
cd ./HgvsGo/
```

### Step2: Download and Prepare the Human Genome
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

### Step5: Run The Program Using docker (optional, only when ./HgvsGo not work)
```
docker pull ubuntu:22.04
docker build -t hgvsgo:1.0.0 .
docker run -v `pwd`:/home/ hgvsgo:1.0.0  ./HgvsGo GRCh37_latest_rna.fna human.genome.fa refseq.select.hg19.parsed.txt demo.input.txt demo.output.txt
```

## Input Format 
Input is a tab-delimited txt file. It must contain at least four columns with chrom,pos,ref,alt as their header. See demo.input.txt as an example.

## Output Format
Output is also a tab-delimited txt file. It will contian all the columns in the input, and will write the annotated hgvs result after each line. See demo.output.txt as an example. Notice that for each line in the input, if the variant in this line is overlaped with multiple transcripts, then multiple hgvs results will be generated, thus for each line in the input file, one or multiple line will be generated in the output file.

## How does HgvsGo work?
For each snv or small indel, HgvsGo will firstly find those transcripts overlapped with this variant. Then, for each overlapped transcript, HgvsGo will calculate "c." and "p." hgvs. The "refseq.select.hg19.parsed.txt" is download from ucsc table browser and parsed, it is used for find overlapping transcripts given a variant. The GRCh37 genome is used for applying 3' rule to "c." hgvs. The "GRCh37_latest_rna.fna" is used for translating dna to amino acids and applying 3' rule to "p." hgvs. 
Specially, for the gene "TERT", those variants located <500 bp before the transcript start position of "TERT" will also be annotated as "TERT", since there is several hot spots on TERT promoter region.

## System Dependency
- HgvsGo is written in C++ and build with gcc version 9.4.0. 
