# HgvsGo

Program For Analyzing Hgvs for SNV and Small INDEL

### Get The repo
```
git clone git@github.com:SoloEdward/HgvsGo.git
cd ./HgvsGo/
```

### Prepare Human Genome
```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
gunzip  GRCh37_latest_genomic.fna.gz 
python parse_genome.py
```

### Download rna Sequence for all transcripts 
```
https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_rna.fna.gz
gunzip  GRCh37_latest_rna.fna.gz
```

### Run The Program
```
./HgvsGo ./GRCh37_latest_rna.fna.gz ./human.genome.fa ./refseq.select.hg19.parsed.txt demo.input.txt demo.output.txt
```
