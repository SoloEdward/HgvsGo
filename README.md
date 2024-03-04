# HgvsGo

HgvsGo is a program designed for analyzing "c." and "p." HGVS (Human Genome Variation Society) notations for single nucleotide variations (SNVs) and small insertions/deletions (indels) after variant calling. It serves as an alternative to tools like snpEff and VEP.

## Why HgvsGo?

HgvsGo was specifically developed for clinical use, making it well-suited for medical applications. It provides accurate annotations by applying the 3' rule to both "c." and "p." HGVS notations. Additionally, HgvsGo offers fast performance, requiring only 20 seconds on a Mac to annotate over 1,300,000 variants (such as those from the ClinVar database downloaded from NCBI).

## How to Use HgvsGo

### Step 1: Download the Repository and Build HgvsGo

```
git clone https://github.com/SoloEdward/HgvsGo.git
cd ./HgvsGo/src/
mkdir build
cd build/
cmake ..
make
cd ../../
```

### Step 2: Download and Prepare the Human Genome

```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
gunzip GRCh37_latest_genomic.fna.gz
python parse_genome.py
```

### Step 3: Download RNA Sequences for All Transcripts

```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_rna.fna.gz
gunzip GRCh37_latest_rna.fna.gz
```

### Step 4: Run the Program

```
./src/build/HgvsGo ./GRCh37_latest_rna.fna.gz ./human.genome.fa ./refseq.select.hg19.parsed.txt demo.input.txt demo.output.txt
```

## Input Format

The input should be a tab-delimited text file with a minimum of four columns: "chrom", "pos", "ref" and "alt" as their headers. Please refer to the provided demo.input.txt file for an example.

## Output Format

The output is also a tab-delimited text file. It includes all the columns from the input file, followed by the annotated HGVS results for each variant. For cases where a variant overlaps with multiple transcripts, multiple lines will be generated in the output file for each input line. Please refer to the provided demo.output.txt file for an example.

## How Does HgvsGo Work?

For each SNV or small indel, HgvsGo first identifies the transcripts that overlap with the variant. It then calculates the "c." and "p." HGVS notations for each overlapping transcript. The "refseq.select.hg19.parsed.txt" file, obtained from the UCSC Table Browser and parsed, is used to find overlapping transcripts based on a variant. The GRCh37 genome is employed to apply the 3' rule to "c." HGVS notations. The "GRCh37_latest_rna.fna" file is used to translate DNA to amino acids and apply the 3' rule to "p." HGVS notations. Notably, variants located within 500 bp before the transcription start position of the "TERT" gene are also annotated as "TERT," as this region contains several hotspots on the TERT promoter.
