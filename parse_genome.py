import re
#les GRCh37_latest_genomic.fna |grep '>'  |grep Primary  |grep chromosome  |grep -v unlocalized

def parse_human_genome_fasta(fasta_file):
    results = []
    seq_name = None
    seqs = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                if seq_name is not None:
                    results.append((seq_name, seqs))
                seq_name = line
                seqs = []
            else:
                seqs.append(line)
        results.append((seq_name, seqs))
    return results

with open('human.genome.fa', 'w') as f:
    for seq_name, seqs in parse_human_genome_fasta('GRCh37_latest_genomic.fna'):
        if 'Primary' in seq_name and 'chromosome' in seq_name and 'unlocalized' not in seq_name:
            match_obj = re.match(r'.+chromosome\s([0-9XY]+).+', seq_name)
            chrom =  'chr' + match_obj.group(1)
            f.write('>{}\n'.format(chrom))
            f.write(''.join(map(lambda x: x.upper(), seqs)))
      
