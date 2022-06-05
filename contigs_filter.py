#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

parser = argparse.ArgumentParser(description='This script allow to filter SPAdes contigs/scaffolds sequences and obtain cotings in specified range of reads coverag and GC percentage')
parser.add_argument('-i', '--infile', help='.fasta file with contigs/scaffolds sequences from SPAdes', required=True)
parser.add_argument('-cv', '--cov', help='reads coverage range between Colon ":" mark', required= True)
parser.add_argument('-gc', '--gc', help='gc percentage range between Colon ":" mark', required= True)
parser.add_argument('-o', '--outfolder', help='outfolder', required= True)

args = parser.parse_args()
infile = args.infile
gc_per = args.gc
gc_min = float(gc_per.split(':')[0])
gc_max = float(gc_per.split(':')[1])
cov = args.cov
cov_min = float(cov.split(':')[0])
cov_max = float(cov.split(':')[1])
outfile= args.outfolder

sed_d ={}
for record in SeqIO.parse(infile, "fasta"):    
    seq_cov = record.id.split('_')[5]
    seq_gc = str(GC(record.seq))
    sed_d.update( {record.id : [seq_cov, seq_gc]} )

select_id=[]
for key, value in sed_d.items():
    if  cov_min <= float(value[0]) <= cov_max and gc_min <= float(value[1]) <= gc_max:
        select_id.append(key)

with open(outfile, "w") as f:
    seqs_DBs = SeqIO.parse(open(infile),'fasta')
    for seq in seqs_DBs:
        if seq.id in select_id:
            SeqIO.write([seq], f, "fasta")
f.close()



