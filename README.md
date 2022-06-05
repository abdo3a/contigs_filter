# contigs_filter

contigs_filter is a python program to to filter SPAdes contigs/scaffolds sequences and obtain cotings in specified range of reads coverage and GC percentages.

### USAGE

./contigs_filter.py -i [contigs/scaffolds sequences file generated by SPAde] -cv [reads coverage range between Colon ":" mark' (e.g. 0:100)] -gc [gc percentage range between Colon ":" mark' (e.g. 1:40)] -o [output file]

* contigs_filter require Python 3 and Biopython and there are no external dependencies required.


