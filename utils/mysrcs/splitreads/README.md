# grab1chr 
C++ code to read and print fasta/fastq file:
Usage: grab1chr  <reads.fq/fa>  [minL] [ctg] [ipos] [epos] [len_from_end]

1.  minl: print only contigs > minl  (default=0)
2.  ctg: print only contig whose name includes the string ctg (default 'all': print every contigs)
3.  ipos/epos only if ctg defined: print contig ctg from position ipos to position epos
4. elen only if ctg defined: print last elen bases of contig ctg
