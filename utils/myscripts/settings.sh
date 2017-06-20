#### Parameters to set: ####

# the draft assembly contigs/scaffolds are shred in pieces of how many base-pairs? default=10000 bp
shred=10000
# alignment < noise base-pairs will be considered noise and not show in ACT. default=30000 bp. (for smaller genomes reduce up to 5000 bp)
noise=30000
minid=80

# lfs jobs parameters:
lfsjobs=1  # lfsjobs=0 not implemented yet
myqueue=normal
myjobmem=5000
myncpus=15
maxjobs=5  #maximum number of jobs to run at a time
########################



myforACT=MYFORACT
fullpathref=MYREF
notshred=MYDRAFT
dir=MYDESTDIR

scriptdir=MYFORACT/utils/myscripts
srcdir=MYFORACT//utils/mysrcs
runalign=$scriptdir/runalign.sh
runprep=$scriptdir/runprep.sh
runreport=$scriptdir/runreport.sh

wdir=whole

folder=$wdir\_$shred
refdir=$dir/ref
ref=ref.fasta   #$(basename $fullpathref) 
outdir=$dir/$folder
aldir=$outdir/aligns
fastadir=$outdir/fasta
workdir=$outdir/inter
splitdir=$fastadir/split
finaldir=$outdir/unique

# output files:
oalign=align.output
oprep=prepfiles.output

firstal=first  # alignment of shreded initial fasta
secal=second   # alignment of shreaded forward fasta
thirdal=third  # alignment of  shreaded forward fasta, fixed shred position in scaffold 

notshredname=$(basename $notshred)
shreddraft=shred$shred\_$notshredname
forwnotshred=forw$notshredname
singlefolder=$fastadir/$(basename $forwnotshred .fasta)_fastas
forwshred=$fastadir/shred$shred\_$forwnotshred 
forwal=$aldir/$secal.al
finalal=third.al

foractfa=$finaldir/foract$noise.fasta
foractal=$finaldir/foract$noise.al

