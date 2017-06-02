myforACT=MYFORACT
fullpathref=MYREF
notshred=MYDRAFT
dir=MYDESTDIR

scriptdir=MYFORACT/utils/myscripts
srcdir=MYFORACT//utils/mysrcs
runalign=$scriptdir/runalign.sh
runprep=$scriptdir/runprep.sh

wdir=whole
shred=10000

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

foractfa=$finaldir/foract.fasta
foractal=$finaldir/foract.al

