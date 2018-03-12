#### Parameters to set: ####

aligner=minimap2

# alignment < noise base-pairs will be considered noise and not show in ACT. default=30000 bp. (for smaller genomes reduce up to 5000 bp)
noise=2   # 1 == 0.1% of contig   #30000
debug=1

# Aligner
mysmalt=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT//utils/mysrcs/mylibs/smalt-0.7.4/smalt_x86_64
myminimap2=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT//utils/mysrcs/mylibs/minimap2/minimap2
mybwa=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT//utils/mysrcs/mylibs/bwa/bwa



# lfs jobs parameters:
lfsjobs=0  # 
myqueue=normal
maxjobs=100  #maximum number of jobs to run at a time



###### Mis-joint analysis using synteny
min_len=300000     # absolute min_lenght to consider as possible misjoint
min_len_max=500000  # min_lenght for a major al  (== ignore scaffold if major alignment block is < min_len_max
splitlen=5000000


##### the draft assembly contigs/scaffolds are shred in pieces of how many base-pairs?
shred=10000000
minid=10
myjobmem=8000
myncpus=1

if [[ $aligner == "smalt" ]]; then
	shred=10000
	minid=70
	myjobmem=5000
	myncpus=15
fi

########################
myforACT=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT
fullpathref=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT/test_ecoli/forACT_testdata/Escherichiacoli-K-12.fasta
notshred=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT/test_ecoli/forACT_testdata/draft.fasta
dir=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT/test_ecoli/forACT

scriptdir=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT/utils/myscripts
srcdir=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/forACT//utils/mysrcs
runalign=$scriptdir/runalign.sh
runprep=$scriptdir/runprep.sh
runreport=$scriptdir/runreport.sh
runmisjoints=$scriptdir/runlocate_misjoint.sh
runcirclize=$scriptdir/runcirclize.sh
wdir=$aligner

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
contigsizes=$fastadir/contig.sizes

name_fornoise=nonoise0.$noise
foractfa=$finaldir/foract$name_fornoise\_minid$minid.fasta
foractal=$finaldir/foract$name_fornoise\_minid$minid.al

### misjoint analysis ###
misals_file=$outdir/report/misjoints_details_$name_fornoise\_minid$minid\_$min_len.txt
chrass_file=$outdir/report/chr_assignment_$name_fornoise\_minid$minid\_$min_len.txt
bedfile=$outdir/report/bedfile_$name_fornoise\_minid$minid\_$min_len.txt
circl_fig=circos_$name_fornoise\_minid$minid\_$min_len
d=$(basename $dir)
title_fig=$d\_$aligner\_minid$minid\_$min_len
unplaced="unplaced"

