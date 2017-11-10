#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

fastasplit=$1
outname=$2
alsplitdir=$aldir/split_$outname
option=$3  #bwa

mkdir -p $aldir
cd $aldir
mkdir -p $alsplitdir 
cd $alsplitdir

ii=0
temponame=split$ii\_$outname
subqueue="s#MYQUEUE#$myqueue#g"
submem="s#MYJOBMEM#$myjobmem#g"
subcpu="s#MYCPUS#$myncpus#g"
waitfor="split"

if [[ ! -f  $refdir/ref.fasta.bwt ]]; then
  $mybwa index $refdir/ref.fasta
fi

ii=0
if [ ! -f $alsplitdir/$temponame.sam ]; then	    
    for fa in $fastasplit/*; do
	temponame=split$ii\_$outname
	
	if [ ! -f $temponame.out ]; then
	    command=`echo "$mybwa mem -t 5 $refdir/ref.fasta $fa | samtools view -q 1 -F 4 -F 256 -F 0x900 | grep -v XA:Z | grep -v SA:Z   - > $temponame.sam"`
	    echo $command > runalign_$ii.sh
	    

	    
	    if [[ $lfsjobs == 1 ]]; then
	        sed $subqueue $scriptdir/gensplital.sh | sed $submem | sed $subcpu > gensplital.sh
        	chmod +x gensplital.sh
	    	njobs=`bjobs | wc -l`	   	    
	    	#while [ $njobs -ge $maxjobs ]; do
       		#    sleep 30
		#    njobs=`bjobs | wc -l`
	    	#done	   	    
	    	#./gensplital.sh  runalign_$ii.sh
	    	#sleep 1
	    else
		chmod +x ./runalign_$ii.sh
		./runalign_$ii.sh
	    fi
	fi
	ii=$(($ii+1))
    done
fi

echo all jobs launched
exit
sleep 30
if [[ $lfsjobs == 1 ]]; then

	njobs=`bjobs | grep $waitfor | wc -l`
	echo $waitfor $njobs
	if [[ $njobs > 0 ]]; then
	   while [ $njobs -ge 1 ]; do
		sleep 30
		njobs=`bjobs | grep $waitfor  | wc -l`
    	    done
    	    sleep 5
	fi
	sleep 60

	t1=`ls $alsplitdir/runalign_* | wc -l`
	t2=`grep Success $alsplitdir/outs/out.* | wc -l`
	if [[ $t1 != $t2 ]]; then
	    echo " Some alignment jobs possibly failed? Check in $alsplitdir/outs/"
	fi
fi
