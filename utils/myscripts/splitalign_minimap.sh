#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

fastasplit=$1
outname=$2
alsplitdir=$aldir/split_$outname
option=$3  #minimap

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

echo index
time $myminimap2 -x asm5 -d ref.mmi $refdir/ref.fasta

if [[ $inparallel == "yes" ]]; then  # waitfor is not working
    echo inparallel
    ii=0
    tt=0
    rm -f runalign_inpar.sh
    if [ ! -f $alsplitdir/$temponame.out ]; then	    
	for fa in $fastasplit/*; do
	    temponame=split$ii\_$outname
	    	 
	    if [ ! -f $temponame.out ]; then
		if [[ $option == "c" ]]; then		
		    command=`echo "$myminimap2" -c -x asm5 ref.mmi  $fa "> $temponame.out &"`
		else
		    command=`echo "$myminimap2" -x asm5 ref.mmi  $fa "> $temponame.out &"`
		fi
		#if [[ $ii == 0 ]]; then 
		echo $command >> runalign_inpar.sh; 
		if [[ $tt == 9 ]]; then
		    echo "wait" >> runalign_inpar.sh
		    tt=0
		fi
                #fi
	    fi

	    tt=$(($tt+1))
	    ii=$(($ii+1))
	done
	echo "wait" >> runalign_inpar.sh
	subcpu="s#MYCPUS#10#g"
	sed $subqueue $scriptdir/parsplit.sh | sed $submem | sed $subcpu > parsplit.sh
	chmod +x parsplit.sh
	
	./parsplit.sh runalign_inpar.sh
	waitfor="inpar"
    fi
   exit
else

    ii=0
    if [ ! -f $alsplitdir/$temponame.out ]; then	    
	for fa in $fastasplit/*; do
	    temponame=split$ii\_$outname
	
	    if [ ! -f $temponame.out ]; then
		if [[ $option == "c" ]]; then	
		    command=`echo "$myminimap2" -c -x asm5  ref.mmi  $fa "> $temponame.out"`
		else
		    command=`echo "$myminimap2" -x asm5  ref.mmi  $fa "> $temponame.out"`
		fi
		echo $command > runalign_$ii.sh
		
		
		if [[ $lfsjobs == 1 ]]; then
	            sed $subqueue $scriptdir/gensplital.sh | sed $submem | sed $subcpu > gensplital.sh
        	    chmod +x gensplital.sh
	    	    njobs=`bjobs | wc -l`	   	    
	    	    while [ $njobs -ge $maxjobs ]; do
       			sleep 30
			njobs=`bjobs | wc -l`
	    	    done	   	    
	    	    ./gensplital.sh  runalign_$ii.sh
	    	    sleep 1
		else
		    chmod +x ./runalign_$ii.sh
		    ./runalign_$ii.sh
		fi
	    fi
	    ii=$(($ii+1))
	done
    fi
fi # inparallel or not
echo all jobs launched

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
