#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

fastasplit=$1
outname=$2
alsplitdir=$aldir/split_$outname


mkdir -p $aldir
cd $aldir

mkdir -p $alsplitdir 
cd $alsplitdir
ii=0
temponame=split$ii\_$outname
subqueue="s#MYQUEUE#$myqueue#g"
submem="s#MYJOBMEM#$myjobmem#g"
subcpu="s#MYCPUS#$myncpus#g"



if [ ! -f $alsplitdir/$temponame.out ]; then	    
    for fa in $fastasplit/*; do
	temponame=split$ii\_$outname
	
	if [ ! -f $temponame.out ]; then	    
	
	    command=`echo "$mysmalt  map -f ssaha -m 100 -n 15 -O " $refdir/smalt_hash $fa "> $temponame.out"`
	    echo $command > runalign_$ii.sh

	
	    if [[ $lfsjobs == 1 ]]; then
	        sed $subqueue $scriptdir/gensplital.sh | sed $submem | sed $subcpu > gensplital.sh
        	chmod +x gensplital.sh
	    	njobs=`bjobs | wc -l`	   	    
	    	while [ $njobs -ge $maxjobs ]; do
       			sleep 60
			njobs=`bjobs | wc -l`
	    	done	   	    
	    	./gensplital.sh  runalign_$ii.sh
	    	#sleep 10
            else
		chmod +x ./runalign_$ii.sh
		./runalign_$ii.sh
	    fi
	fi
	ii=$(($ii+1))
     done
fi

echo all done

if [[ $lfsjobs == 1 ]]; then

	njobs=`bjobs | grep "split" | wc -l`
	if [[ $njobs > 0 ]]; then
	   while [ $njobs -ge 1 ]; do
		sleep 60
		njobs=`bjobs | grep "split" | wc -l`
    	    done
    	    sleep 10
	fi
	sleep 600

	t1=`ls $alsplitdir/runalign_* | wc -l`
	t2=`grep Success $alsplitdir/outs/out.* | wc -l`
	if [[ $t1 != $t2 ]]; then
	    echo " Some alignment jobs possibly failed? Check in $alsplitdir/outs/"
	fi
fi
