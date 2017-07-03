#!/bin/bash
set -o errexit

single=$1
comp_to=$2

thisdir=`pwd`
source $thisdir/mysettings.sh
debug=0



if [[ $single == 1 ]]; then
	$scriptdir/launch_act.sh $refdir/$ref  $foractal  $foractfa &
	

else

	if [ $# -lt 2 ]; then
		echo; echo "  Please provide the forACT-folder for the draft you want to compare to "
		exit
	fi

	#top
	alfile1=$dir/whole_$shred/unique/foract$noise\_minid$minid.al
	afasta1=$dir/whole_$shred/unique/foract$noise\_minid$minid.fasta

	#bottom
	alfile2=$comp_to/whole_$shred/unique/foract$noise\_minid$minid.al
	afasta2=$comp_to/whole_$shred/unique/foract$noise\_minid$minid.fasta


	if [[ ! -f $alfile1 ]] || [[ ! -f $afasta1 ]]; then
		echo " Something is wrong: cannot find the final results for this forACT: $alfile1 or $afasta1 are missing " ; exit
	elif  [[ ! -f $alfile2 ]] || [[ ! -f $afasta2 ]]; then
		echo " Something is wrong: cannot find the final results for the forACT you want to compare to: $alfile2 or $afasta2 are missing " ; exit
	fi

	$scriptdir/launch_act_pair.sh $refdir/$ref $afasta1 $alfile1 $afasta2 $alfile2  
fi

