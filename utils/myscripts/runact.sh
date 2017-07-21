#!/bin/bash
#set -o errexit

single=$1
comp_to=$2

thisdir=`pwd`
source $thisdir/mysettings.sh
debug=0


# check first if assembly sizes are not too big:
file=$refdir/myn50.dat
refsize=`head -1 $file | awk '{print $2}' `
file=$fastadir/myn50.dat
if [[ ! -f $file ]]; then $srcdir/n50/n50 $notshred > $file; fi
draftsize=`head -1 $file | awk '{print $2}' `

maxsize=2000000000  # ACT cannot handle assembly with size > 2.1 Gb
if [[ $refsize -gt $maxsize ]] || [[ $draftsize -gt $maxsize ]]; then
#	echo $refsize, $draftsize, $maxsize
	echo "  Error: your assembly is too large for ACT to handle. "
	echo "         Please launch ACT on single chromosomes, or up to 5 chromosomes using:"
	echo "         ./mypipeline.sh act_select chr1_name chr2_name ...."
	echo
	exit
fi

if [[ $single == 1 ]]; then
	$scriptdir/launch_act.sh $refdir/$ref  $foractal  $foractfa &
else

	if [ $# -lt 2 ]; then
		echo; echo "  Please provide the forACT-folder for the draft you want to compare to "
		exit
	fi

	#top
	alfile1=$dir/whole_$shred/unique/foract$name_fornoise\_minid$minid.al
	afasta1=$dir/whole_$shred/unique/foract$name_fornoise\_minid$minid.fasta

	#bottom
	alfile2=$comp_to/whole_$shred/unique/foract$name_fornoise\_minid$minid.al
	afasta2=$comp_to/whole_$shred/unique/foract$name_fornoise\_minid$minid.fasta


	if [[ ! -f $alfile1 ]] || [[ ! -f $afasta1 ]]; then
		echo " Something is wrong: cannot find the final results for this forACT: $alfile1 or $afasta1 are missing " ; exit
	elif  [[ ! -f $alfile2 ]] || [[ ! -f $afasta2 ]]; then
		echo " Something is wrong: cannot find the final results for the forACT you want to compare to: $alfile2 or $afasta2 are missing " ; exit
	fi

	$scriptdir/launch_act_pair.sh $refdir/$ref $afasta1 $alfile1 $afasta2 $alfile2  
fi

