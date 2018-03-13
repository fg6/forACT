#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

err=0


mkdir -p  $outdir
cd $outdir

jobs_first=`ls $aldir/split_$firstal/*.out | wc -l`
jobs_second=`ls $aldir/split_$secal/*.out | wc -l`


if (( $jobs_first < 1 )); then
    err=$(($err+1))
    echo Error: Too few jobs! $jobs_first
fi

if (( $jobs_first != $jobs_second )); then
    echo " " Error: Number of jobs are different in first and second alignment! $jobs_first  $jobs_second
    err=$(($err+1))
else
    echo " " Good: same number of jobs in first and second alignment $jobs_first  $jobs_second
fi

if [[ $lfsjobs == 1 ]]; then
	sum_jobals_first=`ls $aldir/split_$firstal/outs/out.* | wc -l | awk '{print $1}'`
	sum_jobals_second=`ls $aldir/split_$secal/outs/out.* | wc -l | awk '{print $1}'`
	if (( $sum_jobals_first != $sum_jobals_second )); then
    	echo " " Error: Number of output are different in first and second alignment! $sum_jobals_first $sum_jobals_second
    	err=$(($err+1))
	else
	    echo " " Good: same number of output in first and second alignment  $sum_jobals_first $sum_jobals_second
	fi

	success_first=`grep Success $aldir/split_$firstal/outs/out.* | wc -l | awk '{print $1}'`
	success_second=`grep Success $aldir/split_$secal/outs/out.* | wc -l | awk '{print $1}'`

	if (( $success_first != $jobs_first )); then
	    echo " " Error: Number of Successes in first alignment is wrong! $success_first instead of $jobs_first
	    err=$(($err+1))
	else
	    echo  " " Good: All jobs were Succesful in first alignment $success_first, $jobs_first
	fi

	if (( $success_second != $jobs_second )); then
	    echo " " Error: Number of Successes in second alignment is wrong! $success_second instead of $jobs_second
	    err=$(($err+1))
	else
	    echo  " " Good: All jobs were Succesful in second alignment  $success_second, $jobs_second
	fi
fi


als_sum_jobs_first=`wc -l $aldir/split_$firstal/*out | tail -1 | awk '{print $1}'`
als_sum_jobs_second=`wc -l $aldir/split_$secal/*out | tail -1 | awk '{print $1}'`

global_first=`wc -l $aldir/$firstal.al | awk '{print $1}' `
global_second=`wc -l $aldir/$secal.al| awk '{print $1}'`
global_third=`wc -l $workdir/$thirdal.al| awk '{print $1}'`


if (( $als_sum_jobs_first != $global_first )); then
    echo " " Error: Number of alignments in first alignment is wrong! sum_single_jobs=$als_sum_jobs_first global=$global_first
    err=$(($err+1))
else
    echo  " " Good:  Number of alignments in first alignment: sum_single_jobs=$als_sum_jobs_first global=$global_first
fi
if (( $als_sum_jobs_second != $global_second )); then
    echo " " Error: Number of alignments in second alignment is wrong! sum_single_jobs=$als_sum_jobs_second global=$global_second
    err=$(($err+1))
else
    echo  " " Good:  Number of alignments in second alignment: sum_single_jobs=$als_sum_jobs_second global=$global_second
fi

if (( $global_second != $global_first )); then
    if [[ $aligner == "smalt" ]]; then
	echo " " Error: Number of alignments in first alignment different from second: $global_second != $global_first
        err=$(($err+1))
    else
	echo "  Warning: Number" of alignments in first alignment different from second: $global_second != $global_first
	echo "           Ignore this warning if the difference is small and every other check is ok"
    fi
else
    echo  " " Good:  Same number of alignments in first and second alignments:  $global_second, $global_first
fi

if (( $global_second != $global_third )); then
    echo " " Error: Number of alignments in second alignment different from third: $global_second != $global_third
    err=$(($err+1))
else
    echo  " " Good:  Same number of alignments in third and second alignments:  $global_second, $global_third
fi


initial_bases=`$srcdir/n50/n50 $notshred | awk '{print $2}'`
forward_bases=`$srcdir/n50/n50 $fastadir/$forwnotshred | awk '{print $2}'`

if (( $initial_bases != $forward_bases )); then
    echo " " Error: Number of bases in intial and forward draft are different! $initial_bases, $forward_bases
    err=$(($err+1))
else
    echo  " " Good:  Same number of  bases in intial and forward draft  $initial_bases, $forward_bases
fi


echo; 
if (( $err > 0 )); then
    echo " Some errors occurred!"
else
    echo " Everything looks ok!"
fi
