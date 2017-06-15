#!/bin/bash
set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh
debug=$1


mkdir -p  $outdir
cd $outdir

jobs_first=`ls $aldir/split_$firstal/*.out | wc -l`
jobs_second=`ls $aldir/split_$secal/*.out | wc -l`

echo " number of jobs:"  $jobs_first  $jobs_second

sum_jobals_first=`ls $aldir/split_$firstal/outs/out.* | wc -l | awk '{print $1}'`
sum_jobals_second=`ls $aldir/split_$secal/outs/out.* | wc -l | awk '{print $1}'`

echo " number of out. jobs:" $sum_jobals_first $sum_jobals_second

success_first=`grep Success $aldir/split_$firstal/outs/out.* | wc -l | awk '{print $1}'`
success_second=`grep Success $aldir/split_$secal/outs/out.* | wc -l | awk '{print $1}'`

echo " number of Successful jobs:" $success_first $success_second


als_sum_jobs_first=`wc -l $aldir/split_$firstal/*out | tail -1 | awk '{print $1}'`
als_sum_jobs_second=`wc -l $aldir/split_$secal/*out | tail -1 | awk '{print $1}'`

echo " jobs-sum number of alignments:" $als_sum_jobs_first $als_sum_jobs_second

global_first=`wc -l $aldir/$firstal.al | awk '{print $1}' `
global_second=`wc -l $aldir/$secal.al| awk '{print $1}'`
global_third=`wc -l $workdir/$thirdal.al| awk '{print $1}'`


echo " global number of alignments:" $global_first $global_second $global_third


# check number of jobs, number of success for first and second, can we find why jobs crashed?
# check number of alignments in split_first, first.al and same for second