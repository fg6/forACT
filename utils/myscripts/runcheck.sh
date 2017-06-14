#!/bin/bash
set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh
debug=$1


mkdir -p  $outdir
cd $outdir




# check number of jobs, number of success for first and second, can we find why jobs crashed?
# check number of alignments in split_first, first.al and same for second