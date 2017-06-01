#!/bin/bash
set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

fasta=$1
alfile=$2
outname=$3
if [ ! -f "$dir/$outname.al" ]; then
    mkdir -p $workdir
    cd $workdir
    $srcdir/samectgpos/samectgpos $fasta $alfile
    mv ctgpos_$(basename $alfile) $outname.al
fi
