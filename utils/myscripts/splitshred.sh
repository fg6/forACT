#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

dir=$1
shredfasta=$2

mkdir -p $dir
cd $dir

# split the shreded fasta in multiple fastas, each with 5000 contigs. 
# if the original shreaded fasta has less than 5000 contigs, a single fasta will be produced
if [ ! -f $dir/split0_$shred ]; then
    cd $dir
    $srcdir/splitinfastas/splitinfastas ../$shredfasta
fi

