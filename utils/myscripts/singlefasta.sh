#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

fasta=$1

faname=$(basename $fasta .fasta)
outdir=$fastadir/$faname\_fastas

if [ ! -d "$outdir" ]; then
    mkdir -p $outdir
    cd $outdir
    $srcdir/grabeachchr/grabeachchr $fasta
fi




