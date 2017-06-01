#!/bin/bash
set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

fasta=$1
draft=shred$shred\_$fasta

mkdir -p $fastadir
cd $fastadir

if [ ! -f $fastadir/$draft ]; then
    $srcdir/splitreads/splitreads $fasta $shred
fi

