#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

fasta=$1

faname=$(basename $fasta .fasta)


if [ ! -d "$singlefolder" ]; then
    mkdir -p $singlefolder
    cd $singlefolder
    $srcdir/grabeachchr/grabeachchr $fasta
fi

## write file with contig sizes
rm -f $contigsizes
for file in $singlefolder/*; do 
#	echo $file; 
	contig=`head -1 $file | sed 's#>##g'`; 
	bases=`$srcdir/n50/n50 $file | awk '{print $2}'`; 
#	echo $contig $bases;  
	echo $contig $bases >> $contigsizes
done




