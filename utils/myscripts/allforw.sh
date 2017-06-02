#!/bin/bash
set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

alfile=$1
outname=$2

mkdir -p $fastadir
cd $fastadir

if [ ! -f "$dir/$outname" ]; then
    #echo Reverting complement scaffolds ...
    $srcdir/revertcompl/revertcompl $notshred $alfile
fi

