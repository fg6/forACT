#!/bin/bash
set -o errexit


thisdir=`pwd`

echo $forACTdir

shscripts=( setup.sh test.sh debug.sh pipeline.sh )

for script in "${shscripts[@]}"; do 
    thiscript=$thisdir/utils/$script
    
    if [[ -f $thiscript ]]; then 
	echo Error, missing script: $thiscript
    else
	ls $thisscript
    fi
done