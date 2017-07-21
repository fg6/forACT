#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

file=$1
location=$2
err=0

if [ ! -f $file ]; then 
    echo " Error! file " $file " not found in location " $location    
    err=1
else
    check=`head -1 $file`

    if [ -z "${check}" ]; then 
	echo; echo " Error! file " $file " appears to be empty in location "  $location 
	err=2
    fi
fi

echo $err
