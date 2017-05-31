#!/bin/bash
set -o errexit

thisdir=`pwd`
myforACT=`echo $forACTdir`
myscripts=$forACTdir/utils/myscripts
#genscripts=$forACTdir/utils/genscripts

mkdir -p $myscripts

#shscripts=( setup.sh test.sh debug.sh pipeline.sh )
#for script in "${shscripts[@]}"; do 
#    genscript=$genscripts/$script
#    myscript=$myscripts/$script

#    if [[ ! -f $genscript ]]; then 
#	echo Error! Missing script: $genscript
#	echo Please try 'github pull'
#    else
#	sub="s#=XXXX#=$myforACT#g"
#	sed $sub $genscript > $myscript
#    fi
#done