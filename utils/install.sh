#!/bin/bash
set -o errexit

thisdir=`pwd`


curl https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz > gzstream.tgz

#myscripts=$forACTdir/utils/myscripts
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
