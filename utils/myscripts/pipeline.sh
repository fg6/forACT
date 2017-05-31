#!/bin/bash
set -o errexit

whattodo=$1
if [ $# -lt 1 ] || [ $1 == '-h' ]; then
    echo; echo "  Usage:" $(basename $0) \<command\> 
    echo "     command: command to be run. Options: align, prepfiles"
    echo "     align: shred draft assemblies and align against Reference. A draft contig is re-oriented "
    echo "                if most of the shreded pieces are complements wrt the Reference"
    echo "     prepfiles: the alignment files and the fasta files are prepared to be compatible with the format required by ACT. "
    echo "                Draft contigs are ordered according to the position of their major alignment."

    exit
fi

# shred in:
shred=10000


thisdir=`pwd`
myforACT=MYFORACT
myref=MYREF
notshred=MYDRAFT
dir=MYDESTDIR
scriptdir=$myforACT/utils/myscripts


runalign=$scriptdir/splitrun.sh
runchrs=$scriptdir/runchrs.sh
wdir=whole

# output dir:
workdir=$wdir\_$shred
refdir=$dir/$workdir/ref
ref=$(basename $myref) 

mkdir -p $dir
cd $dir

errfile=err.runpipe
errchr=errchr.runpipe


#######################################################
###############   GLOBAL PIPELINE    ##################
#######################################################
$runpipeline $wdir $myref $notshred $shred $scriptdir 2>&1 | tee -a $errfile
check=`grep Error $errfile | wc -l`
if [ $check -gt 0 ]; then  
    echo; echo " Runpipeline exited with errors!"; echo;
    cat $errfile
    echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
    exit; 
fi

exit

#######################################################
##################   PREPARE FILES  ###################
#######################################################
if [ ! -f $dir/$workdir/inter/third.al ]; then
    echo final alignment file is missing! $dir/$workdir/inter/third.al
    echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
    exit
fi
rm $errchr
$runchrs $workdir $(basename $myref) $notshred $scriptdir 2>&1 | tee -a $errchr

check=`grep Error $errchr | wc -l`
if [ $check -gt 0 ]; then  
        echo; echo " Runpipeline second step exited with errors! check output in" $errchr; echo " here:";
        cat $errchr
        echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
        exit; 
fi
