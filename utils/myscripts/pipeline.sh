#!/bin/bash
set -o errexit

whattodo=$1
debug=0
if [ $# -lt 1 ] || [ $1 == '-h' ]; then
    echo; echo "  Usage:" $(basename $0) \<command\> 
    echo "     command: command to be run. Options: align, prepfiles"
    echo "     align: shred draft assemblies and align against Reference. A draft contig is re-oriented "
    echo "                if most of the shreded pieces are complements wrt the Reference"
    echo "     prepfiles: the alignment files and the fasta files are prepared to be compatible with the format required by ACT. "
    echo "                Draft contigs are ordered according to the position of their major alignment."

    exit
fi


thisdir=`pwd`
source $thisdir/mysettings.sh
mkdir -p $dir
cd $dir


if [ $whattodo == "align" ]; then
    #######################################################
    ###############   ALIGN PIPELINE    ##################
    #######################################################
    cd $dir
    ok=0
    if [ ! -f $fullpathref ]; then echo; echo "Could not find Reference in" $fullpathref;  else ok=$(($ok+1)); fi
    if [ ! -d $scriptdir ]; then echo; echo "Could not find script-dir in" $scriptdir; else ok=$(($ok+1)); fi
    if [ ! -f $notshred ]; then echo; echo "Could not find draft assembly fasta in" $notshred;  else ok=$(($ok+1)); fi
    if [ ! $ok -eq 3 ]; then echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; exit; 
    else echo; echo " All input files found! Proceeding with pipeline.."; fi

    $runalign $debug 2>&1 #| tee -a $oalign
fi


if [ $whattodo == "prepfiles" ]; then
    #######################################################
    ##################   PREPARE FILES  ###################
    #######################################################
    cd $dir
    ok=0
    if [ ! -f $refdir/$ref ]; then echo; echo "Could not find Reference in" $refdir/$ref;  else ok=$(($ok+1)); fi
    if [ ! -d $scriptdir ]; then echo; echo "Could not find script-dir in" $scriptdir; else ok=$(($ok+1)); fi
    if [ ! -f $fastadir/$forwnotshred ]; then echo; echo "Could not find draft assembly fasta in" $fastadir/$forwnotshred;else ok=$(($ok+1)); fi
    if [ ! -f $workdir/$thirdal.al ]; then echo; echo "Could not find final alignment in" $workdir/$thirdal; else ok=$(($ok+1)); fi
    if [ ! $ok -eq 4 ]; then echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; exit; 
    else echo; echo " All needed files found! Proceeding with pipeline.."; fi
    

    $runprep 
#$folder $ref $notshred $scriptdir 2>&1 #| tee -a $oprep
fi