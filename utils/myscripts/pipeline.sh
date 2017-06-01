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
    ok=0
    if [ ! -f $fullpathref ]; then echo; echo "Could not find Reference in" $fullpathref;  else ok=$(($ok+1)); fi
    if [ ! -d $scriptdir ]; then echo; echo "Could not find script-dir in" $scriptdir; else ok=$(($ok+1)); fi
    if [ ! -f $notshred ]; then echo; echo "Could not find draft assembly fasta in" $notshred;  else ok=$(($ok+1)); fi
    if [ ! $ok -eq 3 ]; then echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; exit; 
    else echo; echo " All input files found! Proceeding with pipeline.."; fi

    $runalign $debug 2>&1 | tee -a $oalign
# $folder $fullpathref $notshred $shred $scriptdir 2>&1 #| tee -a $oalign

exit
    # check for errors:
    check=`grep Error $oalign | wc -l`
    if [ $check -gt 0 ]; then  
	echo; echo " Runpipeline exited with errors!"; echo;
	cat $oalign
	echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
	exit; 
    fi
    
fi
exit

#######################################################
##################   PREPARE FILES  ###################
#######################################################
if [ ! -f $dir/$folder/inter/third.al ]; then
    echo final alignment file is missing! $dir/$folder/inter/third.al
    echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
    exit
fi
rm $errchr
$runchrs $folder $ref $notshred $scriptdir 2>&1 | tee -a $oprep

check=`grep Error $oprep | wc -l`
if [ $check -gt 0 ]; then  
        echo; echo " Runpipeline second step exited with errors! check output in" $oprep; echo " here:";
        cat $oprep
        echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
        exit; 
fi
