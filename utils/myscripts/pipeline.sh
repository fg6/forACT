#!/bin/bash
set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh


whattodo=$1
debug=0
if [ $# -lt 1 ] || [ $1 == '-h' ]; then
    echo; echo "  Usage:" $(basename $0) \<command\> 
    echo "     command: command to be run. Options: align, prepfiles, check, report, act, act_compare, act_select, act_select_compare"
    echo "      * align: shred draft assemblies and align against Reference. A draft contig is re-oriented "
    echo "                if most of the shreded pieces are complements wrt the Reference"
    echo "      * prepfiles: the alignment files and the fasta files are prepared to be compatible with the format required by ACT. "
    echo "                Draft contigs are ordered according to the position of their major alignment."
    echo "      * check: check if pipeline ran smoothly"
    echo "      * report: write a report for the draft, reference assemblies and their mapping" 
    echo "      * act:  launch act for the latest forACT launched"
    echo "      * act_compare folder_to_compare_to:  launch act for the latest forACT launched and compare with another forACT (needs additional input the ull path to the forACT-folder to compare to)"
    echo "      * act_select chr1...chr5:  launch act for the latest forACT launched for chromosomes/ctgs chr1 to chr5 (up to 5 chromosomes)"
    echo "      * act_select_compare folder_to_compare_to chr1..chr5:  launch act for chromosomes/ctgs chr1 to chr5 for the latest forACT launched compared with another forACT (needs additional input the full path to the forACT-folder to compare to), up to 5 chromosomes"

    echo; echo  "  Check" https://github.com/fg6/forACT/blob/master/README.md "for detailed instructions"; echo
    exit
fi


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

if [ $whattodo == "report" ]; then
    #######################################################
    ##################  CREATE REPORT  ###################
    #######################################################
    cd $dir
    ok=0
    
    $runreport | tee report_noise$noise\_minid$minid.txt	


fi

if [ $whattodo == "debug" ]; then
  ###################################################
  echo; echo " Looking for possible issues... "
  ###################################################
  #$thisdir/utils/myscripts/debug.sh $myforACT
fi



if [ $whattodo == "check" ]; then
  ###################################################
  echo; echo " Looking for possible issues... "
  ###################################################
  $myforACT/utils/myscripts/runcheck.sh 
fi

if [ $whattodo == "act" ]; then
  ###################################################
  echo; echo " Launching act... "
  ###################################################
  $myforACT/utils/myscripts/runact.sh 1
fi

if [ $whattodo == "act_compare" ]; then
  ###################################################
  echo; echo " Launching act... "
  ###################################################
  if [ $# -lt 2 ]; then
        echo "  Error: Please provide the full path to the forACT-folder for the draft you want to compare to "
        exit
  fi
  echo " Please notice that this will work only if the comparison is done between two forACTs with same noise cut "

  $myforACT/utils/myscripts/runact.sh 2 $2
fi

if [ $whattodo == "act_select" ]; then
  ###################################################
  ###################################################
  $myforACT/utils/myscripts/runact_select.sh 1 $2
fi

if [ $whattodo == "act_select_compare" ]; then
  ###################################################
  ###################################################
  $myforACT/utils/myscripts/runact_select.sh 2 $2 $3 $4 $5 $6 $7
fi




