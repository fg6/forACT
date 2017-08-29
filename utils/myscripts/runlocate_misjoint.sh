#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh
alsfile=$workdir/$name_fornoise\_minid$minid\_$thirdal.al
# old case alsfile=$workdir/nonoise$noise\_$thirdal.al

mkdir -p  $outdir/report
cd $outdir/report

####################### Report possible misjoint for chromosome-assigned scaffolds (Synteny Groups)  ##########################
## writing misfinder.cpp code for this
echo; 
echo "*****************************************************************************************************************"
echo "************** Report possible misjoints for chromosome-assegned scaffolds (Synteny Groups)  ********************" 
echo "*****************************************************************************************************************"
$srcdir/locate_misjoints/locate_misjoints $workdir/$name_fornoise\_minid$minid\_selctg_$forwnotshred $alsfile  $min_len_max $min_len_perc $min_len
mv $outdir/report/misjoints_details.txt $outdir/report/misjoints_details_$name_fornoise\_minid$minid.txt

echo; echo " Misjoint report summary in" $outdir/report/misjoints_details_$name_fornoise\_minid$minid.txt
echo;echo
