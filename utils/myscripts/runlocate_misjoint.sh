#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh
alsfile=$workdir/$name_fornoise\_minid$minid\_$thirdal.al

mkdir -p  $outdir/report
cd $outdir/report

####################### Report possible misjoint for chromosome-assigned scaffolds (Synteny Groups)  ##########################
echo; 
echo "*****************************************************************************************************************"
echo "************** Report possible misjoints for chromosome-assigned scaffolds (Synteny Groups)  ********************" 
echo "*****************************************************************************************************************"
$srcdir/locate_misjoints/locate_misjoints $workdir/$name_fornoise\_minid$minid\_selctg_$forwnotshred $alsfile $refdir/ref.fasta $splitlen $min_len_max $min_len

if [[ -f $outdir/report/alfile.txt ]] &&  [[ -s $outdir/report/alfile.txt ]] ; then
	mv $outdir/report/alfile.txt  $misals_file
	mv $outdir/report/majorfiles.txt   $chrass_file
	mv $outdir/report/bedfile.txt $bedfile
	echo; echo " Misjoint report in" $outdir/report/misjoints_details_$name_fornoise\_minid$minid\_$min_len.txt	
else
	rm -f $outdir/report/alfile.txt $outdir/report/majorfiles.txt $outdir/report/bedfile.txt
	echo; echo " No misjoints found with the required parameters "
fi
echo;echo
