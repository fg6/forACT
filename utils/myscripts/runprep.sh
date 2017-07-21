#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh



file=$workdir/selctg_$forwnotshred;  location="One" 
if [ ! -f $file ]; then 
    cd $workdir 
    #echo Creating fasta including only  aligned ctgs
    python $scriptdir/selctg.py  $workdir/$finalal | awk '{print $1}' > ctglist.al  
    $srcdir/writeselctg/writeselctg $fastadir/$forwnotshred ctglist.al
fi
cd $dir
checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then  echo; echo "   " $checkfile; exit; fi
echo " 1. Draft assembly including only mapped contigs ready"

file1=$workdir/$name_fornoise\_minid$minid\_selctg_$forwnotshred;  location1="Two" 
file2=$workdir/$name_fornoise\_minid$minid\_$finalal;  location2="Three" 

if [ ! -f $file1 ] ||  [ ! -f $file2 ]; then 
    cd $workdir 
    $srcdir/actnoise/actnoise $refdir/$ref $workdir/selctg_$forwnotshred $workdir/$finalal $noise $minid
fi

cd $dir
checkfile=`$scriptdir/checkfile.sh $file1 $location1`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then  echo; echo "   " $checkfile; exit; fi
checkfile=`$scriptdir/checkfile.sh $file2 $location2`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then  echo; echo "   " $checkfile; exit; fi

echo " 2. Contigs ordered by mapping appearance and noise alignments cut"


file=$workdir/foract$name_fornoise\_minid$minid\_$finalal;  location1="Four" 
if  [ ! -f $file ]; then
    #echo Resetting positioning scheme according to ACT ...
    cd $workdir 
    $srcdir/chrpos/chrpos $refdir/$ref $workdir/$name_fornoise\_minid$minid\_selctg_$forwnotshred $workdir/$name_fornoise\_minid$minid\_$finalal
fi
echo " 3. Alignment positions reset according to ACT scheme"


thisfasta=$workdir/$name_fornoise\_minid$minid\_selctg_$forwnotshred
thisal=$workdir/foract$name_fornoise\_minid$minid\_$finalal


####################################
########## Final  ##################
#rm -rf $finaldir/
mkdir -p $finaldir

file1=$foractfa;  location1="Five" 
file2=$foractal;  location2="Six" 
if [ ! -f $file1 ] || [ ! -f $file2 ]; then
    cp $thisal $foractal
    cp $thisfasta $foractfa
fi


if [ ! -f $file1 ] || [ ! -f $file2 ]; then
    echo; echo " Error! final files  not found " in  $finaldir; exit
else
    echo; echo All done! 
    echo; echo  Launch act this way:
    echo; echo  ./mypipeline.sh act 
    echo " or to select a single chromosome to visualize:"
    echo " ./mypipeline.sh act_select chr"
fi


