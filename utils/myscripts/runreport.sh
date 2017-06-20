#!/bin/bash
set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh
alsfile=$workdir/nonoise$noise\_minid$minid\_$thirdal.al
debug=$1


mkdir -p  $outdir/report
cd $outdir/report

######## create a report of info ############
### variations? > inversions
### contigs that join more chr/contig  (misassemblies if one assembly is the reference)

#################################### Assembly reports ###################################
file1=$refdir/myn50.dat
refsize=`head -1 $file1 | awk '{print $2}' | awk '{ printf("%'"'"'d\n",$1); }' `
nrefsize=`head -1 $file1 | awk '{print $2}'`
refchrs=`head -1 $file1 | awk '{print $4}'| awk '{ printf("%'"'"'d\n",$1); }'`
refn50=`head -1 $file1 | awk '{print $10}'| awk '{ printf("%'"'"'d\n",$1); }'`
reflong=`head -1 $file1 | awk '{print $8}'| awk '{ printf("%'"'"'d\n",$1); }'`
file2=$fastadir/myn50.dat
if [[ ! -f $file2 ]]; then $srcdir/n50/n50 $notshred > $file2; fi
drsize=`head -1 $file2 | awk '{print $2}'| awk '{ printf("%'"'"'d\n",$1); }'`
drchrs=`head -1 $file2 | awk '{print $4}'| awk '{ printf("%'"'"'d\n",$1); }'`
drn50=`head -1 $file2 | awk '{print $10}'| awk '{ printf("%'"'"'d\n",$1); }'`
drlong=`head -1 $file2 | awk '{print $8}'| awk '{ printf("%'"'"'d\n",$1); }'`
echo; 
echo "*****************************************"
echo "*********** Assemblies report ***********"
echo "*****************************************"
echo " Reference info:" 
echo "   Size =" $refsize "bases"
echo "   Number of chrs =" $refchrs
echo "   N50 =" $refn50 "bases"
echo "   Longest chr =" $reflong "bases"
echo; echo " Draft assembly info:" 
echo "   Size =" $drsize "bases"
echo "   Number of ctgs =" $drchrs
echo "   N50 =" $drn50 "bases"
echo "   Longest chr =" $drlong "bases"
echo; 

############################## Mapping report ###############################
file3=$fastadir/myn50_shred.dat
if [[ ! -f $file3 ]]; then $srcdir/n50/n50 $fastadir/$shreddraft > $file3; fi
chunknum=`head -1 $file3 | awk '{print $4}'| awk '{ printf("%'"'"'d\n",$1); }'`
inimapped=`wc -l $aldir/$firstal.al  | awk '{print $1}'| awk '{ printf("%'"'"'d\n",$1); }'`
finalmapped=`wc -l $alsfile  | awk '{print $1}'| awk '{ printf("%'"'"'d\n",$1); }'`
i=`wc -l $aldir/$firstal.al  | awk '{print $1}'`
f=`wc -l $alsfile | awk '{print $1}'`
unmapped=`echo $(($i-$f))`
read avgid refcov <<< $(python $scriptdir/avgid.py $alsfile)
rcov=`echo $refcov | awk '{ printf("%'"'"'d\n",$1); }'`
ratio=`echo $refcov*1./$nrefsize`
percov=`awk "BEGIN {printf \"%.2f\", $refcov*100./$nrefsize}"`
echo "**********************************************************************"
echo "**** Mapping report for Draft ctgs shred in chunks of" $shred "bases ****"
echo "***********  isolated noise <" $noise " and min-id =" $minid "**************"
echo "**********************************************************************"
echo " Total number of mapped chunks:" $chunknum
echo " Initial number of chunks mapped:" $inimapped
echo " Unmapped chunks:" $unmapped 
echo " Final number of chunks mapped (after noise/quality selection):" $finalmapped
echo " Final global average identity:" $avgid"%"
# write avg identy per each contig in a file: how to evaluate if there are chunks?
echo " Reference coverage:" $percov\%   "  ("$rcov "bases)"



####################### Variations report ##########################
## writing misfinder.cpp code for this
echo; 
echo "*****************************************"
echo "*********** Variation report ***********"
echo "*****************************************"
$srcdir/misfinder/misfinder $fullpathref $workdir/nonoise$noise\_minid$minid\_selctg_$forwnotshred $alsfile 
#echo " Number of variations:" 
#echo " Number of inversions:"
#echo " Number of contigs that connects to more than 1 Reference contigs:"
echo; echo " Detailed ctg report summary in $outdir/report/ctg_report.txt"
echo " Detailed misassembly report summary in $outdir/report/misassembly_report.txt"


echo;echo
