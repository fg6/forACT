#!/bin/bash
set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh
debug=$1


mkdir -p  $outdir
cd $outdir

######## create a report of info ############
### myn50 info from $refdir/myn50.dat 
### variations? > inversions
### contigs that join more chr/contig  (misassemblies if one assembly is the reference)

## note: read the alignments after noise cut but before changing positions for act, so to report position in each contig or chr


###### Assemblies and alignment reports #####
file1=$refdir/myn50.dat
refsize=`head -1 $file1 | awk '{print $2}' | awk '{ printf("%'"'"'d\n",$1); }' `
refchrs=`head -1 $file1 | awk '{print $4}'| awk '{ printf("%'"'"'d\n",$1); }'`
refn50=`head -1 $file1 | awk '{print $10}'| awk '{ printf("%'"'"'d\n",$1); }'`
reflong=`head -1 $file1 | awk '{print $8}'| awk '{ printf("%'"'"'d\n",$1); }'`

file2=$fastadir/myn50.dat
if [[ ! -f $file2 ]]; then 
    $srcdir/n50/n50 $notshred > $file2
fi
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
echo "   N50 =" $refn50
echo "   Longest chr =" $reflong "bases"
echo; echo " Draft assembly info:" 
echo "   Size =" $drsize "bases"
echo "   Number of chrs =" $drchrs
echo "   N50 =" $drn50
echo "   Longest chr =" $drlong "bases"
echo; 
echo "**********************************************************************"
echo "**** Mapping report for Draft ctgs shred in chunks of" $shred "bases ****"
echo "**********************************************************************"
echo " Total number of chunks: "
echo " Initial number of chunks mapped:"
echo " Final number of chunks mapped (after noise/quality selection):"

echo; 
echo "*****************************************"
echo "*********** Variation report ***********"
echo "*****************************************"
echo " Number of variations:" 
echo " Number of inversions:"
echo " Number of contigs that connects to more than 1 Reference contigs:"
### contigs that join more chr/contig  (misassemblies if one assembly is the reference)


echo;echo