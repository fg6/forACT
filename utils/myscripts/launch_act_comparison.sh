# to run separately, independently from the pipeline
ref=$1  
folder1=$2
folder2=$3
shred=$4
noise=$5
minid=$6

if  [[ $shred == "" ]]; then
    shred=10000
fi
if  [[ $noise == "" ]]; then
    noise=30000
fi
if  [[ $minid == "" ]]; then
    noise=80
fi


#top
alfile1=$folder1/$aligner\_$shred/unique/foract$noise\_minid$minid.al  
afasta1=$folder1/$aligner\_$shred/unique/foract$noise\_minid$minid.fasta

#bottom
alfile2=$folder2/$aligner\_$shred/unique/foract$noise\_minid$minid.al 
afasta2=$folder2/$aligner\_$shred/unique/foract$noise\_minid$minid.fasta



if [ $# -lt 3 ]  || [ $1 == '-h' ] ; then
        echo; echo "  Usage:" $(basename $0) \<ref\> \<folder1\> \<folder2\>  [shred] [noise]
	echo "   folder1: full path to folder for first draft  => top assembly "
        echo "   folder2: full path to folder for first draft  => bottom assembly "
	echo "   shred: chunk size for alignments (in case there are more). Default=10Kb"
	echo "   noise: noise cut (in case there are more). Default=30Kb"
        exit 1
fi


echo Running ACT with files:
echo Reference: $ref 
echo Top Draft Assembly: $afasta1 
echo Top Alignment file: $alfile1
echo Bottom Draft Assembly: $afasta2 
echo Bottom Alignment file: $alfile2




/Artemis/act  $afasta1 $alfile1 $ref $alfile2 $afasta2 & 

