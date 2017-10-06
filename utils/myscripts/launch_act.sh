ref=$1  
alfile=$2  
fasta=$3  



if [ $# -lt 3 ]  || [ $1 == '-h' ] ; then
        echo; echo "  Usage:" $(basename $0) \<ref\> \<alignments\> \<fasta\> 
        exit 1
fi


echo Running ACT with files:
echo Reference: $ref 
echo Draft Assembly: $fasta 
echo Alignment file: $alfile


source $thisdir/mysettings.sh
$myforACT/utils/mysrcs/Artemis/act $ref $alfile $fasta & 

