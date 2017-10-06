thisref=$1  
alfile=$2  
fasta=$3  

thisdir=`pwd`

if [ $# -lt 3 ]  || [ $1 == '-h' ] ; then
        echo; echo "  Usage:" $(basename $0) \<ref\> \<alignments\> \<fasta\> 
        exit 1
fi


echo Running ACT with files:
echo Reference: $thisref 
echo Draft Assembly: $fasta 
echo Alignment file: $alfile


source $thisdir/mysettings.sh
$myforACT/utils/mysrcs/Artemis/act $thisref $alfile $fasta & 

