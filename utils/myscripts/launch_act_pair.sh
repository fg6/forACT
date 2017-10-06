
ref=$1  
tfasta=$2
tals=$3
bfasta=$4
bals=$5

source $thisdir/mysettings.sh


echo Running ACT with files:
echo Reference: $ref 
echo Top Draft Assembly: $tfasta 
echo Top Alignment file: $tals
echo Bottom Draft Assembly: $bfasta 
echo Bottom Alignment file: $bals


$myforACT/utils/mysrcs/Artemis/act  $tfasta $tals $ref $bals $bfasta & 

