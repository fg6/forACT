
ref=$1  
tfasta=$2
tals=$3
bfasta=$4
bals=$5



echo Running ACT with files:
echo Reference: $ref 
echo Top Draft Assembly: $tfasta 
echo Top Alignment file: $tals
echo Bottom Draft Assembly: $bfasta 
echo Bottom Alignment file: $bals

/software/hpag/Artemis/act  $tfasta $tals $ref $bals $bfasta & 

