#!/bin/bash
#set -o errexit

single=$1
chr=$2
comp_to=$3

thisdir=`pwd`
source $thisdir/mysettings.sh
debug=0


create_chr () {
    echo "hello"
    
}



exit

aref=$folder/../ref/ref.fasta
aalfile=$folder/unique/foract.al  
afasta=$folder/unique/foract.fasta


if [[ $chr == 'all' ]]; then
    ref=$aref
    alfile=$aalfile
    fasta=$afasta
elif [[ $chr == 'list' ]]; then
    ~/ana/cpp/foract_listchrs/foract_listchrs $aref
    exit
else
    ref=$folder/../ref/$chr\_*
    fasta=$folder/unique/chrs/for_$chr\.fasta
    alfile=$folder/unique/chrs/for_$chr\.al

  if [[ ! -f $fasta ]] ||  [[ ! -f $alfile ]]; then 
      mkdir -p $folder/unique/chrs
      cd $folder/unique/chrs
      aalfile=$folder/inter/nonoisethird.al  # from here otherwise positions already changed!

      rm -rf temp
      mkdir temp
      cd temp
      grep -P "\t"$chr"\t"   $aalfile  > tempal


      echo Reordering contigs and cutting noise ...
      ~/ana/cpp/actnoise/actnoise $ref $afasta tempal
      echo ~/ana/cpp/actnoise/actnoise $ref $afasta tempal
      echo Resetting positioning scheme according to ACT ...
      ~/ana/cpp/chrpos/chrpos $ref nonoiseforact.fasta  nonoisetempal
      mv foractnonoisetempal $alfile
      mv nonoiseforact.fasta $fasta
     
  fi
fi


echo; echo " " Running ACT with files:
echo " " Reference: $ref 
echo " " Draft Assembly: $fasta 
echo " " Alignment file: $alfile


/software/hpag/Artemis/act $ref $alfile $fasta & 

