
#dir=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/10Xg/Anopheles/act/Supernova_pseudohap/whole_10000

ref=$1  #$dir/ref/Anopheles_gambiae.AgamP4.dna.toplevel.fa
alfile=$2  #$dir/unique/foract.al
fasta=$3  #$dir/unique/foract.fasta



if [ $# -lt 3 ]  || [ $1 == '-h' ] ; then
        echo; echo "  Usage:" $(basename $0) \<ref\> \<alignments\> \<fasta\> 
        exit 1
fi


echo Running ACT with files:
echo Reference: $ref 
echo Draft Assembly: $fasta 
echo Alignment file: $alfile


/software/hpag/Artemis/act $ref $alfile $fasta & 
#~/software//Artemis/act  $ref $alfile $fasta &
