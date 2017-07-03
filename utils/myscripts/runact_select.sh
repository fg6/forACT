#!/bin/bash
set -o errexit

single=$1
starti=0
if [[ $single != 1 ]]; then 
    comp_to=$2; 
    starti=1
fi

ii=0
chrs=()
chrlist=""
for i; do  # loop over input args
    if [[ $ii > $starti ]]; then
	chrlist+="_"$i
	chrs+=($i) 
    fi
    let ii=$(($ii+1))
done


thisdir=`pwd`
source $thisdir/mysettings.sh
debug=0


if  [[ $chrs == "" ]]; then
    chrs=list
fi
if [[ $chrs == 'list' ]]; then
    $srcdir/listchrs/listchrs $refdir/$ref
    exit
fi

echo "  Selected chromosomes: " "${chrs[@]}"

create_chr () {    
    # global file (for the whole assembly:
    globals=$myfolder/whole_$shred/inter/nonoise$noise\_minid$minid\_$finalal
    globfasta=$myfolder/whole_$shred/inter/selctg_forw* #$forwnotshred  
    herefasta=nonoise$noise\_minid$minid\_$(basename $globfasta)
    hereals=nonoise$noise\_minid$minid\_tempal

    if [[ ! -f $thisfasta ]] ||  [[ ! -f $thisals ]]; then 
	rm -rf temp
	mkdir temp
	cd temp
	grep -P "\t"$chr"\t"   $globals  > tempal
    
	echo Reordering contigs and cutting noise ...         
	$srcdir/actnoise/actnoise $thisref $globfasta tempal $noise $minid
	echo Resetting positioning scheme according to ACT ...
	$srcdir/chrpos/chrpos $thisref $herefasta $hereals
	mv foract$hereals  $thisals
	mv $herefasta $thisfasta

    fi # files created 
}

create_multichr () {    
    # global file (for the whole assembly:
    globals=$myfolder/whole_$shred/inter/nonoise$noise\_minid$minid\_$finalal
    globfasta=$myfolder/whole_$shred/inter/selctg_forw* #$forwnotshred  
    herefasta=nonoise$noise\_minid$minid\_$(basename $globfasta)
    hereals=nonoise$noise\_minid$minid\_tempal
    thisref=$chr_folder/ref$chrlist.fasta   #$refdir/$chr\_*   

    if [[ ! -f $thisfasta ]] ||  [[ ! -f $thisals ]]; then 
	rm -rf temp
	mkdir temp
	cd temp
	rm -f tempal
	rm -f $thisref
	for chr in "${chrs[@]}"; do 
	    grep -P "\t"$chr"\t"   $globals  >> tempal
	    cat $refdir/$chr\_*  >> $thisref
	done

	echo Reordering contigs and cutting noise ...         
	$srcdir/actnoise/actnoise $thisref $globfasta tempal $noise $minid
	echo Resetting positioning scheme according to ACT ...
	$srcdir/chrpos/chrpos $thisref $herefasta $hereals
	mv foract$hereals  $thisals
	mv $herefasta $thisfasta

    fi # files created 
}


myfolder=`pwd`
chr_folder=$myfolder/whole_$shred/unique/chrs$noise/
mkdir -p $chr_folder
cd $chr_folder


# files to plot:
thisfasta=$myfolder/whole_$shred/unique/chrs$noise/for$chrlist\.fasta
thisals=$myfolder/whole_$shred/unique/chrs$noise/for$chrlist\.al
create_multichr

if [[ $single == 1 ]]; then
    /software/hpag/Artemis/act $thisref $thisals $thisfasta & 
    exit
fi

onefasta=$thisfasta
oneals=$thisals

# if comparison
myfolder=$comp_to
chr_folder=$myfolder/whole_$shred/unique/chrs$noise/
mkdir -p $chr_folder
cd $chr_folder

# files to plot:
thisfasta=$myfolder/whole_$shred/unique/chrs$noise/for$chrlist\.fasta
thisals=$myfolder/whole_$shred/unique/chrs$noise/for$chrlist\.al
create_multichr


echo Running ACT with files:
echo Reference: $thisref 
echo Top Draft Assembly: $onefasta
echo Top Alignment file: $oneals
echo Bottom Draft Assembly: $thisfasta
echo Bottom Alignment file: $thisals

/software/hpag/Artemis/act $onefasta $oneals $thisref $thisals $thisfasta & 


   
