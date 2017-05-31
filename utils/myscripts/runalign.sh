wdir=$1
# reference name:
myref=$2   
# initial draft fasta:
notshred=$3 
# shred size:
shred=$4

scriptdir=$5

# only test existing files:
test=0
outputname=$wdir\_$shred

### check if files exists:
ok=0
if [ ! -f $myref ]; then echo; echo "Could not find Reference in" $myref;  else ok=$(($ok+1)); fi
if [ ! -d $scriptdir ]; then echo; echo "Could not find script-dir in" $scriptdir; else ok=$(($ok+1)); fi
if [ ! -f $notshred ]; then echo; echo "Could not find initial fasta in" $notshred;  else ok=$(($ok+1)); fi
if [ ! $ok -eq 3 ]; then echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; exit; 
else echo Input files found!; fi




### Do not change
launchdir=`pwd`
dir=$launchdir
notshredname=$(basename $notshred)
draft=shred$shred\_$notshredname

outdir=$dir/$outputname
aldir=$outdir/aligns
fastadir=$outdir/fasta
workdir=$outdir/inter
refdir=$outdir/../ref
ref=$(basename $myref)

firstal=first  # shreded initial fasta
secal=second   # shreaded forward fasta
thirdal=third  #  shreaded forward fasta, fixed shred position in scaffold 

forwnotshred=forw$notshredname
singlefolder=$fastadir/$(basename $forwnotshred .fasta)_fastas
forwshred=$fastadir/shred$shred\_$forwnotshred 
forwal=$aldir/$secal.al


mkdir -p  $outdir
cd $outdir 



#######################################################
################# PREPARE REFERENCE  ##################
#######################################################
if [ ! -d $refdir ]; then
   mkdir -p $refdir
   cd $refdir
   cp $myref .
   ln -sf $myref ref.fasta  

 
   ~/ana/cpp/findchr/findchr $myref > chr.list 
   ~/ana/cpp/grabeachchr/grabeachchr $ref
fi
if [ ! -d $refdir ]; then
    echo; echo " Error! Reference temp location not found in"  $refdir;
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`ls $refdir | wc -l`
    #shouldbe=`$scriptdir/myn50.sh $myref | awk '{print $4}'`  # chromosomes number
    shouldbe=64  #24 # for seabass
    shouldbe=$(($shouldbe+3)) # chrs + ref + chr.list +ref.fasta

    if [ ! $check -eq $shouldbe ]; then 
	echo; echo " Error! too many or too few single fastas in" $refdir
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi
echo " Reference dir ready"


#######################################################
########## FIRST ALIGNMENT: SHRED AND ALIGN  ##########
#######################################################
### create shreded fasta ###
$scriptdir/createshred.sh $fastadir $notshred $shred $test 
# **** this created shred$shred\_$notshredname  = $draft ****
if [ ! -f $fastadir/shred$shred\_$notshredname ]; then 
    echo; echo " Error! Shred fasta not found in" $fastadir/shred$shred\_$notshredname
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`head -1 $fastadir/shred$shred\_$notshredname`
    if [ -z "${check}" ]; then 
	echo; echo " Error! Shred fasta is empty in" $fastadir/shred$shred\_$notshredname
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi

### split shreded fasta in fasta of 1M contigs ###
$scriptdir/splitshred.sh $fastadir/split shred$shred\_$notshredname $test 
if [ ! -f $fastadir/split/split0_shred$shred\_$notshredname ]; then 
    echo; echo " Error! Shred fasta not split in" $fastadir/split
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`head -1 $fastadir/split/split0_shred$shred\_$notshredname`
    if [ -z "${check}" ]; then 
	echo; echo " Error! Shred fasta split0 is empty in" $fastadir/split/split0_shred$shred\_$notshredname
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi

## align  fastas
$scriptdir/splitalign.sh $aldir $fastadir/split $refdir/$ref $firstal $scriptdir $shred $test 
if [ ! -f $aldir/$firstal.al ]; then
    for file in $aldir/split_first/*.out; do
	awk '{print $3"\t"$4"\t"$11"\t"$2"\t"20"\t"1"\t"$5"\t"$6"\t"$7"\t"$8"\t"0.0"\t"$12"\t"$9}' $file  >> $aldir/$firstal.al
    done
fi

# **** this created $aldir/smalt_hash.sma/smi and $aldir/$firstal.out and $aldir/$firstal.al ****
if [ ! -f $aldir/$firstal.al ]; then # || [ !  -f $aldir/$firstal.out ]; then
    echo; echo " Error! First alignment not found in" $aldir/$firstal.al  # or $aldir/$firstal.out
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`head -1 $aldir/$firstal.al`
    if [ -z "${check}" ]; then 
	echo; echo " Error!  First alignment empty in" $aldir/$firstal.al
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi

echo " First alignment done!"

#######################################################
### SECOND ALIGNMENT: ALL FORWARD, SHRED AND ALIGN  ###
#######################################################
### create fasta with all forward scaffold:
$scriptdir/allforw.sh $fastadir $notshred $aldir/$firstal.al $forwnotshred  $test 
# **** this created $fastadir/$forwnotshred ****
if [ ! -f $fastadir/$forwnotshred ]; then 
    echo; echo " Error! Forward fasta not found in" $fastadir/$forwnotshred
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`head -1 $fastadir/$forwnotshred`
    if [ -z "${check}" ]; then 
	echo; echo " Error!  Forward fasta empty in" $fastadir/$forwnotshred
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi


### create shreded fasta from all_forward_fasta
$scriptdir/createshred.sh $fastadir $fastadir/$forwnotshred $shred $test 
# **** this created $forwshred ****
if [ ! -f $forwshred ]; then 
    echo; echo " Error! Shred forward fasta not found in" $forwshred
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`head -1  $forwshred`
    if [ -z "${check}" ]; then 
	echo; echo " Error!   Shred Forward fasta empty in" $forwshred
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi


### split shreded fasta in fasta of 1M contigs ###
$scriptdir/splitshred.sh $fastadir/splitforw  $(basename $forwshred)  $test 
if [ ! -f $fastadir/splitforw/split0_$(basename $forwshred) ]; then 
    echo; echo " Error! Shred fasta not split in" $fastadir/splitforw
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`head -1 $fastadir/splitforw/split0_$(basename $forwshred)`
    if [ -z "${check}" ]; then 
	echo; echo " Error! Shred fasta split0 is empty in" $fastadir/splitforw/split0_$(basename $forwshred)
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi

## align  fastas
echo " Second alignments:"
$scriptdir/splitalign.sh $aldir $fastadir/splitforw $refdir/$ref $secal $scriptdir $shred $test 

if [ ! -f $forwal ]; then
    for file in $aldir/split_$secal/*.out; do
	awk '{print $3"\t"$4"\t"$11"\t"$2"\t"20"\t"1"\t"$5"\t"$6"\t"$7"\t"$8"\t"0.0"\t"$12}' $file >> $forwal
    done
fi

# this created $aldir/$secal.out and $aldir/$secal.al (= $forwal)
if [ ! -f $forwal ]; then
    echo; echo " Error! Second alignment not found in" $forwal 
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`head -1 $forwal`
    if [ -z "${check}" ]; then 
	echo; echo " Error! Second alignment empty in" $forwal 
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi

echo " Second alignment done!"

#######################################################
###### FIX AL POSITIONS and CREATE SINGLE FASTAS  #####
#######################################################
#### create single fasta from $fastadir/shred$shred\_$forwnotshred
$scriptdir/singlefasta.sh $fastadir $fastadir/$forwnotshred  $test 
# **** this created  $singlefolder
if [ ! -d $singlefolder ]; then
    echo; echo " Error! Single fasta folder not found in"  $singlefolder;
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`ls $singlefolder | wc -l`
    shouldbe=`$scriptdir/myn50.sh $fastadir/$forwnotshred | awk '{print $4}'`
    if [ ! $check -eq $shouldbe ]; then 
	echo; echo " Error! too many or too few single fastas in" $singlefolder 
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi

#### fix positions of shred ctg/scaffold:  (only if not using Zemin tabs)
$scriptdir/fixpos.sh $workdir $forwshred $forwal $thirdal $test 
# **** this created  $workdir/$thirdal.al
if [ ! -f $workdir/$thirdal.al ]; then
    echo; echo " Error! Third alignment not found in" $aldir/$third.al 
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    check=`head -1 $workdir/$thirdal.al`
    if [ -z "${check}" ]; then 
	echo; echo " Error! Third alignment empty in" $workdir/$thirdal.al
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi

echo All done!
