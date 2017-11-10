#!/bin/bash
##set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

mkdir -p  $outdir
cd $outdir 

#######################################################
################# PREPARE REFERENCE  ##################
#######################################################
if [ ! -d $refdir ]; then
   mkdir -p $refdir
   cd $refdir
   cp $fullpathref ref.fasta  
    
   $srcdir/grabeachchr/grabeachchr $ref
  
   if [[ $aligner == "smalt" ]]; then	
       if [ ! -f smalt_hash.sma ] || [ ! -f smalt_hash.smi ]; then
       	  thiscom=`echo $mysmalt  index -k 13 -s 6  smalt_hash $ref `
          if [[ $debug == 1 ]]; then
	     $thiscom 
	  
          else
	      $thiscom &> /dev/null
          fi
       fi
   fi
fi
if [ ! -d $refdir ]; then
    echo; echo " Error! Reference temp location not found in"  $refdir;
    echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
else
    if [[ ! -f $refdir/myn50.dat ]]; then
	$srcdir/n50/n50 $fullpathref > $refdir/myn50.dat
	chrnum=`head -1 $refdir/myn50.dat | awk '{print $4}'`  # chromosomes number
	refbases=`head -1 $refdir/myn50.dat | awk '{print $2}'`  #  number of bases
	echo chrnum=$chrnum > $refdir/refinfo.dat
	echo refbases=$refbases >> $refdir/refinfo.dat
    else
	source $refdir/refinfo.dat
    fi

    check=`ls $refdir | wc -l`
    
    shouldbe=$(($chrnum+3)) # chrs + ref.fasta + refinfo.dat + hash.smi + hash.sma + myn50.dat

    if [ ! $check -ge $shouldbe ]; then 
	echo; echo " Error! too many or too few single fastas in" $refdir $shouldbe $check
	echo; echo " ****  Something went wrong! Giving up! **** "; echo; exit
    fi
fi
echo " 1. Reference ready"
echo

######################################################
########## FIRST ALIGNMENT: SHRED AND ALIGN  ##########
#######################################################

### create shreded fasta ###
cd $thisdir
file=$fastadir/shred$shred\_$notshredname;  location="Two" 
if [ ! -f $file ]; then 
    $scriptdir/createshred.sh $notshred
fi
checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then  echo; echo "   " $checkfile; exit; fi
echo " 2. Draft assembly shreded"
echo

### split shreded fasta in fastas of 5K contigs ###
cd $thisdir
file=$fastadir/split/split0_shred$shred\_$notshredname;  location="Three" 
if [ ! -f $file  ]; then 
    if [[ $aligner == "smalt" ]]; then
	$scriptdir/splitshred.sh  $splitdir $fastadir/$shreddraft  #done
    else  #minimap2 and bwa
	$scriptdir/splitshred.sh  $splitdir $fastadir/$shreddraft  $splitlen
    fi
fi


checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then  echo; echo "   " $checkfile; exit; fi
echo " 3. Shreded draft assembly split in multiple fastas for faster alignment in parallel"
echo

## align  fastas
file=$aldir/$firstal.al; location="Four"
if [ ! -f  $aldir/split_$firstal/split0_$firstal.out  ]; then 
    if [[ $aligner == "smalt" ]]; then
	$scriptdir/splitalign.sh  $splitdir $firstal
    elif [[ $aligner == "minimap2" ]]; then
	$scriptdir/splitalign_minimap.sh  $splitdir $firstal $aldir/split_$firstal
    else
	$scriptdir/splitalign_bwa.sh  $splitdir $firstal $aldir/split_$firstal
    fi
fi


if [ ! -f $file ]; then
    for ofile in $aldir/split_$firstal/*.out; do
	if [[ $aligner == "smalt" ]]; then
	    awk '{print $3"\t"$4"\t"$11"\t"$2"\t"20"\t"1"\t"$5"\t"$6"\t"$7"\t"$8"\t"0.0"\t"$12"\t"$9}' $ofile  >> $aldir/$firstal.al
	else
	    cat $ofile | awk '{print $1"\t"$6"\t"$12"\t"0"\t"0"\t"0"\t"$3"\t"$4"\t"$8"\t"$9"\t"0"\t"$11"\t"$5}'  | sed 's#\t+#\tF#g' | sed 's#\t-#\tR#g' >> $aldir/$firstal.al 	
	    
	fi
    done
fi
checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then  echo; echo "   " $checkfile; exit; fi
echo " 4. First alignment done! "
echo


#######################################################
### SECOND ALIGNMENT: ALL FORWARD, SHRED AND ALIGN  ###
#######################################################
### create fasta with all forward scaffold:
file=$fastadir/$forwnotshred; location="Five"
if [ ! -f  $file ]; then
    $scriptdir/allforw.sh $aldir/$firstal.al $forwnotshred 
fi
checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then  echo; echo "   " $checkfile; exit; fi
echo " 5. All contigs/scaffolds forward now "
echo


### create shreded fasta from all_forward_fasta
file=$forwshred; location="Six"
if [ ! -f $file ]; then 
    $scriptdir/createshred.sh  $forwnotshred 
fi
checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then echo; echo "   " $checkfile; exit; fi
echo " 6. Forward shreded draft assembly split in multiple fastas for faster alignment in parallel "
echo

### split shreded fasta in fasta of 1M contigs ###
file=$fastadir/splitforw/split0_$(basename $forwshred); location="Seven"
if [ ! -f  $file ]; then 
    if [[ $aligner == "smalt" ]]; then
	$scriptdir/splitshred.sh $fastadir/splitforw  $forwshred 
	#$scriptdir/splitshred.sh $fastadir/splitforw  $(basename $forwshred)  
    else
	$scriptdir/splitshred.sh  $fastadir/splitforw  $forwshred $splitlen
    fi
fi
checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then echo; echo "   " $checkfile; exit; fi
echo " 7. Forward draft shreded "
echo
## align  fastas
file=$forwal; location="Eigth"
if [ ! -f  $aldir/split_$secal/split0_$secal.out  ]; then 
    if [[ $aligner == "smalt" ]]; then
	$scriptdir/splitalign.sh $fastadir/splitforw $secal 
    else
    	$scriptdir/splitalign_minimap.sh  $fastadir/splitforw $secal  #c not needed
    fi
fi

if [ ! -f $file ]; then
    for ofile in $aldir/split_$secal/*.out; do
	if [[ $aligner == "smalt" ]]; then
	    cat $ofile | sed 's#:# #g' | awk '{print $5"\t"$6"\t"$13"\t"$4"\t"22"\t"1"\t"$7"\t"$8"\t"$9"\t"$10"\t"0.0"\t"$14"\t"$3}' >> $forwal  ## need to check this
	    #awk '{print $3"\t"$4"\t"$11"\t"$2"\t"20"\t"1"\t"$5"\t"$6"\t"$7"\t"$8"\t"0.0"\t"$12}' $ofile >> $forwal
	else
	    cat $ofile | awk '{print $1"\t"$6"\t"$10*100/$11"\t"0"\t"0"\t"0"\t"$3"\t"$4"\t"$8"\t"$9"\t"0"\t"$11"\t"$12}' >>  $forwal
	fi
    done
fi

checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then echo; echo "   " $checkfile; exit; fi
echo " 8. Second alignment done!"
echo
#######################################################
##################  FIX AL POSITIONS  #################
#######################################################


#######################################################
################ CREATE SINGLE FASTAS  ################
#######################################################
cd $thisdir
if [ ! -d $singlefolder ]; then
    $scriptdir/singlefasta.sh $fastadir/$forwnotshred  
fi
if [ ! -d $singlefolder ]; then
    echo; echo " Error! Single fasta folder not found in"  $singlefolder; exit
else
    check=`ls $singlefolder | wc -l`
    shouldbe=`$srcdir/n50/n50 $fastadir/$forwnotshred | awk '{print $4}'`
    if [ ! $check -eq $shouldbe ]; then 
	echo; echo " Error! too many or too few single fastas in" $singlefolder; exit
    fi
    # contig sizes file
    check=`wc -l $contigsizes | awk '{print $1}'`
    if [ ! $check -eq $shouldbe ]; then 
	echo; echo " Error! too many or too few contigs in" $contigsizes; exit
    fi
fi
echo " 9. Prepared draft contigs"
echo


#### fix positions of shred ctg/scaffold:  (only if not using Zemin tabs) [ it changes ctg names from ctg1_1, ctg1_100000, ... in ctg1,...]
file=$workdir/$thirdal.al; location="Ten"
if [ ! -f $workdir/$thirdal.al ]; then
    $scriptdir/fixpos.sh $forwshred $forwal $thirdal 
fi
checkfile=`$scriptdir/checkfile.sh $file $location`
err=`echo $checkfile | tail -1`
if [[ $err > 0 ]]; then echo; echo "   " $checkfile; exit; fi
echo " 10. Prepared draft contigs order and positions"
echo
echo "   "All done! You can now proceed with \" $ ./mypipeline.sh prepfiles \"
echo
