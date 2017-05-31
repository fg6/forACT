set -o nounset
set -o errexit

scriptdir=~/ana/scripts/foract/
# ref fasta:
myref=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/MyRefs/Anopheles/RefGambiae/Anopheles_gambiae.AgamP4.dna.toplevel.fa
# initial draft fasta:
notshred=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/10Xg/Anopheles/act/Supernova_pseudohap/fasta/250M_pseudohap.fasta
# shred in:
shred=10000
# write in dir:
dir=/lustre/scratch117/sciops/team117/hpag/fg6/analysis/10Xg/Anopheles/act/Supernova_pseudohap



### do no change ###
runpipeline=$scriptdir/splitrun.sh
runchrs=$scriptdir/runchrs.sh
merge=$scriptdir/mergechrs.sh
#output dir label:
wdir=whole

errfile=err.runpipe
errchr=errchr.runpipe
showout=1

# output dir:
workdir=$wdir\_$shred
refdir=$dir/$workdir/ref
ref=$(basename $myref) 

mkdir -p $dir
cd $dir



#######################################################
###############   GLOBAL PIPELINE    ##################
#######################################################
rm $errfile
$runpipeline $wdir $myref $notshred $shred $scriptdir > $errfile
check=`grep Error $errfile | wc -l`
if [ $check -gt 0 ]; then  
    echo; echo " Runpipeline exited with errors!"; echo;
    cat $errfile
    echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
    exit; 
fi
if [ $showout -eq 1 ]; then 
    echo; echo "********* OUTPUT FROM run.sh ***********"
    cat $errfile
fi
exit


#######################################################
############   PIPELINE PER SINGLE CHR  ###############
#######################################################
if [ ! -f $dir/$workdir/inter/third.al ]; then
    echo final alignment file is missing! $dir/$workdir/inter/third.al
    echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
    exit
fi


rm $errchr
$runchrs $workdir $(basename $myref) $notshred $scriptdir >> $errchr

check=`grep Error $errchr | wc -l`
if [ $check -gt 0 ]; then  
	echo; echo " Runpipeline second step exited with errors! check output in" $errchr; echo " here:";
	cat $errchr
	echo; echo " ****  Error! Something is wrong! Giving up! **** "; echo; 
	exit; 
fi
if [ $showout -eq 1 ]; then 
    echo; echo "********* OUTPUT FROM runchrs.sh ***********"
    cat $errchr
fi


