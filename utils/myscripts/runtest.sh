#!/bin/bash
#set -o errexit

myforACT=$1
thisdir=`pwd`


cd $myforACT/test_ecoli

#test done with minimap2, noise 0.1 and id 80. shred says 10000, but actually not shred

### download data

if [[ ! -d forACT_testdata ]]; then
    tar -xzf forACT_testdata.tar.gz
    if [[ "$?" != 0 ]]; then echo " Error during un-compressing test data. Exiting now"; exit; 
    else
	rm forACT_testdata.tar.gz
    fi
fi
echo " 1. Data ready"

### setup project
cd $myforACT/
if [[ $debug == 0 ]]; then
	./launchme.sh setup $myforACT/test_ecoli/forACT_testdata/Escherichiacoli-K-12.fasta  $myforACT/test_ecoli/forACT_testdata/draft.fasta $myforACT/test_ecoli/forACT > /dev/null
else
	./launchme.sh setup $myforACT/test_ecoli/forACT_testdata/Escherichiacoli-K-12.fasta  $myforACT/test_ecoli/forACT_testdata/draft.fasta $myforACT/test_ecoli/forACT
fi
#sed -i "s#aligner=minimap2#aligner=smalt#g#" $myforACT/test_ecoli/forACT/mysettings.sh 
echo " 2. Test project setup in " $myforACT/test_ecoli/forACT


### run pipeline
cd $myforACT/test_ecoli/forACT
sed -i 's/lfsjobs=1/lfsjobs=0/g' mysettings.sh
if [[ $debug == 0 ]]; then
	./mypipeline.sh align > /dev/null
	./mypipeline.sh prepfiles > /dev/null
else
        ./mypipeline.sh align 
        ./mypipeline.sh prepfiles 
fi


### check results
 ./mypipeline.sh check
echo

source mysettings.sh

testal=`diff $myforACT/test_ecoli/forACT/$aligner\_$shred/unique/foract$name_fornoise\_minid$minid.al $myforACT/test_ecoli/forACT_testdata/results/minimap2_foractnonoise0.2_minid10.al | wc -l`
testfasta=`diff $myforACT/test_ecoli/forACT/$aligner\_$shred/unique/foract$name_fornoise\_minid$minid.fasta $myforACT/test_ecoli/forACT_testdata/results/minimap2_foractnonoise0.2_minid10.fasta | wc -l`

if [[ $testal == 0 ]] && [[ $testfasta == 0 ]]; then
    echo " 3. Pipeline run: "
    echo  " *****  Test succeded ***** "
    echo; echo " To launch act:"
    echo " cd " $myforACT"/test_ecoli/forACT; "$myforACT/utils/myscripts/launch_act.sh $myforACT/test_ecoli/forACT/ref/ref.fasta $myforACT/test_ecoli/forACT/$aligner\_$shred/unique/foract$name_fornoise\_minid$minid.al  $myforACT/test_ecoli/forACT/$aligner\_$shred/unique/foract$name_fornoise\_minid$minid.fasta  &

else
    echo " Errors! "
fi
