#!/bin/bash
#set -o errexit

myforACT=$1
thisdir=`pwd`

mkdir -p $myforACT/test_ecoli
cd $myforACT/test_ecoli


### download data
file=ftp://ftp.sanger.ac.uk/pub/users/fg6/forACT_testdata.tar.gz
if [[ ! -d forACT_testdata ]]; then
    if [[ ! -f forACT_testdata.tar.gz ]]; then
	if [[ `wget -S --spider $file 2>&1  | grep exists` ]]; then
	    wget -nv -c $file &> /dev/null 
	    if [[ "$?" != 0 ]]; then echo " Error while downloading the test data. Exiting now"; rm forACT_testdata.tar.gz; exit; fi
	else
	    echo " Error: could not find the test data. Exiting now"; exit; 
	fi
    fi
    tar -xzf forACT_testdata.tar.gz
    if [[ "$?" != 0 ]]; then echo " Error during un-compressing test data. Exiting now"; exit; 
    else
	rm forACT_testdata.tar.gz
    fi
fi
echo " 1. Data downloaded "

### setup project
cd $myforACT/
./launchme.sh setup $myforACT/test_ecoli/forACT_testdata/Escherichiacoli-K-12.fasta  $myforACT/test_ecoli/forACT_testdata/draft.fasta $myforACT/test_ecoli/forACT > /dev/null
echo " 2. Test project setup in " $myforACT/test_ecoli/forACT

### run pipeline
cd $myforACT/test_ecoli/forACT
sed -i 's/lfsjobs=1/lfsjobs=0/g' mysettings.sh
./mypipeline.sh align > /dev/null
./mypipeline.sh prepfiles > /dev/null


testal=`diff $myforACT/test_ecoli/forACT/whole_10000/unique/foract.al $myforACT/test_ecoli/forACT_testdata/results/foract.al | wc -l`
testfasta=`diff $myforACT/test_ecoli/forACT/whole_10000/unique/foract.fasta $myforACT/test_ecoli/forACT_testdata/results/foract.fasta | wc -l`

if [[ $testal == 0 ]] && [[ $testfasta == 0 ]]; then
    echo " 3. Pipeline run: "
    echo  " *****  Test succeded ***** "
    echo; echo " To launch act:"
    echo " " $myforACT/utils/myscripts/launch_act.sh $myforACT/test_ecoli/forACT/ref/ref.fasta $myforACT/test_ecoli/forACT/whole_10000/unique/foract.al $myforACT/test_ecoli/forACT/whole_10000/unique/foract.fasta &

else
    echo " Errors! "
fi
