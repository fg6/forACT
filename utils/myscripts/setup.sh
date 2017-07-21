#!/bin/bash
#set -o errexit

myforACT=$1
myref=$2
notshred=$3   # change to draft?
dir=$4

if [ $# -lt 4 ]  || [ $1 == '-h' ]; then
    echo; echo "  Usage:" $(basename $0) \</full/path/to/forACT\> \</full/path/to/reference\>  \</full/path/to/draft\>  \</full/path/to/destdir\>
    echo "   /full/path/to/forACT: location of forACT "
    echo "   /full/path/to/reference: location of Reference (fasta) "
    echo "   /full/path/to/draft: location of draft assembly (fasta) "
    echo "   /full/path/to/destdir: location of your project "
    exit
fi

thisdir=`pwd`
scriptdir=$myforACT/utils/myscripts

mkdir -p $dir

sub1="s#MYFORACT#$myforACT#g"
sub2="s#MYREF#$myref#g"
sub3="s#MYDRAFT#$notshred#g"
sub4="s#MYDESTDIR#$dir#g"

sed $sub1 $myforACT/utils/myscripts/settings.sh | sed $sub2 | sed $sub3 | sed $sub4 > $dir/mysettings.sh
cp  $myforACT/utils/myscripts/pipeline.sh $dir/mypipeline.sh
#echo $0 $1 $2 $3 $4 > $dir/setupas.txt
echo $myforACT/launchme.sh setup  $2 $3 $4  > $dir/setupas.txt
chmod +x $dir/*.sh

echo; echo  " Your new project is set in folder: " $dir;
echo " Modify relevant paramters in " $dir/mysettings.sh
