#!/bin/bash
set -o errexit

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
# shred in:
shred=10000


mkdir -p $dir

sub1="s#MYFORACT#$myforACT#g"
sub2="s#MYREF#$myref#g"
sub3="s#MYDRAFT#$notshred#g"
sub4="s#MYDESTDIR#$dir#g"

sed $sub1 $myforACT/utils/myscripts/pipeline.sh | sed $sub2 | sed $sub3 | sed $sub4 > $dir/mypipeline.sh
chmod +x $dir/mypipeline.sh

$dir/mypipeline.sh
