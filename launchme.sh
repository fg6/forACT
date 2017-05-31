#!/bin/bash
set -o errexit

thisdir=`pwd`
launchdir="$(dirname $0)"

if [[ $launchdir == '.' ]]; then
    myforACT=`pwd`
else
    myforACT=$launchdir
fi

whattodo=$1


if [ $# -lt 1 ] || [ $1 == '-h' ]; then
    echo; echo "  Usage:" $(basename $0) \<command\> 
    echo "  command: command to be run. Options: install, setup, test, run, debug"
    echo "  "
    echo "  "

    exit
fi


if [ $whattodo == "install" ]; then
  ###################################################
  echo; echo " Installing forACT ..."
  ###################################################
  source $thisdir/utils/install.sh  $myforACT
fi


if [ $whattodo == "setup" ]; then
  ###################################################
  # echo; echo " Setup for your project in progress..."
  ###################################################
  if [ $# -lt 4 ]  || [ $2 == '-h' ]; then
      echo; echo "  Usage:" $(basename $0) setup \</full/path/to/reference\>  \</full/path/to/draft\>  \</full/path/to/destdir\>
      echo "   /full/path/to/reference: location of Reference (fasta) "
      echo "   /full/path/to/draft: location of draft assembly (fasta) "
      echo "   /full/path/to/destdir: location of your project "
      echo
      exit
  fi
  myref=$2
  notshred=$3   # change to draft?
  dir=$4
  $thisdir/utils/myscripts/setup.sh  $myforACT $myref $notshred $dir
fi


if [ $whattodo == "test" ]; then
  ###################################################
  echo; echo " Testing with E.coli data"
  ###################################################
  #$thisdir/utils/myscripts/test.sh
fi

if [ $whattodo == "debug" ]; then
  ###################################################
  echo; echo " Looking for possible issues... "
  ###################################################
  #$thisdir/utils/myscripts/debug.sh $myforACT
fi

if [ $whattodo == "run" ]; then
  ###################################################
  echo; echo " Running your project ..."
  ###################################################
  #$thisdir/utils/myscripts/pipeline.sh $myforACT
fi