#!/bin/bash
#this is a test
whattodo=$1

thisdir=`pwd`
launchdir="$(dirname $0)"

if [[ $launchdir == '.' ]]; then
    myforACT=`pwd`
else
    myforACT=$launchdir
fi



if [ $# -lt 1 ] || [ $1 == '-h' ]; then
    echo; echo "  Usage:" $(basename $0) \<command\> 
    echo "  command: command to be run. Options: install, test, setup, suggestions"
    echo "   Check" https://github.com/fg6/forACT/blob/master/README.md "for detailed instructions"
    echo "  "

    exit
fi


if [ $whattodo == "install" ]; then
  ###################################################
  echo; echo " Installing forACT ..."
  ###################################################
  source $myforACT/utils/install.sh  $myforACT
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
  notshred=$3   
  dir=$4
  $myforACT/utils/myscripts/setup.sh  $myforACT $myref $notshred $dir
fi


if [ $whattodo == "test" ]; then
  ###################################################
  echo; echo " Testing with E.coli data"
  ###################################################
  $myforACT/utils/myscripts/runtest.sh  $myforACT
fi



if [ $whattodo == "suggestions" ]; then
  ###################################################
  #echo  Writing some suggestions for parameters  
  ###################################################
  $myforACT/utils/myscripts/runsuggestions.sh  
fi

