#!/bin/bash
set -o errexit


thisdir=`pwd`
whattodo=$1


if [ $# -lt 1 ]  || [ $1 == '-h' ]; then
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
  export forACTdir=$thisdir
  source $thisdir/utils/install.sh 
fi


if [ $whattodo == "setup" ]; then
  ###################################################
  echo; echo " Setup for your project in progress..."
  ###################################################
  #$thisdir/utils/setup.sh 
fi


if [ $whattodo == "test" ]; then
  ###################################################
  echo; echo " Testing with E.coli data"
  ###################################################
  #$thisdir/utils/test.sh
fi

if [ $whattodo == "debug" ]; then
  ###################################################
  echo; echo " Looking for possible issues... "
  ###################################################
  #$thisdir/utils/debug.sh
fi

if [ $whattodo == "run" ]; then
  ###################################################
  echo; echo " Running your project ..."
  ###################################################
  #$thisdir/utils/pipeline.sh
fi