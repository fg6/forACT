#!/bin/bash
#set -o errexit

myforACT=$1
thisdir=`pwd`

cd $myforACT/utils/mysrcs
mkdir -p mylibs

### Intalling gzstream (it needs zlib!)
if [[ ! -d  gzstream ]]  || [[ ! -f gzstram/gzstream.o ]]; then
    
    rm -rf mylibs
    mkdir mylibs
    cd mylibs
    
    curl -s https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz > gzstream.tgz
    if [[ "$?" != 0 ]]; then
	echo "Error downloading gzstream, try again"
	rm -rf gzstream gzstream.tgz 
	exit
    else
	tar -xvzf gzstream.tgz &> /dev/null
	if [[ "$?" != 0 ]]; then echo " Error during gzstream un-compressing. Exiting now"; exit; fi
	cd gzstream
	make &> /dev/null
	
	if [[ "$?" != 0 ]]; then echo " Error during gzstream compilation. Exiting now"; exit; fi
	test=`make test | grep "O.K" | wc -l`

	if [[ $test == 1 ]]; then echo " "1. gzstream installed; rm ../gzstream.tgz 
	else  echo  " Gzstream test failed. Exiting now"; exit; fi
    fi
fi
 
if [[ ! -f gzstream.o ]]; then 
	echo "  !! Error: gzstream not installed properly!"; 
	exit
fi


#OLD_CPLUS_INCLUDE_PATH=`echo $CPLUS_INCLUDE_PATH`
#export CPLUS_INCLUDE_PATH=$myforACT/utils/mysrcs/mylibs/gzstream:${CPLUS_INCLUDE_PATH}



cd $myforACT/utils/mysrcs/

srcs=( listchrs actnoise  chrpos grabeachchr   n50  samectgpos	splitinfastas  splitreads  writeselctg revertcompl misfinder )
#srcs=( misfinder )

for code in "${srcs[@]}"; do 
    echo $code
    cd $myforACT/utils/mysrcs/$code
    make all
done


cd $myforACT/utils/mysrcs/
echo; echo " Checking installations:"
exes=( mylibs/gzstream/gzstream.o  listchrs/listchrs actnoise/actnoise  chrpos/chrpos grabeachchr/grabeachchr   n50/n50  samectgpos/samectgpos  splitinfastas/splitinfastas  splitreads/splitreads  writeselctg/writeselctg  revertcompl/revertcompl  misfinder/misfinder )

errs=0
for exe in "${exes[@]}"; do
    if [[ ! -f $exe ]]; then 
        echo cannot find $exe: Error! 
        errs=$(($errs+1))
    fi
done
if [  $errs -gt 0 ]; then echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi


## I guess I need to define the gzstream path only during compilation?
#export CPLUS_INCLUDE_PATH=$OLD_CPLUS_INCLUDE_PATH
