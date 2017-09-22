#!/bin/bash
#set -o errexit

myforACT=$1
thisdir=`pwd`

cd $myforACT/utils/mysrcs
mkdir -p mylibs

### Intalling gzstream (it needs zlib!)
if [[ ! -d  mylibs/gzstream ]]  || [[ ! -f mylibs/gzstream/gzstream.o ]]; then
    
    rm -rf mylibs
    mkdir mylibs
    cd mylibs
    
    #curl -s https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz > gzstream.tgz
    wget https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz 

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
 
cd $myforACT/utils/mysrcs
if [[ ! -f mylibs/gzstream/gzstream.o ]]; then 
	echo "  !! Error: gzstream not installed properly!"; 
	exit
fi


cd $myforACT/utils/mysrcs
if [[ ! -d  mylibs/smalt-0.7.4 ]]; then
    cd mylibs
    wget ftp://ftp.sanger.ac.uk/pub/resources/software/smalt/smalt-0.7.4.tgz
    tar -xvzf smalt-0.7.4.tgz
    rm -f smalt-0.7.4.tgz
fi
cd $myforACT/utils/mysrcs
if [[ ! -f mylibs/smalt-0.7.4/smalt_x86_64  ]]; then 
	echo "  !! Error: smalt not installed properly!"; 
	exit
fi

cd $myforACT/utils/mysrcs
if [[ ! -f  mylibs/minimap2/minimap2 ]]; then
    cd mylibs
    git clone https://github.com/lh3/minimap2.git
    cd minimap2
    git checkout ea5a0cd17da3752d396473c56f074a3ab5ede613  # release 2.2
    make 
fi

cd $myforACT/utils/mysrcs
if [[ ! -f mylibs/minimap2/minimap2 ]]; then 
        echo "  !! Error: minimap2 not installed properly!"; 
        errs=$(($errs+1))
        exit
fi


cd $myforACT/utils/mysrcs/

srcs=( listchrs actnoise  chrpos grabeachchr   n50  samectgpos	splitinfastas  splitreads  writeselctg revertcompl misfinder locate_misjoints )

for code in "${srcs[@]}"; do 
    cd $myforACT/utils/mysrcs/$code

    if [[ ! -f $code ]] || [[ $code -ot $code.cpp ]] || [[ $code -ot $myforACT/utils/mysrcs/myinc/macro.h ]]; then
	echo "  " $code
        make all
    fi

done


cd $myforACT/utils/mysrcs/
echo; echo " Checking installations:"
#exes=( mylibs/gzstream/gzstream.o  mylibs/smalt-0.7.4/smalt_x86_64  listchrs/listchrs actnoise/actnoise  chrpos/chrpos grabeachchr/grabeachchr   n50/n50  samectgpos/samectgpos  splitinfastas/splitinfastas  splitreads/splitreads  writeselctg/writeselctg  revertcompl/revertcompl  misfinder/misfinder )

errs=0
for code in "${srcs[@]}"; do
    if [[ ! -f $myforACT/utils/mysrcs/$code/$code ]]; then 
        echo cannot find $exe: Error! 
        errs=$(($errs+1))
    fi
done
if [  $errs -gt 0 ]; then echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi




