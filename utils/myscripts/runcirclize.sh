#!/bin/bash
#set -o errexit

thisdir=`pwd`
source $thisdir/mysettings.sh

mkdir -p  $outdir/figs

rscript=`which Rscript | wc -l`

if [[ $rscript == 0 ]]; then
    echo "  I need Rscript to be installed to use circlize! "
    echo "  Please install Rscript and try again"
    exit
fi

circlize_found=`Rscript -e 'if(!suppressWarnings(suppressMessages(require(circlize)))){print("missing")}else{print("found!")}' | grep found | wc -l`


if [[ $circlize_found == 0 ]]; then
    echo "   the R package \"circlize\" not found, please installed it from R using: install.packages(\"circlize\") "
    echo "      before running mypipeline.sh circlize "
    exit
elif [[ $circlize_found != 1 ]]; then
    echo "   Something went wrong, cannot find circlize. Try again "
    exit
fi

cd $outdir/figs


if [[ ! -f chrs.file ]] || [[ ! -f als.file ]]; then 
	grep -v $unplaced $misals_file > als.file	
	grep -v $unplaced $chrass_file > chrs.file
else
	echo " Using existing files, if need to recreate them pls first remove them:"
	echo $outdir/figs/als.file and $outdir/figs/chrs.file
fi
    
if [[ ! -f ./my_plot.png ]]; then
    echo "now running..."
    rm -f err.log 
    Rscript  $scriptdir/circosplot.R $title_fig > err.log 
    exe_error=$?
    error_status=$(( `grep "Error" err.log | wc -l` == 1 ? 1 : $exe_error ))  

    if [[ $error_status == 1 ]]; then
	cat err.log
	exit
    fi

    err=0
    if [[ -f ./my_plot.png ]]; then
	mv  ./my_plot.png   $circl_fig.png
    else
	echo " Error: figure" my_plot.png "is missing..."
	err=$(($err+1))
    fi
    if [[ -f ./noscaff.png ]]; then
	mv  ./noscaff.png   $circl_fig\_noscaff.png
    else
	echo " Error: figure" noscaff.png "is missing..."
 	err=$(($err+1))
   fi
    rm -f ./my_plot.png ./noscaff.png 
fi

if [[ $err == 0 ]]; then
    echo; echo " Circos plots in folder " $outdir/figs
else
    echo; echo " Error while producing the circos plot "
fi

#rm als.file chrs.file
