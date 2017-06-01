#!/bin/bash

# farm settings
ncpus=10 
nmem=2000 
queue=small

odir=.
mybin=$odir/$1
chmod +x $mybin


# farm settings
hh=1 # hosts
outd=$odir/outs
errd=$odir/errs

 
if [ ! -f $mybin ] ; then
 echo “Sorry binary not available”
 exit 1
fi
mkdir -p $outd $errd
bsub -q $queue -o $outd/out.%J -e $errd/errout.%J -n$ncpus -R"span[ptile=$ncpus] select[mem>$nmem] rusage[mem=$nmem]" -M$nmem  $mybin
