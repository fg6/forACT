#!/bin/bash
#set -o errexit


echo; echo "******** Some suggestions for settings and parameters  ********"; echo
echo Once your project is setup in a folder /full/path/to/destdir the default parameters 
echo are defined in your   /full/path/to/destdir/mysettings.sh
echo Before running the pipeline you might want to have a look at the settings and change some of the parameters

echo; echo Some parameters you might want to change:
echo; echo " ###### The SMALT aligner #######"
echo " the aligner executable is defined in your forACT/utils/myscripts/settings.sh file: it automatically points to a x86_64 compiled version, "
echo " more version you might want to try are in your forACT/utils/mysrcs/mylibs/smalt-0.7.4"
echo " otherwise just pont to your compiled smalt executable in forACT/utils/myscripts/settings.sh"

echo; echo " ###### LSF JOBS PARAMETERS #######"
echo " The pipeline will run mapping jobs using lsf jobs, the parameters for lsf jobs in your /full/path/to/destdir/mysettings.sh are: maxjobs, myqueue, myjobmem, myncpus"
echo "    * maxjobs: this parameter define how many lsf jobs you want to run in parallel, and includes other jobs not related to the forACT pipeline."
echo "               Large genomes can potentially be split in hundreds of jobs, setting maxjob fix a limit to how many jobs are running at a time"
echo "    * myqueue, myjobmem, myncpus should generally work, but you might want to change them if your lsf jobs crash."
echo "                                 f.i. increase the myjobmem if the lsf jobs crashed because the jobs exceeded the requested memory."
echo "                                 If you are comparing genomes with low identity (f.i. different species) the jobs will require longer time and memory"
echo " Alternatively, for very small genomes, the pipeline can run locally changing the parameter":
echo "    * lfsjobs = 0"


echo; echo " ###### For Sinteny Groups: locate inter-chromosome rearrengment/misjoint #######"
echo " The 'misjoints' option will make the pipeline to look for blocks of alignments in the same scaffold that map to different chromosomes of the reference"
echo " These 'blocks' can be real rearrengments in the genome under study, or assembly or scaffolding errors (misjoints)"
echo; echo " * min_len_max"
echo "  Only consider scaffolds with the largest block mapped to a chromosome at least min_len_max bp long"
echo; echo " * min_len, min_len_perc"
echo "  Ignore all shorter blocks which are either < min_len or < min_len_perc % of the largest block"




echo; echo " ###### OTHER PARAMETERS #######"
echo; echo " * shred"
echo "  Each draft-assembly contig is shred in chunks before mapping against the reference. By default, these chunks are 10 K bases long. Change this value by varying the parameter 'shred'"
echo "  Reduce to 1 or 2 kb for smaller genomes, or to look at more detailed mapping. Reducing this parameter will increase the number of lsf jobs."
echo "  The "shred" size corresponds to the red-band size visible in ACT"
echo "  **** Warning!!  **** After changing 'shred' you NEED to run the pipeline from the beginning, so from the 'align' step"

echo; echo " * noise"; echo "  Noise level to cut out. The pipeline cut out isolated alignments < 30K bases. This value can be changed by varying the parameter 'noise'"
echo "  **** Warning!!  **** After changing 'noise' you DON'T NEED to rerun the whole pipeline, simply re-run the 'prepfiles' step" 

echo; echo " * minid"; echo "  Minimum identity: the pipeline will cut out all alignments with identity < 80%. Change this value by varying the parameter 'minid'"
echo "  **** Warning!!  **** After changing 'minid' you DON'T NEED to rerun the whole pipeline, simply re-run the 'prepfiles' step" 

echo;echo
