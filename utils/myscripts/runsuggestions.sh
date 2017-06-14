#!/bin/bash
set -o errexit


echo; echo "******** Some suggestions for settings and parameters  ********"; echo
echo Once your project is setup in a folder /full/path/to/destdir the default parameters 
echo are defined in your   /full/path/to/destdir/mysettings.sh
echo Before running the pipeline you might want to have a look at the settings and change some of the parameters

echo; echo Some parameters you might want to change:
echo; echo " ###### LSF JOBS PARAMETERS #######"
echo " The pipeline will run mapping jobs using lsf jobs, the parameters for lsf jobs in your /full/path/to/destdir/mysettings.sh are: maxjobs, myqueue, myjobmem, myncpus"
echo "    * maxjobs: this parameter define how many lsf jobs you want to run in parallel, and includes other jobs not related to the forACT pipeline."
echo "               Large genomes can potentially be split in hundreds of jobs, setting maxjob fix a limit to how many jobs are running at a time"
echo "    * myqueue, myjobmem, myncpus should generally work, but you might want to change them if your lsf jobs crash."
echo "                                 f.i. increase the myjobmem if the lsf jobs crashed because the jobs exceeded the requested memory."
echo "                                 If you are comparing genomes with low identity (f.i. different species) the jobs will require longer time and memory"


echo; echo " ###### OTHER PARAMETERS #######"
echo "   * noise: Noise level to cut out. IN PROGRESS"
echo "   * shred: before mapping each contig is shred in chunks, the parameter \"shred\" define the lengths of the chunks."
echo "            Default is 10Kb, can reduce to 1 or 2 kb for smaller genomes, or to look at more detailed mapping. Reducing this parameter will increase the number of lsf jobs."


echo;echo
