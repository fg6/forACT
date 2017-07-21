## Requirements
You will need a g++ version >=4.7  and the zlib in your path for installation

# forACT
Pipeline to prepare alignments between a Reference fasta and a draft assembly for ACT.

Main Pipeline steps:

* Shreds draft assembly contigs/scaffolds in 10K bases (can vary length of chunks by changing parameter "shred")
* Maps the shredded draft assembly against Reference  (SMALT aligner)
* Checks orientation of each contig/scaffold with respect to Reference and change to complement the contigs/scaffolds for which majority of chunks are complements. This creates an equivalent draft assembly with only forward contigs/scaffolds
* Shreds the forward draft assembly contigs/scaffolds in 10K bases (can vary length of chunks by changing parameter "shred")
* Maps the shredded forward draft assembly against Reference (SMALT aligner)
* Filters out:

	- alignments with identity < 80% (can vary minimum identity by changing parameter "minid")
	- isolated alignments shorter than 1% of the contig (can vary this limit by changing parameter "noise")
* Re-order position of forward contigs/scaffolds to match each contig/scaffold alignment position in the Reference (when one contig/scaffold map to more than one place in the Reference, the major alignment is considered)
* Re-define alignment positions according to ACT positioning scheme

## Instructions
Download repository and install utilities/compile tools: 

	$ git clone https://github.com/fg6/forACT.git
	$ cd forACT
	$ ./launchme.sh install
	$ ./launchme.sh test

Get information about parameters and settings:

	$ ./launchme.sh suggestions
	
### Run the pipeline for a project: 
#### Step 1: setup folder and scripts:
	$ ./launchme.sh setup </full/path/to/reference> </full/path/to/draft>  </full/path/to/destdir>
	    /full/path/to/reference:  Reference fasta (Please provide full path)
	    /full/path/to/draft:  draft assembly fasta (Please provide full path)
	    /full/path/to/destdir: folder where to run the pipeline (Please provide full path)
#### Step 2: run the pipeline:

	$ cd /full/path/to/destdir
	$ ./mypipeline.sh align
Check if alignment ran smoothly and with no errors:

	$ ./mypipeline.sh check 
	
If everything run smoothly then prepare the files for ACT: 

	$ ./mypipeline.sh prepfiles
If no errors, the newly created files will be in /full/path/to/destdir/whole_10000/unique/ folder.

#### Step 3: write a report of the assemblies and alignments, and possible misassemblies/re-arrengements:

	$ ./mypipeline.sh report
	In progress: some variables like "Reference Coverage" or "Structure coverage" definition 
		are still not fixed, do not take them too seriously for now
	
#### Step 4: Launch ACT

If Step 1 and 2 gave no errors, then launch ACT:

	$ ./mypipeline.sh act  
To view this project compared to another assembly:

	$ ./mypipeline.sh act_compare full/path/to/folder/other_assembly_forACT
Please note this only works if for the other assembly forACT has been run with the same reference and
same shred/noise/minid parameters

#### Step 5: Launch ACT for genome size > 2.1 GB
*** Warning *** 
ACT cannot handle genome size > 2.1 GB, for genome of this size run instead single chromosome
at a time in ACT:

List possible chromosomes to view:

	$ ./mypipeline.sh act_select list
View up to 5 chromosomes:

	$ ./mypipeline.sh act_select chr1 chr2 ...
where chr1..chr2 are the name of the chromosomes to visualize.
	
To view this project compared to another assembly (up to 5 chromosomes):

	$ ./mypipeline.sh act_select_compare full/path/to/folder/other_assembly_forACT  chr1 chr2 ...
Please note this only works if for the other assembly forACT has been run with the same reference and
same shred/noise/minid parameters	






