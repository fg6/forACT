# forACT
Pipeline to prepare alignments between a Reference fasta and a draft assembly for ACT.

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
	
To view this project compared to another assembly (up to 5 chromosomes):

	$ ./mypipeline.sh act_select_compare full/path/to/folder/other_assembly_forACT  chr1 chr2 ...
Please note this only works if for the other assembly forACT has been run with the same reference and
same shred/noise/minid parameters	






