# forACT
Pipeline to prepare alignments between a Reference fasta and a draft assembly for ACT.

## Instructions
Download repository: 

	git clone https://github.com/fg6/forACT.git
Usage:
	$ ./launchme.sh <command> 
	  command: command to be run. Options: install, setup, test, run, debug

### Step 1. Install utilities and compile tools: 
	$ ./launchme.sh install
### Step 2. Prepare pipeline for a project: 
	$ ./launchme.sh setup </full/path/to/reference> </full/path/to/draft>  </full/path/to/destdir>

where:
	/full/path/to/reference:  Reference fasta (Please provide full path)
	/full/path/to/draft:  draft assembly fasta (Please provide full path)
	/full/path/to/destdir: folder where to run the pipeline (Please provide full path)





