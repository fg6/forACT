# forACT
Pipeline to prepare alignments between a Reference fasta and a draft assembly for ACT.

## Instructions
Download repository and install utilities/compile tools: 

	git clone https://github.com/fg6/forACT.git
        cd forACT
	./launchme.sh install

### Run the pipeline for a project: 
#### Step 1: setup folder and scripts:
	$ ./launchme.sh setup </full/path/to/reference> </full/path/to/draft>  </full/path/to/destdir>
	/full/path/to/reference:  Reference fasta (Please provide full path)
	/full/path/to/draft:  draft assembly fasta (Please provide full path)
	/full/path/to/destdir: folder where to run the pipeline (Please provide full path)
#### Step 2: run the pipeline:
	cd /full/path/to/destdir
	./mypipeline.sh





