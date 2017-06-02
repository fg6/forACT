# forACT
## nor ready yet! In progress !!
Pipeline to prepare alignments between a Reference fasta and a draft assembly for ACT.

## Instructions
Download repository and install utilities/compile tools: 

	git clone https://github.com/fg6/forACT.git
	cd forACT
	./launchme.sh install
	./launchme.sh test
	
### Run the pipeline for a project: 
#### Step 1: setup folder and scripts:
	$ ./launchme.sh setup </full/path/to/reference> </full/path/to/draft>  </full/path/to/destdir>
	    /full/path/to/reference:  Reference fasta (Please provide full path)
	    /full/path/to/draft:  draft assembly fasta (Please provide full path)
	    /full/path/to/destdir: folder where to run the pipeline (Please provide full path)
#### Step 2: run the pipeline:

	$ cd /full/path/to/destdir
	$ ./mypipeline.sh align
	$ ./mypipeline.sh prepfiles

If everything run smoothly the newly created files will be in /full/path/to/destdir/whole_10000/unique/ folder.
	
#### Step 3: Launch ACT
If Step 2 gave no errors, then launch ACT:

	$ /full/path/to/your/forACT/utils/myscripts/launch_act.sh /full/path/to/destdir/ref/ref.fasta /full/path/to/destdir//whole_10000/unique/foract.al /full/path/to/destdir/whole_10000/unique/foract.fasta

*** Warning *** 

ACT cannot handle genome size > 2.1 GB, for genome of this size run instead single chromosome
at a time in ACT:





