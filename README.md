
# forACT
Pipeline to prepare alignments between a Reference fasta and a draft assembly; then uses the ACT pipeline (http://www.sanger.ac.uk/science/tools/artemis-comparison-tool-act) to visualize genome comparison interactively.
Printout a summary of the assembly analysed, and possible inter-chromosomal misjoints (if comparing with a reference)
For chromosome assigned scaffolds there is also the possibiliy to have circos-like plots for a quick general view.
Please notice that ACT will not work for genomes ~> 2.1 Gbp; in this case you can still plot the comparison but have to select 1 to 5 chromosomes to show at a time. To decide which chromosomes to plot, you can make use of the misassembly report to see if and which chromosomes get mixed-up by the assembly.

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

## Requirements
The new version of Artemis require maven to properly install act, please make sure it's in your PATH.
You will need a g++ version >=4.7 and the zlib in your path for installation

The aligner executable (Smalt from http://www.sanger.ac.uk/science/tools/smalt-0) is defined in your forACT/utils/myscripts/settings.sh file: it automatically points to a 
x86\_64 compiled version, more versions you might want to try are in your forACT/utils/mysrcs/mylibs/smalt-0.7.4
otherwise just pont to your compiled smalt executable in forACT/utils/myscripts/settings.sh

## External packages
forACT downloads and installs the gzstream library to handle gzip input files (https://www.cs.unc.edu/Research/compgeom/gzstream/). It also uses the external software ACT, smalt and minimap2 (https://github.com/lh3/minimap2).

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
If no errors, the newly created files will be in /full/path/to/destdir/[aligner]_10000/unique/ folder.

#### Step 3: write a report of the assemblies and alignments, and possible misassemblies/re-arrengements:

	$ ./mypipeline.sh report
	In progress: some variables like "Reference Coverage" or "Structure coverage" definition 
		are still not fixed, do not take them too seriously for now
#### Step 4: for chromosome-assigned super-scaffolds: look for inter-chromosome rearrengments/misjoints:

	$ ./mypipeline.sh misjoints
	
and then

	$ ./mypipeline.sh circlize
	
for circos type plots. Be aware that in these plots only
and only the alignment to the major chromosome per each scaffold and the possible inter-chromosomal misassemblies are visualized.
	
	
#### Step 5: Launch ACT

If Step 1 and 2 gave no errors, then launch ACT:

	$ ./mypipeline.sh act  
To view this project compared to another assembly:

	$ ./mypipeline.sh act_compare full/path/to/folder/other_assembly_forACT
Please note this only works if the other assembly forACT has been run with the same reference and
same shred/noise/minid parameters

#### Step 6: Launch ACT for genome size > 2.1 GB
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






