# aligner2counts
This repository contains script to process read mapping file (.bam) on the fly with mapping process and produce 3 outputs (i) read abundance of contigs, (ii) fraction of reads shared between contigs if found (split-read mapping) and (iii) alignment file in `.sam` format.
As of now, it is specialized to process only bowtie mapping.

# Installation
	g++ -o aligner2counts aligner2counts.cpp -O3
	export PATH=${PATH}/$(pwd) (save this line in .bashrc or .bash_profile for all time access)
Require gcc version `>=9.4.0`


# Usage
	Usage: aligner command | aligner2counts outdir outputname [--minlength N] [--no-coverage] [--single] [--only-mapids] [--qcov X] [--seq-id Y]
	Options:
  	<outdir>            Directory where output files will be stored.
  	<outputname>        Base name for the output files (without extensions).

 	 --minlength N       Minimum length of contigs to be considered. (Default: N=1000)
      	                Alignments for contigs shorter than this length will be ignored.

  	--no-coverage       Do not output coverage. (Optional)
      	                This flag disables the output of contig coverage.
	
  	--single            Process input as single-end reads. (Optional)
      	                By default, paired-end reads are expected unless this flag is set.

  	--only-mapids       Output only the mapping identifiers. (Optional)
   	                   This flag restricts the output to just the IDs of reads and mapped contigs.

  	--qcov X            Minimum query/read coverage threshold X. (Optional, default=99.0%)
  	                    Specifies the minimum percentage of query sequence that must be aligned.

  	--seq-id Y          Minimum sequence identity percentage Y. (Optional, default=97.0%)
     	                 Filters alignments by requiring at least Y% identity.

  	-h, --help          Display this help message and exit.

# Example
One go mapping and processing,

 	bowtie2 -q -x <ref_index> <reads.fq> | aligner2counts <output_directory> <sample_id>

Process only `.bam` or `.sam` (require samtools preinstalled)

 	samtools view -h <input_alignment> | aligner2counts <output_directory> <sample_id>

If you have `.sam` file, you can also use `cat <input>.sam | aligner2counts <output_directory> <sample_id>
