# aligner2counts
This repository contains script to process read mapping file (.bam) on the fly with mapping process and produce 3 outputs (i) read abundance of contigs, (ii) fraction of reads shared between contigs if found (split-read mapping) and (iii) alignment file in `.sam` format.
As of now, it is specialized to process only bowtie mapping.

# Installation
	g++ -o aligner2counts aligner2counts.cpp -O3
Require g++ version `>=9.4.0`

# Usage
On the fly mapping and processing,
	bowtie2 -q -x <ref_index> <reads.fq> | ./aligner2counts <output_directory> <sample_id>

Process `.bam` or `.sam` file
	bamtools view <input_alignment> | ./aligner2counts <output_directory> <sample_id>
