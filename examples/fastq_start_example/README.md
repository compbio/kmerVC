# Fastq Input Example

In this example, we start with fastq files for the normal and tumor samples as input to our analysis. After following the directions to download the data in the _resources_ directory, proceed with the instructions below:

1. Create the jellfyish count file for the reference genome using your desired kmer size and place it in the resources directory with the command:
	jellyfish count -m 30 -s 100M -t 24 -C -o chrT_30mer.jf ../resources/chrT.fa

2. Take note of all the input arguments to the kmervc.py script used:
	- k : Kmer size used for analysis. Ensure jellyfish file of reference genome from step one matches size used here.
	- t1 : First fastq file for tumor sample
	- t2 : Second fastq file for tumor sample
	- c1 : First fastq file for normal sample
	- c2 : Second fastq file for normal sample
	- b : Variants file in bed format
	- o : Output file and directory base name
	- fi : Reference genome file

Given these specifications and the inputs from the resources directory, call the command:
```
python ../../kmervc.py compare -k 30 -t1 ../resources/tumor-1.fq -t1 ../resources/tumor-2.fq -c1 ../resources/normal-1.fq -c2 ../resources/normal-2.fq -b ../resources/variants.bed -o example_1 -fi ../resources/chrT.fa
```
