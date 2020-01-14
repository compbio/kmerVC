## Example 1 - Fastq Input File Start

In this example, we start with fastq files for the normal and tumor samples as input to our analysis. After downloading the data, follow the instructions below. Run all commands for this example in the 
_examples/fastq\_start\_example_ directory:

1. Create the jellfyish count file for the reference genome using your desired kmer size and place it in the resources directory with the command:
```
	jellyfish count -m 30 -s 100M -t 24 -C -o ../resources/chrT_30mer.jf ../resources/chrT.fa
```

2. For reference, note of all the input arguments to the kmervc.py script used and their significance:
	- k : Kmer size used for analysis. Ensure jellyfish file of reference genome from step one matches size used here.
	- t1 : First fastq file for tumor sample
	- t2 : Second fastq file for tumor sample
	- c1 : First fastq file for normal sample
	- c2 : Second fastq file for normal sample
	- b : Variants file in bed format
	- o : Output file and directory base name
	- fi : Reference genome file

We will use these downloaded files in the resource directory as the respective inputs:
   - t1 : tumor-1.fq
   - t2 : tumor-2.fq
   - c1 : normal-1.fq
   - c2 : normal-2.fq
   - b : variants.bed
   - fi : chrT.fa
   
Finally, we need to specify the desired kmer size for our analysis. In this example, we will use 30:
   - k : 30
    
Call the program with this command from the _examples/fastq\_start\_example_ directory:

```
python ../../kmervc.py compare -k 30 -t1 ../resources/tumor-1.fq -t2 ../resources/tumor-2.fq -c1 ../resources/normal-1.fq -c2 ../resources/normal-2.fq -b ../resources/variants.bed -o example_1 -fi ../resources/chrT.fa
```
