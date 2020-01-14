## Example 2 - Jellyfish Input File Start

In this example, we start with jellyfish files for the normal and tumor samples as input to our analysis. After downloading the data, follow the instructions below. Run all commands for this example in the 
_examples/jellyfish\_start\_example_ directory:

1. Create the jellfyish count file for the reference genome using your desired kmer size and place it in the resources directory with the command:
```
	jellyfish count -m 30 -s 100M -t 24 -C -o ../resources/chrT_30mer.jf ../resources/chrT.fa
```

2. For reference, note of all the input arguments to the kmervc.py script used and their significance:
	- k : Kmer size used for analysis. Ensure jellyfish file of reference genome from step one matches size used here.
	- j1 : Jellyfish file file for tumor sample
	- j2 : Jellyfish file for normal sample
	- b : Variants file in bed format
	- o : Output file and directory base name
	- fi : Reference genome file

We will use these downloaded files in the resource directory as the respective inputs:
   - j1 : tumor.jf
   - j2 : normal.jf
   - b : variants.bed
   - fi : chrT.fa
   
Finally, we need to specify the desired kmer size for our analysis. In this example, we will use 30:
   - k : 30
    
Call the program with this command from the _examples/jellyfish\_start\_example_ directory:

```
python ../../kmervc.py compare -k 30 -j1 ../resources/tumor.jf -j2 ../resources/normal.jf -b ../resources/variants.bed -o example_2 -fi ../resources/chrT.fa
```
