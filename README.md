# Utilizing k-mer counts for somatic mutation calls
This repo contains the kmerVC software scripts that function to evaluate the somatic mutation validity of a set of variants as input. The code here was used to produce the results in the paper "Feasibility of k-mers counts for somatic mutation calls", which is currently in submission. The lab software page for kmerVC can be accessed [here.](https://dna-discovery.stanford.edu/research/software/kmervc/)

## Requirements
This code has been tested with the following software versions:
* Python 2.7.13
* Pandas 0.20.3
* Numpy 1.13.1
* Scipy 0.19.1
* Bedtools 2.25.0
* Jellyfish 2.2.6

Aside from bedtools and Jellyfish, we recommend using the Anaconda Python Distribution to install all the required packages. The standard instructions from the providers of each package can be used. As for bedtools and Jellyfish, the respective installation instructions are outlined below. For the above installation a requirements.txt document has been created to streamline the process. Creation of a virtual environment with package nstallation can be carried out by with the command:

```
conda create --name <envname> --file requirements.txt
```

with the environment name <envname> of your choosing. Activate your newly created virtual environment with:

```
source activate <envname>
```

and proceed with the installation of bedtools and jellyfish in this virtual environment below.

### Bedtools
Using the Anaconda Python Distribution, this can be done with the command:

```
conda install -c bioconda bedtools
```

### Jellyfish
The k-mer counting command line program can be installed from their respective github release at: https://github.com/gmarcais/Jellyfish. Follow their up-to-date installation instructions to download the Jellyfish software. Finally, ensure that the jellyfish program is added to your path with the shell command (replacing X with number of downloaded version):

```
export PATH=/programs/jellyfish-2.2.X/bin:$PATH
```

### Examples

You can test your installation of kmerVC on the otulined below. These can be carried out in the examples directory.

First, download the necessary files from the [lab software page](https://dna-discovery.stanford.edu/publicmaterial/software/kmervc/reference/) and place them in the _examples/resources_ folder of the cloned github project.

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
