# Utilizing k-mer counts for somatic mutation calls
This repo contains the kmerVC software scripts that function to evaluate the somatic mutation validity of a set of variants as input. The code here was used to produce the results in the paper "Feasibility of k-mers counts for somatic mutation calls", which is currently in submission. The lab software page for kmerVC can be accessed [here.](https://dna-discovery.stanford.edu/research/software/kmervc/)

## Installation and testing run
- git clone https://github.com/compbio/kmerVC.git
- Navigate into TEST directory and download test.zip from https://dna-discovery.stanford.edu/publicmaterial/software/kmervc/ and unzip test.zip.
- Install bedtools and Jellfish. See Requirments below.
- You need several Python pacakges inlcuding pandas. See Requirments below.
- Type the following command that starts with jellyfish outcomes at TEST directory:
```
python ../kmervc.py compare -k 30 -j1 tumor_30mer.jf -j2 normal_30mer.jf -b variants.bed -o test -fi chrT.fa
```

Type the following command that starts with fastq files:
```
python ../kmervc.py compare -k 30 -t1 tumor-1.fq -t2 tumor-2.fq -c1 normal-1.fq -c2 normal-2.fq -b variants.bed -o example_1 -fi chrT.fa
```
- These commands will generate 


## Requirements
This code has been tested with the following software versions:
* Python 2.7.13 and Python 3.6
* Pandas 0.20.3
* Numpy 1.13.1
* Scipy 0.19.1
* Bedtools 2.29.2
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
conda install -c bioconda bedtools=2.29.2
```

### Jellyfish
The k-mer counting command line program can be installed from the respective github release at: https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.6. Follow their up-to-date installation instructions to download the Jellyfish software. Finally, ensure that the jellyfish program is added to your path with the shell command (replacing X with number of downloaded version):

```
export PATH=/programs/jellyfish-2.2.6/bin:$PATH
```

## Script Command Line Usage

usage: kmervc.py [-h] (-v VCF_INPUT | -b BED_INPUT) [-t1 TEST\_FASTQ1]
                 [-t2 TEST_FASTQ2] [-c1 CONTROL_FASTQ1] [-c2 CONTROL_FASTQ2]
                 [-j1 JELLYFISH_TEST] [-j2 JELLYFISH_CONTROL] -k KMER_SIZE -o
                 OUTPUT_NAME [-fi REFERENCE_GENOME_FASTA] [-d DELIMITER] [-m]
                 [-r] [-poi] [-a ALPHA]
                 {compare}

required positional arguments:{compare}

required arguments:
  - -k, --kmer\_size KMER\_SIZE : Size of kmer to use for analysis
  - -o, --output\_name OUTPUT\_NAME : Output file directory name
  - -v, --vcf VCF\_INPUT : Input vcf file
    - OR 
  - -b, --bed BED\_INPUT : Input ved file

fastq\_group: fastq input files

  - -t1, --test1 TEST\_FASTQ1 : Fastq file 1 from test sample
  - -t2, --test2 TEST\_FASTQ2 : Fastq file 2 from test sample
  - -c1, --control1 CONTROL\_FASTQ1 : Fastq file 1 from control sample
  - -c2, --control2 CONTROL\_FASTQ2 : Fastq file 2 from control sample

jellyfish\_group: jellyfish input files

  - -j1, --jellyfish\_test JELLYFISH\_TEST : Jellyfish file of test input
  - -j2, --jellyfish\_control JELLYFISH\_CONTROL : Jellyfish file of control input

optional arguments:
  - -h, --help            show this help message and exit
  - -fi, --reference\_genome\_fasta REFERENCE\_GENOME\_FASTA : Reference genome fasta file to use, if different than default
  - -m, --microsatellite  Flag : indicating if doing microsequence analysis with : respective vcf file
  - -r, --rna             Flag : indicating if doing RNA analysis
  - -poi, --poisson       Flag : indicating if using doing poisson distribution : for variant analysis
  - -a, --alpha ALPHA : Alpha value used in hypothesis testing

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

### Example 3 - Jellyfish Input With Reference Genome

In this example, we start with jellyfish files for the normal and tumor samples as input to our analysis. After downloading the data, follow the instructions below. Run all commands for this example in the
_examples/reference\_genome\_example_ directory:

1. Create the jellfyish count file for the reference genome using your desired kmer size and place it in the resources directory with the command:
```
	jellyfish count -m 30 -s 100M -t 24 -C -o ../resources/grch38_canonical_chrs_chrM_30mer.jf ../resources/grch38_canonical_chrs_chrM.fa
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
   - b : genome\_variants.bed
   - fi : chrT.fa

Finally, we need to specify the desired kmer size for our analysis. In this example, we will use 30:
   - k : 30

Call the program with this command from the _examples/jellyfish\_start\_example_ directory:

```
python ../../kmervc.py compare -k 30 -j1 ../resources/tumor.jf -j2 ../resources/normal.jf -b ../resources/genome_variants.bed -o example_3 -fi ../resources/grch38_canonical_chrs_chrM.fa
```
