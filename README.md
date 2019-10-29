# Utilizing k-mer counts for somatic mutation calls
This repo contains the kmerVC software scripts that function to evaluate the somatic mutation validity of a set of variants as input. The code here was used to produce the results in the paper "Feasibility of k-mers counts for somatic mutation calls", which is currently in submission.

## Requirements
This code has been tested with the following software versions:
* Python 2.7.13 
* Pandas 0.20.3 
* Numpy 1.13.1 
* Scipy 0.19.1 
* Bedtools 2.25.0 
* Jellyfish 2.2.6 

Aside from bedtools and Jellyfish, we recommend the using the Anaconda Python Distribution to install all the required packages. The standard instructions from the providers of each package can be used. As for bedtools and Jellyfish, the respective installation instructions are outlined below:

### Bedtools
Using the Anaconda Python Distribution, this can be done with the command:
> conda install -c bioconda bedtools

### Jellyfish
The k-mer counting command line program can be installed from their respective github release at: https://github.com/gmarcais/Jellyfish

