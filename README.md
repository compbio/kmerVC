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

Aside from bedtools and Jellyfish, we recommend the using the Anaconda Python Distribution to install all the required packages. The standard instructions from the providers of each package can be used. As for bedtools and Jellyfish, the respective installation instructions are outlined below. For the above installation a requirements.txt document has been created to streamline the process. Creation of a virtual environment with package nstallation can be carried out by with the command:

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
The k-mer counting command line program can be installed from their respective github release at: https://github.com/gmarcais/Jellyfish

### Testing

You can test your installation of kmerVC on the examples in the test directory. The tests were fashioned as a shell script that calls the kmerVC program and can be run as follows:

```
cd test
bash run_test.sh
```
