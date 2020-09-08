download test.zip to TEST folder from:
https://dna-discovery.stanford.edu/publicmaterial/software/kmervc/
unzip the file without folder structure.

Before running the command, make sure to install pandas (for Python), bedtools and jellyfish

At TEST folder, run the following command:
python ../kmervc.py compare -k 31 -j1 tumor_31mer.jf -j2 normal_31mer.jf -b variants.bed -o test -fi chrT.fa
