download all files from the following URLs to TEST folder:
https://dna-discovery.stanford.edu/publicmaterial/software/kmervc/example/
https://dna-discovery.stanford.edu/publicmaterial/software/kmervc/reference/

At TEST folder, run the following command:
python ../kmervc.py compare -k 30 -o test -b variants.bed -j1 tumor.jf -j2 normal.jf -fi chrT.fa
