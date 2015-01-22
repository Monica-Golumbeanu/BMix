BMix is a toolbox for analysing PAR-CLIP data and detecting T-to-C substitutions induced following RNA-protein cross-linking. Candidate binding sites are reported starting from the identified T-to-C substitutions.

The toolbox is presented in the manuscript:
Monica Golumbeanu, Pejman Mohammadi and Niko Beerenwinkel.(2015) Probabilistic modeling of occurring mutations in PAR-CLIP data

Contact: monica.golumbeanu@bsse.ethz.ch

The present repository is organized as follows:
1. File README.txt - the present file
2. File documentation.pdf - contains instructions of use of the toolbox
3. Folder source/ - contains the source code of all the scripts of the toolbox
4. Folder test/ - contains example data ready to use with BMix, as well as the expected outcome

Dependencies:
- Matlab (at least R2013)
- samtools (http://www.htslib.org/)
- awk (http://www.gnu.org/software/gawk/manual/gawk.html)
- bedtools (http://bedtools.readthedocs.org/en/latest/)
