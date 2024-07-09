# Relative-Rate-Test
Tajima Relative Rate Test Of Molecular Clock

The script tests the hypothesis of a molecular evolutionary clock (i.e., a constant rate of molecular evolution) between two groups of sequences (DNA or Amino Acid) using an outgroup sample. 

Usage: perl RelativeRateTest.pl [1.aln.fasta] [2.aln.fasta] [outgroup.aln.fasta] [output] 

- The input sequences should be aligned and in fasta format (i.e. the *.aln.fasta files should have the same length). The script also does Holm–Bonferroni correction. 

- If you want to run the relative rate test for clusters of species on a tree, please use LINTRE: http://www.kms.ac.jp/~genomelb/takezaki/lintre/index.html

- References: 

- Tajima, F. (1993) Simple methods for testing molecular clock hypothesis. Genetics, 135, 599--607. (Equation 4)
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics. 6 (2): 65–70.

