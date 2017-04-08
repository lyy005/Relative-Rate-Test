# Relative-Rate-Test
Tajima Relative Rate Test Of Molecular Clock

The script tests the hypothesis of a molecular evolutionary clock (i.e., a constant rate of molecular evolution) between two groups of sequences (DNA or Amino Acid) using an outgroup sample. 

Usage: perl RelativeRateTest.pl [1.aln.fasta] [2.aln.fasta] [outgroup.aln.fasta] [output] 

The input sequences should be aligned and in fasta format. The script also do Holm–Bonferroni correction. 

References: 
Tajima, F. (1993) Simple methods for testing molecular clock hypothesis. Genetics, 135, 599--607. (Equation 4)
Holm, S. (1979). "A simple sequentially rejective multiple test procedure". Scandinavian Journal of Statistics. 6 (2): 65–70.

