# codon_usage

The scripts are associated with the analysis of codon usage. The scripts are developed for the use with Python v2.7.


## Analyze codon usage ##
This script parses a given FASTA file with coding sequences (CDS) and splits these sequences into blocks of 3 nucleotides (codons). The appearance of these codons is counted. A normalization per encoded amino acid is performed to get the contribution of each codon to the respective aminio acid.
If gene expression information is available, the codon usage can be weighted by the transcript abundance. Instead of using a "1" for the appearance of a codon in a CDS, the respective transcript abundance value (e.g. TPM) is used as weight.

```
python analyze_codon_usage.py --in <FILE> --out <DIR>

Mandatory:
  --in STR          FASTA file containing CDS
  --out STR         Output file

Optional:
  --exp STR         Expression data file
```

`--in` specifies a multiple FASTA file containing a collection of protein coding sequences (CDS).

`--out` specifies the output file.

`--exp` specifies a data file with expression values. The first row contains the sample names. The sequence names in the first column need to match the CDS names in the input file.



## Get translational bottleneck ##
Just replacing all codons in a CDS by the most frequently used one for the respective amino acid might not work very well. This could be due to the fact that domains need time to fold which requires the ribosome to slow down or even stall at the right positions on the mRNA. Therefore, it might be necessary to identify and preserve such translational bottlenecks when performing codon optimization. This script helps to identify position on a mRNA where codons with low frequency are clustered.


```
python get_translational_bottle_necks.py --in <FILE> --out <DIR>

Mandatory:
  --in    STR         FASTA file containing CDS
  --codon STR         Codon usage table
  --out   STR         Output folder

Optional:
  --win   INT         Window size
```

`--in` specifies a multiple FASTA file containing a collection of protein coding sequences (CDS).

`--codon` specifies the codon usage table produced by the script described above. Other tables could be used, but the amino acid need to be specified in the second column and the usage (0.0-1.0) needs to be in the fourth column of a TAB-separated file.

`--out` specifies the output folder where all result files will be stored. This folder will be created if it does not exist already.

`--win` specifies a window size for the analysis. The frequency of codons in this window is checked to find regions where the translation slows down.




