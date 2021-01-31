# codon_usage
scripts associated with codon usage


## Anazyle codon usage ##
This script parses a given FASTA file with coding sequences (CDS) and splits these sequences into blocks of 3 nucleotides (codons). The frequency of these codons is counted. A normalization per encoded amino acid is performed to get the contribution of each codon to the respective aminio acid.
If gene expression information is available, the codon usage can be weighted by the transcript abundance. Instead of using a "1" for the appearance of a codon in a CDS, the respective transcript abundance values (e.g. TPM) is used as weight.

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




