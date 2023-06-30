# Prioritizing Unresolved New Genotypes (PUNG)

PUNG further analyzes polymorphic killer cell immunoglobulin like receptors (KIR) short-read sequencing data results analyzed with [PING](https://github.com/Hollenbach-lab/PING), and discriminates potential new alleles sites from annotation errors to streamline new allele validation.

We developed an automated pipeline that mines large dataset outputs from PING to determine whether unresolved genotypes are a product of low-quality sequence data, misaligned reads from another KIR gene, or new alleles. The pipeline performs quality control of new variants, determines if the novel variants are specific to the gene being analyzed, and it aligns reads to a local database containing all KIR variants deposited on the ImmunoPolymorphism Database (IPD-KIR). Finally, PUNG verifies whether the resulting sequences have already been deposited in IPD-KIR or are suggested new alleles.

From our validation dataset, PUNG was able to completely resolve 35.6% of the unresolved genotypes in the dataset. Using Sanger sequencing, we confirmed 81.1% of PUNG suggestions, and 21 newly discovered KIR3DL2 alleles were described while minimizing computational requirements. Thus, we suggest that PING be used in conjunction with PUNG to enable state-of-the-art study of KIR polymorphism.

PUNG currently supports 9 of the 13 _KIR_ loci (KIR3DL3, KIR2DS2, KIR2DL1, KIR2DP1, KIR3DP1, KIR2DL4, KIR2DS1, KIR2DS4, KIR3DL2).

Author: Luciana Vargas

Last updated October 6, 2022

## Installation

Requirements: Linux

Download the PUNG program to your machine.

```bash
git clone https://github.com/lucivargas/PUNG
```

## Usage

In the PUNG_master.R files, modify the input and output variables (lines 8 and 9), and source the script.

## Output files

PUNG produces the following output files:

The *finalAlleleCalls_PUNG.csv* contains an updated genotype data frame.

The *NewAllelesList.csv* contains a data frame with the following columns:
  ind: Sample name
  locus: Locus name
  genotype: Genotype found in PING results
  copy_number: Copy number of the locus in that sample
  genotype_string: Final genotype after analysis with PUNG
  snp_list: Final list of suggested new alleles after analysis with PUNG, if any

The log file (*.log*) contains a summary of the analysis.

The *snp_output* folder provides updated tables of SNP calls for each locus after analysis with PUNG.

## License
[MIT](https://choosealicense.com/licenses/mit/)
