# vcfFilter
Python script to filter VCF files.

`python3 vcfFilter.py -h`

usage: `vcfFilter.py [-h] [--prefix PREFIX] [--minmeanDP MINMEANDP] [--minAD MINAD] [--maf MAF] [--minQ MINQ] [--minGQ MINGQ] [--completeness COMPLETENESS] [-v]`

Alleles are automatically restricted to biallelic.

Non-variants are removed after applying all other filters.

optional arguments:

  -h, --help
    
    show this help message and exit

  --prefix PREFIX
  
    Prefix of the VCF file.

  --minmeanDP MINMEANDP

    The minimum mean read depth supporting a SNP variant.

    Default = 5

  --minAD MINAD
  
    The minimum allele depth for a REF or ALT allele at a given position.

    Default = 0

  --maf MAF
  
    Minimum minor allele frequency. SNPs with a lower minimal minor allele frequency are discarded.

    Default = 0.0

  --minQ MINQ
  
    The minimum quality of the variant site.

    Default = 20

  --minGQ MINGQ
  
    The minimum genotype quality of the call.

    Default = 0.0

  --completeness COMPLETENESS

    SNPs with a lower completeness level are discarded.

    Default = 0.0
