# Tempus Bioinformatics Challenge

1. Input VCF file: Challenge_data.vcf
2. Output Annotated file: tempus_challange_variant_annotation.csv
3. Python script for generating annotation: tempus_create_annotation.py

USAGE: python tempus_create_annotation.py -i Challenge_data.vcf > tempus_variant_annotation.csv

#Description of columns  tempus_challange_variant_annotation.csv
1. Variant: variant id represented in the format
	chr-genomic coordinate-Reference Allele-ALT Allele
2. Type of Variation: 
	represented as ts/tv/ins/deletion/snp/complex/mnp
3. site_coverage/depth: 
	Read depth/coverage for the site/locus
4. read supporting variants: 
	Number of reads supporting the variant
5. percentage_variant_read_coverage|percentage_reference_read_coverage:
	Percentage of reads supporting the variant versus those supporting reference reads
6. Allele_Frequency:
	Allele frequency of variant from Broad Institute ExAC Project API
8.Ensemble_gene_id:
	Ensemble gene identifier
9.SNP_id:
	dbSNP identifier
