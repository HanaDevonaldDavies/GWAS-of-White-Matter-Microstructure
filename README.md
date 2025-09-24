# GWAS-of-White-Matter-Microstructure

## ENIGMA White Matter GWAS in the Cardiff Cohort

This repository documents the workflow for processing diffusion MRI and genetic data from the MBBRAINS Cardiff cohort, harmonised with the ENIGMA consortium
 standards.
The primary goal is to generate quality-controlled, imputed genotype data and harmonised diffusion tensor imaging (DTI) phenotypes (FA, MD, AD, RD) for contribution to ENIGMA’s international GWAS meta-analyses of white matter microstructure.

## Project Overview

1. Imaging Side
   - Preprocessing of diffusion MRI using FSL and the ENIGMA-DTI TBSS protocol
   - Extraction of four DTI metrics:
   - FA (Fractional Anisotropy)
   - MD (Mean Diffusivity)
   - AD (Axial Diffusivity)
   - RD (Radial Diffusivity)
   - Quality control through visual inspection and exclusion of artefactual scans.

2. Genetics Side
   - Pre-imputation QC using PLINK.
   - Note: The official ENIGMA Pre-Imputation QC protocol was attempted but not compatible with the Cardiff dataset. A custom QC workflow was developed (details below).
   - Conversion to VCF and imputation on the Michigan Imputation Server
   - Post-imputation QC and filtering.

3. GWAS
   - Association analysis between imputed variants and DTI phenotypes.
   - Covariates include age, sex, and genetic principal components.
   - Output formatted for ENIGMA contribution.

4. Post-GWAS
   - LD Score Regression (heritability and genetic correlations).
   - Fine-mapping of associated loci.
   - Regional visualisation with LocusZoom.

## Pipeline

1. Imaging Pipeline
   - Data acquisition: 30+2 directions, b=1000 s/mm².
   - Preprocessing: BET → Eddy → DTIFIT → FA/MD/AD/RD maps.
   - TBSS: skeletonisation and atlas-based extraction (25 tracts).
   - QC: visual inspection, exclusion of poor alignments.

2. Genetics Pipeline
1. Pre-Imputation QC (Custom Workflow)
   - Start: MBBrains_R123_3oct2016_step5
   - SNPs: 233,075
   - Individuals: 385

Steps applied:
   - Remove SNPs with missing chromosome/position (0 removed).
   - Remove indels/missing alleles (97 removed).
   - Rename SNPs to HRC1.1 nomenclature.
   - Convert build hg18 → hg19.
   - Filter:
     - Allele frequency outliers (>10 SD vs reference, 13 removed).
     - SNP missingness >2% (341 removed).
     - HWE p < 1e-6 (0 removed).
     - Individual missingness >2% (0 removed).
     - Heterozygosity >4 SD (2 removed).
     - Relatedness (IBD > .01, 0 removed).

Finish: MBBRAINS_R123-HumanCoreExome-12v1
    - SNPs: 232,624
    - Individuals: 385

2. Imputation
    - Reference panel: 1000 Genomes Phase 3 v5 (hg19).
    - Phasing: Eagle v2.4.
    - Rsq filter = 0.
    - EUR allele frequency check.
    - Mode: QC & Imputation.

3. Post-Imputation QC
    - Example (chr1):
    - Input = 18,367 markers.
    - Post-imputation = 3,738,261 markers.
    - After INFO (R² < 0.3 removed): 977,805 markers.
    - After MAF > 0.01: 706,237 markers.
