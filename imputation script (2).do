cd D:\Sandbox
global dataset MBBRAINS_chrall_rsq3_final_with_rsids_QC3
display "$dataset"

// changing ids to hrc
bim2dta, bim(D:/Sandbox/$dataset)
use "D:\OneDrive\Obsidian\Software\STATA\data\bim2dta\MBBRAINS_chrall_rsq3_final_with_rsids_QC3\MBBRAINS_chrall_rsq3_final_with_rsids_QC3-bim2dta.dta", clear
keep snpid locname_b37 
rename snpid original
merge 1:1 locname_b37 using D:/OneDrive/Obsidian/Software/STATA/data/.dependencies/bim2hrc/hrc_r1-1_b37.dta
drop if _merge == 2

// identify marker to exclude that are not in hrc1.1
outsheet original if _merge == 1 using $dataset.exclude, non noq replace
!$plink --bfile $dataset --exclude $dataset.exclude --make-bed --out temp

quietly { // run HRC-1000G-check-bim.pl
*!bash -c "wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip"
*!bash -c "wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
local perl  HRC-1000G-check-bim.pl
local hrc   HRC.r1-1.GRCh37.wgs.mac5.sites.tab 
!$plink --freq --bfile temp --out temp

// move to hawk 
perl HRC-1000G-check-bim.pl -b temp.bim -f temp.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -p EUR -v

// back from hawk
!$plink --bfile temp --exclude Exclude-temp-HRC.txt --make-bed --out TEMP1
!$plink --bfile TEMP1 --update-chr Chromosome-temp-HRC.txt  --make-bed --out TEMP2
!$plink --bfile TEMP2 --update-map Position-temp-HRC.txt --make-bed --out TEMP3
!$plink --bfile TEMP3 --flip Strand-Flip-temp-HRC.txt --make-bed --out TEMP4
!$plink --bfile TEMP4 --reference-allele Force-Allele1-temp-HRC.txt --make-bed --out temp-updated
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 1 --out temp-updated-chr1
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 2 --out temp-updated-chr2
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 3 --out temp-updated-chr3
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 4 --out temp-updated-chr4
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 5 --out temp-updated-chr5
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 6 --out temp-updated-chr6
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 7 --out temp-updated-chr7
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 8 --out temp-updated-chr8
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 9 --out temp-updated-chr9
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 10 --out temp-updated-chr10
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 11 --out temp-updated-chr11
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 12 --out temp-updated-chr12
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 13 --out temp-updated-chr13
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 14 --out temp-updated-chr14
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 15 --out temp-updated-chr15
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 16 --out temp-updated-chr16
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 17 --out temp-updated-chr17
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 18 --out temp-updated-chr18
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 19 --out temp-updated-chr19
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 20 --out temp-updated-chr20
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 21 --out temp-updated-chr21
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 22 --out temp-updated-chr22
!$plink --bfile temp-updated --reference-allele Force-Allele1-temp-HRC.txt --make-bed --chr 23 --out temp-updated-chr23

// convert to vcf
foreach chr of numlist 1/23 {
	import delim using temp-updated-chr`chr'.bim, delim(whitespace,collapse) clear
	keep v2 v6 
	outsheet using temp-updated-chr`chr'.reference-allele, non noq replace
	!$plink ///
		--bfile temp-updated-chr`chr' ///
		--chr `chr' ///
		--reference-allele temp-updated-chr`chr'.reference-allele /// Changes in Allele Swaps Handling in Imputation Server 2.0
		--recode vcf ///
		--out temp-updated-chr`chr'
	}
// back to hawk
cat temp-updated-chr23.vcf | sed -e 's/^23/X/' > TEMP.vcf
cat  TEMP.vcf  | sed -e 's/ID=23/ID=X/' > temp-updated-chrX.vcf
			
// archive all the vcf files
/*
bcftools sort temp-updated-chr1.vcf -Oz -o temp-updated-chr1.vcf.gz
bcftools sort temp-updated-chr2.vcf -Oz -o temp-updated-chr2.vcf.gz
bcftools sort temp-updated-chr3.vcf -Oz -o temp-updated-chr3.vcf.gz
bcftools sort temp-updated-chr4.vcf -Oz -o temp-updated-chr4.vcf.gz
bcftools sort temp-updated-chr5.vcf -Oz -o temp-updated-chr5.vcf.gz
bcftools sort temp-updated-chr6.vcf -Oz -o temp-updated-chr6.vcf.gz
bcftools sort temp-updated-chr7.vcf -Oz -o temp-updated-chr7.vcf.gz
bcftools sort temp-updated-chr8.vcf -Oz -o temp-updated-chr8.vcf.gz
bcftools sort temp-updated-chr9.vcf -Oz -o temp-updated-chr9.vcf.gz
bcftools sort temp-updated-chr10.vcf -Oz -o temp-updated-chr10.vcf.gz
bcftools sort temp-updated-chr11.vcf -Oz -o temp-updated-chr11.vcf.gz
bcftools sort temp-updated-chr12.vcf -Oz -o temp-updated-chr12.vcf.gz
bcftools sort temp-updated-chr13.vcf -Oz -o temp-updated-chr13.vcf.gz
bcftools sort temp-updated-chr14.vcf -Oz -o temp-updated-chr14.vcf.gz
bcftools sort temp-updated-chr15.vcf -Oz -o temp-updated-chr15.vcf.gz
bcftools sort temp-updated-chr16.vcf -Oz -o temp-updated-chr16.vcf.gz
bcftools sort temp-updated-chr17.vcf -Oz -o temp-updated-chr17.vcf.gz
bcftools sort temp-updated-chr18.vcf -Oz -o temp-updated-chr18.vcf.gz
bcftools sort temp-updated-chr19.vcf -Oz -o temp-updated-chr19.vcf.gz
bcftools sort temp-updated-chr20.vcf -Oz -o temp-updated-chr20.vcf.gz
bcftools sort temp-updated-chr21.vcf -Oz -o temp-updated-chr21.vcf.gz
bcftools sort temp-updated-chr22.vcf -Oz -o temp-updated-chr22.vcf.gz
bcftools sort temp-updated-chrX.vcf -Oz -o temp-updated-chrX.vcf.gz
/*

** module load bcftools
** dos2unix convert2vcf.gz.sh 
** bash convert2vcf.gz.sh 

** upload to michigan imputation server
/** job parameter
Reference Panel: 1000g-phase-3-v5 (hg19)
Array Build: hg19
rsq Filter: 0
Phasing: eagle v2.4
Allele Frequency Check: EUR
Mode: QC & Imputation		
			
