# COBI-8302-Statistical-Genetics
This is the course materials for Session 19: building polygenic risk scores for course COBI 8302 in Fall 2026 at the University of Virginia.

## Table of Contents
- [Calculate PRS using pre-trained weights](#calculate-prs-using-pre-trained-weights)
  - [Download Summary Statistics from PGSCatalog](#download-summary-statistics-from-pgscatalog)
  - [Data cleaning](#data-cleaning)
    - [Ambiguous SNPs](#ambiguous-snps)
    - [Standard GWAS QC (will do this step when calculating the PRS)](#standard-gwas-qc-will-do-this-step-when-calculating-the-prs)
    - [Mismatching SNPs and/or Duplicate SNPs](#mismatching-snps-andor-duplicate-snps)
    - [Save the score file for PLINK](#save-the-score-file-for-plink)
  - [Calculation of PRS](#calculation-of-prs)
  - [Visualize the PRS](#visualize-the-prs)
- [Train PRS models using GWAS summary statistics](#train-prs-models-using-gwas-summary-statistics)
  - [Prepare the summary statistics for PRScs](#prepare-the-summary-statistics-for-prscs)
  - [Calculate the effect size of GWAS for case-control study](#calculate-the-effect-size-of-gwas-for-case-control-study)
  - [Run PRScs](#run-prscs)

---

# Calculate PRS using pre-trained weights

## Download Summary Statistics from PGSCatalog

First, we look for the available pre-trained PRS weights from PGSCatalog, see <https://www.pgscatalog.org/>.
Note that to use this, we need to make sure that there is no overlap between our data and the original GWAS, as well as the training sample in PGSCatalog. For downstream analyses, sample relatedness is another issue. A high degree of relatedness between individuals between the base and target data can also generate inflation of the association between the PRS and target phenotype. 


Pay attention to the Genome build (HG19 or HG38) of your own data. As an example, let's use pre-trained PRS weights for autism spectrum disorder: <https://www.pgscatalog.org/score/PGS000327/>. Ensure that the base and target data SNPs have genomic positions assigned on the same genome build. LiftOver is an excellent tool for standardizing genome build across different data sets.


## Data cleaning

### Ambiguous SNPs
Here, we use the HGDP project data. Our data version is HG38. If the base and target data were generated using different genotyping chips, and the chromosome strand (+/−) that was used for either is unknown, then it is not possible to pair up the alleles of ambiguous SNPs (i.e., those with complementary alleles, either C/G or A/T SNPs) across the data sets, because it will be unknown whether the base and target data are referring to the same allele or not. While allele frequencies could be used to infer which alleles are on the same strand, the accuracy of this could be low for SNPs with MAF close to $0.5$ or when the base and target data are from different populations. Therefore, we recommend removing all ambiguous/palindromic SNPs to avoid introducing this potential source of systematic error when populations do not match between the training GWAS and your target data, or when the effect allele frequency is not provided in the downloaded score file.
```{r,eval=F}
library(data.table)
in.dir = "/standard/cobi8302/PRS_calculation/"
out.dir = "your-own-directory"  #Please set up your own directory for analysis results
pgs_weights <- fread(paste0(in.dir,"PGS000327_hmPOS_GRCh38.txt.gz"), header = T)
head(pgs_weights)

# Change all alleles to upper case for easy comparison
pgs_weights$effect_allele <- toupper(pgs_weights$effect_allele)
pgs_weights$other_allele <- toupper(pgs_weights$other_allele)

#remove ambiguous snps due to the effect allele frequencies not provided in the score file
rm_id=which( (pgs_weights$effect_allele=="A" & pgs_weights$other_allele=="T") |
               (pgs_weights$effect_allele=="T" & pgs_weights$other_allele=="A") |
               (pgs_weights$effect_allele=="C" & pgs_weights$other_allele=="G") |
               (pgs_weights$effect_allele=="G" & pgs_weights$other_allele=="C") )
length(rm_id)
pgs_weights=pgs_weights[-rm_id,]
#Note that this score file has some NA values in the HG38 positions, we remove those entries
pgs_weights=pgs_weights[-which(is.na(pgs_weights$hm_pos)),]

```

### Standard GWAS QC (will do this step when calculating the PRS)
We recommend the following QC criteria for standard analyses: genotyping rate $>0.99$, sample missingness $<0.02$, Hardy-Weinberg Equilibrium P $>1 × 10^{−6}$, heterozygosity within 3 standard deviations of the mean, minor allele frequency (MAF) $>0.01$ (MAF $>0.05$ if target sample $N <1000$) and imputation ‘info score’ $>0.8$. If both the base and target data are large (e.g., $N >50\,000$), then SNPs with MAF $<0.01$ may be included, in which case we recommend a minor allele count $>100$ in the base and target data to ensure the integrity of normality assumptions implicit in association testing and LD calculation. 


### Mismatching SNPs and/or Duplicate SNPs
SNPs that have mismatching alleles reported in the base and target data are either resolvable by strand-flipping the alleles to their complementary alleles in, for example, the target data, such as for a SNP with A/C in the base data and G/T in the target, or non-resolvable, such as for a SNP with C/G in the base and C/T in the target. Most polygenic score software programs perform strand-flipping automatically for SNPs that are resolvable and remove non-resolvable mismatching SNPs. Ensure that there are no duplicated SNPs in either the base or target data.
```{R,eval=F}
#Identify SNPs that require strand flipping
bim <- read.table(paste0(in.dir,"data/merged.bim"))
colnames(bim) <- c("hm_chr", "snp", "CM", "hm_pos", "B.A1", "B.A2")

# Change all alleles to upper case for easy comparison
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)

# Merge summary statistic with target
library(dplyr)
snp_clean <- pgs_weights %>% 
  distinct(.keep_all = TRUE)

info <- merge(bim, pgs_weights, by = c("hm_chr", "hm_pos"))

id_keep=which( ((info$B.A1 == info$effect_allele) | (info$B.A1 == info$other_allele)) & ((info$B.A2 == info$effect_allele) | (info$B.A2 == info$other_allele)) )
length(id_keep)

#let's do some strand flipping
#Function for finding the complementary allele
complement <- function(x) {
  switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  )
}

# Identify SNPs that are complementary between base and target
info$effect_allele.complement <- sapply(info$effect_allele, complement)
info$other_allele.complement <- sapply(info$other_allele, complement)
id_complement = which( ((info$B.A1 == info$effect_allele.complement) | (info$B.A1 == info$other_allele.complement)) & ((info$B.A2 == info$effect_allele.complement) | (info$B.A2 == info$other_allele.complement)) )
# Update the complementary alleles in the info file
# This allow us to match the allele in subsequent analysis
info$effect_allele[id_complement] <- info$effect_allele.complement[id_complement]
info$other_allele[id_complement] <- info$other_allele.complement[id_complement]

info_final = info[c(id_keep,id_complement),]
info_score=info_final[,c("snp","effect_allele","effect_weight")]
dim(info_score)

```


### Save the score file for PLINK
```{R,eval=F}
#save this file for calculating scores in plink
write.table(info_score,file=paste0(out.dir,"asd_score.txt"),row.names = F,col.names = F,quote = F)
```
Other software tools for data harmonization: [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/articles/introduction.html)

### Additional note when palindromic SNPs are able to be sorted out
Suppose the EAF information is available and the training and target populations match, we can sort out and include palindromic SNPs in our PRS calculations
```
#Additional note: Suppose we have the EAF information in the PGS weights file for palindromic SNPs
info <- merge(bim, pgs_weights, by = c("hm_chr", "hm_pos"))
# Identify palindromic allele pairs
info$ambiguous <- (info$effect_allele == "A" & info$other_allele == "T") | 
  (info$effect_allele == "T" & info$other_allele == "A") |
  (info$effect_allele == "C" & info$other_allele == "G") | 
  (info$effect_allele == "G" & info$other_allele == "C")

# 2. Check that alleles belong to the same pair (prevents multi-allelic A/T vs C/G collisions)
same_pair <- ((info$B.A1 == info$effect_allele) | (info$B.A1 == info$other_allele)) & 
  ((info$B.A2 == info$effect_allele) | (info$B.A2 == info$other_allele))

# 3. Palindromic SNPs: frequency check + allele pair check, EAF is loaded from the PGS weight file downloaded online, ALT_FREQ is loaded from allele frequency counts for B.A1
id_match <- which(info$ambiguous & same_pair & 
                    ((info$EAF < 0.35 & info$ALT_FREQ < 0.35) | (info$EAF > 0.65 & info$ALT_FREQ > 0.65)))

id_flip  <- which(info$ambiguous & same_pair & 
                    ((info$EAF < 0.35 & info$ALT_FREQ > 0.65) | (info$EAF > 0.65 & info$ALT_FREQ < 0.35)))

# Harmonize palindromic effect_allele to match genotype file coding
info$effect_allele[id_match] <- info$B.A1[id_match]
info$effect_allele[id_flip]  <- info$B.A2[id_flip]

# 4. Non-ambiguous SNPs: direct match vs. strand complement
id_keep <- which(!info$ambiguous & same_pair)

id_complement <- which(!info$ambiguous & 
                         ((info$B.A1 == info$effect_allele.complement) | (info$B.A1 == info$other_allele.complement)) & 
                         ((info$B.A2 == info$effect_allele.complement) | (info$B.A2 == info$other_allele.complement)))

# Update complementary alleles
info$effect_allele[id_complement] <- info$effect_allele.complement[id_complement]
info$other_allele[id_complement]  <- info$other_allele.complement[id_complement]

# 5. final output (sorted to preserve chromosome order)
keep_idx <- sort(unique(c(id_keep, id_complement, id_match, id_flip)))
info_final <- info[keep_idx, ]
info_score=info_final[,c("snp","effect_allele","effect_weight")]
dim(info_score)
```


## Calculation of PRS
This will be done in PLINK2.0
```{R, echo=T,eval=F}
####################
/standard/cobi8302/PRS_calculation/tools/plink2 \
--threads 5 \
--rm-dup exclude-all \
--mach-r2-filter \
--maf 0.01 \
--hwe 1e-6 \
--geno 0.01 \
--mind 0.05 \ 
--bfile /standard/cobi8302/PRS_calculation/data/merged  \
--score ${RESULTS_DIR}/asd_score.txt 1 2 3 list-variants cols=+scoresums,-scoreavgs \
--out ${RESULTS_DIR}/ASD_score_final

```
Please set up your own directory for results, need to fill in ${RESULTS_DIR}. Note that --mind 0.05 is a relatively lenient threshold, as we are only using 663 individuals, which is a small set of data

## Visualize the PRS
```{R}
library(data.table)
prs_asd <- fread(paste0(out.dir,"ASD_score_final.sscore"), header = T)
hist(prs_asd$SCORE1_SUM, col = "blue", breaks = 20, xlab = "PRS value", main = "Histogram of Autism Spectrum Disorder PRS in 1000G+HGDP EUR", freq = F)

```


# Train PRS models using GWAS summary statistics
First, download GWAS summary statistics from [GWASCatalog](https://www.ebi.ac.uk/gwas/). We use the example dataset for ADHD (available for download [download](https://livejohnshopkins-my.sharepoint.com/:t:/g/personal/zwang389_jh_edu/EaNT0Lh2rj9FjIsPmjlEhbYBr58WxfokseOIP9uqegZu3g?e=WaJBfi)).
Second, install PRScs into your computer or computing server [Instructions](https://github.com/getian107/PRScs).

## Prepare the summary statistics for PRScs
```{R,eval = F}
library(data.table)
adhd = fread("ADHD2022_iPSYCH_deCODE_PGC.meta.gz") #This was downloaded from GWAS Catalog

#create input for PRScs
sum_stats = adhd[,c("SNP","A1","A2","OR","P")]
write.table(sum_stats,file="sum_stats.txt",col.names = T,row.names = F,quote = F)
```

## Calculate the effect size of GWAS for case-control study
https://rdrr.io/github/neurogenomics/MungeSumstats/man/compute_sample_size.html

4/(1/N_case+1/N_control)

In the original ADHD paper, N_case = 38691, N_control = 186843
```{R,eval = F}
4/(1/38691+1/186843)
#128213.8
```

## Run PRScs
First download the reference LD panels. Specify your paths of the data, genotype data, and output file.
```{R, eval=F}
python PRScs.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --out_dir=OUTPUT_DIR [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --chrom=CHROM --beta_std=BETA_STD --write_psi=WRITE_PSI --write_pst=WRITE_POSTERIOR_SAMPLES --seed=SEED]


python PRScs.py --ref_dir=/standard/cobi8302/PRS_calculation/PRScs/data/ldblk_ukbb_eur --bim_prefix=/standard/cobi8302/PRS_calculation/PRScs/data/mergedplink --sst_file=/standard/cobi8302/PRS_calculation/PRScs/data/sum_stats.txt --n_gwas=128214 --out_dir=/standard/cobi8302/PRS_calculation/PRScs/results/eur
```



## Then prepare data for input in PLINK
Combine the scoring files from 22 chromosomes to 1 file, save to rsid, effect allele, effect size.
```{R, eval=F}
score_file = fread("/standard/cobi8302/PRS_calculation/PRScs/results/eur_pst_eff_a1_b0.5_phiauto_chr1.txt",header = F)
score_file = data.frame(score_file)
score_file = score_file[,c(2,4,6)]
for(i in 2:22){
  tmp = fread(paste0("/standard/cobi8302/PRS_calculation/PRScs/results/eur_pst_eff_a1_b0.5_phiauto_chr",i,".txt"),header = F)
  tmp = data.frame(tmp)
  tmp = tmp[,c(2,4,6)]
  score_file = rbind(score_file,tmp)
}
write.table(score_file,file="/standard/cobi8302/PRS_calculation/PRScs/results/info_score.txt",col.names = F,row.names = F,quote = F)

```

Finally, use the same steps in previous PLINK code to calculate PRS for ADHD.

Note that for tuning parameters, you will need to split the target data into validation set and then assess the performance of PRS in testing set (example, 1:1 validation/testing or 1:2 validation/testing). Here, we directly used the default parameters provided by PRScs.

# References

Choi, S.W., Mak, T.SH. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc 15, 2759–2772 (2020). https://doi.org/10.1038/s41596-020-0353-1
