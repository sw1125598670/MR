set.seed(123)
source("Get_MR2.0.r")

# ---- 01 outcome data ----
outcomefile_fp="iPSYCH-PGC_ASD_Nov2017.gz"
outcomeINFO=read.table(outcomefile_fp,header=TRUE)
outcome_exp_dat_origin<-data.frame(
    SNP=outcomeINFO$SNP,
    chr_col=outcomeINFO$CHR,
    snp_col=outcomeINFO$SNP,
    beta_col=log(outcomeINFO$OR),
    se_col=outcomeINFO$SE,
    effect_allele_col=outcomeINFO$A1,
    other_allele_col=outcomeINFO$A2,
    pval_col = outcomeINFO$P
)

outcome_exp_dat<-get_eaf_from_1000G(outcome_exp_dat_origin,"1kg.v3", type = "outcome")
outcome_exp_dat <- format_data(
 outcome_exp_dat_origin,
 type='outcome',
 snp_col = "snp_col",
 beta_col = "beta_col",
 se_col = "se_col",
 effect_allele_col ="effect_allele_col",
 other_allele_col = "other_allele_col",
 eaf_col = "eaf.outcome",
 pval_col = "pval_col"

)

# ---- 02 exposure data ----
exposurefileDIR="MiBioGen"
exposurefiles=dir(exposurefileDIR)

for(exposurefile in exposurefiles){
exposurefile_fp=paste(exposurefileDIR,exposurefile,sep="/")

MiBioID=sub("\\.summary.*$", "", exposurefile)
p_fn=paste(MiBioID,'mr.pdf',sep="_")
p_fp=paste("results",p_fn,sep="/")

exposureINFO=fread(exposurefile_fp)

exposure_exp_dat_origin <-data.frame(
    SNP=exposureINFO$rsID,
    chr_col=exposureINFO$chr,
    snp_col=exposureINFO$rsID,
    beta_col=exposureINFO$beta,
    se_col=exposureINFO$SE,
    effect_allele_col=exposureINFO$eff.allele,
    other_allele_col=exposureINFO$ref.allele,
    pval_col=exposureINFO$P.weightedSumZ,
    pos_col=exposureINFO$bp,
    samplesize_col=exposureINFO$N
)
exposure_exp_dat_origin=get_eaf_from_1000G(exposure_exp_dat_origin,"1kg.v3",type="exposure")

exposure_exp_dat <- format_data(
 exposure_exp_dat_origin,
 type='exposure',
 snp_col = "snp_col",
 beta_col = "beta_col",
 se_col = "se_col",
 effect_allele_col ="effect_allele_col",
 other_allele_col = "other_allele_col",
 eaf_col = "eaf.exposure",
 pval_col = "pval_col",
 samplesize_col="samplesize_col"
)

##  filter snp ---- 
exposure_exp_dat=exposure_exp_dat[exposure_exp_dat$pval.exposure<1e-5,]

exposure_exp_dat=clump_data(exposure_exp_dat,clump_r2=0.01,pop = "EUR",
                            plink_bin="~/.local/bin/plink",
                            bfile="1kg.v3/EUR")

exposure_exp_dat=get_f(exposure_exp_dat,F_value=10)

##  harmonise data ----
dat <- harmonise_data(
    exposure_dat = exposure_exp_dat, 
    outcome_dat = outcome_exp_dat)
    
## MR ----
res <- mr(dat)
res_fn=paste(MiBioID,"mr_res.csv",sep="_")
res_fp=paste("results",res_fn ,sep = "/")
write.csv(res, file = res_fp, quote = TRUE)

##  mr_heterogeneity
mr_hetero_res=mr_heterogeneity(dat)
mr_hetero_res_fn=paste(MiBioID,"mr_heterogeneity_res.csv",sep="_")
mr_hetero_res_fp=paste("results",mr_hetero_res_fn ,sep = "/")
write.csv(mr_hetero_res, file = mr_hetero_res_fp, quote = TRUE)

##  mr_pleiotropy_test ----

mr_pleio_res=mr_pleiotropy_test(dat)
mr_pleio_res_fn=paste(MiBioID,"mr_pleiotropy_res.csv",sep="_")
mr_pleio_res_fp=paste("results",mr_pleio_res_fn ,sep = "/")
write.csv(mr_pleio_res, file = mr_pleio_res_fp, quote = TRUE)

## single snp ----
res_single <- mr_singlesnp(dat)
res_single_fn=paste(MiBioID,"mr_res_single_snp.csv",sep="_")
res_single_fp=paste("results",res_single_fn ,sep = "/")
write.csv(res_single, file = res_single_fp, quote = TRUE)

## leave one out ----
res_loo <- mr_leaveoneout(dat)
res_loo_fn=paste(MiBioID,"mr_res_loo.csv",sep="_")
res_loo_fp=paste("results",res_loo_fn ,sep = "/")
write.csv(res_loo, file = res_loo_fp, quote = TRUE)

## mpresso out ----
res_mpresso <- try(run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05))
res_mpresso_fn=paste(MiBioID,"mr_res_mpresso.txt",sep="_")
res_mpresso_fp=paste("results",res_mpresso_fn ,sep = "/")
write.csv(res_mpresso, file = res_mpresso_fp, quote = TRUE)

pdf(paste("results",p_fn,sep="/"))
### mr_scatter_plot ----
p1 <-mr_scatter_plot(res, dat)
print(p1)
### mr_singlesnp----
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
print(p2)
### leave-one-out
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
print(p3)

p4 <- mr_funnel_plot(res_single)
print(p4)
dev.off()
}
