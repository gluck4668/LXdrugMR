
# VCF data were downloaed from https://gwas.mrcieu.ac.uk/datasets/
exp_vcf_data <- function(tart_gene,chr_exp,pos_exp_start,pos_exp_end,exp_dat,exp_p){

exp_df <- VariantAnnotation::readVcf(exp_dat)  # 读取vcf数据
exp_df <- gwasglue::gwasvcf_to_TwoSampleMR(vcf =exp_df,type = "exposure") #格式化数据
exp_df <- distinct(exp_df,SNP,.keep_all = T) %>%
            subset(pval.exposure<exp_p) %>%
            subset(grepl("rs",exp_df$SNP))

if(nrow(exp_df)>0)
  print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_df), " SNPs.")) else
    stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))

#-------Linkage Disequilibrium (LD) test--------
ld_clum <- exp_ld_clum(exp_df,clump_kb,clump_r2)

#------------------------------
exp_df <- ld_clum$exp_clum
target_gene_snp <- subset(exp_df,chr.exposure==chr_exp & pos.exposure>pos_exp_start-100000 & pos.exposure<pos_exp_end+100000)

exp_data <- target_gene_snp
exp_file <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),".csv")
write.csv(exp_data,exp_file)

exp_snp <- data.frame(exp_data$SNP)
names(exp_snp) <- "SNP"

result <- list(exp_data=exp_data,exp_snp=exp_snp)
return(result)

}
