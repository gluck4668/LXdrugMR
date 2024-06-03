
# The gwas data ID were obtained from https://gwas.mrcieu.ac.uk/datasets/

exp_gwas_id <- function(tart_gene,chr_exp,pos_exp_start,pos_exp_end,
                        exp_dat,exp_p,clump_r2,clump_kb){


print("Begin extract_instruments...")
exp_df <- extract_instruments(outcomes=exp_dat,clump=F)
exp_df <- subset(exp_df,pval.exposure<exp_p)

if(nrow(exp_df)>0)
  print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_df), " SNPs.")) else
    stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))

#-------Linkage Disequilibrium (LD) test--------
print("Begin the linkage disequilibrium (LD) test......")
ld_clum <- exp_ld_clum(exp_df,clump_kb,clump_r2)

#----------------------------------
exp_df <- ld_clum$exp_clum
target_gene_snp <- subset(exp_df,chr.exposure==chr_exp & pos.exposure>pos_exp_start-100000 & pos.exposure<pos_exp_end+100000)

exp_file <- paste0("analysis results/",tart_gene," snp.csv")
write.csv(target_gene_snp,exp_file)

exp_data <- target_gene_snp
exp_snp <- data.frame(target_gene_snp$SNP)
names(exp_snp) <- "SNP"

result <- list(exp_data=exp_data,exp_snp=exp_snp)
return(result)

}
