
exp_other_format <- function(tart_gene,chr_exp,pos_exp_start,pos_exp_end,
                             exp_dat,snp_exp,pval_exp,beta_exp,
                             se_exp,effect_allele_exp,other_allele_exp,
                             eaf_exp,clump_kb,clump_r2){

#-----read data file-----
exp_file_type <- str_extract(exp_dat,"(?<=\\.)[^\\.]+$")

if(tolower(exp_file_type)=="gz")
  exp_df <- fread(exp_dat) else
      {if(tolower(exp_file_type)=="txt")
        exp_df <- read_table(exp_dat) else
          exp_df <- eval(str2expression(paste0("read.",exp_file_type,"(exp_dat)")))
        }

#---------------
exp_10 <- head(exp_df)

 for (i in 1:ncol(exp_10)) {
      col_val <- eval(str2expression(paste0("exp_10[,",i,"]")))
      if(any(grepl("rs",col_val))){
        snp_exp_sit <- i
        snp_exp <- names(exp_10)[i]}
    }

exp_df <- eval(str2expression(paste0("distinct(exp_df,",snp_exp,",.keep_all = T)")))
exp_df <-eval(str2expression(paste0("subset(exp_df,",pval_exp,"<",exp_p,")")))

if(nrow(exp_df)>0)
  print(paste0("Note: At the condition of exp_p<",exp_p,", there are ",nrow(exp_df), " SNPs.")) else
  stop(paste0("Note: At the condition of exp_p<",exp_p,", there is no SNP."))

exp_file <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),".csv")
write.csv(exp_df,exp_file)

colnames(exp_df)[snp_exp_sit] <- "SNP"
exp_file02 <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),"02.csv")
write.csv(exp_df,exp_file02)

chr_sit <- grep("chr",colnames(exp_df),ignore.case = T)
pos_sit <- grep("pos",colnames(exp_df),ignore.case = T)
SNP_sit <- grep("SNP",colnames(exp_df),ignore.case = T)
df_chr <- c(chr_sit,pos_sit,SNP_sit) %>% as.numeric()
paste0(df_chr,collapse = ",")
exp_df_chr <- eval(str2expression(paste0("exp_df[,","c(",paste0(df_chr,collapse = ","),")","]")))
colnames(exp_df_chr) <- c("chr.exposure","pos.exposure","SNP")

exp_df <- read_exposure_data(filename = exp_file02,
                                   sep = ",",
                                   snp_col = "SNP",
                                   beta_col = beta_exp,
                                   se_col = se_exp,
                                   effect_allele_col = effect_allele_exp,
                                   other_allele_col = other_allele_exp,
                                   eaf_col = eaf_exp,
                                   pval_col = pval_exp)

exp_df <- subset(exp_df,grepl("rs",exp_df$SNP))

if(file.exists(exp_file02))
  file.remove(exp_file02)

if(trimws(toupper(beta_exp)) == "OR")
exp_df$beta.exposure <- log(exp_df$beta.exposure)

#-------Linkage Disequilibrium (LD) test--------
print("Begin the linkage disequilibrium (LD) test......")
ld_clum <- exp_ld_clum(exp_df,clump_kb,clump_r2)

#------------------------------
exp_df <- ld_clum$exp_clum

if(!"chr.exposure" %in% colnames(exp_df))
  exp_df <- inner_join(exp_df_chr[,c("chr.exposure","SNP")],exp_df,"SNP")

if(!"pos.exposure" %in% colnames(exp_df))
  exp_df <- inner_join(exp_df_chr[,c("pos.exposure","SNP")],exp_df,"SNP")

target_gene_snp <- subset(exp_df,chr.exposure==chr_exp & pos.exposure>pos_exp_start-100000 & pos.exposure<pos_exp_end+100000)

if(nrow(target_gene_snp)==0)
  stop("Note: there is no target gene SNP, please check the chr.exposure and pos.exposure in the exposure data file.")

exp_file <- paste0("analysis results/",tart_gene," snp.csv")
write.csv(target_gene_snp,exp_file)

exp_data <- target_gene_snp
exp_snp <- data.frame(target_gene_snp$SNP)
names(exp_snp) <- "SNP"

result <- list(exp_data=exp_data,exp_snp=exp_snp)
return(result)


}
