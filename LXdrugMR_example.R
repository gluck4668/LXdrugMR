
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXdrugMR")

library(LXdrugMR)


#---Drug target MR -----------------------------------------------------------

# 暴露(exp_dat)和结局(out_dat)数据类型可以是 GWAS ID (如ieu-a-300,ebi-a-GCST003156),
# 也可以是本地GWAS Summary文件(格式包括gz;vcf.gz;csv;txt)
# 暴露和结局数据可以来自不同的数据库

rm(list=ls()) #*每次都要运行，以清除内存中旧的数据

{
#---target gene---
tart_gene ="HMGCR"
chr_exp= 5
pos_exp_start = 74632993
pos_exp_end = 74657941


#---exposure data---
exp_dat="ieu-a-300" # 暴露，必填

exp_p=1e-6        #* 关联性，必填
clump_kb = 100  #* 连锁不平衡kb, 必填
clump_r2 = 0.3  #* 连锁不平衡r2,必填

exposure_name="sphingolipid" # 暴露，一般是疾病名称，选填

# 以下6项，如果GWAS data不是来自https://gwas.mrcieu.ac.uk/datasets/，需要提供。如果是，则不需要填
beta_exp = "beta" #如果有beta,用beta;如果没有，则用OR或者logOR
se_exp = "se"
effect_allele_exp = "alt"
other_allele_exp = "ref"
eaf_exp = "AF"
pval_exp = "pval"


#---outcome data---
out_dat="ebi-a-GCST003156" #* 结局，必填

outcome_name="NAFLD"  # 结局，一般是疾病名称，选填

# 以下6项，如果GWAS data不是来自https://gwas.mrcieu.ac.uk/datasets/，需要提供。如果是，则不需要填
beta_out = "beta" #如果有beta,用beta;如果没有，则用OR或者logOR
se_out = "sebeta"
effect_allele_out = "alt"
other_allele_out = "ref"
eaf_out = "af_alt"
pval_out = "pval"
}

#---run the main function---
devtools::load_all()

LXdrugMR()


#---The end--------------------------------------------------------------------






