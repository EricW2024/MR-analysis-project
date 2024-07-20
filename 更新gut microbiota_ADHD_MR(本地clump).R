#library
library(TwoSampleMR)
library(ggplot2)
library(data.table)
# 
# install.packages("pacman", update=F, ask=F)
library(pacman)
p_load(data.table, dplyr, tidyr)
# remotes::install_github("yulab-smu/yulab.utils")
library(yulab.utils)
# install_zip_gh("MRCIEU/ieugwasr")
library(ieugwasr)
rm(list=ls())
#set up your work directory
setwd("F:\\肠道菌群MR\\细胞因子 and 神经性疾病\\Mic and 细胞因子 and ADHD(多动症)\\肠道菌群与ADHD")
Mic_total <- read.csv("肠道微生物_tophits_2024-06-05_10_08_57.csv", header = T, sep=",")  ###211种肠道菌群;已过滤p < 1e-5
dim(Mic_total)
colnames(Mic_total)

library(TwoSampleMR)
library(data.table)

# 初始化一个空列表用于存储所有数据
all_data <- list()

# 获取 Mic_total 数据表中所有唯一的细菌名
unique_bac <- unique(Mic_total$bac)

# 循环遍历每个细菌名，并将其子集分配给变量 M1 到 M211
for (i in seq_along(unique_bac)) {
  # 获取当前细菌名
  bac_name <- unique_bac[i]
  
  # 动态创建变量名
  var_name <- paste0("M", i)
  
  # 使用子集数据框
  data <- subset(Mic_total, bac == bac_name & P < 1e-5)
  
  # 检查筛选后的数据是否为空
  if (nrow(data) == 0) {
    next  # 如果筛选后数据为空，跳过当前细菌名
  }
  
  # 使用 TwoSampleMR::format_data 函数格式化数据
  formatted_data <- TwoSampleMR::format_data(
    data,
    phenotype_col = "bac",
    type = "exposure",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "eff.allele",
    other_allele_col = "re.allele",
    pval_col = "P",
    samplesize_col = "N"
  )
  
  # 本地 clump
  biof_iv <- formatted_data[, c("SNP", "pval.exposure")]
  
  # 将 biof_iv 的数据框中的列名修改为 rsid 和 pval
  colnames(biof_iv) <- c("rsid", "pval")
  
  # 打印调试信息
  print(paste0("Processing ", bac_name, " with ", nrow(formatted_data), " rows after format_data"))
  
  # 检查 biof_iv 数据
  print(head(biof_iv))
  
  # 使用本地 clump 函数
  clump_dat <- ld_clump_local(
    dat = biof_iv,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p = 1,
    bfile = "F:/肠道菌群MR/细胞因子 and 自闭症/肠道菌群 and 细胞因子 and 自闭症/细胞因子 and 自闭症/本地clump包/data_maf0.01_rs_ref/data_maf0.01_rs_ref",
    plink_bin = "F:/肠道菌群MR/细胞因子 and 自闭症/肠道菌群 and 细胞因子 and 自闭症/细胞因子 and 自闭症/本地clump包/plink_win64_20231211/plink.exe"
  )
  
  # 检查 clump_dat 数据
  print(paste0("Clumped data for ", bac_name, ": ", nrow(clump_dat), " rows"))
  print(head(clump_dat))
  
  # 保留 formatted_data 中在 clump_dat$rsid 中的 SNP
  formatted_data <- formatted_data[formatted_data$SNP %in% clump_dat$rsid, ]
  
  # 将筛选后的数据添加到列表中
  all_data[[i]] <- formatted_data
  
  # 打印调试信息
  print(paste0("Processed ", bac_name, " with ", nrow(formatted_data), " rows after clumping"))
}


save(all_data, file="all_data_Mic211_clump.Rdata")
# 将所有列表中的数据框合并成一个数据框
combined_data <- bind_rows(all_data)
save(combined_data, file="combined_data_Mic211_clump.Rdata")
# 保存合并后的数据框为文件
write.table(combined_data, "combined_expose_Mic211_clump.txt", row.names = FALSE, sep = "\t", quote = FALSE)

setwd("F:\\肠道菌群MR\\细胞因子 and 神经性疾病\\Mic and 细胞因子 and ADHD(多动症)\\肠道菌群与ADHD")
GWAS_1 <- fread("F:/肠道菌群MR/细胞因子 and 自闭症/Mic and 细胞因子 and ADHD(多动症)/源数据/ADHD2022_iPSYCH_deCODE_PGC.meta.gz", header = T) ##读取疾病GWAS数据；其他疾病可直接替换数据
save(GWAS_1 , file="ADHD.Rdata")
load("ADHD.Rdata")
dim(GWAS_1)
GWAS_1 <- as.data.frame(GWAS_1)
# GWAS_1$beta <- log(GWAS_1$OR)  ##如果部分疾病只有OR；则通过该步转换成beta

results <- list()

for (i in c(1:length(all_data))) {
  # 对应all_data[i]中的SNP进行筛选
  GWAS_2 <- subset(GWAS_1, GWAS_1$SNP %in% all_data[[i]]$SNP)   ##merge出共同的交集SNP
  
  # 添加表型
  GWAS_2$PHENO <- "ADHD"
  GWAS_2$beta <- log(GWAS_2$OR)   ##如果部分疾病只有OR；则通过该步转换成beta
  
  # 格式转换为outcome
  out_data <- TwoSampleMR::format_data(
    GWAS_2,
    type = "outcome",
    phenotype_col = "PHENO",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "A2", ##注意观察原文的effet；该原始数据中A2是指的effect allele
    other_allele_col = "A1",
    pval_col = "P",
    eaf_col = "FRQ_A_38691"
  )
  
  out_dat <- list()
  out_dat[[1]] <- out_data
  
  # 定义outcome表型向量
  names(out_dat) <- c("ADHD")
  
  for (j in c(1:length(out_dat))) {
    # 调和数据
    dat <- harmonise_data(
      exposure_dat = all_data[[i]],
      outcome_dat = out_dat[[j]],
      action = 1   ### 注意：暴露和结局没有eaf时用action=1
    )
    
    res <- mr(dat)
    
    res$exposure = names(all_data)[i]
    res$outcome = names(out_dat)[j]
    print(paste0("......", names(all_data)[i], "&", names(out_dat)[j], "....."))
    print(generate_odds_ratios(res))
    
    results[[length(out_dat)*(i-1)+j]] <- generate_odds_ratios(res)
  }
}

## 删除GWAS_1释放内存
rm(GWAS_1)


# ################更新代码###################################
# results <- list()
# 
# for (i in c(1:length(all_data))) {
#   for (j in c(1:length(out_dat))) {
#     
#     #调和数据
#     dat <- harmonise_data(
#       exposure_dat = all_data[[i]],
#       outcome_dat = out_dat[[j]],
#       action = 1
#     )
#     res <- mr(dat)
#     
#     # Extract the names of exposure and outcome for printing
#     exposure_name <- names(all_data)[i]
#     outcome_name <- names(out_dat)[j]
#     
#     print(paste0("...... ", exposure_name, " & ", outcome_name, " ....."))
#     
#     odds_ratios <- generate_odds_ratios(res)
#     print(odds_ratios)
#     
#     # Add exposure and outcome names to the results
#     res$exposure <- exposure_name
#     res$outcome <- outcome_name
#     
#     # Store the results in the list
#     results[[length(out_dat) * (i - 1) + j]] <- odds_ratios
#   }
# }
##########################################################################################
save(results, file="MR_Mic211_ADHD.Rdata")
#绑定所有结果
results_allIV <- do.call(rbind, results)
save(results_allIV, file="results_Mic211_ADHD_clump.Rdata")

write.csv(results_allIV, file="result_Mic211_ADHD_clump.csv")


##格式化估计值列
results_allIV$estimate <- paste0(round(results_allIV$or, 2), "(", 
                                 round(results_allIV$or_lci95, 2), "-", 
                                 round(results_allIV$or_uci95, 2), ")" )

##提取p值大于0.05的行号
row_x <- rownames(results_allIV[which(results_allIV$pval > 0.05), ])

##格式化p值列
results_allIV$pvalue <- format(results_allIV$pval, scientific = TRUE, digits = 2)

######输出##############
#输出文件路径
write.table(results_allIV[, c(1:ncol(results_allIV))],
            file="results_MR_Mic211_ADHD_clump.csv",
            sep=",",
            row.names = FALSE,
            col.names = TRUE, 
            quote = TRUE)





