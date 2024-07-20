####储存harmonise_data分析后的dat结果
setwd("F:\\肠道菌群MR\\细胞因子 and 神经性疾病\\Mic and 细胞因子 and ADHD(多动症)\\肠道菌群与ADHD\\dat_harmo数据")
all_data <- all_data[c(22, 40, 50, 51, 77, 96, 137, 138, 140, 181, 191)] ##提取阳性肠道菌群进行下一步敏感性分析
dat_results <- list()  # 用于存储harmonise_data后的dat结果

for (i in c(1:length(all_data))) {
  # 对应all_data[i]中的SNP进行筛选
  GWAS_2 <- subset(GWAS_1, GWAS_1$SNP %in% all_data[[i]]$SNP)
  
  # 添加表型
  GWAS_2$PHENO <- "ADHD"
  GWAS_2$beta <- log(GWAS_2$OR)
  
  # 格式转换为outcome
  out_data <- TwoSampleMR::format_data(
    GWAS_2,
    type = "outcome",
    phenotype_col = "PHENO",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "A2",
    other_allele_col = "A1",
    pval_col = "P"
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
    
    # 储存调和后的数据
    dat_results[[length(out_dat)*(i-1)+j]] <- dat
    
    # 如果仍需计算并打印 res 的结果，可以保留以下代码
    # res <- mr(dat)
    
    # res$exposure = names(all_data)[i]
    # res$outcome = names(out_dat)[j]
    # print(paste0("......", names(all_data)[i], "&", names(out_dat)[j], "....."))
    # print(generate_odds_ratios(res))
    
    # results[[length(out_dat)*(i-1)+j]] <- generate_odds_ratios(res)
  }
}

# 返回或处理 dat_results
dat_results

save(dat_results, file="dat_results_ADHD.Rdata")
results_dat_all <- do.call(rbind, dat_results)
save(results_dat_all, file="results_dat_ADHD.Rdata")
write.csv(results_dat_all, file="results_dat_ADHD.csv")
