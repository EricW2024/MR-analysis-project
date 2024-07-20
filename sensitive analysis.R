rm(list=ls())
setwd("F:\\肠道菌群MR\\细胞因子 and 神经性疾病\\Mic and 细胞因子 and ADHD(多动症)\\肠道菌群与ADHD\\dat_harmo数据")
load("dat_results_ADHD.Rdata")
setwd("F:\\肠道菌群MR\\细胞因子 and 神经性疾病\\Mic and 细胞因子 and ADHD(多动症)\\肠道菌群与ADHD\\dat_harmo数据\\敏感性分析")

#heterogeneity test
harm_rt <- dat_results[[11]]
mr_heterogeneity(harm_rt)
mr_het <- mr_heterogeneity(harm_rt)
write.csv(mr_het, file="mr_heterogeneity_order.Coriobacteriales.id.810.csv")
#outlier test
run_mr_presso(harm_rt,NbDistribution = 5000)

#pleiotropy test
mr_pleiotropy_test(harm_rt)

#obtain the beta for each SNP
singlesnp_res <- mr_singlesnp(harm_rt)
# View(singlesnp_res)
singlesnpOR=generate_odds_ratios(singlesnp_res)
write.table(singlesnpOR,"singlesnpOR.txt",row.names = F,sep = "\t",quote = F)
mr_funnel_plot(singlesnp_res)
##留一法分析
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harm_rt))

p2 <- mr_forest_plot(singlesnp_res)
p2[[1]]

#sensitivity analysis
sen_res<- mr_leaveoneout(harm_rt)
View(sen_res)

#Scatter plots of several statistical methods
mr_result<- mr(harm_rt)
p1 <- mr_scatter_plot(mr_result, harm_rt)
p1[[1]]
ggsave(p1[[1]], file="scatter.pdf", width=8, height=8)
