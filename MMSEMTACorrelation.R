data <- read.csv("E:/VANS/data/total_mmse_genecounts_normalized.csv")

coeff <- c()
pv <- c()
mmse <- as.numeric(data[1,2:length(data)])

for (i in 2:length(data$X))
{
  test <- cor.test(mmse, as.numeric(data[i,2:length(data)]), 
                   method="spearman", 
                   alternative="two.sided")
  coeff <- append(coeff, test$estimate)
  pv <- append(pv, test$p.value)
}

df<-data.frame(data$X[-1],coeff,abs(coeff),pv)
##take index
idx <- which(df$pv < 0.05)
df_sub <- df[c(idx),]

ord <- order(df_sub$abs.coeff.)
df_plot <- df_sub[ord,]
df_plot_ <- setNames(df_plot, 
                         c("gene",
                           "coeff_rho",
                           "abs_coeff_rho",
                           "p_value"
                           ))
row.names(df_plot_) <- NULL #reset index

#idx_plot <- as.numeric(tail(rownames(df_plot))[-1]) + 1 

write.csv(df_plot_, file='E:/VANS/data/total_gene_corr_mmse.csv')

library("ggpubr")


