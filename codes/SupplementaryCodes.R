################################################################################
#                                                                              #  
#                         Main Make Supplementary Figures                      #      
#                         Benjamin Hivert - 23/01/2026                         #
#                                                                              #
################################################################################

# This files contains codes to generate all the main figures of the paper. 
# All results were pre-imported and prepared in 20250417_PrepareResults.R script

# Library 

library(dplyr)
library(ggplot2)
library(patchwork)
library(latex2exp)
library(mvtnorm)
library(readxl)

theme_set(theme_classic())

#------------------------------------------------------------------------------#
#                 Sup. Figure X: QQ-Plot of pvalues when 
#                         correlated data generation            
#------------------------------------------------------------------------------#

neg_bin_cor_res <- read.csv("results/CorrelationAndRelativeBiais.csv")

qqplot_cor <- neg_bin_cor_res %>% 
  filter(Error == 1) %>%
  filter(Variable %in% c("Var1", "Var50")) %>%
  filter(Estimation == "Oracle") %>%
  mutate(VarName = ifelse(Variable == "Var1", "X[1]", "X[50]")) %>%
  mutate(MethodName = ifelse(Method == "NB", "Negative~Binomial~Thinning", "Gaussian~Fisson")) %>%
  ggplot() +
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pval, 
              colour = factor(Rho)),
          distribution = qunif, 
          size = 1,
          alpha = .5) +
  facet_grid(VarName~MethodName, labeller = label_parsed) +
  scale_colour_manual(name = TeX(r'(Correlation ($\rho$))'),
                      values = c(
                        "#FFE5B4",  # abricot clair
                        "#FFBFA0",  # corail doux
                        "#FF8CA0",  # rose saumon
                        "#C478B8",  # violet chaud
                        "#8050A0",  # violet saturé
                        "#50286F",  # prune foncé
                        "#2B1240"   # aubergine très foncé
                      )) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic() +
  theme(text = element_text(size = 14))

ggsave(filename = "Supplementary Figures/SuppFigure_QQPlotCor.pdf", 
       plot = qqplot_cor,
       width = 200,
       height = 150,
       units = "mm")

#------------------------------------------------------------------------------#
#                 Sup. Figure X: QQ-Plot of pvalues when 
#                   correlated data generation as a
#                       function of relative biais 
#------------------------------------------------------------------------------#

neg_bin_cor_res <- read.csv("results/CorrelationAndRelativeBiais.csv")

qqplot_cor_biais <- neg_bin_cor_res %>% 
  mutate(Theta_hat = Theta*Error) %>% 
  filter(Rho == 0) %>%
  filter(Variable %in% c("Var1", "Var50")) %>%
  filter(Estimation == "Wrong") %>%
  filter(Error %in% c(0.001, 0.5, 1, 1.5, 5)) %>%
  group_by(Error) %>%
  mutate(RelativeBiais = (Theta_hat-Theta)/Theta) %>% 
  mutate(VarName = ifelse(Variable == "Var1", "X[1]", "X[50]")) %>%
  mutate(MethodName = ifelse(Method == "NB", "Negative~Binomial~Thinning", "Gaussian~Fisson")) %>%
  ggplot() +
  geom_abline(slope=1, 
              intercept=0, 
              col="red", 
              size = 1.2, 
              alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pval, 
              colour = factor(RelativeBiais)),
          distribution = qunif, 
          size = 1, 
          alpha = .5) +
  facet_grid(VarName~MethodName, labeller = label_parsed) +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  guides(colour = guide_legend("Relative Biais")) +
  # scale_colour_manual(name = TeX(r'(Correlation ($\rho$))'),
  #                     values = c(
  #                       "#FFE5B4",  # abricot clair
  #                       "#FFBFA0",  # corail doux
  #                       "#FF8CA0",  # rose saumon
  #                       "#C478B8",  # violet chaud
  #                       "#8050A0",  # violet saturé
  #                       "#50286F",  # prune foncé
  #                       "#2B1240"   # aubergine très foncé
  #                     )) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic() +
  theme(text = element_text(size = 14))

ggsave(filename = "Supplementary Figures/SuppFigure_QQPlotCorBiais.pdf", 
       plot = qqplot_cor_biais,
       width = 200,
       height = 150,
       units = "mm")
#------------------------------------------------------------------------------#
#                               Sup. Figure X                                  #
#                         Applications on Bonne marrow                         #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                     Gene-specific overdispersion parameters                  #

# Parameter
cell_pop_to_test <- c("neutrophil", 
                      "macrophage",
                      "monocyte",
                      "granulocyte",
                      "CD4-positive, alpha-beta T cell",
                      "memory B cell")

cell_theta <- read.csv(file = "results/Application_CellPopulationOverdispersion.csv")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

pair_cellPop <- combn(2:(length(cell_pop_to_test)+1),2)
allPairplot <- lapply(1:ncol(pair_cellPop), function(p){
  if (colnames(cell_theta)[pair_cellPop[1,p, drop = T]] == "memory.b"){
    lab_x <- "memory b cells"
  }
  else{
    lab_x <- firstup(colnames(cell_theta)[pair_cellPop[1,p, drop = T]])
  }
  
  if (colnames(cell_theta)[pair_cellPop[2,p, drop = T]] == "memory.b"){
    lab_y <- "memory b cells"
  }
  else{
    lab_y <- firstup(colnames(cell_theta)[pair_cellPop[2,p, drop = T]])
  }
  df_temp <- data.frame(Gene = cell_theta$Gene, 
                        Pop1 = cell_theta[,pair_cellPop[1,p, drop = T]],
                        Pop2 = cell_theta[,pair_cellPop[2,p, drop = T]])
  
  rmse <- round(sqrt(mean((df_temp$Pop1 - df_temp$Pop2)^2)), 2)
  df_temp$rmse <- paste0("RMSE=", rmse)
  ggplot(df_temp) + aes(x=Pop1, y = Pop2) +
    # geom_point(alpha = .5) +
    scattermore::geom_scattermore(pointsize = 4, alpha = .3) +
    geom_abline(slope = 1, intercept = 0, colour = "darkred", linetype = "dashed", linewidth = 1) +
    xlab(lab_x) +
    ylab(lab_y) +
    # geom_label(aes(x = Inf, y = Inf, label = paste("RMSE =", round(rmse, 2))), 
    #            hjust = 1, vjust = 1, size = 6, color = "white", fill = "darkred") +
    facet_grid(~rmse) +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10), 
                  labels = c(0.01, 0.1, 1, 10)) +
    scale_y_log10() +
    annotation_logticks(side = "bl") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
})

# plt_application <- ((allPairplot[[1]]) + plot_spacer()  + allPairplot[[2]] + plot_spacer()  + allPairplot[[3]] + plot_spacer()  + allPairplot[[4]]) + plot_layout(nrow = 1, widths = c(4,.5,4,.5,4,.5,4))
plt_application <- ((allPairplot[[1]]) + allPairplot[[2]] + allPairplot[[3]] + allPairplot[[4]]) +
  plot_layout(nrow = 1) & 
  theme_classic() +
  theme(text = element_text(size = 18))

#------------------------------------------------------------------------------#
#                 Sup. Figure X: Applications on real world data sets          #
#                         of a two cell populations mixture                    #
#------------------------------------------------------------------------------#

results_application <- read.csv("results/SupplementaryResultsApplicationOnScRNASeq.csv") %>%
  mutate(ClusteringName = paste(Clustering, "Clustering", sep = " "),
         ThinningName = paste(Thinning, "Thinning", sep = " "))

tagName <- results_application %>% 
  group_by(ThinningName, ClusteringName) %>%
  arrange(pvalue) %>%
  slice_head(n = 5)

plt_scatter_pval <-results_application  %>%
  ggplot() +
  aes(x=-log10(TruePvalue),
      y=-log10(pvalue), 
      colour = ARI) +
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) +
  geom_point(alpha = .5) +
  ggrepel::geom_label_repel(data = tagName, aes(label = GeneSymbol)) +
  # geom_hline(yintercept = -log10(0.05)) +
  # geom_vline(xintercept = -log10(0.05)) +
  facet_grid(ClusteringName~ThinningName) +
  xlab("-log10(True p-value)") +
  ylab("-log10(p-value)") +
  scale_colour_viridis_c(option = "plasma",
                         direction = -1, limits = c(min(results_application$ARI),1))

ggsave(filename = "Supplementary Figures/SuppFigure_ApplicationTwoCellPops.pdf",
       width = 150, 
       height = 100, 
       units = "mm")


# venn_data <- results_application %>% 
#   mutate(Name = as.factor(paste(Clustering, Thinning, sep = "-")))
# 
# res <- lapply(levels(venn_data$Name), function(l){
#   venn_data %>% filter(Name == l) %>%
#     filter(pvalue < 0.05) %>%
#     pull(Gene)
# })
# names(res) <- levels(venn_data$Name)
# res$True <- venn_data %>% filter(TruePvalue<0.05) %>% pull(Gene) %>% unique()
# UpSetR::upset(fromList(res), nintersects = NA)

indicator_performances <- results_application %>% 
  group_by(Thinning, Clustering) %>%
  summarise(Power = mean(pvalue < 0.05 & TruePvalue<0.05), 
            typeI = mean(pvalue < 0.05 & TruePvalue > 0.05))
