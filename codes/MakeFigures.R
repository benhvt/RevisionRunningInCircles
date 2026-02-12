################################################################################
#                                                                              #  
#                         Main Make Figures Generations                        #      
#                         Benjamin Hivert - 17/04/2025                         #
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
library(scales)
theme_set(theme_classic())

#------------------------------------------------------------------------------#
#                           Internals functions                                #
#------------------------------------------------------------------------------#

# Marginal variance of a Gaussian Mixture
mixture_variance <- function(pi, mu_list, sigma_list) {
  # Function that computes the marginal variance of a Gaussian mixtures models of
  # parameters :
  #pi: a vector of length G that contains mixing proportion 
  # mu_list: a list of length G that contains vector of length p that are the means of each components
  #sigma_list: a list of G that contains matrix of dimension p x p that are the component-specific covariance
  #
  #return a matrix of dimension p x p that contains the marginal covariance of the mixture 
  
  G <- length(pi)  #Number of component in the mixture
  d <- length(mu_list[[1]])  # Dimension
  mu_Z <- Reduce("+", Map("*", pi, mu_list))  # Marginal mean of the mixture
  
  # Initialisation of the covariance matrix
  Sigma <- matrix(0, nrow = d, ncol = d)
  
  # Marginal variance computation
  for (g in 1:G) {
    Sigma <- Sigma + pi[g] * (sigma_list[[g]] + tcrossprod(mu_list[[g]]))
  }
  Sigma <- Sigma - tcrossprod(mu_Z)
  
  return(Sigma)
}

# Conditional covariance in marginal fission 
margFission_condCovariance <- function(g,
                                       pi,
                                       list_mu,
                                       list_Sigma){
  # Function that compute the conditional covariance between X1 and X2 when marginal
  # data fission is applied. This covariance is equale to Sigma_g - Sigma
  #
  # parameters :
  #g: The number of the component of the mixture where the conditional covariance must be computed
  #pi: a vector of length G that contains mixing proportion 
  # mu_list: a list of length G that contains vector of length p that are the means of each components
  #sigma_list: a list of G that contains matrix of dimension p x p that are the component-specific covariance
  #
  #return a matrix of dimension p x p that contains the conditional covariance 
  Sigma <- mixture_variance(pi, list_mu, list_Sigma)
  margFission <- Sigma - list_Sigma[[g]]
}


# Marginal covariance in conditional fission
condFission_margCovariance <- function(pi,
                                       list_mu,
                                       list_Sigma){
  # Function that compute the marginal covariance between X1 and X2 when conditional
  # data fission is applied.
  #
  # parameters :
  #pi: a vector of length G that contains mixing proportion 
  # mu_list: a list of length G that contains vector of length p that are the means of each components
  #sigma_list: a list of G that contains matrix of dimension p x p that are the component-specific covariance
  #
  #return a matrix of dimension p x p that contains the marginal covariance 
  
  # Marginal covariance in conditional fission
  list_pi <- list(pi[1], pi[2])
  margMu <- Reduce(`+`, Map(`*`,list_pi, list_mu))
  
  margSigma <- Reduce(`+`, Map(function(pi_k, mu_k) {
    diff <- matrix(mu_k - margMu, ncol = 1)
    pi_k * (diff %*% t(diff))  # Produit matriciel
  }, list_pi, list_mu))
  
}


compute_typeI <- function(n, tau, sigma, sigma_hat, alpha){
  # Function that computes theoritical type I error of data fission
  qu_t <- qt(alpha/2, df = n-2, lower.tail = FALSE)
  qu_n <- qnorm(alpha/2, lower.tail = F)
  cor_fg <- (sigma^2-sigma_hat^2)/sqrt((sigma^2 + (tau^2)*sigma_hat^2)*(sigma^2 + (1/tau^2)*sigma_hat^2))  
  mean_t <- sqrt(n)*sqrt((2/pi)*cor_fg^2)/sqrt(1-(2/pi)*cor_fg^2)
  
  return(pnorm(qu_n, mean = mean_t, sd= 1, lower.tail = FALSE) + pnorm(-qu_n, mean = mean_t, sd= 1, lower.tail = TRUE))
}

#------------------------------------------------------------------------------#
#                                  Figures                                     #
#               Behaviour of marginal and conditional data fission             # 
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                           Simulations setting                                #


# Parameters of the mixtures (2 components)
pi <- c(0.5, 0.5)
mu <- list(c(0,5),
           c(0,0))
Sigma <- list(diag(1,2),
              cbind(c(1, 0.5),
                    c(0.5, 1)))

# Generation of the illustrative data example 
set.seed(20250404)

G <- length(pi)
sample_size <- c(50, 100, 250, 500, 1000, 5000, 10000)

Z <- sample(1:G, sample_size[4], replace = TRUE, prob = pi)
X <- matrix(NA, nrow = sample_size[4], ncol = 2)
for (g in 1:G){
  X[Z==g,] <- rmvnorm(n = sum(Z==g), 
                      mean = mu[[g]], 
                      sigma = Sigma[[g]])
}

#------------------------------------------------------------------------------#
#                        First panel of the Figure                             #
#                  Statistical power in the ideal scenario                     #


# Figure generation 
plt_illu_power <- data.frame(X1 = X[,1],
                             X2 = X[,2],
                             Cluster = as.factor(kmeans(X, centers = 2, nstart = 100)$cluster)) %>%
  ggplot() +
  aes(x=X1, 
      y = X2, 
      colour = Cluster) +
  geom_point(size = 2,
             alpha = .8) +
  scale_colour_manual(name = "Estimated clusters",
                      values = c("#294122", "#EB3D00")) +
  theme_classic() + 
  theme(legend.position = "top") +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)'))

# Importation of results
ideal_sc <- read.csv(file = "results/IdealScenario_GaussianFission.csv")

# Derivation: compute Statistical power and ARI

powerResults <- ideal_sc %>%
  mutate(Fission_lab = paste(Fission, "fission", sep = " ")) %>%
  group_by(Fission_lab, Variable, tau, n) %>%
  summarise(Power = mean(pvalues < 0.05)) %>%
  mutate(Variable_lab = ifelse(Variable == "X1", "X[1]", "X[2]")) 

ariResults <- ideal_sc %>%
  mutate(Fission_lab = paste(Fission, "fission", sep = " ")) %>%
  group_by(Fission_lab, tau, n) %>%
  summarise(ARI_m = mean(ARI),
            ARI_sd = sd(ARI))

plt_power <- powerResults %>%
  mutate(Fission_lab = gsub(" ", "~", Fission_lab)) %>%
  ggplot() +
  aes(x=tau, y = Power, 
      colour = as.factor(n)) +
  geom_line(linewidth = .75,
            alpha = .8) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  scale_linetype_manual(name = "Data fission",
                        values = c(1, 6)) +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal level"),
             linewidth = .75, 
             linetype = "dashed",
             show.legend = FALSE) +
  scale_colour_manual(name = "",
                      values = "#DB2763") +
  facet_grid(Variable_lab~Fission_lab, 
             labeller = label_parsed) +
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Statistical Power \n 5% level")

plt_ari <- ggplot(ariResults) +
  aes(x=tau, y = ARI_m, 
      colour = as.factor(n)) +
  geom_line(linewidth = .75, 
            alpha = .8) +
  facet_grid(~Fission_lab) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Adjusted Rand \n Index")

figure2A <- ((plt_illu_power + theme(legend.position = "bottom")) +
               (plt_ari/plt_power + plot_layout(guides = "collect",
                                                heights = c(1,4)))) +
  plot_layout(widths = c(1.5,1.25)) +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(text = element_text(size = 9)) &
  theme(plot.tag = element_text(face = "bold"))

ggsave(figure2A,
       filename = "Figures/figure2A.pdf",
       width = 180,
       height = 100,
       units = "mm")
#------------------------------------------------------------------------------#
#                      Second panel of the Figure                              #
#                 Type I error rate in the worst scenario                      #


# Clustering: Only the first component of the mixture is clusterized into two spurious clusters
clust_comp <- kmeans(X[Z==1,], centers = 2, nstart = 100)$cluster
cluster <- rep(NA, sample_size[4])
cluster[Z==1] <- clust_comp
cluster[Z==2] <- 3

# Figure generation 
plt_illu_typeI <- data.frame(X1 = X[,1],
                             X2 = X[,2],
                             Cluster = as.factor(cluster)) %>%
  ggplot() +
  aes(x=X1, 
      y = X2, 
      colour = Cluster) +
  geom_point(size = 2,
             alpha = .8) +
  scale_colour_manual(name = "Estimated clusters",
                      values = c("#294122", "#EB3D00", "#FFBBA6")) +
  theme_classic() + 
  theme(legend.position = "bottom") +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)'))

# Import files of results
adverse_sc <- read.csv(file = "results/AdverseScenario_GaussianFission.csv")

plt_typeI <- adverse_sc %>%
  mutate(Fission_lab = paste(Fission, "fission", sep = "~")) %>%
  group_by(Fission_lab,Variable, tau, n) %>%
  summarise(TypeI = mean(pvalues < 0.05)) %>%
  mutate(Variable_lab = ifelse(Variable == "X1", "X[1]", "X[2]")) %>%
  ggplot() +
  aes(x=tau, y = TypeI, 
      colour = as.factor(n)) +
  geom_line(linewidth = .75, 
            alpha = .8) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal level"),
             linewidth = .75, 
             linetype = "dashed",
             show.legend = FALSE) +
  scale_colour_manual(name = "",
                      values = "#DB2763") +
  facet_grid(Variable_lab~Fission_lab, 
             labeller = label_parsed) +
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Type I error rate \n 5% nominal levels") +
  theme(legend.box = "vertical")


figure2B <- plt_illu_typeI + plt_typeI +
  plot_layout(widths = c(1.5,1.25)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(text = element_text(size = 8)) &
  theme(plot.tag = element_text(face = "bold"))

ggsave(figure2B,
       filename = "Figures/figure2B.pdf",
       width = 180,
       height = 100,
       units = "mm")


#------------------------------------------------------------------------------#
#                                  Figures                                     #
#               Impact of variance estimation on Type I error rate             # 
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                           Simulations setting                                #
rm("pi")
n <- 100
n_grid <- c(50, 100, 200, 500, 1000)
sigma <- c(0.1, 0.5, 1, 2)
sigma_grid <- sort(c(seq(0, 4, length.out = 50), 2, 0.1))
sigma_grid2 <- seq(0.8, 1.8, length = 50)
tau <- .4

#------------------------------------------------------------------------------#
#             First panel: Function of the original variance                   #


sigma_grid_plot <- seq(0, 6*max(sigma), length.out = 1000)

typeI_theo <- lapply(sigma_grid_plot, function(s){
  temp <- compute_typeI(n, tau , sigma, sigma_hat = s, alpha = .05)
  return(data.frame(TypeItheo = temp,
                    sigma_hat = s,
                    sigma = sigma))
})

typeI_df <- do.call(rbind.data.frame, typeI_theo) %>% mutate(sigma_name = paste0("sigma^2==", sigma^2))


var_bias <- read.csv(file = "results/VarianceEstimationAndTypeIError.csv")

var_bias_res <- var_bias %>%
  group_by(sigma, sigma_hat, n) %>%
  summarise(EmpTypeI = mean(pvalues < 0.05)) %>%
  ungroup() %>%
  # mutate(TheTypeI = compute_typeI(n = n, tau = tau, sigma = sigma, sigma_hat = sigma_hat, alpha = .05)) %>%
  mutate(sigma_name = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = (sigma^2-sigma_hat^2)/sigma^2) %>%
  reshape2::melt(id.vars = c("Ratio", "sigma_name", "sigma", "sigma_hat"), measure_vars = c("TheTypeI", "EmpTypeI"))

plot_sigma <-  ggplot(var_bias_res) + aes(x=Ratio, y = value, colour = sigma_name) +
  geom_point(data = subset(var_bias_res, variable == "EmpTypeI"), aes(shape = "Empirical"), size = 4) +
  scale_shape_manual(name = "Empirical", values = 2, labels = '', guide = "legend") +
  # geom_line(data = subset(df_res_sigma, variable == "TheTypeI"), aes(group = sigma_name, linetype = "Theoritical"), size = 1) +
  geom_line(data = typeI_df, aes(x= (sigma^2-sigma_hat^2)/sigma^2, y = TypeItheo, group = sigma_name, colour = sigma_name, linetype = "Theoritical"), size = 1) +
  scale_linetype_manual(name = "Theoritical", values = 1, labels = "", guide = "legend") +
  scale_colour_manual(name = '',
                      values = c("#93B5C6", "#DBC2CF", "#998bc0", "#BD4F6C"),
                      labels = c(TeX(r'($\sigma^2 = 0.01$)'),
                                 TeX(r'($\sigma^2 = 0.25$)'),
                                 TeX(r'($\sigma^2 = 1$)'),
                                 TeX(r'($\sigma^2 = 4$)'))) +
  xlab(TeX(r'($(\sigma^2 - \widehat{\sigma^2})/\sigma^2$)')) +
  xlim(c(-5,2)) +
  ylab("Type I error rate") +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal levels"),
             linetype = 2,
             size = 1.2) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  guides(
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    color = guide_legend(order = 3)
  ) +
  NULL


#------------------------------------------------------------------------------#
#               Second panel: Function of the sample size                      #


var_bias_sampsize <- read.csv(file = "results/VarianceEstimationAndTypeIErrorSampleSize.csv")

var_bias_sampsize_res <-  var_bias_sampsize %>%
  group_by(n, sigma_hat, sigma) %>%
  rename(SampSize = n) %>%
  summarise(EmpTypeI = mean(pvalues < 0.05)) %>%
  ungroup() %>%
  mutate(TheTypeI = compute_typeI(n = SampSize, 
                                  tau = tau, 
                                  sigma = sigma,
                                  sigma_hat = sigma_hat, 
                                  alpha = .05)) %>%
  mutate(sigma_name = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = (sigma^2-sigma_hat^2)/sigma^2) %>%
  reshape2::melt(id.vars = c("Ratio", "sigma_name", "SampSize", "sigma", "sigma_hat"), measure_vars = c("TheTypeI", "EmpTypeI"))

plot_sigma_n <- ggplot(var_bias_sampsize_res) + aes(x=Ratio, y = value, colour = as.factor(SampSize)) +
  geom_point(data = subset(var_bias_sampsize_res, variable == "EmpTypeI"),
             aes(shape = "Empirical"), 
             size = 4) +
  scale_shape_manual(name = "Empirical",
                     values = 2, labels = '',
                     guide = "legend") +
  geom_line(data = subset(var_bias_sampsize_res, variable == "TheTypeI"), 
            aes(group = SampSize, linetype = "Theoritical"), size = 1) +
  scale_linetype_manual(name = "Theoritical", 
                        values = 1,
                        labels = "", 
                        guide = "legend") +
  scale_colour_manual(name = "", 
                      values = colorRampPalette(c("#008154", "#0092a4", "#2a2956"))(length(n_grid)),
                      labels = paste0("n=", n_grid)) +
  xlab(TeX(r'($(\sigma^2 - \widehat{\sigma^2})/\sigma^2$)')) +
  ylab("Type I error rate") +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal levels"),
             linetype = 2,
             size = 1.2) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  guides(
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    color = guide_legend(order = 3)
  ) +
  NULL

(plot_sigma + plot_sigma_n)  +
  plot_annotation(tag_levels = "A") &
  theme_classic() &
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.tag = element_text(face = "bold", size = 24),
        legend.spacing = unit(0.01, "cm"))

ggsave(filename = "Figures/figure3.pdf",
       width = 350, 
       height = 120, 
       units = "mm",
       dpi = 600)

#------------------------------------------------------------------------------#
#                                  Figure 4                                    #
#         Data fission for scRNA-seq (Binomial Negatives Simulations)          # 
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#           Panel A and B: Simulations results under Mixture Settings          #

#-- Parameters 
n <- 100
Z <- rep(1:2, each = n/2)
probs <- c(0.5, 0.4)
size <- c(5, 40)
set.seed(09022024)

x1 <- rnbinom(n=n/2, prob = probs[1], size = size[1])
x2 <- rnbinom(n=n/2, prob = probs[2], size = size[2])
y <- rnbinom(n=n, prob = .5, size = 5)
X <- cbind.data.frame(X1=c(x1,x2), 
                      X2=y, 
                      TrueClasses=as.factor(Z))
cl_X <- kmeans(cbind(c(x1,x2), y), centers = 3)$cluster

plt_illu_negbin <- ggplot(X) + aes(x=X1, 
                                   y=X2, 
                                   colour = as.factor(cl_X), 
                                   shape = as.factor(Z)) + 
  geom_point(size = 3) +
  scale_shape_manual(name = "True classes",
                     values = c(15,16)) +
  scale_colour_manual(name = "Clusters", 
                      values = c("#294122", "#EB3D00", "#FFBBA6"),
                      labels = c(TeX(r'($C_1$)'),
                                 TeX(r'($C_2$)'),
                                 TeX(r'($C_3)'))) +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)')) +
  theme_classic() +
  # theme(legend.position = "bottom") +
  NULL

neg_bin_typeI <- read.csv("results/TypeIThinningNegBin.csv")

plt_typeI_negbin <- neg_bin_typeI %>%
  ggplot() + 
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pvalues, colour = factor(Fission)),
          distribution = qunif, size = 1.5) +
  scale_colour_manual(name = "", 
                      values = c("#334EAC", "#BAD6EB")) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic()  +
  theme(legend.position = "bottom")

plt_negbin <- (plt_illu_negbin + ggtitle("A") + theme(plot.title = element_text(face = "bold", size = 20))) + 
  (plt_typeI_negbin + ggtitle("B") + theme(plot.title = element_text(face = "bold", size = 20)))


#------------------------------------------------------------------------------#
#            Panel C: Simulations results under Correlated setting             #


neg_bin_cor_res <- read.csv("results/CorrelationAndRelativeBiais.csv")

plt_neg_bin_cor <- neg_bin_cor_res %>% 
  mutate(Theta_hat = Theta*Error) %>% 
  mutate(RelativeBiais = (Theta_hat-Theta)/Theta) %>% 
  group_by(Method, Estimation, Rho, RelativeBiais, Error) %>% 
  filter(Rho != 0.01) %>%
  filter(!(RelativeBiais<0 & Method == "NB")) %>%
  filter(RelativeBiais < 5) %>%
  summarise(FDR = mean(p.adjust(pval, method = "BH") < 0.05),
            typeI = mean(pval < 0.05)) %>%
  mutate(Name = paste(Method, Estimation, collapse = "_")) %>%
  mutate(
    LabelName = case_when(
      Name == "Gauss Oracle" ~ "Gaussian~(theta)",
      Name == "Gauss Wrong"  ~ "Gaussian~(widehat(theta))",
      Name == "NB Oracle"    ~ "Negative~Binomial~(theta)",
      Name == "NB Wrong"     ~ "Negative~Binomial~(widehat(theta))",
      TRUE ~ Name
    )
  ) %>%
  ggplot() +
  aes(x=RelativeBiais,
      y = typeI, 
      colour = as.factor(Rho)) +
  geom_point(size = 2) +
  geom_line(linewidth = .9) +
  facet_grid(~LabelName, labeller = label_parsed, scales = "free_x") +
  geom_hline(yintercept = 0.05,
             color = "grey20",
             linewidth = 1.2,
             linetype = "dashed") +
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
  xlab(TeX(r'($\frac{\widehat{\theta} - \theta}{\theta}$)')) +
  ylab("Type I error") + 
  theme_classic() +
  # scale_x_continuous(breaks = c(round(seq(-0.99, 9, length.out = 5), 0) ,0)) +
  theme(legend.position = 'bottom') +
  guides(
    colour = guide_legend(nrow = 1))

((plt_illu_negbin + ggtitle("A")) + 
  (plt_typeI_negbin + ggtitle("B"))) /
  (plt_neg_bin_cor + ggtitle("C")) +
  plot_layout(heights = c(1.5, 2)) &
  theme(plot.title = element_text(face = "bold", size = 22),
        text = element_text(size = 20))



ggsave("Figures/figure4.pdf", width = 315, height = 200, units = "mm")


#------------------------------------------------------------------------------#
#                                 Figure 5                                     #
#                         Applications on Bonne marrow                         #
#------------------------------------------------------------------------------#


Xtrain <- read_xlsx("results/ResultsOnScRNASeq.xlsx", sheet = "Xtrain")
Xtest <- read_xlsx("results/ResultsOnScRNASeq.xlsx", sheet = "Xtest")
cluster <- read_xlsx("results/ResultsOnScRNASeq.xlsx", sheet = "Cluster")
results <- read_xlsx("results/ResultsOnScRNASeq.xlsx", sheet = "Results")


#------------------------------------------------------------------------------#
#                 Panel A Gene-specific overdispersion parameters              #

pval_toQQ <- results %>% dplyr::select(Gene,
                                PvaluesAfterUnivariateClustering,
                                PValuesAfterMultivariateClustering) %>%
  tidyr::pivot_longer(
    cols = c(PvaluesAfterUnivariateClustering, PValuesAfterMultivariateClustering),
    names_to = "Method",
    values_to = "Pvalue"
  ) %>%
  mutate(Method = case_when(
    Method == "PvaluesAfterUnivariateClustering" ~ "Univariate \n Clustering",
    Method == "PValuesAfterMultivariateClustering" ~ "Multivariate \n Clustering"
  ))

plt_qqApplication <- ggplot(pval_toQQ) +
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + 
  stat_qq(aes(sample = Pvalue, colour = factor(Method)),
          distribution = qunif, size = 1.5) +
  scale_colour_manual(name = "Method",
                      values = c("#8C33FF",  # bleu
                                 "#2CA02C")) +
  xlab("Theoretical Quantiles") + 
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic() +
  theme(text = element_text(size = 14))

#------------------------------------------------------------------------------#
#         Panel B-C Correlation between genes lead to problems even on         #
#                         homogeneous sub-population                           #

# Correlation plots 
plt_cor_orig_thin <- ggplot(results) +
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + 
  
  aes(x=OriginalCorWithGene1, 
      y = ThinningCorWithGene1TrainAndTest) +
  geom_point(size = 1,
             alpha = .5) +
  ylab(TeX(r'(Cor$\left(X^{(1)}_1, X^{(2)}_j\right)$)')) +
  xlab(TeX(r'(Cor$\left(X_1, X_j\right)$)')) +
  xlim(c(-1,1)) +
  ylim(c(-1,1)) +
  theme_classic()


## Sélection du gène
GeneToSel <- sort(
  results$PValuesAfterMultivariateClustering,
  index.return = TRUE
)$ix[1]

## Clustering univarié (Train)
km_uni <- kmeans(
  log2(Xtrain[, GeneToSel] + 1),
  centers = 2,
  nstart  = 100
)

## Réplication Train + Test
cluster_uni <- c(km_uni$cluster, km_uni$cluster)

Xtrain1 <- tibble(
  Gene1 = c(
    Xtrain[, GeneToSel, drop = TRUE],
    Xtest[, GeneToSel, drop = TRUE],
    Xtrain[, GeneToSel, drop = TRUE],
    Xtest[, GeneToSel, drop = TRUE]
  ),
  Cluster = factor(c(
    rep(cluster$Cluster, 2), # Multivarié
    cluster_uni # Univarié
  )),
  Method = factor(
    rep(c("Multivariate", "Univariate"),
        each = 2 * nrow(Xtrain))
  ),
  Thinning = factor(
    rep(rep(c("Train", "Test"), each = nrow(Xtrain)), 2),
    levels = c("Train", "Test"))
) %>%
  mutate(
    MethodName = paste(Method, "Clustering", sep = " ")
  )


pvals <- Xtrain1 %>%
  group_by(MethodName, Thinning) %>%
  summarise(
    p_value = wilcox.test(log2(Gene1 + 1) ~ Cluster)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Wilcoxon~italic(p)~\"=\"~\"",
      scales::pvalue(p_value),
      "\""
    ),
    y_pos = max(log2(Xtrain1$Gene1 + 1), na.rm = TRUE) * 1.25
  )

plt_illuApplication <- 
  Xtrain1%>%
  ggplot(
    aes(
      x    = Thinning,
      y    = log2(Gene1 + 1),
      fill = Cluster
    )
  ) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE) +
  geom_boxplot(width = .2, alpha = 1, show.legend = FALSE) +
  
  ## P-values (solution robuste)
  geom_text(
    data = pvals,
    aes(
      x     = Thinning,
      y     = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    parse = TRUE,
    size  = 4
  ) +
  
  scale_fill_manual(
    name   = "PGK1",
    values = rev(colorspace::lighten(c("#294122", "#EB3D00"), 0.25)),
    labels = c(TeX(r'($C_1$)'), TeX(r'($C_2$)'))
  ) +
  scale_colour_manual(
    name   = "PGK1",
    values = rev(colorspace::lighten(c("#294122", "#EB3D00"), 0.25)),
    labels = c(TeX(r'($C_1$)'), TeX(r'($C_2$)'))
  ) +
  facet_grid(~MethodName) +
  ylab("log2(counts + 1)") +
  xlab("") +
  scale_x_discrete(
    labels = c(
      "Train" = expression(X^{(1)}),
      "Test"  = expression(X^{(2)})
    )
  )

plt_illuApplication

figure5 <- ((plt_cor_orig_thin + plt_qqApplication ) / (plt_illuApplication)) +
  plot_annotation(tag_levels = "A") +
  plot_layout(height = c(2, 4)) &
  theme_classic() &
  theme(text = element_text(size = 14),
        plot.tag = element_text(face = "bold", size = 16))

ggsave("Figures/figure5.pdf", plot = figure5, width = 250, height = 200, units = "mm")


# Table pvalues 
typeI_table <- data.frame(Method = c("Multivariate Clustering",
                                     "Univariate Clustering"),
                          `Type I` = c(mean(results$PValuesAfterMultivariateClustering < 0.05),
                                       mean(results$PvaluesAfterUnivariateClustering < 0.05))) 
xt <- xtable::xtable(typeI_table, 
             caption = "Type I error rates from Wilcoxon tests after Negative Binomial data thinning. 
  In the multivariate setting, k-means clustering is performed on all 500 genes in X^{(1)}, 
  and Wilcoxon tests are applied gene-wise on X^{(2)}. 
  In the univariate setting, clustering is performed individually for each gene in X^{(1)}."
)

xtable::print.xtable(xt, 
      type = "latex",           # génère du LaTeX
      file = "results/TypeIErrorRateNBOnScRNASeq.txt",    # nom du fichier
      include.rownames = FALSE) 


