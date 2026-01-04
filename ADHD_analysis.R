library(tidyverse)
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(MendelianRandomization)
setwd("D://RHDD/")
sourceCpp("gibbssexbiasedmr.cpp")

res <- matrix(nrow = 0, ncol = 16)
colnames(res) = c("exposure", "samples", "numIVs", "est_beta", "est_beta_SX", "est_beta_mix", "pval_beta", "pval_beta_SX", "pval_beta_mix", "beta_init", "beta_SX_init", "beta_mix_init", "p_thr", "n_iter", "clump_kb", "clump_r2")

setwd("D://RHDD/clumpeddata")
for (i in 1:length(list.files())) # length(list.files())
  { 
    print(i)
    input <- list.files()[i]
    inputinfo <- str_match(input, "^clumpedGWAS_(\\d+)_(\\d\\.\\d+)_(\\de-\\d+)_imp_imputed_(.+).txt.gz.txt")
    df <- fread(input)
    numIV <- length(df$variant_id)
    n_iter = numIV*50
    
    if (numIV <= 1000)
      {
        mr_ivw_single1 <- mr_ivw(mr_input(bx = df$exposure.beta, bxse = df$exposure.se, by = df$outcome.male.beta, byse = df$outcome.male.se)) # male is 1
        mr_ivw_single2 <- mr_ivw(mr_input(bx = df$exposure.beta, bxse = df$exposure.se, by = df$outcome.female.beta, byse = df$outcome.female.se)) # female is 2
        mr_ivw_mix <- mr_ivw(mr_input(bx = df$exposure.beta, bxse = df$exposure.se, by = df$outcome.combined.beta, byse = df$outcome.combined.se))
        
        beta_init <- mr_ivw_single1$Estimate # Single is male
        beta_SX_init <- mr_ivw_single2$Estimate - mr_ivw_single1$Estimate # Mixed is female
    
        # 2 sample
        post_samples <- gibbs_sampling_cpp2mix(n_iter = n_iter, n_thin = 5, p = numIV,
                                               Shat_Gamma_mix_sq = (df$outcome.female.se)^2, Shat_Gamma_single_sq = (df$outcome.male.se)^2, Shat_gamma_sq = (df$exposure.se)^2,
                                               hat_Gamma_mix = df$outcome.female.beta, hat_Gamma_single = df$outcome.male.beta, hat_gamma = df$exposure.beta,
                                               1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1, 1, 1, 1,
                                               prob = 1, beta_init = 0, beta_SX_init = 0)
        
        selected_samples_beta<-post_samples$samples_beta[seq((round(n_iter/5)+1), n_iter, 1)]
        selected_samples_beta_SX<-post_samples$samples_beta_SX[seq((round(n_iter/5)+1), n_iter, 1)]
        est_beta <- mean(selected_samples_beta)
        est_beta_SX <- mean(selected_samples_beta_SX)
        est_beta_combined <- est_beta + 1 * est_beta_SX
        pval_beta<-2*(1-pnorm(abs(mean(selected_samples_beta))/sd(selected_samples_beta)))
        pval_beta_SX<-2*(1-pnorm(abs(mean(selected_samples_beta_SX))/sd(selected_samples_beta_SX)))
        pval_beta_combined <- 2*(1-pnorm(abs(est_beta + 1 * est_beta_SX)/sd(selected_samples_beta_SX+selected_samples_beta)))
        
        res <- rbind(res, c(inputinfo[5], "2", numIV, est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, "", inputinfo[4], n_iter, inputinfo[2], inputinfo[3]))
        # print(c(inputinfo[5], "IVW", "2", length(df$variant_id), est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, "", inputinfo[4], n_iter, inputinfo[2], inputinfo[3]))
        
        # 3 sample
        post_samples <- gibbs_sampling_cpp(n_iter = n_iter, n_thin = 5, p = numIV, Shat_Gamma_mix1_sq = (df$outcome.combined.se)^2,
                                           Shat_Gamma_mix2_sq = (df$outcome.male.se)^2, Shat_Gamma_mix3_sq = (df$outcome.female.se)^2, Shat_gamma_sq = (df$exposure.se)^2,
                                           hat_Gamma_mix1 = df$outcome.combined.beta, hat_Gamma_mix2 = df$outcome.male.beta, hat_Gamma_mix3 = df$outcome.female.beta, hat_gamma = df$exposure.beta,
                                           1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1, 1, 1, 1, 1, 1,
                                           prob1 = 111897/225534, prob2 = 0, prob3 = 1, beta_init = 0, beta_SX_init = 0)
        
        selected_samples_beta<-post_samples$samples_beta[seq((round(n_iter/5)+1), n_iter, 1)]
        selected_samples_beta_SX<-post_samples$samples_beta_SX[seq((round(n_iter/5)+1), n_iter, 1)]
        est_beta <- mean(selected_samples_beta)
        est_beta_SX <- mean(selected_samples_beta_SX)
        est_beta_combined <- est_beta + 1 * est_beta_SX
        pval_beta<-2*(1-pnorm(abs(mean(selected_samples_beta))/sd(selected_samples_beta)))
        pval_beta_SX<-2*(1-pnorm(abs(mean(selected_samples_beta_SX))/sd(selected_samples_beta_SX)))
        pval_beta_combined <- 2*(1-pnorm(abs(est_beta + 1 * est_beta_SX)/sd(selected_samples_beta_SX+selected_samples_beta)))
        
        res <- rbind(res, c(inputinfo[5], "3", numIV, est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, mr_ivw_mix$Estimate, inputinfo[4], n_iter, inputinfo[2], inputinfo[3]))
        # print(c(inputinfo[5], "IVW", "3", length(df$variant_id), est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, mr_ivw_mix$Estimate, inputinfo[4], n_iter, inputinfo[2], inputinfo[3]))
    }
}

resdf <- as.data.frame(res)
nrow(resdf) 
nrow(filter(resdf, pval_beta_SX < 0.001)) 
nrow(filter(resdf, samples == 2, pval_beta_SX < 0.001)) 
nrow(filter(resdf, samples == 3, pval_beta_SX < 0.001)) 

save(res, file="res814.RData")
