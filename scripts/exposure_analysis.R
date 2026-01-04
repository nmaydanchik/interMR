library(tidyverse)
library(data.table)
library(ieugwasr)
library(Rcpp)
library(RcppArmadillo)
library(MendelianRandomization)
library(mr.raps)
setwd("D://RHDD/")
sourceCpp("gibbssexbiasedmr.cpp")

exposure = "CUD"
set.seed(1234)

# 2 sample InterMR --------------------------------------------------------

res <- matrix(nrow = 0, ncol = 18)
colnames(res) = c("exposure", "study", "year", "init_method", "samples", "numIVs", "est_beta", "est_beta_SX", "est_beta_combined", "pval_beta", "pval_beta_SX", "pval_beta_combined", "beta_init", "beta_SX_init", "p_thr", "n_iter", "clump_kb", "clump_r2")

y1 <- fread("ADHD/processed_ADHD2018_Martin_male.gz", select=c(1:5)); print("Loaded Y1") # Male is 1
y2 <- fread("ADHD/processed_ADHD2018_Martin_female.gz", select=c(1:5)); print("Loaded Y2") # Female is 2

for (p_thr in c(5e-6, 5e-8)) # c(5e-6, 5e-8)
  {
  for (n_iter in c(5000))
    {
    for (kb in c(100, 250, 500)) # c(100, 250, 500)
      {
      for (r2 in c(0.1, 0.5)) # c(0.1, 0.5)
        {
        for (i in 1:length(str_subset(list.files(paste0(exposure,"/")), "processed")))
          {
            input <- str_subset(list.files(paste0(exposure,"/")), "processed")[i]
            x <- fread(paste0(exposure,"/",input)); print("Loaded X")

            df <- inner_join(x, y1, by = "rsid"); rm(x); print("Combined df 1")
            df <- inner_join(df, y2, by = "rsid"); print("Combined df 2")
            df <- filter(df, eff.y == eff & noneff.y == noneff) # Ensure alleles of single-sex ADHD match
            
            dfp <- filter(df, pval < p_thr); print("Filtered pval")
            
            # Match alleles of exposures and outcome. Inverse beta of exposure if alleles are switched
            for (j in 1:nrow(dfp)) {
              if (dfp$eff.x[j] == dfp$noneff.y[j] && dfp$noneff.x[j] == dfp$eff.y[j]) {
                dfp$eff.x[j] = dfp$eff.y[j]
                dfp$noneff.x[j] = dfp$noneff.y[j]
                dfp$beta.x[j] = -1*dfp$beta.x[j]
              } 
            }
            dfp <- filter(dfp, eff.x == eff.y & noneff.x == noneff.y) # Ensure all alleles are now matched
            dfp <- select(dfp, c(1,4,5,6,9,10,13,14))
            colnames(dfp) = c("rsid", "gammah", "seg", "pval", "Gammahs1", "se1", "Gammahs2", "se2")  
            print("Matched alleles")
            
            dfp <- filter(dfp, se1 != 0) %>% filter(se2 != 0) %>% filter(seg != 0) # Ensure no standard errors are 0 (IVW will break)
            
            dfl <- ld_clump(
              dfp,
              plink_bin = "plink/plink.exe",
              bfile = "1kg.v3/EUR",
              clump_kb = kb,
              clump_r2 = r2
            ) %>% select(c(1:8))
            print("LD clumped")
            
            #IVW prior estimate
            mr_ivw_single1 <- mr_ivw(mr_input(bx = dfl$gammah, bxse = dfl$seg, by=dfl$Gammahs1, byse=dfl$se1))
            mr_ivw_single2 <- mr_ivw(mr_input(bx = dfl$gammah, bxse = dfl$seg, by=dfl$Gammahs2, byse=dfl$se2))
            
            beta_init <- mr_ivw_single1$Estimate # Single is male
            beta_SX_init <- mr_ivw_single2$Estimate - mr_ivw_single1$Estimate # Mixed is female
            print("Initialized priors (IVW)")
            
            # Gibbs Sampler
            post_samples <- gibbs_sampling_cpp2mix(n_iter = n_iter, n_thin = 5, p = length(dfl$rsid),
                                                   Shat_Gamma_mix_sq = (dfl$se2)^2, Shat_Gamma_single_sq = (dfl$se1)^2, Shat_gamma_sq = (dfl$seg)^2,
                                                   hat_Gamma_mix = dfl$Gammahs2, hat_Gamma_single = dfl$Gammahs1, hat_gamma = dfl$gammah,
                                                   1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1, 1, 1, 1,
                                                   prob = 1, beta_init = beta_init, beta_SX_init = beta_SX_init)
            print("Ran Gibbs sampler")
            
            selected_samples_beta<-post_samples$samples_beta[seq((round(n_iter/5)+1), n_iter, 1)]
            selected_samples_beta_SX<-post_samples$samples_beta_SX[seq((round(n_iter/5)+1), n_iter, 1)]
            est_beta <- mean(selected_samples_beta)
            est_beta_SX <- mean(selected_samples_beta_SX)
            est_beta_combined <- est_beta + 1 * est_beta_SX
            pval_beta<-2*(1-pnorm(abs(mean(selected_samples_beta))/sd(selected_samples_beta)))
            pval_beta_SX<-2*(1-pnorm(abs(mean(selected_samples_beta_SX))/sd(selected_samples_beta_SX)))
            pval_beta_combined <- 2*(1-pnorm(abs(est_beta + 1 * est_beta_SX)/sd(selected_samples_beta_SX+selected_samples_beta)))
            print("Calculated result")
            
            res <- rbind(res, c(exposure, str_match(input, "^\\w+_\\w+_(\\w+)")[2], str_match(input, "^\\w+_[:alpha:]+([:digit:]+)")[2], "IVW", "2", length(dfl$rsid), 
                                est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, p_thr, n_iter, kb, r2))
            print(c(exposure, str_match(input, "^\\w+_\\w+_(\\w+)")[2], str_match(input, "^\\w+_[:alpha:]+([:digit:]+)")[2], "IVW", "2", length(dfl$rsid),
                    est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, p_thr, n_iter, kb, r2))
            

            #MR-RAPS prior estimate
            mr_raps_single1 <- mr.raps(b_exp = dfl$gammah, b_out = dfl$Gammahs1, se_exp = dfl$seg, se_out = dfl$se1)
            mr_raps_single2 <- mr.raps(b_exp = dfl$gammah, b_out = dfl$Gammahs2, se_exp = dfl$seg, se_out = dfl$se2)

            beta_init <- mr_raps_single1$beta.hat # Single is male
            beta_SX_init <- mr_raps_single2$beta.hat - mr_raps_single1$beta.hat # Mixed is female
            print("Initialized priors (MR-RAPS)")
            
            # Gibbs Sampler
            post_samples <- gibbs_sampling_cpp2mix(n_iter = n_iter, n_thin = 5, p = length(dfl$rsid),
                                                   Shat_Gamma_mix_sq = (dfl$se2)^2, Shat_Gamma_single_sq = (dfl$se1)^2, Shat_gamma_sq = (dfl$seg)^2,
                                                   hat_Gamma_mix = dfl$Gammahs2, hat_Gamma_single = dfl$Gammahs1, hat_gamma = dfl$gammah,
                                                   1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1, 1, 1, 1,
                                                   prob = 1, beta_init = beta_init, beta_SX_init = beta_SX_init)
            print("Ran Gibbs sampler")
            selected_samples_beta<-post_samples$samples_beta[seq((round(n_iter/5)+1), n_iter, 1)]
            selected_samples_beta_SX<-post_samples$samples_beta_SX[seq((round(n_iter/5)+1), n_iter, 1)]
            est_beta <- mean(selected_samples_beta)
            est_beta_SX <- mean(selected_samples_beta_SX)
            est_beta_combined <- est_beta + 1 * est_beta_SX
            pval_beta<-2*(1-pnorm(abs(mean(selected_samples_beta))/sd(selected_samples_beta)))
            pval_beta_SX<-2*(1-pnorm(abs(mean(selected_samples_beta_SX))/sd(selected_samples_beta_SX)))
            pval_beta_combined <- 2*(1-pnorm(abs(est_beta + 1 * est_beta_SX)/sd(selected_samples_beta_SX+selected_samples_beta)))
            print("Calculated result")

            res <- rbind(res, c(exposure, str_match(input, "^\\w+_\\w+_(\\w+)")[2], str_match(input, "^\\w+_[:alpha:]+([:digit:]+)")[2], "MR-RAPS", "2", length(dfl$rsid),
                                est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, p_thr, n_iter, kb, r2))
            print(c(exposure, str_match(input, "^\\w+_\\w+_(\\w+)")[2], str_match(input, "^\\w+_[:alpha:]+([:digit:]+)")[2], "MR-RAPS", "2", length(dfl$rsid),
                    est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, p_thr, n_iter, kb, r2))
            
        }
      }
    }
  }
}

save(res, file = "res2sampleCUD.RData")
write_csv(as.data.frame(res), "res2sampleCUD.csv")

# 3 Sample InterMR --------------------------------------------------------

res <- matrix(nrow = 0, ncol = 19)
colnames(res) = c("exposure", "study", "year", "init_method", "samples", "numIVs", "est_beta", "est_beta_SX", "est_beta_combined", "pval_beta", "pval_beta_SX", "pval_beta_combined", "beta_init", "beta_SX_init", "beta_mix", "p_thr", "n_iter", "clump_kb", "clump_r2")

ymix <- fread("ADHD/processed_ADHD2022_Demontis.gz", select=c(1:5)); print("Loaded Ymix")
y1 <- fread("ADHD/processed_ADHD2018_Martin_male.gz", select=c(1:5)); print("Loaded Y1") # Male is 1
y2 <- fread("ADHD/processed_ADHD2018_Martin_female.gz", select=c(1:5)); print("Loaded Y2") # Female is 2

for (p_thr in c(5e-6, 5e-8))
{
  for (n_iter in c(5000))
  {
    for (kb in c(100, 250, 500))
    {
      for (r2 in c(0.1, 0.5))
      {
        for (i in 1:length(str_subset(list.files(paste0(exposure,"/")), "processed"))) {
          input <- str_subset(list.files(paste0(exposure,"/")), "processed")[i]
          x <- fread(paste0(exposure,"/", input)); print("Loaded X")
          
          df <- inner_join(x, ymix, by = "rsid"); rm(x); print("Combined X Ymix")
          df <- inner_join(df, y1, by = "rsid") %>% inner_join(y2, by = "rsid"); print("Combined all")
          
          dfp <- filter(df, pval < p_thr); print("Filtered pval")
          
          dfp <- filter(dfp, eff.x.x == eff.y.y & noneff.x.x == noneff.y.y) # Ensure alleles of single-sex ADHD match
          
          # match alleles of mix and single sex
          for (j in 1:nrow(dfp)) {
            if (dfp$eff.y[j] == dfp$noneff.y.y[j] && dfp$noneff.y[j] == dfp$eff.y.y[j]) {
              dfp$eff.y[j] = dfp$eff.y.y[j]
              dfp$noneff.y[j] = dfp$noneff.y.y[j]
              dfp$beta.y[j] = -1*dfp$beta.y[j]
            } 
          }
          dfp <- filter(dfp, eff.y == eff.y.y & noneff.y == noneff.y.y) # Ensure all Y alleles are now matched
          print("Matched Y alleles")
          
          # match alleles of x and y
          for (j in 1:nrow(dfp)) {
            if (dfp$eff.x[j] == dfp$noneff.y.y[j] && dfp$noneff.x[j] == dfp$eff.y.y[j]) {
              dfp$eff.x[j] = dfp$eff.y.y[j]
              dfp$noneff.x[j] = dfp$noneff.y.y[j]
              dfp$beta.x[j] = -1*dfp$beta.x[j]
            } 
          }
          print("Matched all alleles")
          dfp <- filter(dfp, eff.x == eff.y.y & noneff.x == noneff.y.y) # Ensure all alleles are now matched
          dfp <- select(dfp, c(1,4,5,6,9,10,13,14,17,18))
          colnames(dfp) = c("rsid", "gammah", "seg", "pval", "Gammah", "seG", "Gammahs1", "se1", "Gammahs2", "se2")  
          
          dfp <- filter(dfp, se1 != 0) %>% filter(se2 != 0) %>% filter(seg != 0) %>% filter(seG != 0)
        
          dfl <- ld_clump(
            dfp,
            plink_bin = "plink/plink.exe",
            bfile = "1kg.v3/EUR",
            clump_kb = kb,
            clump_r2 = r2
          ) %>% select(c(1:10))
          print("LD Clumped")
          
          mr_ivw_mix <- mr_ivw(mr_input(bx = dfl$gammah, bxse = dfl$seg, by = dfl$Gammah, byse = dfl$seG))
          mr_ivw_single1 <- mr_ivw(mr_input(bx = dfl$gammah, bxse = dfl$seg, by=dfl$Gammahs1, byse=dfl$se1))
          mr_ivw_single2 <- mr_ivw(mr_input(bx = dfl$gammah, bxse = dfl$seg, by=dfl$Gammahs2, byse=dfl$se2))
          
          beta_init <- mr_ivw_single1$Estimate # Single is male
          beta_SX_init <- mr_ivw_single2$Estimate - mr_ivw_single1$Estimate # Mixed is female
          print("Initialized priors (IVW)")
          
          post_samples <- gibbs_sampling_cpp(n_iter = n_iter, n_thin = 5, p = length(dfl$rsid), Shat_Gamma_mix1_sq = (dfl$seG)^2,
                                             Shat_Gamma_mix2_sq = (dfl$se1)^2, Shat_Gamma_mix3_sq = (dfl$se2)^2, Shat_gamma_sq = (dfl$seg)^2,
                                             hat_Gamma_mix1 = dfl$Gammah, hat_Gamma_mix2 = dfl$Gammahs1, hat_Gamma_mix3 = dfl$Gammahs2, hat_gamma = dfl$gammah,
                                             1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1, 1, 1, 1, 1, 1,
                                             prob1 = 111897/225534, prob2 = 0, prob3 = 1, beta_init = beta_init, beta_SX_init = beta_SX_init)
          print("Ran Gibbs sampler")
          
          selected_samples_beta<-post_samples$samples_beta[seq((round(n_iter/5)+1), n_iter, 1)]
          selected_samples_beta_SX<-post_samples$samples_beta_SX[seq((round(n_iter/5)+1), n_iter, 1)]
          est_beta <- mean(selected_samples_beta)
          est_beta_SX <- mean(selected_samples_beta_SX)
          est_beta_combined <- est_beta + 1 * est_beta_SX
          pval_beta<-2*(1-pnorm(abs(mean(selected_samples_beta))/sd(selected_samples_beta)))
          pval_beta_SX<-2*(1-pnorm(abs(mean(selected_samples_beta_SX))/sd(selected_samples_beta_SX)))
          pval_beta_combined <- 2*(1-pnorm(abs(est_beta + 1 * est_beta_SX)/sd(selected_samples_beta_SX+selected_samples_beta)))
          print("Calculated result")
          
          res <- rbind(res, c(exposure, str_match(input, "^\\w+_\\w+_(\\w+)")[2], str_match(input, "^\\w+_[:alpha:]+([:digit:]+)")[2], "IVW", "3", length(dfl$rsid), 
                              est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, mr_ivw_mix$Estimate, p_thr, n_iter, kb, r2))
          print(c(exposure, str_match(input, "^\\w+_\\w+_(\\w+)")[2], str_match(input, "^\\w+_[:alpha:]+([:digit:]+)")[2], "IVW", "3", length(dfl$rsid),
                  est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, mr_ivw_mix$Estimate, p_thr, n_iter, kb, r2))
          
          mr_raps_mix <- mr.raps(b_exp = dfl$gammah, b_out = dfl$Gammah, se_exp = dfl$seg, se_out = dfl$seG)
          mr_raps_single1 <- mr.raps(b_exp = dfl$gammah, b_out = dfl$Gammahs1, se_exp = dfl$seg, se_out = dfl$se1)
          mr_raps_single2 <- mr.raps(b_exp = dfl$gammah, b_out = dfl$Gammahs2, se_exp = dfl$seg, se_out = dfl$se2)
          
          beta_init <- mr_raps_single1$beta.hat # Single is male
          beta_SX_init <- mr_raps_single2$beta.hat - mr_raps_single1$beta.hat # Mixed is female
          print("Initialized priors (MR-RAPS)")
          
          post_samples <- gibbs_sampling_cpp(n_iter = n_iter, n_thin = 5, p = length(dfl$rsid), Shat_Gamma_mix1_sq = (dfl$seG)^2,
                                             Shat_Gamma_mix2_sq = (dfl$se1)^2, Shat_Gamma_mix3_sq = (dfl$se2)^2, Shat_gamma_sq = (dfl$seg)^2,
                                             hat_Gamma_mix1 = dfl$Gammah, hat_Gamma_mix2 = dfl$Gammahs1, hat_Gamma_mix3 = dfl$Gammahs2, hat_gamma = dfl$gammah,
                                             1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1, 1, 1, 1, 1, 1,
                                             prob1 = 111897/225534, prob2 = 0, prob3 = 1, beta_init = beta_init, beta_SX_init = beta_SX_init)
          print("Ran Gibbs sampler")
          
          selected_samples_beta<-post_samples$samples_beta[seq((round(n_iter/5)+1), n_iter, 1)]
          selected_samples_beta_SX<-post_samples$samples_beta_SX[seq((round(n_iter/5)+1), n_iter, 1)]
          est_beta <- mean(selected_samples_beta)
          est_beta_SX <- mean(selected_samples_beta_SX)
          est_beta_combined <- est_beta + 1 * est_beta_SX
          pval_beta<-2*(1-pnorm(abs(mean(selected_samples_beta))/sd(selected_samples_beta)))
          pval_beta_SX<-2*(1-pnorm(abs(mean(selected_samples_beta_SX))/sd(selected_samples_beta_SX)))
          pval_beta_combined <- 2*(1-pnorm(abs(est_beta + 1 * est_beta_SX)/sd(selected_samples_beta_SX+selected_samples_beta)))
          print("Calculated result")
          
          res <- rbind(res, c(exposure, str_match(input, "^\\w+_\\w+_(\\w+)")[2], str_match(input, "^\\w+_[:alpha:]+([:digit:]+)")[2], "MR-RAPS", "3", length(dfl$rsid),
                              est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, mr_raps_mix$beta.hat, p_thr, n_iter, kb, r2))
          print(c(exposure, str_match(input, "^\\w+_\\w+_(\\w+)")[2], str_match(input, "^\\w+_[:alpha:]+([:digit:]+)")[2], "MR-RAPS", "3", length(dfl$rsid),
                  est_beta, est_beta_SX, est_beta_combined, pval_beta, pval_beta_SX, pval_beta_combined, beta_init, beta_SX_init, mr_raps_mix$beta.hat, p_thr, n_iter, kb, r2))
        }
      }
    }
  }
}

save(res, file = "res3sampleCUD.RData")
write_csv(as.data.frame(res), "res3sampleCUD.csv")

