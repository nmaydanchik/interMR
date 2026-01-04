library(tidyverse)
library(data.table)
setwd("D://RHDD/")


setwd("D://RHDD/BW/")
prepro <- fread("BW_EGG3.txt.gz")
postpro <- prepro %>% select(c(3,4,5,7,8,9)) %>% 
  filter(str_detect(prepro$rsid, "^rs")) %>% 
  drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_BW2016_Horikoshi.gz", compress = "gzip")

prepro <- fread("BW_EGG2.txt.gz")
postpro <- prepro %>% filter(str_detect(prepro$SNP, "^rs")) %>% drop_na()
postpro$EFFECT_ALLELE <- toupper(postpro$EFFECT_ALLELE)
postpro$OTHER_ALLELE <- toupper(postpro$OTHER_ALLELE)
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_BW2013_Horikoshi.gz", compress = "gzip")

prepro <- fread("BW2017Neale.gz")
postpro <- prepro %>% select(c(2,3,4,5,6,7)) %>% relocate(variant_id, effect_allele)
postpro <- postpro %>% filter(str_detect(postpro$variant_id, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_BW2017_Neale.gz", compress = "gzip")

prepro <- fread("BW2018Elsworth.gz")
postpro <- prepro %>% select(c(2,3,4,5,6,7)) %>% relocate(variant_id, effect_allele)
postpro <- postpro %>% filter(str_detect(postpro$variant_id, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_BW2018_Elsworth.gz", compress = "gzip")

prepro <- read_tsv("BW/BW_UKB.bgz")

setwd("D://RHDD/ADHD/")
prepro <- fread("ADHD2018_male.gz")
postpro <- prepro %>% select(c(1,4,5,6,7,8)) %>% filter(str_detect(prepro$SNP, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_ADHD2018_Martin_male.gz", compress = "gzip")

prepro <- fread("ADHD2018_female.gz")
postpro <- prepro %>% select(c(1,4,5,6,7,8)) %>% filter(str_detect(prepro$SNP, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_ADHD2018_Martin_female.gz", compress = "gzip")

prepro <- fread("ADHD2017_IEU.gz")
postpro <- prepro %>% select(c(2,3,4,5,6,7)) %>% relocate(variant_id, effect_allele)
postpro <- postpro %>% filter(str_detect(postpro$variant_id, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_ADHD2017_DemontisIEU.gz", compress = "gzip")

prepro <- fread("ADHD2017_iPSYCH.gz")
postpro <- prepro %>% select(c(2,4,5))

prepro <- fread("ADHD2022.gz")
postpro <- prepro %>% select(2,4,5,9,10,11) %>% filter(str_detect(SNP, "^rs")) %>% drop_na()
postpro <- mutate(postpro, "beta" = log(OR)) %>% mutate("se" = abs(beta/qnorm(P/2))) %>% select(c(1,2,3,6,7,8))
postpro <- relocate(postpro, P, .after = se) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_ADHD2022_Demontis.gz", compress = "gzip")

prepro <- fread("GCST90275136.tsv")
postpro <- prepro %>% select(c(3,4,5,6,8,10)) %>% relocate(rs_id)
postpro <- postpro %>% filter(str_detect(postpro$rs_id, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_ADHD2023_Pedersen.gz", compress = "gzip")

setwd("D://RHDD/ASD/")
prepro <- fread("Grove2017.gz")
postpro <- prepro %>% select(2,4,5,7,8,9) %>% filter(str_detect(SNP, "^rs")) %>% drop_na()
postpro <- mutate(postpro, "beta" = log(OR))
postpro <- mutate(postpro, "se" = abs(beta/qnorm(P/2)))
postpro <- select(postpro, c(1,2,3,6,7,8))
postpro <- relocate(postpro, P, .after = se)
postpro <- drop_na(postpro)
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_ASD2017_Grove.gz", compress = "gzip")

prepro <- fread("PedersenASD.tsv")
postpro <- prepro %>% select(c(3,4,5,6,8,10)) %>% relocate(rs_id)
postpro <- postpro %>% filter(str_detect(postpro$rs_id, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_ASD2023_Pedersen.gz", compress = "gzip")

prepro <- fread("PedersenASD_ADuLT.tsv")
postpro <- prepro %>% select(c(3,4,5,6,8,10)) %>% relocate(rs_id)
postpro <- postpro %>% filter(str_detect(postpro$rs_id, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_ASD2023_PedersenADuLT.gz", compress = "gzip")

prepro <- fread("ZhouTopMed.gz")
prepro <- select(prepro, c(3,4,5,8,9,10)) %>% drop_na()
prepro <- filter(prepro, rsids != "")

prepro <- fread("ZhouSaige.gz")

setwd("D://RHDD/CUD/")
prepro <- fread("CUDDemontis2019.gz")
postpro <- prepro %>% select(2,4,5,7,8,9) %>% filter(str_detect(SNP, "^rs")) %>% drop_na()
postpro <- mutate(postpro, "beta" = log(OR)) %>% mutate("se" = abs(beta/qnorm(P/2))) %>% select(c(1,2,3,6,7,8))
postpro <- relocate(postpro, P, .after = se) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_CUD2019_Demontis.gz", compress = "gzip")

prepro <- fread("CUDJohnson2020.gz")
postpro <- prepro %>% select(2,4,5,6,7,8) %>% filter(str_detect(SNP, "^rs")) %>% drop_na()
colnames(postpro) <- c("rsid", "eff", "noneff", "beta", "se", "pval")
fwrite(postpro, file = "processed_CUD2020_Johnson.gz", compress = "gzip")
