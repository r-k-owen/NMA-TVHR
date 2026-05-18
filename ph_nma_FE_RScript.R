# Author: Rhiannon K Owen #
# Date: 18/05/2026 #
# R version: 4.5.0 #

#############################################
# Clear environment & set working directory #
#############################################
rm(list = ls())

setwd("wd")

# ══════════════════════════════════════════════════════════════════════════════
# JAGS MODEL - PHHR NMA Fixed Effect
# ══════════════════════════════════════════════════════════════════════════════

ph_nma_model <- "
model {
  ##########################################################################
  # LIKELIHOOD — 2-arm studies
  ##########################################################################
  for(i in 1:NS2){
    y[i] ~ dnorm(delta2[i], prec[i])
    delta2[i] <- d[t[i,2]] - d[t[i,1]]
  }

  ##########################################################################
  # LIKELIHOOD — 3-arm study (CheckMate067, study NS)
  ##########################################################################
  y_joint3[1:2] ~ dmnorm(delta3[1:2], prec_joint3[1:2, 1:2])
  delta3[1] <- d[t[NS,2]] - d[t[NS,1]]
  delta3[2] <- d[t[NS,3]] - d[t[NS,1]]

  ##########################################################################
  # DEVIANCE — 2-arm studies
  ##########################################################################
  for(i in 1:NS2){
    dev2[i] <- (y[i] - delta2[i])^2 * prec[i]
  }

  ##########################################################################
  # DEVIANCE — 3-arm study
  ##########################################################################
  dev_mat3[1:2] <- y_joint3[1:2] - delta3[1:2]
  Pdev3[1] <- prec_joint3[1,1] * dev_mat3[1] + prec_joint3[1,2] * dev_mat3[2]
  Pdev3[2] <- prec_joint3[2,1] * dev_mat3[1] + prec_joint3[2,2] * dev_mat3[2]
  dev3 <- inprod(dev_mat3[], Pdev3[])

  ##########################################################################
  # TOTAL RESIDUAL DEVIANCE
  ##########################################################################
  totresdev <- sum(dev2[1:NS2]) + dev3

  ##########################################################################
  # PRIORS
  ##########################################################################
  d[1] <- 0
  for(k in 2:NT){
    d[k] ~ dnorm(0, 0.001)
  }

  ##########################################################################
  # TRANSFORMATIONS — single HR per treatment (no time grid needed)
  ##########################################################################
  for(k in 1:NT){
    loghr[k] <- d[k]
    hr[k]    <- exp(d[k])
  }

  # Treatment rankings
  for(k in 1:NT){
    for(j in 1:NT){
      better[j, k] <- step(loghr[k] - loghr[j] - 0.0001)
    }
    rk[k]   <- 1 + sum(better[, k])
    best[k] <- equals(rk[k], 1)
  }
}
"

# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════
library(dplyr)
library(rjags)
library(coda)
library(surrosurv)
library(coxme)
library(survival)
library(MASS)
library(tidyverse)
library(ggplot2)
library(coda)

# ══════════════════════════════════════════════════════════════════════════════
# DATA PREP
# ══════════════════════════════════════════════════════════════════════════════
dat <- read.csv("melanoma_ipd.csv")

# Re-order studies to loop through 2-arm studies sequentially
dat$studyCode[dat$studyCode==4]<- 14
dat$studyCode[dat$studyCode==13]<- 4
dat$studyCode[dat$studyCode==14]<- 13
dat$trt<- dat$arm-1

# Assign treatment indicators
dat$trt1 <- ifelse(dat$trt == 0, 1, 0)
dat$trt2 <- ifelse(dat$trt == 1, 1, 0)
dat$trt3 <- ifelse(dat$trt == 2, 1, 0)

# ══════════════════════════════════════════════════════════════════════════════
# STUDY STRUCTURE
# ══════════════════════════════════════════════════════════════════════════════
# NS - number of studies
# NS2 - number of 2 arm studies
NS  <- 13
NS2 <- 12 

# Treatment codes 
# Treatment 1 - Dacarbazine
# Treatment 2 - Tremelimumab
# Treatment 3 - Ipilimumab + Dacarbazine
# Treatment 4 - Dabrafenib
# Treatment 5 - Dabrafenib + Trametinib
# Treatment 6 - Vemurafenib
# Treatment 7 - Vemurafenib + Cobimetinib
# Treatment 8 - Ipilimumab
# Treatment 9 - Nivolumab
# Treatment 10 - Nivolumab + Ipilimumab
# Treatment 11 - Pembrolizumab
# Treatment 12 - Ipilimumab + Sargramostin
# Treatment 13 - Selumetinib + Dacarbazine

study_trts <- list(
  c(1,  4),    # Study 1  - BREAK3
  c(1,  6),    # Study 2  - BRIM3
  c(1,  9),    # Study 3  - CheckMate066
  c(1, 13),    # Study 4  - Robert 2013
  c(8, 10),    # Study 5  - CheckMate069
  c(6,  7),    # Study 6  - coBRIM
  c(4,  5),    # Study 7  - COMBI-d
  c(5,  6),    # Study 8  - COMBI-v
  c(8, 12),    # Study 9  - Hodi 2014
  c(8, 11),    # Study 10 - Keynote 006
  c(1,  2),    # Study 11 - Ribas 2013
  c(1,  3),    # Study 12 - Robert 2011
  c(8,  9, 10) # Study 13 - CheckMate067 (3-arm)
)

# Number of arms
na       <- sapply(study_trts, length)
max_arms <- max(na)
# Number of treatments
NT       <- max(unlist(study_trts))

# Study specific treatment matrix
t_mat <- matrix(0, nrow = NS, ncol = max_arms)
for(i in 1:NS){
  t_mat[i, 1:na[i]] <- study_trts[[i]]
}

# ══════════════════════════════════════════════════════════════════════════════
# FIT COX MODELS - STANDARD PH, NO TIME INTERACTION
# ══════════════════════════════════════════════════════════════════════════════
fit.cox.ph <- vector("list", NS)

for(i in 1:NS){
  if(i != 13){   # 2-arm trials
    fit.cox.ph[[i]] <- coxph(
      Surv(time, event) ~ trt,
      data   = dat,
      subset = studyCode == i
    )
  } else {       # 3-arm trial (study 13)
    fit.cox.ph[[i]] <- coxph(
      Surv(time, event) ~ trt2 + trt3,
      data   = dat,
      subset = studyCode == i
    )
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# EXTRACT COEFFICIENTS AND BUILD JOINT ARRAYS FOR JAGS
# ══════════════════════════════════════════════════════════════════════════════

# y[i, c]    = log-HR for contrast c in study i
# prec[i, c] = precision (1/var) for 2-arm studies
y    <- matrix(NA, nrow = NS, ncol = max_arms - 1)
vmat <- matrix(NA, nrow = NS, ncol = max_arms - 1)

# Covariance matrix for 3-arm study (2x2), filled below
cov3 <- matrix(NA, nrow = 2, ncol = 2)

for(i in 1:NS){
  coefs <- coef(fit.cox.ph[[i]])
  vc    <- vcov(fit.cox.ph[[i]])
  
  if(i != 13){   # 2-arm trial — single contrast
    y[i, 1]    <- coefs["trt"]
    vmat[i, 1] <- vc["trt", "trt"]
    
  } else {       # 3-arm trial — two contrasts
    y[i, 1] <- coefs["trt2"]
    y[i, 2] <- coefs["trt3"]
    
    # Full 2x2 covariance matrix (off-diagonal captures shared reference arm)
    cov3[1, 1] <- vc["trt2", "trt2"]
    cov3[2, 2] <- vc["trt3", "trt3"]
    cov3[1, 2] <- vc["trt2", "trt3"]
    cov3[2, 1] <- vc["trt3", "trt2"]
  }
}

# Precision for 2-arm studies (scalar per study)
prec_2arm <- 1 / vmat[1:NS2, 1]

# Precision matrix for 3-arm study
prec_joint3 <- solve(cov3)

# ══════════════════════════════════════════════════════════════════════════════
# JAGS DATA LIST
# ══════════════════════════════════════════════════════════════════════════════
jags_data <- list(
  NS2         = NS2,
  NS          = NS,
  NT          = NT,
  t           = t_mat,
  y           = y[1:NS2, 1],       # log-HRs for 2-arm studies
  prec        = prec_2arm,          # precisions for 2-arm studies
  y_joint3    = y[NS, 1:2],         # log-HRs for 3-arm study
  prec_joint3 = prec_joint3         # 2x2 precision matrix for 3-arm study
)

# ══════════════════════════════════════════════════════════════════════════════
# INITIAL VALUES - 3 CHAINS
# ══════════════════════════════════════════════════════════════════════════════
inits1 <- 
  list(
    d = c(NA, rep(0, NT-1)))


inits2 <- 
  list(
    d = c(NA, rep(1, NT-1)))

inits3 <- 
  list(
    d = c(NA, rep(1, NT-1)))


# ══════════════════════════════════════════════════════════════════════════════
# SET PARAMETERS TO MONITOR
# ══════════════════════════════════════════════════════════════════════════════
params <- c(
  "hr",       
  "rk",        
  "best"      
)

# ══════════════════════════════════════════════════════════════════════════════
# COMPILE AND RUN JAGS MODEL
# ══════════════════════════════════════════════════════════════════════════════
fit_ph <- jags.model(file=textConnection(ph_nma_model),
                       data=jags_data,
                       n.chains=3,
                       inits = list(inits1,inits2,inits3),
                       n.adapt=500)

# Burn-in
update(fit_ph, 50000)

# Sample
params <- c("hr", "rk")
samples_ph <- coda.samples(fit_ph, variable.names=params, n.iter=50000, thin=5)

# ══════════════════════════════════════════════════════════════════════════════
# RESULTS
# ══════════════════════════════════════════════════════════════════════════════
summary(samples_ph)

samples_mat_ph <- as.matrix(samples_ph)

hr_cols <- grep("^hr\\[", colnames(samples_mat_ph), value = TRUE)

hr_summary_ph <- data.frame(
  treatment = 1:NT,
  mean      = apply(samples_mat_ph[, hr_cols], 2, mean),
  lower     = apply(samples_mat_ph[, hr_cols], 2, quantile, 0.025),
  upper     = apply(samples_mat_ph[, hr_cols], 2, quantile, 0.975),
  median    = apply(samples_mat_ph[, hr_cols], 2, median)
)

rk_cols <- grep("^rk\\[", colnames(samples_mat_ph), value = TRUE)

rk_summary_ph <- data.frame(
  treatment = 1:NT,
  mean      = apply(samples_mat_ph[, rk_cols], 2, mean),
  lower     = apply(samples_mat_ph[, rk_cols], 2, quantile, 0.025),
  upper     = apply(samples_mat_ph[, rk_cols], 2, quantile, 0.975),
  median    = apply(samples_mat_ph[, rk_cols], 2, median)
)

trt_labels <- c("1"="DTIC", 
                "2"="TR", 
                "3"="IPI+DTIC", 
                "4"="DB", 
                "5"="DB+TR", 
                "6"="VM",
                "7"="VM+COB", 
                "8"="IPI", 
                "9"="NIV", 
                "10"="NIV+IPI", 
                "11"="PEM",
                "12"="IPI+SRG", 
                "13"="SEL+DTIC")

hr_summary_ph$trt_label <- trt_labels[as.character(hr_summary_ph$treatment)]

# Exclude reference treatment 
hr_plot <- hr_summary_ph %>% filter(treatment != 1)

# Define a time range (for visual only — HR is constant under PH)
time_seq <- seq(0, 6, by = 1)   

# Expand so each treatment has a row per time point
hr_time <- hr_plot %>%
  tidyr::crossing(time = time_seq)

# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════

ggplot(hr_time, aes(x = time, y = mean)) +
  geom_line(colour = "steelblue", linewidth = 0.9) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "steelblue", alpha = 0.2) +
  geom_hline(yintercept = 1, 
             linetype = "dashed", colour = "red", linewidth = 0.5) +
  facet_wrap(~ trt_label, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6)) +
  labs(
    title    = "Time-Varying Hazard Ratio vs DTIC",
    subtitle = "Posterior mean with 95% credible interval",
    x        = "Time (years)",
    y        = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Re-labelled for combined plot
hr_time$treatment <- hr_time$trt_label

###############
# End of file #
###############

