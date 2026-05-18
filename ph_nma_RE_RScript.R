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

ph_nma_re_model <- "
model {
  ##########################################################################
  # LIKELIHOOD & RANDOM EFFECTS — all studies
  ##########################################################################
  for(i in 1:NS){
    # Arm k=2 (all studies)
    y[i,1] ~ dnorm(delta[i,2], prec[i,1])
    delta[i,2] ~ dnorm(md[i,2], precd[i,2])
    md[i,2]    <- d[t[i,2]] - d[t[i,1]]
    precd[i,2] <- 2 / tau2           # k=2: 2*(k-1)/k * tau^2 = 1*tau^2, precision = 1/tau^2
                                      # Note: JAGS uses precision so 2(k-1)/k*tau^2 -> k/(2(k-1)) * 1/tau^2
    w[i,2]     <- delta[i,2] - (d[t[i,2]] - d[t[i,1]])

    # Arm k=3 (3-arm study only; for 2-arm studies na[i]<3 so loop does not execute)
    for(k in 3:na[i]){
      y[i,k-1] ~ dnorm(delta[i,k], prec[i,k-1])
      delta[i,k] ~ dnorm(md[i,k], precd[i,k])
      md[i,k]    <- d[t[i,k]] - d[t[i,1]] + sw[i,k]
      precd[i,k] <- 2*(k-1) / (k * tau2)   # precision = k / (2*(k-1)*tau^2)
      w[i,k]     <- delta[i,k] - (d[t[i,k]] - d[t[i,1]])
      sw[i,k]    <- sum(w[i,2:(k-1)]) / (k-1)
    }
  }

  ##########################################################################
  # RESIDUAL DEVIANCE
  ##########################################################################
  for(i in 1:NS){
    for(k in 1:(na[i]-1)){
      dev[i,k] <- (y[i,k] - delta[i,k+1])^2 * prec[i,k]
    }
    resdev[i] <- sum(dev[i, 1:(na[i]-1)])
  }
  totresdev <- sum(resdev[])

  ##########################################################################
  # PRIORS
  ##########################################################################
  d[1] <- 0
  for(k in 2:NT){
    d[k] ~ dnorm(0, 0.01)    # N(0, 100) on variance scale = precision 0.01
  }

  tau  ~ dunif(0, 1)
  tau2 <- tau * tau

  ##########################################################################
  # TRANSFORMATIONS
  ##########################################################################
  for(k in 1:NT){
    loghr[k] <- d[k]
    hr[k]    <- exp(d[k])
  }

  for(k in 1:NT){
    for(j in 1:NT){
      better[j,k] <- step(loghr[k] - loghr[j] - 0.0001)
    }
    rk[k]   <- 1 + sum(better[,k])
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

# y is NS * (max_arms-1): log-HR for each contrast within each study
# prec is NS * (max_arms-1): precision for each contrast

y_re    <- matrix(NA, nrow = NS, ncol = max_arms - 1)
prec_re <- matrix(NA, nrow = NS, ncol = max_arms - 1)

for(i in 1:NS){
  coefs <- coef(fit.cox.ph[[i]])
  vc    <- vcov(fit.cox.ph[[i]])
  
  if(i != 13){
    y_re[i, 1]    <- coefs["trt"]
    prec_re[i, 1] <- 1 / vc["trt", "trt"]
  } else {
    y_re[i, 1]    <- coefs["trt2"]
    y_re[i, 2]    <- coefs["trt3"]
    prec_re[i, 1] <- 1 / vc["trt2", "trt2"]
    prec_re[i, 2] <- 1 / vc["trt3", "trt3"]
    # Off-diagonal correlation handled via the w[i,k] correction
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# JAGS DATA LIST
# ══════════════════════════════════════════════════════════════════════════════

jags_data_re <- list(
  NS   = NS,
  NT   = NT,
  t    = t_mat,
  na   = na,
  y    = y_re,
  prec = prec_re
)

# ══════════════════════════════════════════════════════════════════════════════
# INITIAL VALUES - 3 CHAINS
# ══════════════════════════════════════════════════════════════════════════════
delta_init <- matrix(NA, nrow = NS, ncol = max_arms)
delta_init[, 2] <- 0        # arm 2 for all studies
delta_init[NS, 3] <- 0      # arm 3 for 3-arm study only

inits1 <- list(
  d     = c(NA, rep(0, NT - 1)),
  delta = delta_init,
  tau   = 0.5
)

inits2 <- list(
  d     = c(NA, rep(0.1, NT - 1)),
  delta = replace(delta_init, !is.na(delta_init), 0.1),
  tau   = 0.25
)

inits3 <- list(
  d     = c(NA, rep(-0.1, NT - 1)),
  delta = replace(delta_init, !is.na(delta_init), -0.1),
  tau   = 0.1
)

inits_re <- list(inits1, inits2, inits3)

# ══════════════════════════════════════════════════════════════════════════════
# SET PARAMETERS TO MONITOR
# ══════════════════════════════════════════════════════════════════════════════
params_re <- c("hr", "rk", "best", "tau")

# ══════════════════════════════════════════════════════════════════════════════
# COMPILE AND RUN JAGS MODEL
# ══════════════════════════════════════════════════════════════════════════════
fit_ph_re <- jags.model(
  textConnection(ph_nma_re_model),
  data    = jags_data_re,
  inits   = inits_re,
  n.chains = 3,
  n.adapt  = 10000
)

# Burn-in
update(fit_ph_re, 50000)

# Sample
samples_re <- coda.samples(
  fit_ph_re,
  variable.names = params_re,
  n.iter = 50000,
  thin   = 5
)

# ══════════════════════════════════════════════════════════════════════════════
# RESULTS
# ══════════════════════════════════════════════════════════════════════════════

summary(samples_re)

###############
# End of file #
###############