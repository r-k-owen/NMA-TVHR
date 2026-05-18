# Author: Rhiannon K Owen #
# Date: 18/05/2026 #
# R version: 4.5.0 #

#############################################
# Clear environment & set working directory #
#############################################
rm(list = ls())

setwd("wd")

# ══════════════════════════════════════════════════════════════════════════════
# JAGS MODEL - TVHR NMA Fixed Effect
# ══════════════════════════════════════════════════════════════════════════════
tvhr_nma_fe_model <- "
model {
  ##########################################################################
  # LIKELIHOOD — 2-arm studies
  ##########################################################################
  for(i in 1:NS2){
    y_joint[i, 1:2] ~ dmnorm(delta2[i, 1:2], prec_joint[i, 1:2, 1:2])

    delta2[i, 1] <- d1[t[i,2]] - d1[t[i,1]]
    delta2[i, 2] <- d2[t[i,2]] - d2[t[i,1]]
  }

  ##########################################################################
  # LIKELIHOOD — 3-arm study (CheckMate067, study 13)
  ##########################################################################
  y_joint[NS, 1:4] ~ dmnorm(delta3[1:4], prec_joint[NS, 1:4, 1:4])

  delta3[1] <- d1[t[NS,2]] - d1[t[NS,1]]
  delta3[2] <- d2[t[NS,2]] - d2[t[NS,1]]
  delta3[3] <- d1[t[NS,3]] - d1[t[NS,1]]
  delta3[4] <- d2[t[NS,3]] - d2[t[NS,1]]

  ##########################################################################
  # DEVIANCE — 2-arm studies
  ##########################################################################
  for(i in 1:NS2){
    dev_mat2[i, 1:2] <- y_joint[i, 1:2] - delta2[i, 1:2]
    dev2[i] <- dev_mat2[i,1] * (prec_joint[i,1,1] * dev_mat2[i,1] +
                                  prec_joint[i,1,2] * dev_mat2[i,2]) +
               dev_mat2[i,2] * (prec_joint[i,2,1] * dev_mat2[i,1] +
                                  prec_joint[i,2,2] * dev_mat2[i,2])
  }

  ##########################################################################
  # DEVIANCE — 3-arm study
  ##########################################################################
  dev_mat3[1:4] <- y_joint[NS, 1:4] - delta3[1:4]
  Pdev3[1] <- prec_joint[NS,1,1] * dev_mat3[1] + prec_joint[NS,1,2] * dev_mat3[2] +
              prec_joint[NS,1,3] * dev_mat3[3] + prec_joint[NS,1,4] * dev_mat3[4]
  Pdev3[2] <- prec_joint[NS,2,1] * dev_mat3[1] + prec_joint[NS,2,2] * dev_mat3[2] +
              prec_joint[NS,2,3] * dev_mat3[3] + prec_joint[NS,2,4] * dev_mat3[4]
  Pdev3[3] <- prec_joint[NS,3,1] * dev_mat3[1] + prec_joint[NS,3,2] * dev_mat3[2] +
              prec_joint[NS,3,3] * dev_mat3[3] + prec_joint[NS,3,4] * dev_mat3[4]
  Pdev3[4] <- prec_joint[NS,4,1] * dev_mat3[1] + prec_joint[NS,4,2] * dev_mat3[2] +
              prec_joint[NS,4,3] * dev_mat3[3] + prec_joint[NS,4,4] * dev_mat3[4]
  dev3 <- inprod(dev_mat3[], Pdev3[])

  ##########################################################################
  # TOTAL RESIDUAL DEVIANCE
  ##########################################################################
  totresdev <- sum(dev2[1:NS2]) + dev3

  ##########################################################################
  # PRIORS
  ##########################################################################
  d1[1] <- 0
  d2[1] <- 0
  for(k in 2:NT){
    d1[k] ~ dnorm(0, 0.001)
    d2[k] ~ dnorm(0, 0.001)
  }

  ##########################################################################
  # TRANSFORMATIONS
  ##########################################################################
  for(tg in 1:200){
    for(k in 1:NT){
      loghr[tg, k] <- d1[k] + d2[k] * lt[tg]
      hr[tg, k]    <- exp(loghr[tg, k])
    }
    for(k in 1:NT){
      for(j in 1:NT){
        better[tg, j, k] <- step(loghr[tg,k] - loghr[tg,j] - 0.0001)
      }
      rk[tg, k]   <- 1 + sum(better[tg, ,k])
      best[tg, k] <- equals(rk[tg, k], 1)
    }
  }
}
"

# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════
library(survival)
library(rjags)
library(coda)
library(dplyr)
library(tidyverse)

# ══════════════════════════════════════════════════════════════════════════════
# DATA PREP
# ══════════════════════════════════════════════════════════════════════════════
dat <- read.csv("melanoma_ipd.csv")

# Re-order studies to loop through 2-arm studies sequentially
dat$studyCode[dat$studyCode == 4]  <- 14
dat$studyCode[dat$studyCode == 13] <- 4
dat$studyCode[dat$studyCode == 14] <- 13

# Assign treatment indicators
dat$trt  <- dat$arm - 1
dat$trt1 <- ifelse(dat$trt == 0, 1, 0)
dat$trt2 <- ifelse(dat$trt == 1, 1, 0)
dat$trt3 <- ifelse(dat$trt == 2, 1, 0)

# Time scale in years
dat$time_years <- dat$time / 12

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
# FIT COX MODELS WITH TIME INTERACTION
# ══════════════════════════════════════════════════════════════════════════════
fit.cox.tt <- vector("list", NS)

# Cox models for two arm studies
for(i in 1:NS){
  if(i != 13){
    fit.cox.tt[[i]] <- coxph(
      Surv(time_years, event) ~ trt + tt(trt),
      data   = dat,
      subset = studyCode == i,
      tt     = function(x, t, ...) x * log(t)
    )
  } 
# Cox models for three arm studies
  else {
    fit.cox.tt[[i]] <- coxph(
      Surv(time_years, event) ~ trt2 + trt3 + tt(trt2) + tt(trt3),
      data   = dat,
      subset = studyCode == i,
      tt     = function(x, t, ...) x * log(t)
    )
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# EXTRACT COEFFICIENTS AND BUILD JOINT ARRAYS FOR JAGS
# ══════════════════════════════════════════════════════════════════════════════
n_param_max <- 4   # 2 params x 2 contrasts (i.e., 3-arm trials), 2 params for 2 arm trials

y_joint        <- matrix(NA, nrow = NS, ncol = n_param_max)
prec_joint_arr <- array(NA,  dim  = c(NS, n_param_max, n_param_max))

for(i in 1:NS){
  coefs <- coef(fit.cox.tt[[i]])
  vc    <- vcov(fit.cox.tt[[i]])
  
  if(i != 13){
    # 2-arm: 2 parameters in slots 1:2
    y_joint[i, 1:2]             <- coefs[c("trt", "tt(trt)")]
    prec_joint_arr[i, 1:2, 1:2] <- solve(vc[c("trt","tt(trt)"), c("trt","tt(trt)")])
    
  } else {
    # 3-arm: 4 parameters in slots 1:4
    params_s13               <- c("trt2", "tt(trt2)", "trt3", "tt(trt3)")
    vc_full                  <- vc[params_s13, params_s13]
    y_joint[i, 1:4]             <- coefs[params_s13]
    prec_joint_arr[i, 1:4, 1:4] <- solve(vc_full)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# TIME GRID IN YEARS
# ══════════════════════════════════════════════════════════════════════════════
time_seq <- seq(0.1, 6, length.out = 200)
lt       <- log(time_seq)

# ══════════════════════════════════════════════════════════════════════════════
# JAGS DATA LIST
# ══════════════════════════════════════════════════════════════════════════════
jags_data <- list(
  NS         = NS,
  NS2        = NS2,
  NT         = NT,
  t          = t_mat,
  y_joint    = y_joint,
  prec_joint = prec_joint_arr,
  lt         = lt
)

# ══════════════════════════════════════════════════════════════════════════════
# INITIAL VALUES - 3 CHAINS
# ══════════════════════════════════════════════════════════════════════════════
inits1 <- list(
  d1 = c(NA, rep( 0.0, NT - 1)),
  d2 = c(NA, rep( 0.0, NT - 1))
)

inits2 <- list(
  d1 = c(NA, rep( 0.1, NT - 1)),
  d2 = c(NA, rep( 0.1, NT - 1))
)

inits3 <- list(
  d1 = c(NA, rep(-0.1, NT - 1)),
  d2 = c(NA, rep(-0.1, NT - 1))
)

inits <- list(inits1, inits2, inits3)

# ══════════════════════════════════════════════════════════════════════════════
# SET PARAMETERS TO MONITOR
# ══════════════════════════════════════════════════════════════════════════════
params <- c(
  "hr", "loghr",
  "rk", "best",
  "totresdev"
  #"d1", "d2"
)

# ══════════════════════════════════════════════════════════════════════════════
# COMPILE AND RUN JAGS MODEL
# ══════════════════════════════════════════════════════════════════════════════
fit_tvhr_fe <- jags.model(
  textConnection(tvhr_nma_fe_model),
  data     = jags_data,
  inits    = inits,
  n.chains = 3,
  n.adapt  = 10000
)

update(fit_tvhr_fe, n.iter = 50000)

samples_tvhr_fe <- coda.samples(
  fit_tvhr_fe,
  variable.names = params,
  n.iter = 50000,
  thin   = 5
)

# ══════════════════════════════════════════════════════════════════════════════
# RESULTS
# ══════════════════════════════════════════════════════════════════════════════
summary(samples_tvhr_fe)

samples_mat <- as.matrix(samples_tvhr_fe)

# Total residual deviance
cat("Total residual deviance:", round(mean(samples_mat[, "totresdev"]), 2), "\n")
cat("Number of data points:  ", NS2 * 2 + 4, "\n")   # 12*2 two-arm + 1*4 three-arm

# ══════════════════════════════════════════════════════════════════════════════
# CONVERGENCE DIAGNOSTICS
# ══════════════════════════════════════════════════════════════════════════════
 gelman.diag(samples_tvhr_fe, multivariate = FALSE)
 effectiveSize(samples_tvhr_fe)

# Extract individual parameters at specific time points
# t = 166, year 5
 params_of_interest <- grep("hr\\[166,", varnames(samples_tvhr_fe), value = TRUE)
 #params_of_interest <- grep("rk\\[166,", varnames(samples_tvhr_fe), value = TRUE)
 #params_of_interest <- grep("best\\[166,", varnames(samples_tvhr_fe), value = TRUE)
 
 samples_subset <- samples_tvhr_fe[, params_of_interest]
 autocorr.plot(samples_subset)
 traceplot(samples_subset)
 densplot(samples_subset)

# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════

d.mcmc <- as.matrix(samples_tvhr_fe)

T  <- 200
NT <- 13

# P(best) 
p_best <- matrix(NA, nrow = T, ncol = NT)

for(tg in 1:T){
  for(k in 1:NT){
    col_name      <- paste0("best[", tg, ",", k, "]")
    p_best[tg, k] <- mean(d.mcmc[, col_name])
  }
}

# Mean rank over time 
mean_rank <- matrix(NA, nrow = T, ncol = NT)

for(tg in 1:T){
  for(k in 1:NT){
    col_name          <- paste0("rk[", tg, ",", k, "]")
    mean_rank[tg, k]  <- mean(d.mcmc[, col_name])
  }
}

# SUCRA
rank_probs <- array(0, dim = c(T, NT, NT))

for(tg in 1:T){
  for(k in 1:NT){
    col_name <- paste0("rk[", tg, ",", k, "]")
    rk_samples <- d.mcmc[, col_name]
    for(r in 1:NT){
      rank_probs[tg, k, r] <- mean(rk_samples == r)
    }
  }
}

sucra <- matrix(NA, nrow = T, ncol = NT)
for(tg in 1:T){
  for(k in 1:NT){
    cum_probs      <- cumsum(rank_probs[tg, k, ])
    sucra[tg, k]   <- sum(cum_probs[1:(NT-1)]) / (NT - 1)
  }
}

# Extract HR samples into array [iterations, time, treatment]
n_iter <- nrow(d.mcmc)

hr_mean <- matrix(NA, nrow = T, ncol = NT)
hr_low  <- matrix(NA, nrow = T, ncol = NT)
hr_high <- matrix(NA, nrow = T, ncol = NT)

for(tg in 1:T){
  for(k in 1:NT){
    col_name        <- paste0("hr[", tg, ",", k, "]")
    hr_samples      <- d.mcmc[, col_name]
    hr_mean[tg, k]  <- mean(hr_samples)
    hr_low[tg, k]   <- quantile(hr_samples, 0.025)
    hr_high[tg, k]  <- quantile(hr_samples, 0.975)
  }
}

trt_labels <- c("DTIC", "TR", "IPI+DTIC", "DB", "DB+TR", "VM",
                "VM+COB", "IPI", "NIV", "NIV+IPI", "PEM",
                "IPI+SRG", "SEL+DTIC")

n_time<-200
time_vec <- seq(0.1, 6, length.out = 200)


hr_df <- expand.grid(time = time_vec, treatment = trt_labels) %>%
  arrange(treatment, time) %>%
  mutate(
    tg   = rep(1:n_time, times = NT),
    k    = rep(1:NT, each  = n_time),
    mean = as.vector(hr_mean),
    low  = as.vector(hr_low),
    high = as.vector(hr_high)
  )

# Remove reference treatment
hr_df_noref <- hr_df %>% filter(treatment != "DTIC")

# Subset at landmark times (years 1, 2 and 5)
hr_yr5<-subset(hr_df_noref, hr_df_noref$tg==166)
hr_yr2<-subset(hr_df_noref, hr_df_noref$tg==65)
hr_yr1<-subset(hr_df_noref, hr_df_noref$tg==31)
hr_yr0.5<-subset(hr_df_noref, hr_df_noref$tg==14)


#pdf(file = "TVHR_NMA_FE.pdf")
ggplot(hr_df_noref, aes(x = time, y = mean)) +
  geom_line(colour = "steelblue", linewidth = 0.9) +
  geom_ribbon(aes(ymin = low, ymax = high),
              fill = "steelblue", alpha = 0.2) +
  geom_hline(yintercept = 1, 
             linetype = "dashed", colour = "red", linewidth = 0.5) +
  facet_wrap(~ treatment, ncol = 4, scales = "free_y") +
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
#dev.off()

# Key treatments plot
key_trts <- c("DB+TR", "VM+COB", "NIV", "NIV+IPI", "PEM")

hr_df_noref %>%
  filter(treatment %in% key_trts) %>%
  ggplot(aes(x = time, y = mean, 
             colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = low, ymax = high),
              alpha = 0.15, colour = NA) +
  geom_hline(yintercept = 1, 
             linetype = "dashed", colour = "black") +
  scale_colour_manual(values = c(
    "DB+TR"   = "#E69F00",
    "VM+COB"  = "#D55E00",
    "NIV"     = "#0072B2",
    "NIV+IPI" = "#009E73",
    "PEM"     = "#CC79A7"
  )) +
  scale_fill_manual(values = c(
    "DB+TR"   = "#E69F00",
    "VM+COB"  = "#D55E00",
    "NIV"     = "#0072B2",
    "NIV+IPI" = "#009E73",
    "PEM"     = "#CC79A7"
  )) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(limits = c(0, 1.5)) +
  labs(
    title    = "Time-Varying Hazard Ratio — Key Treatments vs DTIC",
    subtitle = "Posterior mean with 95% credible interval",
    x        = "Time (years)",
    y        = "Hazard Ratio",
    colour   = "Treatment",
    fill     = "Treatment"
  ) +
  theme_minimal(base_size = 12)


trt_labels <- c("DTIC", "TR", "IPI+DTIC", "DB", "DB+TR", "VM",
                "VM+COB", "IPI", "NIV", "NIV+IPI", "PEM",
                "IPI+SRG", "SEL+DTIC")

time <- seq(0.1, 6, len = 200)
T    <- 200
NT   <- 13

# P(best) 
pbest_df <- as.data.frame(p_best)
colnames(pbest_df) <- trt_labels
pbest_df$time <- time

pbest_long <- pbest_df %>%
  pivot_longer(-time, names_to = "treatment", values_to = "p_best") %>%
  mutate(treatment = factor(treatment, levels = trt_labels))

# Mean Rank 
meanrk_df <- as.data.frame(mean_rank)
colnames(meanrk_df) <- trt_labels
meanrk_df$time <- time

meanrk_long <- meanrk_df %>%
  pivot_longer(-time, names_to = "treatment", values_to = "mean_rank") %>%
  mutate(treatment = factor(treatment, levels = trt_labels))

# SUCRA
sucra_df <- as.data.frame(sucra)
colnames(sucra_df) <- trt_labels
sucra_df$time <- time

sucra_long <- sucra_df %>%
  pivot_longer(-time, names_to = "treatment", values_to = "sucra") %>%
  mutate(treatment = factor(treatment, levels = trt_labels))

# Colour palette 
trt_colours <- c(
  "DTIC"     = "#999999", 
  "TR"       = "#BBBBBB",  
  "IPI+DTIC" = "#DDDDDD",   
  "DB"       = "#E69F00",    
  "DB+TR"    = "#D55E00",   
  "VM"       = "#F0E442",    
  "VM+COB"   = "#AE123E",    
  "IPI"      = "#56B4E9",   
  "NIV"      = "#0072B2",    
  "NIV+IPI"  = "#009E73",    
  "PEM"      = "#006040",    
  "IPI+SRG"  = "#2F4B7C",    
  "SEL+DTIC" = "#D45087"     
)

pbest_nma<- ggplot(pbest_long %>% filter(treatment != "DTIC"),
                   aes(x = time, y = p_best, 
                       colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = trt_colours) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Probability of Being Best Treatment Over Time",
    x        = "Time (years)",
    y        = "P(best)",
    colour   = "Treatment"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

 ggsave("pbest_NMA_FE.pdf",
        pbest_nma,
        width = 12,
        height = 10,
        dpi = 300)

ggplot(pbest_long,
       aes(x = time, y = p_best)) +
  geom_line(aes(colour = treatment), linewidth = 0.9) +
  scale_colour_manual(values = trt_colours) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ treatment, ncol = 4) +
  labs(
    title    = "Probability of Being Best Treatment Over Time",
    x        = "Time (years)",
    y        = "P(best)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text      = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggplot(meanrk_long,
       aes(x = time, y = mean_rank,
           colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = trt_colours) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(breaks = 1:NT,
                     limits = c(1, NT),
                     trans  = "reverse") +   # rank 1 at top
  geom_hline(yintercept = (NT + 1) / 2, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  annotate("text", x = 5.5, y = (NT + 1) / 2 - 0.3,
           label = "Average\nrank", size = 3, colour = "grey50") +
  labs(
    title    = "Mean Rank Over Time",
    subtitle = "Reversed y-axis - best treatments appear at top (Rank 1 = best)",
    x        = "Time (years)",
    y        = "Mean Rank",
    colour   = "Treatment"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggplot(meanrk_long,
       aes(x = time, y = mean_rank)) +
  geom_line(aes(colour = treatment), linewidth = 0.9) +
  geom_hline(yintercept = (NT + 1) / 2, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(values = trt_colours) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6)) +
  scale_y_continuous(breaks = c(1, 4, 7, 10, 13),
                     limits = c(1, NT),
                     trans  = "reverse") +
  facet_wrap(~ treatment, ncol = 4) +
  labs(
    title    = "Mean Rank Over Time",
    subtitle = "(Rank 1 = best)",
    x        = "Time (years)",
    y        = "Mean Rank"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

ggplot(sucra_long %>% filter(treatment != "DTIC"),
       aes(x = time, y = sucra,
           colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = trt_colours) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  labs(
    title    = "SUCRA Over Time",
    subtitle = "Values above 0.5 indicate better than average rank",
    x        = "Time (years)",
    y        = "SUCRA",
    colour   = "Treatment"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

landmark_years <- c(0.5, 1, 2, 5) 
landmark_idx <- sapply(landmark_years, function(yr) which.min(abs(time - yr)))
data.frame(year = landmark_years, idx = landmark_idx, time_at_idx = time[landmark_idx])

# Extract landmark values
landmark_df <- map_dfr(seq_along(landmark_years), function(j){
  tg <- landmark_idx[j]
  data.frame(
    year      = landmark_years[j],
    treatment = trt_labels,
    p_best    = p_best[tg, ],
    mean_rank = mean_rank[tg, ],
    sucra     = sucra[tg, ]
  )
}) %>%
  mutate(
    treatment = factor(treatment, levels = trt_labels),
    year      = factor(paste("Year", year))
  ) %>%
  filter(treatment != "DTIC")

# Heatmap of p(best) at landmark times
heatmap<- ggplot(landmark_df,
       aes(x = year, y = reorder(treatment, p_best),
           fill = p_best)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = scales::percent(p_best, accuracy = 1)),
            size = 3.5, colour = "black") +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = 1/NT,
    limits   = c(0, 1),
    labels   = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title    = "P(Best) at Landmark Time Points",
    subtitle = "Treatments ordered by average P(best)",
    x        = "Time Point",
    y        = "Treatment",
    fill     = "P(best)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_text(face = "bold")
  )

ggsave("heatmap_pbest_NMA_FE.pdf",
       heatmap,
       width = 12,
       height = 10,
       dpi = 300)

# Compute rank probability distribution over time for key treatments
key_trts <- c("VM+COB", "NIV+IPI", "NIV", "PEM", "DB+TR")

rank_ribbon_df <- map_dfr(key_trts, function(trt){
  k <- which(trt_labels == trt)
  map_dfr(1:T, function(tg){
    rk_samp <- d.mcmc[, paste0("rk[", tg, ",", k, "]")]
    data.frame(
      time      = time[tg],
      treatment = trt,
      q25       = quantile(rk_samp, 0.25),
      q50       = quantile(rk_samp, 0.50),
      q75       = quantile(rk_samp, 0.75),
      q025      = quantile(rk_samp, 0.025),
      q975      = quantile(rk_samp, 0.975),
      mean_rank = mean(rk_samp)
    )
  })
})

ggplot(rank_ribbon_df,
       aes(x = time, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = q025, ymax = q975),
              alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              alpha = 0.25, colour = NA) +
  geom_line(aes(y = mean_rank), linewidth = 1) +
  scale_colour_manual(values = trt_colours) +
  scale_fill_manual(values = trt_colours) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(breaks = 1:NT,
                     limits = c(1, NT),
                     trans  = "reverse") +
  facet_wrap(~ treatment, ncol = 3) +
  labs(
    title    = "Rank Distribution Over Time — Select Treatments",
    subtitle = "Line = mean rank, inner ribbon = IQR, outer ribbon = 95% CrI",
    x        = "Time (years)",
    y        = "Rank (1 = best)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

# Rank data frame

rk_mean <- matrix(NA, nrow = T, ncol = NT)
rk_low  <- matrix(NA, nrow = T, ncol = NT)
rk_high <- matrix(NA, nrow = T, ncol = NT)

for(tg in 1:T){
  for(k in 1:NT){
    col_name        <- paste0("rk[", tg, ",", k, "]")
    rk_samples      <- d.mcmc[, col_name]
    rk_mean[tg, k]  <- mean(rk_samples)
    rk_low[tg, k]   <- quantile(rk_samples, 0.025)
    rk_high[tg, k]  <- quantile(rk_samples, 0.975)
  }
}

rk_df <- expand.grid(time = time_vec, treatment = trt_labels) %>%
  arrange(treatment, time) %>%
  mutate(
    tg   = rep(1:n_time, times = NT),
    k    = rep(1:NT, each  = n_time),
    mean = as.vector(rk_mean),
    low  = as.vector(rk_low),
    high = as.vector(rk_high)
  )

rk_yr0.5 <- subset(rk_df, rk_df$tg==14)
rk_yr1 <- subset(rk_df, rk_df$tg==31)
rk_yr2 <- subset(rk_df, rk_df$tg==65)
rk_yr5 <- subset(rk_df, rk_df$tg==166)

# Combined plot with constant HR
pdf(file = "TVHR_NMA_FE_comb.pdf")
ggplot(hr_df_noref, aes(x = time, y = mean)) +
  geom_line(aes(colour = "TVHR + 95% CrI"), linewidth = 0.9) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = "TVHR + 95% CrI"), alpha = 0.2) +
  geom_line(data = hr_time,
            aes(x = time, y = mean, colour = "Constant HR + 95% CrI"),
            linewidth = 0.9) +
  geom_ribbon(data = hr_time,
              aes(x = time, ymin = lower, ymax = upper, fill = "Constant HR + 95% CrI"),
              alpha = 0.2,
              inherit.aes = FALSE) +
  geom_hline(yintercept = 1,
             linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_colour_manual(values = c("TVHR + 95% CrI" = "steelblue",
                                 "Constant HR + 95% CrI"     = "plum"),
                      name = NULL) +
  scale_fill_manual(values = c("TVHR + 95% CrI" = "steelblue",
                               "Constant HR + 95% CrI"     = "plum"),
                    name = NULL) +
  facet_wrap(~ treatment, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6)) +
  labs(
    title    = "Time-Varying and Constant Hazard Ratio vs DTIC",
    subtitle = "Posterior mean with 95% credible interval",
    x        = "Time (years)",
    y        = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )
dev.off()

###############
# End of file #
###############