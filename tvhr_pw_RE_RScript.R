# IPD 2-stage PWMA with TVHR

# arXiv - 2-stage using JAGS & Spherical parameterisation for PFS (YEARS)

setwd("wd")

install.packages("surrosurv")

library(surrosurv)
library(meta)
library(rjags)

data("gastadv")
names(gastadv)

gastadv$tx <- gastadv$trt + 0.5

#################################
#################################

# 1-stage Analyses

library(coxme)

# Random effects Cox model — random trial-specific baseline hazard
# and random treatment effect
fit <- coxme(Surv(timeS/365, statusS) ~ tx + (tx | trialref), data = gastadv)
summary(fit)

#################################
#################################

# 2-Stage Analysis

# PH model 
fit.cox <- list()
loghr <- matrix(NA,ncol=1,nrow=20)
covar <- matrix(NA,ncol=1,nrow=20)
precLHR <- matrix(NA,ncol=1,nrow=20)

for (i in 1:20){
fit.cox[[i]] <- coxph(Surv(timeS/365,statusS)~tx,data=gastadv,subset=trialref==i) 
loghr[i,1] <- coef(fit.cox[[i]])[1]
covar[i,1] <- vcov(fit.cox[[i]])[1] 
precLHR[i,1] <- 1/covar[i,1] }

# MA model 
model.ma.lhr <- metagen(loghr,sqrt(covar))
model.ma.lhr
pdf("forest_figure_PMA.pdf",
    width = 14,
    height = 10)
plot(model.ma.lhr) # Forest plot 
dev.off()

exp(model.ma.lhr$TE.random)
exp(model.ma.lhr$lower.random)
exp(model.ma.lhr$upper.random)
model.ma.lhr$method.random

# JAGS MA model for LHR

# JAGS data
lhr.jags.data <- list()
lhr.jags.data$y <- loghr[,1]
lhr.jags.data$prec <- precLHR[,1]
lhr.jags.data$ns <- length(loghr[,1])

# JAGS inits 
lhr.jags.inits1 <- list()
lhr.jags.inits2 <- list()
lhr.jags.inits3 <- list()

lhr.jags.inits1$d <- 0
lhr.jags.inits1$mu <- rep(0,len=lhr.jags.data$ns)
lhr.jags.inits1$sd <- 0.5

lhr.jags.inits2$d <- -1
lhr.jags.inits2$mu <- rep(1,len=lhr.jags.data$ns)
lhr.jags.inits2$sd <- 0.9

lhr.jags.inits3$d <- 1
lhr.jags.inits3$mu <- rep(-1,len=lhr.jags.data$ns)
lhr.jags.inits3$sd <- 0.1

lhr.jags.inits <- list(lhr.jags.inits1,lhr.jags.inits2,lhr.jags.inits3)

lhr.jags.model <- "model {
  #################################################################
  # PMA 
  # 
  # Data expected:
  #  - ns                : number of studies
  #  - y[i]              : LogHR for study i 
  #  - prec              : Precision for study i 
  #
  #################################################################
  
  for (i in 1:ns){
  y[i] ~ dnorm(mu[i], prec[i]) # likelihood
  mu[i] ~ dnorm(d,tau) }       # Random effects distribution 
  
  # Priors
  d ~ dnorm(0, 0.001)
  tau <- pow(sd,-2)
  sd ~ dunif(0,1)
  
  HR <- exp(d)  # Calculate pooled HR
}"

# Compile JAGS model 
model.lhr <- jags.model(textConnection(lhr.jags.model),
                         data=lhr.jags.data,
                         inits=lhr.jags.inits,n.chains=3,n.adapt=10000)

# Running JAGS model 

update(model.lhr, 50000)
summary.lhr <- coda.samples(model.lhr, c("d","HR","sd"),
                             n.iter=50000,thin=5)

save(summary.lhr, file = "summarypmalhr.RData")

summary(summary.lhr)
traceplot(summary.lhr)
autocorr.plot(summary.lhr)
gelman.diag(summary.lhr)
effectiveSize(summary.lhr)

#################################################################
#################################################################

# log(time) interaction model 

fit.cox.tt <- list()
loghr.tt <- matrix(NA,ncol=2,nrow=20)
covar.tt <- matrix(NA,ncol=3,nrow=20)
loghr.tv <- matrix(NA,ncol=100,nrow=20)
loghr.tv.low <- matrix(NA,ncol=100,nrow=20)
loghr.tv.up <- matrix(NA,ncol=100,nrow=20)
PrecW.tt <- array(NA, dim = c(2, 2, 20))

for (i in 1:20){
  fit.cox.tt[[i]] <- coxph(Surv(timeS/365,statusS)~tx+tt(tx),data=gastadv,subset=trialref==i,
                           tt      = function(x, t, ...) x * log(t))  
  loghr.tt[i,1:2] <- coef(fit.cox.tt[[i]])
  covar.tt[i,1:3] <- c(vcov(fit.cox.tt[[i]])[1,1],vcov(fit.cox.tt[[i]])[1,2],vcov(fit.cox.tt[[i]])[2,2])
  # Precision matrix per study 
  PrecW.tt[1:2,1:2,i] <- solve(vcov(fit.cox.tt[[i]]))
  loghr.tv[i,1:100] <- loghr.tt[i,1] + loghr.tt[i,2]*log(seq(0,1200,len=100))
  loghr.tv.low[i,1:100] <- loghr.tv[i,1:100] - 1.96*sqrt(covar.tt[1] + log(seq(0,1200,len=100))^2*covar.tt[3] + 
                                                           2*log(seq(0,1200,len=100))*covar.tt[2])
  loghr.tv.up[i,1:100] <- loghr.tv[i,1:100] + 1.96*sqrt(covar.tt[1] + log(seq(0,1200,len=100))^2*covar.tt[3] + 
                                                           2*log(seq(0,1200,len=100))*covar.tt[2])
  }

# MA model 
library(mvmeta)

model.re.tt <- mvmeta(loghr.tt, covar.tt, method = "reml")
summary(model.re.tt) 

#################################
#################################

# Plot of study-specific TVHR & Constant HR 

# ── Helper function to fit Cox and Cox-tt models ───────────────────────
fit_study <- function(data, trial) {
  sub <- data[data$trialref == trial, ]
  
  fit.cox <- coxph(Surv(timeS/365, statusS) ~ tx, data = sub)
  b       <- coef(fit.cox)
  vc      <- vcov(fit.cox)
  
  fit.coxtt <- coxph(Surv(timeS/365, statusS) ~ tx + tt(tx), data = sub,
                     tt = function(x, t, ...) x * log(t))
  b.tt  <- coef(fit.coxtt)
  vc.tt <- vcov(fit.coxtt)
  
  time_grid <- seq(0, 3.5, length.out = 200)
  lp    <- b.tt["tx"] + b.tt["tt(tx)"] * log(time_grid)
  se_lp <- sqrt(vc.tt["tx","tx"] +
                  2 * vc.tt["tx","tt(tx)"] * log(time_grid) +
                  vc.tt["tt(tx)","tt(tx)"] * log(time_grid)^2)
  
  list(
    df       = data.frame(time  = time_grid,
                          mean  = exp(lp),
                          lower = exp(lp - 1.96 * se_lp),
                          upper = exp(lp + 1.96 * se_lp)),
    HR.mean  = exp(b),
    HR.lower = exp(b - 1.96 * sqrt(vc)),
    HR.upper = exp(b + 1.96 * sqrt(vc))
  )
}

# ── Helper function to build one TVHR plot ─────────────────────────────
make_tvhr_plot <- function(data, HR.mean, HR.lower, HR.upper, title) {
  ggplot(data, aes(x = time, y = mean)) +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = HR.lower, ymax = HR.upper,
             fill = "plum", alpha = 0.2) +
    geom_hline(yintercept = HR.mean,
               linetype = "solid", colour = "plum", linewidth = 0.9) +
    geom_ribbon(aes(ymin = lower, ymax = pmin(upper, 5)),
                fill = "steelblue", alpha = 0.2) +
    geom_line(colour = "steelblue", linewidth = 0.9) +
    geom_hline(yintercept = 1,
               linetype = "dashed", colour = "red", linewidth = 0.5) +
    scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
    scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
    coord_cartesian(ylim = c(0, 5), clip = "off") +
    labs(title = title, x = "Time (years)", y = "Hazard Ratio") +
    theme_minimal(base_size = 11) +
    theme(strip.text       = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

# ── Fit and plot each study ────────────────────────────────────────────
s1  <- fit_study(gastadv, 1)
s5  <- fit_study(gastadv, 5)
s9  <- fit_study(gastadv, 9)
s12 <- fit_study(gastadv, 12)
s18 <- fit_study(gastadv, 18)

p1  <- make_tvhr_plot(s1$df,  s1$HR.mean,  s1$HR.lower,  s1$HR.upper,  "Study 1")
p5  <- make_tvhr_plot(s5$df,  s5$HR.mean,  s5$HR.lower,  s5$HR.upper,  "Study 5")
p9  <- make_tvhr_plot(s9$df,  s9$HR.mean,  s9$HR.lower,  s9$HR.upper,  "Study 9")
p12 <- make_tvhr_plot(s12$df, s12$HR.mean, s12$HR.lower, s12$HR.upper, "Study 12")
p18 <- make_tvhr_plot(s18$df, s18$HR.mean, s18$HR.lower, s18$HR.upper, "Study 18")

# ── Legend panel (vertical) ────────────────────────────────────────────
swatch_w <- 0.25
gap      <- 0.1
swatch_h <- 0.06
lx       <- 0.1
row1     <- 0.75
row2     <- 0.50
row3     <- 0.25

p_legend <- ggplot() +
  xlim(0, 2) +
  ylim(0, 1) +
  
  # Row 1: TVHR steelblue shading + line
  annotate("rect",
           xmin = lx, xmax = lx + swatch_w,
           ymin = row1 - swatch_h, ymax = row1 + swatch_h,
           fill = "steelblue", alpha = 0.2, colour = NA) +
  annotate("segment",
           x = lx, xend = lx + swatch_w, y = row1, yend = row1,
           colour = "steelblue", linewidth = 0.9) +
  annotate("text",
           x = lx + swatch_w + gap, y = row1,
           label = "TVHR + 95% CI",
           hjust = 0, vjust = 0.5, size = 3, colour = "grey20") +
  
  # Row 2: Cox HR plum shading + line
  annotate("rect",
           xmin = lx, xmax = lx + swatch_w,
           ymin = row2 - swatch_h, ymax = row2 + swatch_h,
           fill = "plum", alpha = 0.2, colour = NA) +
  annotate("segment",
           x = lx, xend = lx + swatch_w, y = row2, yend = row2,
           colour = "plum", linewidth = 0.9) +
  annotate("text",
           x = lx + swatch_w + gap, y = row2,
           label = "Constant HR + 95% CI",
           hjust = 0, vjust = 0.5, size = 3, colour = "grey20") +
  
  # Row 3: HR = 1 red dashed
  annotate("segment",
           x = lx, xend = lx + swatch_w, y = row3, yend = row3,
           colour = "red", linetype = "dashed", linewidth = 0.5) +
  annotate("text",
           x = lx + swatch_w + gap, y = row3,
           label = "HR = 1",
           hjust = 0, vjust = 0.5, size = 3, colour = "grey20") +
  
  theme_void()

# ── Assemble and save ──────────────────────────────────────────────────
extract_ggcoxzph <- function(x) {
  if (inherits(x, "ggplot")) return(x)
  if (is.list(x) && !is.null(x[[1]]) && inherits(x[[1]], "ggplot")) return(x[[1]])
  stop("Not a ggplot-containing object")
}

plots       <- list(p1, p5, p9, p12, p18, p_legend)
study_names <- c("Study 1", "Study 5", "Study 9", "Study 12", "Study 18", " ")
stopifnot(length(plots) == length(study_names))

clean_plots <- lapply(seq_along(plots), function(i) {
  p <- extract_ggcoxzph(plots[[i]])
  
  if (i == length(plots)) return(p)
  
  p +
    labs(title = study_names[i]) +
    scale_y_continuous(limits = c(0, 5), breaks = c(0, 1, 2, 3, 4, 5)) +
    theme_classic(base_size = 10) +
    theme(
      axis.title      = element_text(size = 9),
      axis.text       = element_text(size = 8),
      plot.title      = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
})

study_tvhr_fig <- wrap_plots(clean_plots, ncol = 3)
study_tvhr_fig

ggsave("study_tvhr_PMA.pdf",
       study_tvhr_fig,
       width = 14,
       height = 10,
       dpi = 300)

#################################
#################################

# 2-stage JAGS 

# JAGS data
tvhr.jags.data <- list()
tvhr.jags.data$y <- loghr.tt
tvhr.jags.data$PrecW <- PrecW.tt
tvhr.jags.data$times <- c(250/365,500/365,750/365,1000/365)
tvhr.jags.data$ntimes <- length(tvhr.jags.data$times)
tvhr.jags.data$time_grid <- seq(0.001,3.5,len=200)
tvhr.jags.data$ntime_grid <- length(tvhr.jags.data$time_grid)
tvhr.jags.data$NS <- length(tvhr.jags.data$y[,1])

# JAGS inits
tvhr.jags.inits1 <- list()
tvhr.jags.inits2 <- list()
tvhr.jags.inits3 <- list()

tvhr.jags.inits1$d <- c(0,0)
tvhr.jags.inits1$u <- matrix(0,ncol=2,nrow=20)
tvhr.jags.inits1$sd <- c(0.2,0.2)
tvhr.jags.inits1$a12 <- 0.5

tvhr.jags.inits2$d <- c(1,1)
tvhr.jags.inits2$u <- matrix(0.1,ncol=2,nrow=20)
tvhr.jags.inits2$sd <- c(0.8,0.8)
tvhr.jags.inits2$a12 <- 1

tvhr.jags.inits3$d <- c(-1,-1)
tvhr.jags.inits3$u <- matrix(0.2,ncol=2,nrow=20)
tvhr.jags.inits3$sd <- c(0.5,0.5)
tvhr.jags.inits3$a12 <- 2

tvhr.jags.inits <- list(tvhr.jags.inits1,tvhr.jags.inits2,tvhr.jags.inits3)

# JAGS model 

tvhr.jags.model <- "model {
  #################################################################
  # Bivariate PMA 
  # Outcomes:
  #   1 = log(HR)
  #   2 = log(interaction with log(time))
  #
  # Data expected:
  #  - NS                : number of studies
  #  - y[i,c,1:2]        : observed contrast (logHR, logInt) for study i 
  #  - PrecW[i,1:2,1:2]  : within-study precision matrix (inverse covariance) for study i 
  #  - times             : specific time points at which to calculate HR
  #  - ntimes            : number of specific time points at which to calculate HR 
  #  - time_grid         : time points at which to calculate HR for plotting
  #  - ntime_grid        : number of time points at which to calculate HR for plotting 
  #
  #################################################################
  
  pi <- 3.141592653589793
  
  # small zero vector used below
  zero2[1] <- 0
  zero2[2] <- 0
  
  ###############
  # Within-study likelihood 
  ###############
  for (i in 1:NS) {
    # Random effects model 
      u[i,1:2] ~ dmnorm(zero2[], PrecB[1:2,1:2])
    
    # 
      for (o in 1:2) {
        mu[i,o] <- d[o] + u[i, o ] 
      }
      # observed bivariate data with known within-study precision
      y[i,1:2] ~ dmnorm(mu[i,1:2], PrecW[1:2,1:2,i])
    }
  # Prior on pooled log(HR) & interaction 
      d[1] ~ dnorm(0, 0.001)
      d[2] ~ dnorm(0, 0.001)
      
  ###############
  # Between-study heterogeneity (spherical parameterisation)
  # Build SigmaB (2x2) from outcome SDs and correlation = cos(angle)
  ###############
  # priors for heterogeneity SDs (weakly informative)
  sd[1] ~ dunif(0, 1)    # SD for outcome 1 (logHR)
  sd[2] ~ dunif(0, 1)    # SD for outcome 2 (log(interaction))
  
  # spherical / angular parameter for correlation (maps angle in (0, pi) to rho in (-1,1))
  a12 ~ dunif(0, pi)
  
  rho12 <- cos(a12)
  
  # covariance matrix SigmaB
  SigmaB[1,1] <- pow(sd[1], 2)
  SigmaB[2,2] <- pow(sd[2], 2)
  SigmaB[1,2] <- rho12 * sd[1] * sd[2]
  SigmaB[2,1] <- SigmaB[1,2]
  
  # precision for the between-study random effects (used by u ~ dmnorm)
  PrecB[1:2,1:2] <- inverse(SigmaB[,])
  
  # Calculate HR at specific time points for table 
  for (i in 1:ntimes){
    HR[i] <- exp(d[1] + log(times[i])*d[2])
}

 # Calculate HR at specific time points for plotting  
  for (i in 1:ntime_grid){
    HRp[i] <- exp(d[1] + log(time_grid[i])*d[2])
}
  
}"

# Compile JAGS model 
model.tvhr <- jags.model(textConnection(tvhr.jags.model),
                          data=tvhr.jags.data,
                          inits=tvhr.jags.inits,n.chains=3,n.adapt=10000)

# Running JAGS model 

update(model.tvhr, 50000)
summary.tvhr <- coda.samples(model.tvhr, c("d","HR","HRp","sd","rho12"),
                                 n.iter=50000,thin=5)

save(summary.tvhr, file = "summarypmatvhr.RData")

summary(summary.tvhr)
traceplot(summary.tvhr)
autocorr.plot(summary.tvhr)
gelman.diag(summary.tvhr)
effectiveSize(summary.tvhr)

# Summary of HR at specific time points
params_of_interest <- "^HR\\["

samples_HR <- summary.tvhr[, grep(paste(params_of_interest, collapse = "|"),
                              varnames(summary.tvhr)), drop = FALSE]
summary(samples_HR)
plot(samples_HR)
plot(samples_HR[, 1, drop = FALSE])
plot(samples_HR[, 2, drop = FALSE])
plot(samples_HR[, 3, drop = FALSE])
plot(samples_HR[, 4, drop = FALSE])

# Extract HR (+ low/high) at time_grid 

summary_fit <- summary(summary.tvhr)

subset_summary <- summary(summary.tvhr[, grep("^HRp", varnames(summary.tvhr))])

HR_df <- data.frame(time = seq(0.001,3.5,len=200),
             mean = subset_summary$statistics[, "Mean"],
lower = subset_summary$quantiles[, "2.5%"],
upper = subset_summary$quantiles[, "97.5%"])

# Overall TVHR plot: 

# ── Main plot ──────────────────────────────────────────────────────────
tvhr_fig <- ggplot(HR_df, aes(x = time, y = mean)) +
  geom_line(colour = "steelblue", linewidth = 0.9) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "steelblue", alpha = 0.2) +
  geom_hline(yintercept = 1,
             linetype = "dashed", colour = "red", linewidth = 0.6) +
  geom_hline(yintercept = 0.78,
             linetype = "solid", colour = "plum", linewidth = 0.9) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0.70, ymax = 0.86,
           fill = "plum", alpha = 0.2) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1, 1.25)) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "Time-Varying Hazard Ratio",
    subtitle = "Posterior mean with 95% credible interval",
    x        = "Time (years)",
    y        = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ── Legend panel ───────────────────────────────────────────────────────
# Total legend width = swatch (0.25) + gap (0.08) + label (~0.7) + spacer (0.2) + swatch + gap + label
# Centred on x = 2 (midpoint of xlim 0–4)
swatch_w <- 0.25
gap      <- 0.08
spacer   <- 0.25
ly       <- 0.5
swatch_h <- 0.5

# Entry 1 starts at lx1, entry 2 starts at lx2
lx1 <- 0.9
lx2 <- lx1 + swatch_w + gap + 0.75 + spacer   # ~2.13

p_legend <- ggplot() +
  xlim(0, 4) + ylim(0, 1) +
  
  # TVHR: steelblue shading + line
  annotate("rect",
           xmin = lx1, xmax = lx1 + swatch_w,
           ymin = ly - swatch_h, ymax = ly + swatch_h,
           fill = "steelblue", alpha = 0.2, colour = NA) +
  annotate("segment",
           x = lx1, xend = lx1 + swatch_w, y = ly, yend = ly,
           colour = "steelblue", linewidth = 0.9) +
  annotate("text",
           x = lx1 + swatch_w + gap, y = ly,
           label = "TVHR + 95% CrI",
           hjust = 0, vjust = 0.5, size = 3, colour = "grey20") +
  
  # REML HR: plum shading + line
  annotate("rect",
           xmin = lx2, xmax = lx2 + swatch_w,
           ymin = ly - swatch_h, ymax = ly + swatch_h,
           fill = "plum", alpha = 0.2, colour = NA) +
  annotate("segment",
           x = lx2, xend = lx2 + swatch_w, y = ly, yend = ly,
           colour = "plum", linetype = "solid", linewidth = 0.9) +
  annotate("text",
           x = lx2 + swatch_w + gap, y = ly,
           label = "Constant HR + 95% CrI",
           hjust = 0, vjust = 0.5, size = 3, colour = "grey20") +
  theme_void()

# ── Combine ────────────────────────────────────────────────────────────
tvhr_fig_comb <- tvhr_fig / p_legend +
  plot_layout(heights = c(10, 0.6))

tvhr_fig_comb

ggsave("TVHR_figure_PMA.pdf",
       tvhr_fig_comb,
       width = 14,
       height = 10)

#################################
# END OF FILE
#################################



