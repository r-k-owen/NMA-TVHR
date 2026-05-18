# Clear environment
rm(list=ls())

# Plot KM curves

# Load libraries
library(survival)
library(survminer)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(grid)

# Set working directory
wd <- "wd"
setwd(wd)

# Load data
data <- read.csv("melanoma_ipd.csv")
df1 <- data[data$studyCode==1,]
df2 <- data[data$studyCode==2,]
df3 <- data[data$studyCode==3,]
df4 <- data[data$studyCode==4,]
df5 <- data[data$studyCode==5,]
df6 <- data[data$studyCode==6,]
df7 <- data[data$studyCode==7,]
df8 <- data[data$studyCode==8,]
df9 <- data[data$studyCode==9,]
df10 <- data[data$studyCode==10,]
df11 <- data[data$studyCode==11,]
df12 <- data[data$studyCode==12,]
df13 <- data[data$studyCode==13,]

# Kaplan-Meier estimates
fit1 <- survfit(Surv(time, event) ~ arm, data = df1)
fit2 <- survfit(Surv(time, event) ~ arm, data = df2)
fit3 <- survfit(Surv(time, event) ~ arm, data = df3)
fit4 <- survfit(Surv(time, event) ~ arm, data = df4)
fit5 <- survfit(Surv(time, event) ~ arm, data = df5)
fit6 <- survfit(Surv(time, event) ~ arm, data = df6)
fit7 <- survfit(Surv(time, event) ~ arm, data = df7)
fit8 <- survfit(Surv(time, event) ~ arm, data = df8)
fit9 <- survfit(Surv(time, event) ~ arm, data = df9)
fit10 <- survfit(Surv(time, event) ~ arm, data = df10)
fit11 <- survfit(Surv(time, event) ~ arm, data = df11)
fit12 <- survfit(Surv(time, event) ~ arm, data = df12)
fit13 <- survfit(Surv(time, event) ~ arm, data = df13)

# Combine fits
fit <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit13)
df_list <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13)

# Plot survival curves for each study
ggsurvplot_list(fit=fit, data=df_list, censor=F)

# Colour palette
trt_colours <- c("#999999", 
                 "#BBBBBB",
                 "#DDDDDD",   
                 "#E69F00",   
                 "#D55E00",   
                 "#F0E442",   
                 "#AE123E",   
                 "#56B4E9",   
                 "#0072B2",   
                 "#009E73",   
                 "#006040",   
                 "#2F4B7C",   
                 "#D45087"    
)

p1 <- ggsurvplot(fit=fit1, data = df1, censor=F, palette=c(trt_colours[1], trt_colours[4]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="BREAK-3",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Dabrafenib"),
                 font.legend=18)
p1$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p2 <- ggsurvplot(fit=fit2, data = df2, censor=F, palette=c(trt_colours[1], trt_colours[6]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="BRIM-3",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Vemurafenib"),
                 font.legend=18)
p2$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p3 <- ggsurvplot(fit=fit3, data = df3, censor=F, palette=c(trt_colours[1], trt_colours[9]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="CheckMate 066",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Nivolumab"),
                 font.legend=18)
p3$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p4 <- ggsurvplot(fit=fit4, data = df4, censor=F, palette=c(trt_colours[8], trt_colours[9], trt_colours[10]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="CheckMate 067",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Ipilimumab", "Nivolumab", "Nivolumab + Ipilimumab"),
                 font.legend=18)
p4$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p5 <- ggsurvplot(fit=fit5, data = df5, censor=F, palette=c(trt_colours[8], trt_colours[10]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="CheckMate 069",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Ipilimumab", "Nivolumab + Ipilimumab"),
                 font.legend=18)
p5$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p6 <- ggsurvplot(fit=fit6, data = df6, censor=F, palette=c(trt_colours[6], trt_colours[7]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="COBRIM",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Vemurafenib", "Vemurafenib + Cobimetinib"),
                 font.legend=18)
p6$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p7 <- ggsurvplot(fit=fit7, data = df7, censor=F, palette=c(trt_colours[4], trt_colours[5]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="COMBI-d",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dabrafenib", "Dabrafenib + Trametinib"),
                 font.legend=18)
p7$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p8 <- ggsurvplot(fit=fit8, data = df8, censor=F, palette=c(trt_colours[5], trt_colours[6]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="COMBI-v",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dabrafenib + Trametinib", "Vemurafenib"),
                 font.legend=18)
p8$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p9 <- ggsurvplot(fit=fit9, data = df9, censor=F, palette=c(trt_colours[8], trt_colours[12]),
                 xlab="Time (months)", ylab="Survival Probability", 
                 xlim=c(0, 72), ylim=c(0, 1), title="Hodi 2014",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Ipilimumab", "Ipilimumab + Sargramostin"),
                 font.legend=18)
p9$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p10 <- ggsurvplot(fit=fit10, data = df10, censor=F, palette=c(trt_colours[8], trt_colours[11]),
                  xlab="Time (months)", ylab="Survival Probability", 
                  xlim=c(0, 72), ylim=c(0, 1), title="Keynote 006",
                  legend=c(0.5,0.9), legend.title="", legend.labs=c("Ipilimumab", "Pembrolizumab"),
                  font.legend=18)
p10$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p11 <- ggsurvplot(fit=fit11, data = df11, censor=F, palette=c(trt_colours[1], trt_colours[2]),
                  xlab="Time (months)", ylab="Survival Probability", 
                  xlim=c(0, 72), ylim=c(0, 1), title="Ribas 2013",
                  legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Tremelimumab"),
                  font.legend=18)
p11$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p12 <- ggsurvplot(fit=fit12, data = df12, censor=F, palette=c(trt_colours[1], trt_colours[3]),
                  xlab="Time (months)", ylab="Survival Probability", 
                  xlim=c(0, 72), ylim=c(0, 1), title="Robert 2011",
                  legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Ipilimumab + Dacarbazine"),
                  font.legend=18)
p12$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p13 <- ggsurvplot(fit=fit13, data = df13, censor=F, palette=c(trt_colours[1], trt_colours[13]),
                  xlab="Time (months)", ylab="Survival Probability", 
                  xlim=c(0, 72), ylim=c(0, 1), title="Robert 2013",
                  legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Selumetinib + Dacarbazine"),
                  font.legend=18)
p13$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))



plots <- list(p1, p2, p3, p4, p5, p6,
              p7, p8, p9, p10, p11, p12, p13)

uniform_plots <- lapply(plots, function(x) {
  x$plot + scale_x_continuous(breaks=sort(c(0, 12, 24, 36, 48, 60, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1))) +
    theme_classic(base_size = 10) +
    theme(
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = c(0.8, 0.8),  
      legend.background = element_rect(fill = "white", colour = NA),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.2, "cm")
    ) 
})

km_fig <- wrap_plots(uniform_plots, ncol = 4)
km_fig

ggsave("KM_figure_NMA.pdf",
       km_fig,
       width = 14,
       height = 10)

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


# Fit the Cox PH model
cox1 <- coxph(formula = Surv(time, event) ~ arm, data=df1)
cox2 <- coxph(formula = Surv(time, event) ~ arm, data=df2)
cox3 <- coxph(formula = Surv(time, event) ~ arm, data=df3)
cox4 <- coxph(formula = Surv(time, event) ~ arm, data=df4)
cox5 <- coxph(formula = Surv(time, event) ~ arm, data=df5)
cox6 <- coxph(formula = Surv(time, event) ~ arm, data=df6)
cox7 <- coxph(formula = Surv(time, event) ~ arm, data=df7)
cox8 <- coxph(formula = Surv(time, event) ~ arm, data=df8)
cox9 <- coxph(formula = Surv(time, event) ~ arm, data=df9)
cox10 <- coxph(formula = Surv(time, event) ~ arm, data=df10)
cox11 <- coxph(formula = Surv(time, event) ~ arm, data=df11)
cox12 <- coxph(formula = Surv(time, event) ~ arm, data=df12)
cox13 <- coxph(formula = Surv(time, event) ~ arm, data=df13)

# Test the PH assumption
test1 <- cox.zph(cox1)
test2 <- cox.zph(cox2)
test3 <- cox.zph(cox3)
test4 <- cox.zph(cox4)
test5 <- cox.zph(cox5)
test6 <- cox.zph(cox6)
test7 <- cox.zph(cox7)
test8 <- cox.zph(cox8)
test9 <- cox.zph(cox9)
test10 <- cox.zph(cox10)
test11 <- cox.zph(cox11)
test12 <- cox.zph(cox12)
test13 <- cox.zph(cox13)

# Graphical test of PH
ph1<-ggcoxzph(test1)
ph2<-ggcoxzph(test2)
ph3<-ggcoxzph(test3)
ph4<-ggcoxzph(test4)
ph5<-ggcoxzph(test5)
ph6<-ggcoxzph(test6)
ph7<-ggcoxzph(test7)
ph8<-ggcoxzph(test8)
ph9<-ggcoxzph(test9)
ph10<-ggcoxzph(test10)
ph11<-ggcoxzph(test11)
ph12<-ggcoxzph(test12)
ph13<-ggcoxzph(test13)

# Visual evidence of proportional hazards violation in many studies
# Studies 2, 4, and 10 are statistically significant

plots <- list(ph1, ph2, ph3, ph4, ph5, ph6,
              ph7, ph8, ph9, ph10, ph11, ph12, ph13)

lapply(plots, class)

study_names <- c(
  "BREAK-3", "BRIM-3", "CheckMate 066", "CheckMate 067", "CheckMate 069", "COBRIM",
  "COMBI-d", "COMBI-v", "Hodi 2014", "Keynote 006", "Ribas 2013", "Robert 2011", "Robert 2013"
)

extract_ggcoxzph <- function(x) {
  if (inherits(x, "ggplot")) return(x)
  if (is.list(x) && !is.null(x[[1]]) && inherits(x[[1]], "ggplot")) return(x[[1]])
  stop("Not a ggplot-containing object")
}

uniform_plots <- lapply(plots, extract_ggcoxzph)

uniform_plots <- lapply(seq_along(uniform_plots), function(i) {
  uniform_plots[[i]] +
    labs(title = study_names[i])
})

library(patchwork)

wrap_plots(uniform_plots, ncol = 4)

uniform_plots <- lapply(uniform_plots, function(p) {
  p +
    theme_classic(base_size = 10) +
    theme(
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "none"  
    )
})

cox_ph_fig<-wrap_plots(uniform_plots, ncol = 4)

ggsave("cox_ph_diagnostics_NMA.pdf",
       cox_ph_fig,
       width = 12,
       height = 10,
       dpi = 300)

# Log-log plots
s1 <- Surv(df1$time, df1$event)
s2 <- Surv(df2$time, df2$event)
s3 <- Surv(df3$time, df3$event)
s4 <- Surv(df4$time, df4$event)
s5 <- Surv(df5$time, df5$event)
s6 <- Surv(df6$time, df6$event)
s7 <- Surv(df7$time, df7$event)
s8 <- Surv(df8$time, df8$event)
s9 <- Surv(df9$time, df9$event)
s10 <- Surv(df10$time, df10$event)
s11 <- Surv(df11$time, df11$event)
s12 <- Surv(df12$time, df12$event)
s13 <- Surv(df13$time, df13$event)

plot(survfit(s1 ~ df1$txCode), col=c("blue", "red"), fun="cloglog", 
     xlab="Log(time)", ylab="log(-log(survival probability))")
plot(survfit(s2 ~ df2$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s3 ~ df3$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s4 ~ df4$txCode), col=c("blue", "red", "green"), fun="cloglog")
plot(survfit(s5 ~ df5$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s6 ~ df6$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s7 ~ df7$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s8 ~ df8$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s9 ~ df9$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s10 ~ df10$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s11 ~ df11$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s12 ~ df12$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s13 ~ df13$txCode), col=c("blue", "red"), fun="cloglog")

###############
# End of file #
###############