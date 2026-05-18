# IPD PMA with TVHR

# KM + Residual Plots - arXiv 

# Uses YEARS as time scale 

setwd("wd")

# https://search.r-project.org/CRAN/refmans/surrosurv/html/gastadv.html 

install.packages("surrosurv")
library(surrosurv)
library(coxme)
library(survival)

data("gastadv")
names(gastadv)

gastadv$tx <- gastadv$trt + 0.5

##############################
##############################

install.packages("ggplot2")
install.packages("patchwork")
install.packages("grid")
install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("RColorBrewer")

library(ggplot2)
library(patchwork)
library(grid)
library(survival)
library(survminer)
library(dplyr)
library(RColorBrewer)

##############################

# Check which studies fail PH
cox.all <- list()
coxzph.all <- list()
for (i in 1:20){
  cox.all[[i]] <- coxph(formula = Surv(timeS/365, statusS) ~ tx, data=gastadv,
                     subset=trialref==i) 
  coxzph.all[[i]] <- cox.zph(cox.all[[i]])
}

# KM Plots 

# Only plot 1,5,9,12, 18

df1 <- gastadv[gastadv$trialref==1,]
df2 <- gastadv[gastadv$trialref==5,]
df3 <- gastadv[gastadv$trialref==9,]
df4 <- gastadv[gastadv$trialref==12,]
df5 <- gastadv[gastadv$trialref==18,]

fit1 <- survfit(Surv(timeS/365, statusS) ~ tx, data = df1)
fit2 <- survfit(Surv(timeS/365, statusS) ~ tx, data = df2)
fit3 <- survfit(Surv(timeS/365, statusS) ~ tx, data = df3)
fit4 <- survfit(Surv(timeS/365, statusS) ~ tx, data = df4)
fit5 <- survfit(Surv(timeS/365, statusS) ~ tx, data = df5)


trt_colours <- c("#009E73",   # SoC — green
                 "#D45087"    # Chemo - pink
)

p1 <- ggsurvplot(fit=fit1, data = df1, censor=F, palette=c(trt_colours[1], trt_colours[2]),
                  xlab="Time (years)", ylab="Survival Probability", 
                  xlim=c(0, 2.5), ylim=c(0, 1), title="Study 1",
                  legend=c(0.5,0.9), legend.title="", legend.labs=c("SoC", "Chemotherapy"),
                  font.legend=18)
p1$plot + scale_x_continuous(breaks=sort(c(0, 0.5, 1, 1.5, 2, 2.5))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

p2 <- ggsurvplot(fit=fit2, data = df2, censor=F, palette=c(trt_colours[1], trt_colours[2]),
                 xlab="Time (years)", ylab="Survival Probability", 
                 xlim=c(0, 2.5), ylim=c(0, 1), title="Study 5",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("SoC", "Chemotherapy"),
                 font.legend=18)
p2$plot + scale_x_continuous(breaks=sort(c(0, 0.5, 1, 1.5, 2, 2.5))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

p3 <- ggsurvplot(fit=fit3, data = df3, censor=F, palette=c(trt_colours[1], trt_colours[2]),
                 xlab="Time (years)", ylab="Survival Probability", 
                 xlim=c(0, 2.5), ylim=c(0, 1), title="Study 9",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("SoC", "Chemotherapy"),
                 font.legend=18)
p3$plot + scale_x_continuous(breaks=sort(c(0, 0.5, 1, 1.5, 2, 2.5))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

p4 <- ggsurvplot(fit=fit4, data = df4, censor=F, palette=c(trt_colours[1], trt_colours[2]),
                 xlab="Time (years)", ylab="Survival Probability", 
                 xlim=c(0, 2.5), ylim=c(0, 1), title="Study 12",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("SoC", "Chemotherapy"),
                 font.legend=18)
p4$plot + scale_x_continuous(breaks=sort(c(0, 0.5, 1, 1.5, 2, 2.5))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))


p5 <- ggsurvplot(fit=fit5, data = df5, censor=F, palette=c(trt_colours[1], trt_colours[2]),
                 xlab="Time (years)", ylab="Survival Probability", 
                 xlim=c(0, 2.5), ylim=c(0, 1), title="Study 18",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("SoC", "Chemotherapy"),
                 font.legend=18)
p6$plot + scale_x_continuous(breaks=sort(c(0, 0.5, 1, 1.5, 2, 2.5))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

plots <- list(p1, p2, p3, p4, p5)

extract_plot <- function(x) {
  if (inherits(x, "ggplot"))                                    return(x)
  if (!is.null(x$plot) && inherits(x$plot, "ggplot"))           return(x$plot)
  if (inherits(x, "ggsurvplot"))                                return(x$plot)
  if (is.list(x) && length(x) > 0 && inherits(x[[1]], "ggplot")) return(x[[1]])
  stop(paste("Could not extract ggplot from object of class:", paste(class(x), collapse = ", ")))
}

clean_plots <- lapply(plots, function(x) {
  p <- extract_plot(x)
  
  p +
    scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       limits = c(0, 1)) +
    theme_classic(base_size = 10) +
    theme(
      axis.title        = element_text(size = 9),
      axis.text         = element_text(size = 8),
      plot.title        = element_text(size = 10, face = "bold"),
      legend.position   = c(0.8, 0.8),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.text       = element_text(size = 6),
      legend.title      = element_text(size = 7),
      legend.key.size   = unit(0.2, "cm")
    )
})

km_fig <- wrap_plots(clean_plots, ncol = 3)
km_fig

ggsave("KM_figure_PMA.pdf",
       km_fig,
       width = 14,
       height = 10)

##############################

# Residual plots 

# Fit the Cox PH model
cox1 <- coxph(formula = Surv(timeS/365, statusS) ~ tx, data=df1)
cox2 <- coxph(formula = Surv(timeS/365, statusS) ~ tx, data=df2)
cox3 <- coxph(formula = Surv(timeS/365, statusS) ~ tx, data=df3)
cox4 <- coxph(formula = Surv(timeS/365, statusS) ~ tx, data=df4)
cox5 <- coxph(formula = Surv(timeS/365, statusS) ~ tx, data=df5)

# Test the PH assumption
test1 <- cox.zph(cox1)
test2 <- cox.zph(cox2)
test3 <- cox.zph(cox3)
test4 <- cox.zph(cox4)
test5 <- cox.zph(cox5)

# Graphical test of PH
ph1<-ggcoxzph(test1)
ph2<-ggcoxzph(test2)
ph3<-ggcoxzph(test3)
ph4<-ggcoxzph(test4)
ph5<-ggcoxzph(test5)

plots <- list(ph1, ph2, ph3, ph4, ph5)

study_names <- c("Study 1", "Study 5", "Study 9", "Study 12", "Study 18")

stopifnot(length(plots) == length(study_names))

extract_ggcoxzph <- function(x) {
  if (inherits(x, "ggplot")) return(x)
  if (is.list(x) && !is.null(x[[1]]) && inherits(x[[1]], "ggplot")) return(x[[1]])
  stop("Not a ggplot-containing object")
}

clean_plots <- lapply(seq_along(plots), function(i) {
  p <- extract_ggcoxzph(plots[[i]])
  
  p +
    labs(title = study_names[i]) +
    scale_y_continuous(limits = c(-4, 2), breaks = c(-4, -3, -2, -1, 0, 1, 2)) +
    theme_classic(base_size = 10) +
    theme(
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
})

cox_ph_fig <- wrap_plots(clean_plots, ncol = 3)

cox_ph_fig

ggsave("cox_ph_diagnostics_PMA.pdf",
       cox_ph_fig,
       width = 12,
       height = 10,
       dpi = 300)

##############################
##############################


