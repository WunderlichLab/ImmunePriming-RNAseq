library(survival)
library(survminer)
library(tidyverse)
library(relsurv)

# Useful References
## https://www.datacamp.com/community/tutorials/survival-analysis-R
## http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
## https://cran.r-project.org/web/packages/survminer/readme/README.html
## https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

#### Modeling Survival from Your Data ####
### Read in your data ###
surv_data <- read.delim("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/FINAL_SurvCurv.txt")


### After reading in your data, filter out just the conditions you're interested in analyzing ###
surv_data <- dplyr::filter(surv_data, Genotype == "p38b") %>%
  dplyr::filter(., Treatment == "PBS" | Treatment == "Efae OD 0.05" | Treatment == "Efae OD 0.5")

### Treating your conditions as factors, you can re-order them so that when you plot them they show up in the order you want ###
surv_data$Treatment <- factor(surv_data$Treatment, levels = c("PBS","Efae OD 0.05","Efae OD 0.5"))

### Create an object that specifically models survival ###
surv_surv <- Surv(time = surv_data$Day, event = surv_data$Event)

### Create a linear fit of your survival data so you can do statistical analysis on it ###
surv_fit <- survfit(surv_surv ~ Treatment, data = surv_data)

### Create a tabular object you can use to look at the specific stats you did ###
surv_diff <- survdiff(surv_surv ~ Treatment, data = surv_data)

#### Plotting the Data ####
png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/SurvCurve_p38b_Single.png",width = 1800, height=1800, res=500)

ggsurvplot(surv_fit, 
           data = surv_data, 
           pval = F, 
           fun = "pct",
           risk.table = FALSE, 
           palette = c("#93c1c9","#a7f542", "#4b6e1e"),
           break.time.by = 1,
           #title ="Eater Mutant Single Injections",
           legend = "none", 
           conf.int = T,
           surv.median.line = "v",
           size = 2.5,
           font.x = 18,
           font.y = 18,
           font.tickslab = 15,
           xlab = " ",
           ylab = " ",
           xlim = c(0,7)) 

dev.off()