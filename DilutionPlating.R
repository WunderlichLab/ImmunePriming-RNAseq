library(tidyverse)
library(extrafont)

setwd("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/")

dilpl8_efae <- read_tsv("/Users/cabweebwub/Documents/Wunderlich Lab Docs/RData/FINAL_DilPl8.txt") %>% 
  dplyr::filter(., .$Pathogen=="Efae")

#### Box Plot - Single Injections ####
png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/DilPl8_Single.png",width = 1750, height=1350, res=300)
dilpl8_efae %>%
  dplyr::filter(., .$Conc=="0.05" | .$Conc=="0.5") %>%
  ggplot(., aes(x = as.character(Days), y = CFU, fill =as.factor(Conc))) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(cex=2, alpha = 0.4, position = position_dodge(0.75)) +
  theme_bw(base_size = 16) +
  scale_y_log10(labels=scales::comma, name = " ", limits = c(0.5,10000000)) +
  xlab(" ") +
  scale_fill_manual(name="Condition", 
                     breaks=c("0.05","0.5"),
                     values=c("#a7f542","#4b6e1e")) +
  theme(text = element_text(family = "Arial"), legend.position = "none")
dev.off()

#### Violin Plot - Single Injections ####
png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/DilPl8_Violin_Single.png",width = 4000, height=1500, res=500)
dilpl8_efae %>%
  dplyr::filter(., .$Conc=="0.05" | .$Conc=="0.5") %>%
  ggplot(., aes(x = as.character(Days), y = CFU, fill =as.factor(Conc))) +
  geom_violin(scale = "width", position = position_dodge(0.89), adjust = 0.9) +
  geom_jitter(cex=1, alpha = 0.25, position = position_dodge(0.89)) +
  theme_bw(base_size = 16) +
  scale_y_log10(labels=scales::comma, name = " ", limits = c(0.5,10000000)) +
  xlab(" ") +
  scale_fill_manual(name="Condition", 
                    breaks=c("0.05","0.5"),
                    values=c("#a7f542","#4b6e1e")) +
  theme(text = element_text(family = "Arial"), legend.position = "none")
dev.off()

#### Box Plot - Double Injections ####
png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/DilPl8_Double.png",width = 1750, height=1350, res=300)
dilpl8_efae %>%
  dplyr::filter(., .$Conc=="0.5" | .$Conc=="Trained" | .$Conc=="Mock-Trained") %>%
  dplyr::filter(., .$Days=="0" | .$Days=="1" | .$Days=="2"| .$Days=="3"| .$Days=="4") %>%
  ggplot(., aes(x = as.character(Days), y = CFU, fill =as.factor(Conc))) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(cex=1.5, alpha = 0.4, position = position_dodge(0.75)) +
  theme_bw(base_size = 18) +
  scale_y_log10(labels=scales::comma, name = "CFU/Fly", limits = c(0.5,10000000)) +
  #scale_x_discrete(name=" ", limits = seq(0:6)) +
  xlab("Days Post Injection") +
  scale_fill_manual(name="Condition", 
                     breaks=c("0.5","Mock-Trained","Trained"),
                     labels=c("EfaeHi","Mock-Trained","E.fae-Trained"),
                     values=c("#4b6e1e","#355a78","#1cd9b6")) +
  theme(text = element_text(family = "Arial"), legend.position = "none")
dev.off()

#### Violin Plot - Double Injections ####
png(filename = "/Users/cabweebwub/Desktop/Efae RNA-seq Paper/Figures/DilPl8_Violin_Double.png",width = 3500, height=1500, res=500)
dilpl8_efae %>%
  dplyr::filter(., .$Conc=="Trained" | .$Conc=="Mock-Trained") %>%
  dplyr::filter(., .$Days=="0" | .$Days=="1" | .$Days=="2"| .$Days=="3"| .$Days=="4"| .$Days=="5") %>%
  ggplot(., aes(x = as.character(Days), y = CFU, fill =as.factor(Conc))) +
  geom_violin(scale = "area", position = position_dodge(0.89), adjust = 0.9) +
  geom_jitter(cex=1.5, alpha = 0.4, position = position_dodge(0.89)) +
  theme_bw(base_size = 18) +
  scale_y_log10(labels=scales::comma, name = " ", limits = c(0.5,10000000)) +
  #scale_x_continuous(name=" ", limits = c(0,6)) +
  xlab(" ") +
  scale_fill_manual(name="Condition", 
                    breaks=c("Mock-Trained","Trained"),
                    labels=c("Mock-Trained","E.fae-Trained"),
                    values=c("#355a78","#1cd9b6")) +
  theme(text = element_text(family = "Arial"), legend.position = "none")
dev.off()

dilpl8_efae %>%
  dplyr::filter(., .$Conc=="Trained" | .$Conc=="Mock-Trained") %>%
  kruskal.test(CFU ~ Days, data = .)