# Script created December 3, 2025 by Nichole Ginnan, Assistant Project Scientist (nginn001@ucr.edu)
# Last updated March 26, 2026
# Amplicon sequencing microbiome data analysis
# ECDRE Project, Lead PI: Caroline Roper, funded by the USDA
# Lindcove Research and Extension Center, Exeter, California | 91C plot experiment
# Field experiment testing the impacts of mulch, glyphosate, and humic acid citrus tree microbiome
# 16S bacterial and ITS fungal communities
###############################################################################
# Load Libraries ####
library(tidyverse)
library(writexl)
library(readxl)
library(phyloseq)
library(lme4)
library(lmerTest)
library(scales)
library(emmeans)
library(multcompView)
library(multcomp)
library(RColorBrewer)
library(vegan)
library(microbiome)
library(ALDEx2)
###############################################################################
# Load data objects; centered log-ratio transformed phyloseq objects ####
setwd("/Users/nicholeginnan/Documents/UCR- Current/CA.citrus.paper") # set working directory
bac.clr.ps<-readRDS("data_tables/microbiome/CA.CLR.16S.greengenes.phyloseq.obj.RDS")
fun.clr.ps<-readRDS("data_tables/microbiome/CA.CLR.ITS.UNITE.phyloseq.obj.RDS")
###############################################################################
# subset data #####
# tissue type #
bac.r.ps<- subset_samples(bac.clr.ps, Tissue_type == "Roots")
bac.z.ps<- subset_samples(bac.clr.ps, Tissue_type == "Rhizosphere")
fun.r.ps<- subset_samples(fun.clr.ps, Tissue_type == "Roots")
fun.z.ps<- subset_samples(fun.clr.ps, Tissue_type == "Rhizosphere")
# remove year zero (pre-treatment year and only use years 1-3) #
bac.r.ps.treat <- subset_samples(bac.r.ps, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
bac.z.ps.treat <- subset_samples(bac.z.ps, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.r.ps.treat <- subset_samples(fun.r.ps, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.z.ps.treat <- subset_samples(fun.z.ps, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
###############################################################################
####### Alpha Diversity ###################################################
###############################################################################
# Plot settings ####
plot.theme<- theme(legend.title = element_blank(), 
                   legend.text=element_text(size=10,family="sans"),
                   legend.position="right",
                   legend.box.spacing = unit(2.0, "pt"),
                   legend.spacing.x = unit(10.0, 'pt'),
                   axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
                   strip.text.x= element_text(size=10),  
                   axis.title.y=element_text(size=10, family="sans",vjust = 1), 
                   strip.background=element_blank(),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.spacing=unit(0.5, "lines"), 
                   axis.text.x= element_text(colour="black", size=8, family="sans"), 
                   axis.text.y= element_text(colour="black", size=8, family="sans"))
# Shannon Diversity - Bacteria roots ####
md.bac.r<-data.frame(sample_data(bac.r.ps.treat)) # extract metadata from the phyloseq object
md.bac.r$z_logObs<-scale(md.bac.r$logObs) # create new col with logObs z-score to control for sequencing depth variability
# check for outliers
ggplot(md.bac.r.v2, aes(y=Shannon, x=Mulch, color=Glyphosate, shape=Humic, label = rownames(md.bac.r.v2))) +
  geom_boxplot(outliers = FALSE) +
  facet_wrap(.~Tissue_type)+
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("R.T2.CCC1","R.T2.CGH1","R.T3.MGC4","R.T3.MCH3")) 
md.bac.r.v2<-md.bac.r %>% filter(!rownames(md.bac.r) %in% outliers)
# Mixed Effects Model
model<-lmer(Shannon~Mulch*Glyphosate+Humic + z_logObs + (1|Timepoint), data=data.frame(md.bac.r.v2)) # model 1 
plot(resid(model)~fitted(model)) # looks good
text(x = fitted(model),y = resid(model),labels = rownames(model.frame(model)),pos = 4,cex = 0.6)
anova(model, type = "III") # 
#Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# Mulch            0.3063  0.3063     1 77.409  1.5856    0.2117    
#Glyphosate       0.0928  0.0928     1 78.382  0.4801    0.4904    
#Humic            0.0055  0.0055     1 78.405  0.0283    0.8668    
#z_logObs         8.0621  8.0621     1 59.193 41.7308 2.178e-08 ***
# Mulch:Glyphosate 0.7946  0.7946     1 77.291  4.1131    0.0460 *  
emm<-emmeans(model, specs=~Humic:Mulch:Glyphosate, type="response")
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# plotting emmeans #
cld %>% as.data.frame %>%
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic, group = interaction(Glyphosate, Humic))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.65), color = "black") +
  #geom_point(data=md.bac.r.v2,mapping=aes(x=Mulch, y= Shannon, fill = Glyphosate, color = Glyphosate,shape=Humic), alpha=0.2, position=position_jitterdodge(0.03), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.65)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.65)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 4, vjust = -1, hjust = 0.5,
            position = position_dodge(width = 0.65)) +
  ylab("Shannon Diversity") +
  xlab("Treatment") +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "bac.shannon.rt.mulch.humic.glyphosate.zlogObs.nopoints.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="plots/")

# Shannon Diversity - Bacteria rhizosphere ####
md.bac.z<-data.frame(sample_data(bac.z.ps.treat)) # extract metadata from the phyloseq object
md.bac.z$z_logObs<-scale(md.bac.z$logObs) # create new col with logObs z-score to control for sequencing depth variability
# check for outliers
ggplot(md.bac.z.v2, aes(y=Shannon, x=Mulch, color=Glyphosate, shape=Humic, label = rownames(md.bac.z.v2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("Z.T3.MCH1","Z.T2.CGH2","Z.T3.CCC3","Z.T3.MGC1","Z.T3.MGC3","Z.T3.MGH3","Z.T3.MGH2","Z.T3.CCC4","Z.T3.CCC2","Z.T1.CGC4")) 
md.bac.z.v2<-md.bac.z %>% filter(!rownames(md.bac.z) %in% outliers)
# Mixed Effects Model
model<-lmer(Shannon~Mulch+Humic*Glyphosate + z_logObs + (1|Timepoint), data=data.frame(md.bac.z.v2)) # model 1 
plot(resid(model)~fitted(model)) # looks good
anova(model, type = "III") # 
#Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# Mulch            0.5427  0.5427     1 77.142  8.4934  0.004662 ** 
# Humic            0.2121  0.2121     1 77.032  3.3199  0.072329 .  
# Glyphosate       0.0104  0.0104     1 77.656  0.1624  0.688075    
# z_logObs         4.1462  4.1462     1 78.140 64.8837 7.282e-12 ***
# Humic:Glyphosate 0.2532  0.2532     1 77.529  3.9618  0.050069 .  
emm<-emmeans(model, specs=~Humic:Mulch:Glyphosate, type="response")
cld<-cld(emm, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# plotting emmeans #
cld %>% as.data.frame %>%
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic, group = interaction(Glyphosate, Humic))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.65), color = "black") +
  #geom_point(data=md.bac.z.v2,mapping=aes(x=Mulch, y= Shannon, fill = Glyphosate, color = Glyphosate,shape=Humic), alpha=0.2, position=position_jitterdodge(0.03), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.65)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.65)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 4, vjust = -1, hjust = 0.5,
            position = position_dodge(width = 0.65)) +
  ylab("Shannon Diversity") +
  xlab("Treatment") +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "bac.shannon.rhizo.mulch.humic.glyphosate.zlogObs.nopoints.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="plots/")

# Shannon Diversity - Fungi roots ####
md.fun.r<-data.frame(sample_data(fun.r.ps.treat)) # extract metadata from the phyloseq object
md.fun.r$z_logObs<-scale(md.fun.r$logObs) # create new col with logObs z-score to control for sequencing depth variability
# check for outliers
ggplot(md.fun.r.v2, aes(y=Shannon, x=Mulch, color=Glyphosate, shape=Humic_acid, label = rownames(md.fun.r.v2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("R.T2.CCH3", "R.T1.MCH1")) 
md.fun.r.v2<-md.fun.r %>% filter(!rownames(md.fun.r) %in% outliers)
# Mixed Effects Model
model<-lmer(Shannon~Mulch+Glyphosate*Humic_acid + z_logObs + (1|Timepoint), data=data.frame(md.fun.r.v2)) # model 1 
plot(resid(model)~fitted(model)) # looks good
anova(model, type = "III") # 
#Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# Mulch          0.3521  0.3521     1 82.003  0.6639   0.41755    
#Glyphosate     1.8669  1.8669     1 82.040  3.5206   0.06417 .  
#Humic_acid      0.7915  0.7915     1 82.170  1.4927   0.22530    
#z_logObs       14.1164 14.1164     1 82.823 26.6201 1.665e-06 ***
# Glyphosate:Humic_acid  1.5299  1.5299     1 82.036  2.8851   0.09319 .  
emm<-emmeans(model, specs=~Humic_acid:Mulch:Glyphosate, type="response")
cld<-cld(emm, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# plotting emmeans #
cld %>% as.data.frame %>%
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic_acid, group = interaction(Glyphosate, Humic_acid))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.65), color = "black") +
  #geom_point(data=md.fun.r.v2,mapping=aes(x=Mulch, y= Shannon, fill = Glyphosate, color = Glyphosate,shape=Humic_acid), alpha=0.2, position=position_jitterdodge(0.03), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.65)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.65)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 4, vjust = -1, hjust = 0.5,
            position = position_dodge(width = 0.65)) +
  ylab("Shannon Diversity") +
  xlab("Treatment") +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "fun.shannon.rt.mulch.Humic_acid.glyphosate.zlogObs.nopoints.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="plots/")

# Shannon Diversity - Fungi rhizosphere ####
md.fun.z<-data.frame(sample_data(fun.z.ps.treat)) # extract metadata from the phyloseq object
md.fun.z$z_logObs<-scale(md.fun.z$logObs) # create new col with logObs z-score to control for sequencing depth variability
# check for outliers
ggplot(md.fun.z.v2, aes(y=Shannon, x=Mulch, color=Glyphosate, shape=Humic_acid, label = rownames(md.fun.z.v2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("Z.T2.MCC2")) 
md.fun.z.v2<-md.fun.z %>% filter(!rownames(md.fun.z) %in% outliers)
# Mixed Effects Model
model<-lmer(Shannon~Mulch+Humic_acid+Glyphosate + z_logObs + (1|Timepoint), data=data.frame(md.fun.z.v2)) # model 1 
plot(resid(model)~fitted(model)) # looks good
anova(model, type = "III") # 
#Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# Mulch       0.723   0.723     1 87.071  0.9877   0.32305    
# Humic_acid  0.274   0.274     1 87.303  0.3745   0.54217    
# Glyphosate  2.659   2.659     1 86.965  3.6350   0.05988 .  
# z_logObs   37.414  37.414     1 88.746 51.1467 2.328e-10 ***
emm<-emmeans(model, specs=~Humic_acid:Mulch:Glyphosate, type="response")
cld<-cld(emm, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# plotting emmeans #
cld %>% as.data.frame %>%
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic_acid, group = interaction(Glyphosate, Humic_acid))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.65), color = "black") +
  #geom_point(data=md.fun.z.v2,mapping=aes(x=Mulch, y= Shannon, fill = Glyphosate, color = Glyphosate,shape=Humic_acid), alpha=0.2, position=position_jitterdodge(0.03), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.65)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.65)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 4, vjust = -1, hjust = 0.5,
            position = position_dodge(width = 0.65)) +
  ylab("Shannon Diversity") +
  xlab("Treatment") +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "fun.shannon.rhizo.mulch.humic.glyphosate.zlogObs.nopoints.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="plots/")

# Inverse Simpson Index - Bacteria roots ####
md.bac.r<-data.frame(sample_data(bac.r.ps.treat)) # extract metadata from the phyloseq object
md.bac.r$z_logObs<-scale(md.bac.r$logObs) # create new col with logObs z-score to control for sequencing depth variability
# check for outliers
ggplot(md.bac.r.v2, aes(y=InvSimpson, x=Mulch, color=Glyphosate, shape=Humic, label = rownames(md.bac.r.v2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("R.T1.CGH4","R.T3.MCH2","R.T1.CGC2","R.T1.CGC3","R.T1.MGC4")) 
md.bac.r.v2<-md.bac.r %>% filter(!rownames(md.bac.r) %in% outliers)
# Mixed Effects Model
#model<-lmer((InvSimpson)~Mulch*Glyphosate+Humic + z_logObs + (1|Timepoint), data=data.frame(md.bac.r.v2)) # model 1 
model<-lm(log10(InvSimpson)~Mulch*Glyphosate+Humic + z_logObs, data=data.frame(md.bac.r.v2)) # model 1 
plot(resid(model)~fitted(model)) # looks good
text(x = fitted(model),y = resid(model),labels = rownames(model.frame(model)),pos = 4,cex = 0.6)
anova(model) # 
#Response: log10(InvSimpson)
#                   Df Sum Sq Mean Sq F value    Pr(>F)    
# Mulch             1 0.0669 0.06691  1.3608 0.2469529    
# Glyphosate        1 0.0031 0.00308  0.0626 0.8030281    
# Humic             1 0.0010 0.00100  0.0203 0.8871418    
# z_logObs          1 0.7533 0.75332 15.3213 0.0001928 ***
# Mulch:Glyphosate  1 0.1563 0.15627  3.1782 0.0785190 .  
# Residuals        78 3.8351 0.04917                       
emm<-emmeans(model, specs=~Humic:Mulch:Glyphosate, type="response")
cld<-cld(emm, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# plotting emmeans #
cld %>% as.data.frame %>%
  ggplot(aes(x = Mulch, y = response, fill = Glyphosate, color = Glyphosate,shape=Humic, group = interaction(Glyphosate, Humic))) +
  geom_pointrange(aes(ymin = response-SE, ymax = response+SE),size = 0.8, position = position_dodge(width = 0.65), color = "black") +
  #geom_point(data=md.bac.r.v2,mapping=aes(x=Mulch, y= InvSimpson, fill = Glyphosate, color = Glyphosate,shape=Humic), alpha=0.2, position=position_jitterdodge(0.03), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.65)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.65)) +
  geom_text(aes(label = .group, y = response + SE),
            color = "black", size = 4, vjust = -1, hjust = 0.5,
            position = position_dodge(width = 0.65)) +
  ylab("Inverse Simpson Index") +
  xlab("Treatment") +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "bac.InvSimpson.rt.mulch.humic.glyphosate.zlogObs.nopoints.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="plots/")

# Inverse Simpson Index  - Bacteria rhizosphere ####
md.bac.z<-data.frame(sample_data(bac.z.ps.treat)) # extract metadata from the phyloseq object
md.bac.z$z_logObs<-scale(md.bac.z$logObs) # create new col with logObs z-score to control for sequencing depth variability
# check for outliers
ggplot(md.bac.z.v2, aes(y=InvSimpson, x=Mulch, color=Glyphosate, shape=Humic, label = rownames(md.bac.z.v2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("Z.T2.MCH3","Z.T2.MGH4","Z.T3.MGH1","Z.T3.MCH1")) 
md.bac.z.v2<-md.bac.z %>% filter(!rownames(md.bac.z) %in% outliers)
# Mixed Effects Model
model<-lmer(InvSimpson~Mulch+Humic*Glyphosate + z_logObs + (1|Timepoint), data=data.frame(md.bac.z.v2)) # model 1 
plot(resid(model)~fitted(model)) # looks good
anova(model, type = "III") # 
ranova(model)
#Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# Mulch              5172    5172     1 82.975  2.0643   0.15454    
# Humic              3414    3414     1 82.934  1.3627   0.24642    
# Glyphosate          176     176     1 83.693  0.0701   0.79188    
# z_logObs          63925   63925     1 83.447 25.5146 2.545e-06 ***
# Humic:Glyphosate  17373   17373     1 83.491  6.9342   0.01007 *  
emm<-emmeans(model, specs=~Humic:Mulch:Glyphosate, type="response")
cld<-cld(emm, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# plotting emmeans #
cld %>% as.data.frame %>%
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic, group = interaction(Glyphosate, Humic))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.65), color = "black") +
  #geom_point(data=md.bac.z.v2,mapping=aes(x=Mulch, y= InvSimpson, fill = Glyphosate, color = Glyphosate,shape=Humic), alpha=0.2, position=position_jitterdodge(0.03), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.65)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.65)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 4, vjust = -1, hjust = 0.5,
            position = position_dodge(width = 0.65)) +
  ylab("Inverse Simpson Index") +
  xlab("Treatment") +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "bac.InvSimpson.rhizo.mulch.humic.glyphosate.zlogObs.nopoints.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="plots/")

# Inverse Simpson Index - Fungi roots ####
md.fun.r<-data.frame(sample_data(fun.r.ps.treat)) # extract metadata from the phyloseq object
md.fun.r$z_logObs<-scale(md.fun.r$logObs) # create new col with logObs z-score to control for sequencing depth variability
# check for outliers
ggplot(md.fun.r.v2, aes(y=InvSimpson, x=Mulch, color=Glyphosate, shape=Humic_acid, label = rownames(md.fun.r.v2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("R.T2.MGH1","R.T1.MCH2","R.T2.MCH1","R.T1.CGC2","R.T2.CCH3")) 
md.fun.r.v2<-md.fun.r %>% filter(!rownames(md.fun.r) %in% outliers)
# Mixed Effects Model
model<-lmer(log(InvSimpson)~Mulch+Glyphosate+Humic_acid + z_logObs + (1|Timepoint), data=data.frame(md.fun.r.v2)) # model 1 
plot(resid(model)~fitted(model)) # looks good
anova(model, type = "III") # 
#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# Mulch      0.9838  0.9838     1 79.992  2.5005 0.1177531    
# Glyphosate 0.6263  0.6263     1 80.030  1.5919 0.2107129    
# Humic_acid 0.0177  0.0177     1 80.220  0.0449 0.8326552    
# z_logObs   6.4420  6.4420     1 80.801 16.3742 0.0001183 ***
emm<-emmeans(model, specs=~Humic_acid:Mulch:Glyphosate, type="response")
cld<-cld(emm, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# plotting emmeans #
cld %>% as.data.frame %>%
  ggplot(aes(x = Mulch, y = response, fill = Glyphosate, color = Glyphosate,shape=Humic_acid, group = interaction(Glyphosate, Humic_acid))) +
  geom_pointrange(aes(ymin = response-SE, ymax = response+SE),size = 0.8, position = position_dodge(width = 0.65), color = "black") +
  #geom_point(data=md.fun.r.v2,mapping=aes(x=Mulch, y= InvSimpson, fill = Glyphosate, color = Glyphosate,shape=Humic_acid), alpha=0.2, position=position_jitterdodge(0.03), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.65)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.65)) +
  geom_text(aes(label = .group, y = response + SE),
            color = "black", size = 4, vjust = -1, hjust = 0.5,
            position = position_dodge(width = 0.65)) +
  ylab("Inverse Simpson Index") +
  xlab("Treatment") +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "fun.InvSimpson.rt.mulch.Humic_acid.glyphosate.zlogObs.nopoints.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="plots/")

# Inverse Simpson Index - Fungi rhizosphere ####
md.fun.z<-data.frame(sample_data(fun.z.ps.treat)) # extract metadata from the phyloseq object
md.fun.z$z_logObs<-scale(md.fun.z$logObs) # create new col with logObs z-score to control for sequencing depth variability
# check for outliers
ggplot(md.fun.z.v2, aes(y=InvSimpson, x=Mulch, color=Glyphosate, shape=Humic_acid, label = rownames(md.fun.z.v2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("Z.T2.CCC1","Z.T1.CGC1","Z.T1.CCH3","Z.T3.CCC4","Z.T2.CGH4","Z.T2.MCH4","Z.T1.MGC3","Z.T2.CGC2","Z.T1.CGH2","Z.T2.CGH1")) 
md.fun.z.v2<-md.fun.z %>% filter(!rownames(md.fun.z) %in% outliers)
# Mixed Effects Model
model<-lmer(log(InvSimpson)~Mulch+Humic_acid+Glyphosate + z_logObs + (1|Timepoint), data=data.frame(md.fun.z.v2)) # model 1 
plot(resid(model)~fitted(model)) # looks good
anova(model, type = "III") # 
#Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# Mulch      2.8772  2.8772     1 78.006  5.8707   0.01771 *  
# Humic_acid 0.6119  0.6119     1 78.181  1.2484   0.26727    
# Glyphosate 1.7228  1.7228     1 77.978  3.5151   0.06455 .  
# z_logObs   9.5529  9.5529     1 78.846 19.4918 3.175e-05 ***
emm<-emmeans(model, specs=~Mulch:Glyphosate|Humic_acid, type="response")
cld<-cld(emm, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# plotting emmeans #
cld %>% as.data.frame %>%
  ggplot(aes(x = Mulch, y = response, fill = Glyphosate, color = Glyphosate,shape=Humic_acid, group = interaction(Glyphosate, Humic_acid))) +
  geom_pointrange(aes(ymin = response-SE, ymax = response+SE),size = 0.8, position = position_dodge(width = 0.65), color = "black") +
  #geom_point(data=md.fun.z.v2,mapping=aes(x=Mulch, y= InvSimpson, fill = Glyphosate, color = Glyphosate,shape=Humic_acid), alpha=0.2, position=position_jitterdodge(0.03), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.65)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.65)) +
  geom_text(aes(label = .group, y = response + SE),
            color = "black", size = 4, vjust = -1, hjust = 0.5,
            position = position_dodge(width = 0.65)) +
  ylab("Inverse Simpson Index") +
  xlab("Treatment") +
  facet_wrap(.~Humic_acid)+
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "fun.InvSimpson.rhizo.mulch.humic.glyphosate.split.zlogObs.nopoints.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="plots/")

###############################################################################
####### Beta Diversity - constrained #############
###############################################################################
# Bacteria roots ####
cap.r <- ordinate(bac.r.ps.treat, method='CAP',distance='euclidean',formula=~Mulch*Glyphosate+Humic+Condition(scale(logObs))+Condition(Timepoint))
anova.cca(cap.r) #significant
#         Df Variance      F Pr(>F)    
# Model     4   208.13 1.5803  0.001 ***
# Residual 81  2667.10                  
anova.cca(cap.r,by="term")
#             Df  Variance    F   Pr(>F)    
# Mulch             1    80.03 2.4306  0.001 ***
# Glyphosate        1    39.99 1.2146  0.085 .  
# Humic             1    49.04 1.4892  0.012 *  
# Mulch:Glyphosate  1    39.07 1.1866  0.089 .  
# Residual         81  2667.10                        
anova.cca(cap.r,by="axis") #  CAP1 and CAP2 are significant
summary(cap.r)
# Partitioning of mean squared Euclidean distance:
#               Inertia Proportion
#Total           3508.1    1.00000
#Conditioned     632.9    0.18040
#Constrained     208.1    0.05933 <--
#Unconstrained  2667.1    0.76027
#                          CAP1      CAP2 
#Eigenvalue            81.71426 53.06250
#Proportion Explained    0.02842  0.01846
#Cumulative Proportion   0.02842  0.04688
RDAscores.g <- data.frame(scores(cap.r, c(1,2,3,4),display='sites'))
cap.r.df<- cbind(sample_data(bac.r.ps.treat),RDAscores.g)
# Plot
ggplot(cap.r.df,aes(x=CAP1,y=CAP2,color=Mulch,fill = Mulch, shape=Humic))+
  geom_point(size = 3) +
  stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05,aes(fill = Mulch)) +
  scale_shape_manual(values=c(1, 17),limits=c("humic","no"),labels = c("humic acid","no humic acid")) +
  scale_color_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  xlab("CAP 1 [0.028%]") +
  facet_wrap(.~Glyphosate) +
  ylab("CAP 2 [0.018%]") +
  theme_classic() +
  theme(legend.title = element_text(size=10,family="sans"), 
        legend.text=element_text(size=10,family="sans"),
        legend.position="right",
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.x = unit(5.0, 'pt'),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"), 
        axis.text.x= element_text(colour="black", size=8, family="sans"), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "CAP.bac.root.humic.gly.mulch.classic.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="plots/")

# Bacteria rhizosphere ####
cap.z <- ordinate(bac.z.ps.treat, method='CAP',distance='euclidean',formula=~Mulch*Glyphosate+Mulch*Humic+Condition(scale(logObs))+Condition(Timepoint))
anova.cca(cap.z) #significant
#          Df Variance   F    Pr(>F)    
#Model     5    585.1 1.9817  0.001 ***
#Residual 86   5078.6                  
anova.cca(cap.z,by="term")
#                  Df  Variance    F   Pr(>F)    
# Mulch             1    226.8 3.8411  0.001 ***
# Glyphosate        1     70.2 1.1892  0.099 .  
# Humic             1    101.4 1.7170  0.002 ** 
# Mulch:Glyphosate  1     99.1 1.6775  0.001 ***
# Mulch:Humic       1     87.6 1.4839  0.014 *  
# Residual         86   5078.6                                       
anova.cca(cap.z,by="axis") # CAP1 and CAP2 are significant
summary(cap.z)
# Partitioning of mean squared Euclidean distance:
#               Inertia Proportion
#Total           7132.3    1.00000
#Conditioned    1468.5    0.20590
#Constrained     585.1    0.08204 <--
#Unconstrained  5078.6    0.71206
#                          CAP1      CAP2 
#Eigenvalue            234.18586 142.35597
#Proportion Explained    0.04135   0.02513
#Cumulative Proportion   0.04135   0.06648
RDAscores.g <- data.frame(scores(cap.z, c(1,2,3,4),display='sites'))
cap.z.df<- cbind(sample_data(bac.z.ps.treat),RDAscores.g)
# Plot
ggplot(cap.z.df,aes(x=CAP1,y=CAP2,color=Mulch,fill = Mulch, shape=Humic))+
  geom_point(size = 3) +
  stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05,aes(fill = Mulch)) +
  scale_shape_manual(values=c(1, 17),limits=c("humic","no"),labels = c("humic acid","no humic acid")) +
  scale_color_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  xlab("CAP 1 [0.041%]") +
  facet_wrap(.~Glyphosate) +
  ylab("CAP 2 [0.025%]") +
  theme_classic() +
  theme(legend.title = element_text(size=10,family="sans"), 
        legend.text=element_text(size=10,family="sans"),
        legend.position="right",
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.x = unit(5.0, 'pt'),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"), 
        axis.text.x= element_text(colour="black", size=8, family="sans"), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "CAP.bac.rhizosphere.humic.gly.mulch.classic.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="plots/")

# Fungi roots ####
cap.r <- ordinate(fun.r.ps.treat, method='CAP',distance='euclidean',formula=~Mulch+Glyphosate+Humic_acid+Condition(scale(logObs))+Condition(Timepoint))
anova.cca(cap.r) #significant
#         Df Variance      F Pr(>F)    
#Model     3    55.37 1.5893  0.001 ***
#Residual 85   987.15                  
anova.cca(cap.r,by="term")
#             Df  Variance    F   Pr(>F)    
# Mulch       1   25.41 2.1879   0.004 **
# Glyphosate  1    17.60 1.5152  0.042 * 
# Humic_acid  1    12.37 1.0649  0.277   
# Residual   85   987.15                                 
anova.cca(cap.r,by="axis") # only CAP1 is significant
summary(cap.r)
# Partitioning of mean squared Euclidean distance:
#               Inertia Proportion
#Total           1293.20    1.00000
#Conditioned    250.68    0.19385
#Constrained     55.37    0.04282 <--
#Unconstrained  987.15    0.76334
#                          CAP1      CAP2 
#Eigenvalue            27.24852 17.14005
#Proportion Explained    0.02614  0.01644
#Cumulative Proportion   0.02614  0.04258
RDAscores.g <- data.frame(scores(cap.r, c(1,2,3,4),display='sites'))
cap.r.df<- cbind(sample_data(fun.r.ps.treat),RDAscores.g)
# CAP1 plot
ggplot(cap.r.df,aes(x=Mulch,y=CAP1,color=Glyphosate,fill = Glyphosate, shape=Humic_acid))+
  geom_boxplot(alpha=0.5) +
  geom_point(position = position_jitterdodge(0.1)) +
  theme_bw()
# Plot
ggplot(cap.r.df,aes(x=CAP1,y=CAP2,color=Mulch,fill = Mulch, shape=Humic_acid))+
  stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05,aes(fill = Mulch)) +
  geom_point(size = 3) +
  #scale_shape_manual(values=c(17, 2,16,3),limits=c("humic.Glyphosate","no.Glyphosate","humic.no","no.no"),labels = c("humic acid & glyphosate","glyphosate","humic acid","control")) +
  scale_shape_manual(values=c(1, 17),limits=c("humic","no"),labels = c("humic acid","no humic acid")) +
  scale_color_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  xlab("CAP 1 [0.026%]") +
  facet_wrap(.~Glyphosate) +
  ylab("CAP 2 [0.016%]") +
  theme_classic() +
  theme(legend.title = element_text(size=10,family="sans"), 
        legend.text=element_text(size=10,family="sans"),
        legend.position="right",
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.x = unit(5.0, 'pt'),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        #strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"), 
        axis.text.x= element_text(colour="black", size=8, family="sans"), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "CAP.fun.root.humic.gly.mulch.classic.box.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="plots/")

# Fungi rhizosphere ####
cap.z <- ordinate(fun.z.ps.treat, method='CAP',distance='euclidean',formula=~Mulch+Glyphosate+Humic_acid+Condition(scale(logObs))+Condition(Timepoint))
anova.cca(cap.z) #significant
#         Df Variance      F Pr(>F)    
# Model     3   176.26 2.4285  0.001 ***
# Residual 88  2128.95                  
anova.cca(cap.z,by="term")
#             Df  Variance    F   Pr(>F)    
# Mulch       1   119.21 4.9277  0.001 ***
# Glyphosate  1    35.87 1.4826  0.061 .  
# Humic_acid  1    21.18 0.8753  0.591    
# Residual   88  2128.95                                   
anova.cca(cap.z,by="axis") # only CAP1 is significant
summary(cap.z)
# Partitioning of mean squared Euclidean distance:
#               Inertia Proportion
#Total           3052.1    1.00000
#Conditioned     746.9    0.24471
#Constrained     176.3    0.05775 <--
#Unconstrained  2128.9    0.69754
#                          CAP1      CAP2 
#Eigenvalue            119.24568 35.94806
#Proportion Explained    0.05173  0.01559
#Cumulative Proportion   0.05173  0.06732
RDAscores.g <- data.frame(scores(cap.z, c(1,2,3,4),display='sites'))
cap.z.df<- cbind(sample_data(fun.z.ps.treat),RDAscores.g)
# Plot
ggplot(cap.z.df,aes(x=CAP1,y=CAP2,color=Mulch,fill = Mulch, shape=Humic_acid))+
  geom_point(size = 3) +
  stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05,aes(fill = Mulch)) +
  scale_shape_manual(values=c(1, 17),limits=c("humic","no"),labels = c("humic acid","no humic acid")) +
  scale_color_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  xlab("CAP 1 [0.052%]") +
  facet_wrap(.~Glyphosate) +
  ylab("CAP 2 [0.016%]") +
  theme_classic() +
  theme(legend.title = element_text(size=10,family="sans"), 
        legend.text=element_text(size=10,family="sans"),
        legend.position="right",
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.x = unit(5.0, 'pt'),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"), 
        axis.text.x= element_text(colour="black", size=8, family="sans"), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "CAP.fun.rhizosphere.humic.gly.mulch.classic.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="plots/")

##########################################################################################
# Taxa enrichment | ALDeX2 | Mulch vs. no mulch, Gly vs. no gly, & humic vs. no humic ###########
##########################################################################################
# Load data objects; raw counts phyloseq objects ####
bac.raw.ps<-readRDS("data_tables/microbiome/CA.raw.counts.16S.greengenes.phyloseq.obj.RDS")
fun.raw.ps<-readRDS("data_tables/microbiome/CA.raw.counts.ITS.UNITE.phyloseq.obj.RDS")
# subset data #####
# tissue type #
bac.r.ps.raw<- subset_samples(bac.raw.ps, Tissue_type == "Roots")
bac.z.ps.raw<- subset_samples(bac.raw.ps, Tissue_type == "Rhizosphere")
fun.r.ps.raw<- subset_samples(fun.raw.ps, Tissue_type == "Roots")
fun.z.ps.raw<- subset_samples(fun.raw.ps, Tissue_type == "Rhizosphere")
# remove year zero (pre-treatment) #
bac.r.ps.raw.treat <- subset_samples(bac.r.ps.raw, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
bac.z.ps.raw.treat <- subset_samples(bac.z.ps.raw, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.r.ps.raw.treat <- subset_samples(fun.r.ps.raw, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.z.ps.raw.treat <- subset_samples(fun.z.ps.raw, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 

#############################################
#### ASV level enrichment #################
###########################################
#### Set df for factors being tested ####
# bac root #
mulch.bac.r.df <- sample_data(bac.r.ps.raw.treat)$Mulch %>% as.character() # 2 levels
humic.bac.r.df <- sample_data(bac.r.ps.raw.treat)$Humic %>% as.character() # 2 levels
gly.bac.r.df <- sample_data(bac.r.ps.raw.treat)$Glyphosate %>% as.character() # 2 levels
# bac rhiz #
mulch.bac.z.df <- sample_data(bac.z.ps.raw.treat)$Mulch %>% as.character() # 2 levels
humic.bac.z.df <- sample_data(bac.z.ps.raw.treat)$Humic %>% as.character() # 2 levels
gly.bac.z.df <- sample_data(bac.z.ps.raw.treat)$Glyphosate %>% as.character() # 2 levels
# fun root #
mulch.fun.r.df <- sample_data(fun.r.ps.raw.treat)$Mulch %>% as.character() # 2 levels
humic.fun.r.df <- sample_data(fun.r.ps.raw.treat)$Humic %>% as.character() # 2 levels
gly.fun.r.df <- sample_data(fun.r.ps.raw.treat)$Glyphosate %>% as.character() # 2 levels
# fun rhiz #
mulch.fun.z.df <- sample_data(fun.z.ps.raw.treat)$Mulch %>% as.character() # 2 levels
humic.fun.z.df <- sample_data(fun.z.ps.raw.treat)$Humic %>% as.character() # 2 levels
gly.fun.z.df <- sample_data(fun.z.ps.raw.treat)$Glyphosate %>% as.character() # 2 levels
#### Extract taxa tables ####
ASV.bac.r <- as(otu_table(bac.r.ps.raw.treat, taxa_are_rows = TRUE), 'matrix') 
ASV.bac.z <- as(otu_table(bac.z.ps.raw.treat, taxa_are_rows = TRUE), 'matrix') 
ASV.fun.r <- as(otu_table(fun.r.ps.raw.treat, taxa_are_rows = TRUE), 'matrix') 
ASV.fun.z <- as(otu_table(fun.z.ps.raw.treat, taxa_are_rows = TRUE), 'matrix') 
ASV.bac.r <- t(ASV.bac.r) # run this if aldex.clr() gives error "mismatch between number of samples and condition vector"
ASV.bac.z <- t(ASV.bac.z)
ASV.fun.r <- t(ASV.fun.r)
ASV.fun.z <- t(ASV.fun.z)
#### 2-level factors aldex CLR/MC sampling ####
# bac root
mulch.bac.r.clr <- aldex.clr(ASV.bac.r, conds=mulch.bac.r.df, mc.samples=200) # 
humic.bac.r.clr <- aldex.clr(ASV.bac.r, conds=humic.bac.r.df, mc.samples=200) # 
gly.bac.r.clr <- aldex.clr(ASV.bac.r, conds=gly.bac.r.df, mc.samples=200) # 
# bac rhizo
mulch.bac.z.clr <- aldex.clr(ASV.bac.z, conds=mulch.bac.z.df, mc.samples=200) # 
humic.bac.z.clr <- aldex.clr(ASV.bac.z, conds=humic.bac.z.df, mc.samples=200) # 
gly.bac.z.clr <- aldex.clr(ASV.bac.z, conds=gly.bac.z.df, mc.samples=200) # 
# fun root
mulch.fun.r.clr <- aldex.clr(ASV.fun.r, conds=mulch.fun.r.df, mc.samples=200) # 
humic.fun.r.clr <- aldex.clr(ASV.fun.r, conds=humic.fun.r.df, mc.samples=200) # 
gly.fun.r.clr <- aldex.clr(ASV.fun.r, conds=gly.fun.r.df, mc.samples=200) # 
# fun rhizo
mulch.fun.z.clr <- aldex.clr(ASV.fun.z, conds=mulch.fun.z.df, mc.samples=200) # 
humic.fun.z.clr <- aldex.clr(ASV.fun.z, conds=humic.fun.z.df, mc.samples=200) # 
gly.fun.z.clr <- aldex.clr(ASV.fun.z, conds=gly.fun.z.df, mc.samples=200) # 
#### 2 level factor statistical test ####
#kw.gly.bac.z <- aldex.kw(gly.bac.z.clr) # other test option, but it takes a long time to run and gave a similar result
# bac root
ttest.mulch.bac.r <- aldex.ttest(mulch.bac.r.clr)
ttest.humic.bac.r <- aldex.ttest(humic.bac.r.clr)
ttest.gly.bac.r <- aldex.ttest(gly.bac.r.clr)
# bac rhizo
ttest.mulch.bac.z <- aldex.ttest(mulch.bac.z.clr)
ttest.humic.bac.z <- aldex.ttest(humic.bac.z.clr)
ttest.gly.bac.z <- aldex.ttest(gly.bac.z.clr)
# fun root
ttest.mulch.fun.r <- aldex.ttest(mulch.fun.r.clr)
ttest.humic.fun.r <- aldex.ttest(humic.fun.r.clr)
ttest.gly.fun.r <- aldex.ttest(gly.fun.r.clr)
# fun rhizo
ttest.mulch.fun.z <- aldex.ttest(mulch.fun.z.clr)
ttest.humic.fun.z <- aldex.ttest(humic.fun.z.clr)
ttest.gly.fun.z <- aldex.ttest(gly.fun.z.clr)

#### check for significant ASVs; BH = adjusted p-value ####
# bac root
sum(ttest.mulch.bac.r$wi.eBH < 0.10) # Wilcoxon test #2
sum(ttest.mulch.bac.r$we.eBH < 0.10) # Welch’s t test #1
sum(ttest.humic.bac.r$wi.eBH < 0.10) # #0
sum(ttest.humic.bac.r$we.eBH < 0.10) # #0
sum(ttest.gly.bac.r$wi.eBH < 0.10) # #0
sum(ttest.gly.bac.r$we.eBH < 0.10) # #0
which(ttest.mulch.bac.r$wi.eBH < 0.10) # 
ttest.mulch.bac.r[751,] #                 we.ep     we.eBH       wi.ep    wi.eBH
#                         bASV_454 1.821362e-05 0.01459005 3.030338e-06 0.006861652
#                         bASV_1248 0.0005080577 0.1324709 4.544337e-05 0.04883881
# bac rhizo
sum(ttest.mulch.bac.z$wi.eBH < 0.051) # Wilcoxon test #65
sum(ttest.mulch.bac.z$we.eBH < 0.051) # Welch’s t test #38
sum(ttest.humic.bac.z$wi.eBH < 0.1) # #0
sum(ttest.humic.bac.z$we.eBH < 0.1) # #0
sum(ttest.gly.bac.z$wi.eBH < 0.1) # #0
sum(ttest.gly.bac.z$we.eBH < 0.1) # #0
#which(kw.mulch.bac.z$kw.eBH < 0.1) # 
# fun root
sum(ttest.mulch.fun.r$wi.eBH < 0.10) # Wilcoxon test #2
sum(ttest.mulch.fun.r$we.eBH < 0.10) # Welch’s t test #2
sum(ttest.humic.fun.r$wi.eBH < 0.10) # #0
sum(ttest.humic.fun.r$we.eBH < 0.10) # #0
sum(ttest.gly.fun.r$wi.eBH < 0.10) # #0
sum(ttest.gly.fun.r$we.eBH < 0.10) # #0
which(ttest.mulch.fun.r$wi.eBH < 0.10) # line 57 75
ttest.mulch.fun.r[75,] #                 we.ep     we.eBH       wi.ep    wi.eBH
#                         fASV_82 5.829072e-07 0.0004250923 2.972437e-06 0.002655061 candolleomyces
#                         fASV_110 3.096898e-06 0.002065975 3.973622e-05 0.02239465  candolleomyces
# fun rhizo
sum(ttest.mulch.fun.z$wi.eBH < 0.051) # Wilcoxon test #30
sum(ttest.mulch.fun.z$we.eBH < 0.051) # Welch’s t test #30
sum(ttest.humic.fun.z$wi.eBH < 0.10) # #0
sum(ttest.humic.fun.z$we.eBH < 0.10) # #0
sum(ttest.gly.fun.z$wi.eBH < 0.10) # #0
sum(ttest.gly.fun.z$we.eBH < 0.10) # #0
###########################################
# merge results with sample data and taxa table and write tables
# bac root mulch
tax <- data.frame(tax_table(bac.r.ps.raw.treat)) %>% rownames_to_column("ASV")
ttest.mulch.bac.r.df <- ttest.mulch.bac.r %>% data.frame() %>% rownames_to_column("ASV")
results.mulch.bac.r <- dplyr::left_join(tax,ttest.mulch.bac.r.df,by = "ASV")
sig.mulch.bac.r <- results.mulch.bac.r %>% filter(wi.eBH < 0.051)
write.table(sig.mulch.bac.r, file="data_tables/taxa_enrichments/ASV.sig.mulch.bac.r.csv", sep=",")
# bac rhizosphere mulch
tax <- data.frame(tax_table(bac.z.ps.raw.treat)) %>% rownames_to_column("ASV")
ttest.mulch.bac.z.df <- ttest.mulch.bac.z %>% data.frame() %>% rownames_to_column("ASV")
results.mulch.bac.z <- dplyr::left_join(tax,ttest.mulch.bac.z.df,by = "ASV")
sig.mulch.bac.z <- results.mulch.bac.z %>% filter(wi.eBH < 0.051)
write.table(sig.mulch.bac.z, file="data_tables/taxa_enrichments/ASV.sig.mulch.bac.z.csv", sep=",")
# fun root mulch
tax <- data.frame(tax_table(fun.r.ps.raw.treat)) %>% rownames_to_column("ASV")
ttest.mulch.fun.r.df <- ttest.mulch.fun.r %>% data.frame() %>% rownames_to_column("ASV")
results.mulch.fun.r <- dplyr::left_join(tax,ttest.mulch.fun.r.df,by = "ASV")
sig.mulch.fun.r <- results.mulch.fun.r %>% filter(wi.eBH < 0.051)
write.table(sig.mulch.fun.r, file="data_tables/taxa_enrichments/ASV.sig.mulch.fun.r.csv", sep=",")
# fun rhizosphere mulch
tax <- data.frame(tax_table(fun.z.ps.raw.treat)) %>% rownames_to_column("ASV")
ttest.mulch.fun.z.df <- ttest.mulch.fun.z %>% data.frame() %>% rownames_to_column("ASV")
results.mulch.fun.z <- dplyr::left_join(tax,ttest.mulch.fun.z.df,by = "ASV")
sig.mulch.fun.z <- results.mulch.fun.z %>% filter(wi.eBH < 0.051)
write.table(sig.mulch.fun.z, file="data_tables/taxa_enrichments/ASV.sig.mulch.fun.z.csv", sep=",")

#############################################
#### Genus level enrichment #################
###########################################
#### agglomerate to Genus level ####
# bac root #
bac.r.ps.raw.treat.genus <- tax_glom(bac.r.ps.raw.treat, taxrank = "Genus", NArm = FALSE)
tax <- tax_table(bac.r.ps.raw.treat.genus)
taxa_names(bac.r.ps.raw.treat.genus) <- make.unique(as.character(tax[, "Genus"]))
# bac rhzo #
bac.z.ps.raw.treat.genus <- tax_glom(bac.z.ps.raw.treat, taxrank = "Genus", NArm = FALSE)
tax <- tax_table(bac.z.ps.raw.treat.genus)
taxa_names(bac.z.ps.raw.treat.genus) <- make.unique(as.character(tax[, "Genus"]))
# fun root #
fun.r.ps.raw.treat.genus <- tax_glom(fun.r.ps.raw.treat, taxrank = "Genus", NArm = FALSE)
tax <- tax_table(fun.r.ps.raw.treat.genus)
taxa_names(fun.r.ps.raw.treat.genus) <- make.unique(as.character(tax[, "Genus"]))
# fun rhzo #
fun.z.ps.raw.treat.genus <- tax_glom(fun.z.ps.raw.treat, taxrank = "Genus", NArm = FALSE)
tax <- tax_table(fun.z.ps.raw.treat.genus)
taxa_names(fun.z.ps.raw.treat.genus) <- make.unique(as.character(tax[, "Genus"]))
#### Set df for factors being tested ####
# bac root #
mulch.bac.r.df <- sample_data(bac.r.ps.raw.treat.genus)$Mulch %>% as.character() # 2 levels
humic.bac.r.df <- sample_data(bac.r.ps.raw.treat.genus)$Humic %>% as.character() # 2 levels
gly.bac.r.df <- sample_data(bac.r.ps.raw.treat.genus)$Glyphosate %>% as.character() # 2 levels
# bac rhiz #
mulch.bac.z.df <- sample_data(bac.z.ps.raw.treat.genus)$Mulch %>% as.character() # 2 levels
humic.bac.z.df <- sample_data(bac.z.ps.raw.treat.genus)$Humic %>% as.character() # 2 levels
gly.bac.z.df <- sample_data(bac.z.ps.raw.treat.genus)$Glyphosate %>% as.character() # 2 levels
# fun root #
mulch.fun.r.df <- sample_data(fun.r.ps.raw.treat.genus)$Mulch %>% as.character() # 2 levels
humic.fun.r.df <- sample_data(fun.r.ps.raw.treat.genus)$Humic %>% as.character() # 2 levels
gly.fun.r.df <- sample_data(fun.r.ps.raw.treat.genus)$Glyphosate %>% as.character() # 2 levels
# fun rhiz #
mulch.fun.z.df <- sample_data(fun.z.ps.raw.treat.genus)$Mulch %>% as.character() # 2 levels
humic.fun.z.df <- sample_data(fun.z.ps.raw.treat.genus)$Humic %>% as.character() # 2 levels
gly.fun.z.df <- sample_data(fun.z.ps.raw.treat.genus)$Glyphosate %>% as.character() # 2 levels
#### Extract taxa tables ####
ASV.bac.r <- as(otu_table(bac.r.ps.raw.treat.genus, taxa_are_rows = TRUE), 'matrix') 
ASV.bac.z <- as(otu_table(bac.z.ps.raw.treat.genus, taxa_are_rows = TRUE), 'matrix') 
ASV.fun.r <- as(otu_table(fun.r.ps.raw.treat.genus, taxa_are_rows = TRUE), 'matrix') 
ASV.fun.z <- as(otu_table(fun.z.ps.raw.treat.genus, taxa_are_rows = TRUE), 'matrix') 
ASV.bac.r <- t(ASV.bac.r) # run this if aldex.clr() gives error "mismatch between number of samples and condition vector"
ASV.bac.z <- t(ASV.bac.z)
ASV.fun.r <- t(ASV.fun.r)
ASV.fun.z <- t(ASV.fun.z)
#### 2-level factors aldex CLR/MC sampling ####
# bac root
mulch.bac.r.clr <- aldex.clr(ASV.bac.r, conds=mulch.bac.r.df, mc.samples=200) # 
humic.bac.r.clr <- aldex.clr(ASV.bac.r, conds=humic.bac.r.df, mc.samples=200) # 
gly.bac.r.clr <- aldex.clr(ASV.bac.r, conds=gly.bac.r.df, mc.samples=200) # 
# bac rhizo
mulch.bac.z.clr <- aldex.clr(ASV.bac.z, conds=mulch.bac.z.df, mc.samples=200) # 
humic.bac.z.clr <- aldex.clr(ASV.bac.z, conds=humic.bac.z.df, mc.samples=200) # 
gly.bac.z.clr <- aldex.clr(ASV.bac.z, conds=gly.bac.z.df, mc.samples=200) # 
# fun root
mulch.fun.r.clr <- aldex.clr(ASV.fun.r, conds=mulch.fun.r.df, mc.samples=200) # 
humic.fun.r.clr <- aldex.clr(ASV.fun.r, conds=humic.fun.r.df, mc.samples=200) # 
gly.fun.r.clr <- aldex.clr(ASV.fun.r, conds=gly.fun.r.df, mc.samples=200) # 
# fun rhizo
mulch.fun.z.clr <- aldex.clr(ASV.fun.z, conds=mulch.fun.z.df, mc.samples=200) # 
humic.fun.z.clr <- aldex.clr(ASV.fun.z, conds=humic.fun.z.df, mc.samples=200) # 
gly.fun.z.clr <- aldex.clr(ASV.fun.z, conds=gly.fun.z.df, mc.samples=200) # 
#### 2 level factor statistical test ####
#kw.gly.bac.z <- aldex.kw(gly.bac.z.clr) # other test option, but it takes a long time to run and gave a similar result
# bac root
ttest.mulch.bac.r <- aldex.ttest(mulch.bac.r.clr)
ttest.humic.bac.r <- aldex.ttest(humic.bac.r.clr)
ttest.gly.bac.r <- aldex.ttest(gly.bac.r.clr)
# bac rhizo
ttest.mulch.bac.z <- aldex.ttest(mulch.bac.z.clr)
ttest.humic.bac.z <- aldex.ttest(humic.bac.z.clr)
ttest.gly.bac.z <- aldex.ttest(gly.bac.z.clr)
# fun root
ttest.mulch.fun.r <- aldex.ttest(mulch.fun.r.clr)
ttest.humic.fun.r <- aldex.ttest(humic.fun.r.clr)
ttest.gly.fun.r <- aldex.ttest(gly.fun.r.clr)
# fun rhizo
ttest.mulch.fun.z <- aldex.ttest(mulch.fun.z.clr)
ttest.humic.fun.z <- aldex.ttest(humic.fun.z.clr)
ttest.gly.fun.z <- aldex.ttest(gly.fun.z.clr)

#### check for significant ASVs; BH = adjusted p-value ####
# bac root
sum(ttest.mulch.bac.r$wi.eBH < 0.051) # Wilcoxon test #2
sum(ttest.mulch.bac.r$we.eBH < 0.051) # Welch’s t test #1
sum(ttest.humic.bac.r$wi.eBH < 0.10) # #2
sum(ttest.humic.bac.r$we.eBH < 0.10) # #0
sum(ttest.gly.bac.r$wi.eBH < 0.10) # #0
sum(ttest.gly.bac.r$we.eBH < 0.10) # #0
which(ttest.mulch.bac.r$wi.eBH < 0.10) # 
ttest.mulch.bac.r[751,] #           
# bac rhizo
sum(ttest.mulch.bac.z$wi.eBH < 0.051) # Wilcoxon test #59
sum(ttest.mulch.bac.z$we.eBH < 0.051) # Welch’s t test #22
sum(ttest.humic.bac.z$wi.eBH < 0.1) # #0
sum(ttest.humic.bac.z$we.eBH < 0.1) # #0
sum(ttest.gly.bac.z$wi.eBH < 0.1) # #0
sum(ttest.gly.bac.z$we.eBH < 0.1) # #0
#which(kw.mulch.bac.z$kw.eBH < 0.1) # 
# fun root
sum(ttest.mulch.fun.r$wi.eBH < 0.051) # Wilcoxon test #2
sum(ttest.mulch.fun.r$we.eBH < 0.051) # Welch’s t test #3
sum(ttest.humic.fun.r$wi.eBH < 0.10) # #0
sum(ttest.humic.fun.r$we.eBH < 0.10) # #0
sum(ttest.gly.fun.r$wi.eBH < 0.10) # #0
sum(ttest.gly.fun.r$we.eBH < 0.10) # #0
which(ttest.mulch.fun.r$wi.eBH < 0.10) # 57 75
ttest.mulch.fun.r[57,] #                        we.ep     we.eBH       wi.ep    wi.eBH
#                         Candolleomyces 5.611166e-08 1.560816e-05 7.260413e-07 0.0001236018
#                         Coniochaeta 0.0003533539 0.02010949 0.0002239974 0.01740064
# fun rhizo
sum(ttest.mulch.fun.z$wi.eBH < 0.051) # Wilcoxon test #27
sum(ttest.mulch.fun.z$we.eBH < 0.051) # Welch’s t test #23
sum(ttest.humic.fun.z$wi.eBH < 0.10) # #0
sum(ttest.humic.fun.z$we.eBH < 0.10) # #0
sum(ttest.gly.fun.z$wi.eBH < 0.10) # #0
sum(ttest.gly.fun.z$we.eBH < 0.10) # #0
###########################################
# merge results with sample data and taxa table and write tables
# bac root mulch
tax <- data.frame(tax_table(bac.r.ps.raw.treat.genus)) %>% rownames_to_column("ASV")
ttest.mulch.bac.r.df <- ttest.mulch.bac.r %>% data.frame() %>% rownames_to_column("ASV")
results.mulch.bac.r <- dplyr::left_join(tax,ttest.mulch.bac.r.df,by = "ASV")
sig.mulch.bac.r <- results.mulch.bac.r %>% filter(wi.eBH < 0.051)
write.table(sig.mulch.bac.r, file="data_tables/taxa_enrichments/genus.sig.mulch.bac.r.csv", sep=",")
write.table(results.mulch.bac.r, file="data_tables/taxa_enrichments/genus.all.mulch.bac.r.taxa.csv", sep=",")

# bac rhizosphere mulch
tax <- data.frame(tax_table(bac.z.ps.raw.treat.genus)) %>% rownames_to_column("ASV")
ttest.mulch.bac.z.df <- ttest.mulch.bac.z %>% data.frame() %>% rownames_to_column("ASV")
results.mulch.bac.z <- dplyr::left_join(tax,ttest.mulch.bac.z.df,by = "ASV")
sig.mulch.bac.z <- results.mulch.bac.z %>% filter(wi.eBH < 0.051)
write.table(sig.mulch.bac.z, file="data_tables/taxa_enrichments/genus.sig.mulch.bac.z.csv", sep=",")
write.table(results.mulch.bac.z, file="data_tables/taxa_enrichments/genus.all.mulch.bac.z.taxa.csv", sep=",")

# fun root mulch
tax <- data.frame(tax_table(fun.r.ps.raw.treat.genus)) %>% rownames_to_column("ASV")
ttest.mulch.fun.r.df <- ttest.mulch.fun.r %>% data.frame() %>% rownames_to_column("ASV")
results.mulch.fun.r <- dplyr::left_join(tax,ttest.mulch.fun.r.df,by = "ASV")
sig.mulch.fun.r <- results.mulch.fun.r %>% filter(wi.eBH < 0.051)
write.table(sig.mulch.fun.r, file="data_tables/taxa_enrichments/genus.sig.mulch.fun.r.csv", sep=",")
write.table(results.mulch.fun.r, file="data_tables/taxa_enrichments/genus.all.mulch.fun.r.taxa.csv", sep=",")

# fun rhizosphere mulch
tax <- data.frame(tax_table(fun.z.ps.raw.treat.genus)) %>% rownames_to_column("ASV")
ttest.mulch.fun.z.df <- ttest.mulch.fun.z %>% data.frame() %>% rownames_to_column("ASV")
results.mulch.fun.z <- dplyr::left_join(tax,ttest.mulch.fun.z.df,by = "ASV")
sig.mulch.fun.z <- results.mulch.fun.z %>% filter(wi.eBH < 0.051)
write.table(sig.mulch.fun.z, file="data_tables/taxa_enrichments/genus.sig.mulch.fun.z.csv", sep=",")
write.table(results.mulch.fun.z, file="data_tables/taxa_enrichments/genus.all.mulch.fun.z.taxa.csv", sep=",")
##### Calculate effect size and create a volcano plot #####
## bac root ####
# ALDEx2 effect sizes
eff <- aldex.effect(mulch.bac.r.clr)
# Build plotting table
volc <- eff %>% 
  data.frame() %>% 
  rownames_to_column("Feature") %>%
  left_join(ttest.mulch.bac.r %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(wi.eBH),sig = wi.eBH < 0.05, direction = ifelse(rab.win.mulch > rab.win.no, "mulch", "no"),
         color_group = case_when(
           sig & direction == "mulch" ~ "mulch",
           sig & direction == "no"    ~ "no",
           TRUE                       ~ "ns"))
write.table(volc, file="data_tables/taxa_enrichments/genus.mulch.bac.r.eff.csv", sep=",")

ggplot(volc, aes(x = effect, y = neglog10_q)) +
  geom_point(aes(color = color_group), size = 3.5,alpha=0.9) +
  scale_color_manual(values = c(mulch = "#0065a2", no = "#00b0ba", ns = "gray30"),
    breaks = c("mulch", "no"),   # legend shows only these
    labels = c("Mulch-enriched", "Mulch–depleted")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color="gray50") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2,color="gray50") +
  ggrepel::geom_text_repel(
    data = subset(volc, sig),
    aes(label = Feature),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.25,
    point.padding = 1,
    min.segment.length = 0) +
  labs(x = "Effect size (CLR difference / dispersion)",
    y = expression(-log[10]("FDR q-value"),color = "Enriched in")) + theme_classic()+
  theme(legend.position="none",
      axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
      strip.text.x = element_text(size = 10),
      axis.title.y=element_text(size=10, family="sans",vjust = 1), 
      axis.text.x= element_text(colour="black", size=8, family="sans"), 
      axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "bac.r.mulch.volcano.labeled.pdf", device="pdf", units="mm", dpi=300, width=75, height=100, path="plots/")
## bac rhizo ####
# ALDEx2 effect sizes
eff <- aldex.effect(mulch.bac.z.clr)
# Build plotting table
volc <- eff %>% 
  data.frame() %>% 
  rownames_to_column("Feature") %>%
  left_join(ttest.mulch.bac.z %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(wi.eBH),sig = wi.eBH < 0.05, direction = ifelse(rab.win.mulch > rab.win.no, "mulch", "no"),
         color_group = case_when(
           sig & direction == "mulch" ~ "mulch",
           sig & direction == "no"    ~ "no",
           TRUE                       ~ "ns"))
write.table(volc, file="data_tables/taxa_enrichments/genus.mulch.bac.z.eff.csv", sep=",")

ggplot(volc, aes(x = effect, y = neglog10_q)) +
  geom_point(aes(color = color_group), size = 3.5,alpha=0.9) +
  scale_color_manual(values = c(mulch = "#0065a2", no = "#00b0ba", ns = "gray30"),
                     breaks = c("mulch", "no"),   # legend shows only these
                     labels = c("Mulch-enriched", "Mulch–depleted")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color="gray50") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2,color="gray50") +
  ggrepel::geom_text_repel(data = subset(volc, sig),aes(label = Feature),size = 3,max.overlaps = 30,box.padding = 0.25,point.padding = 1,min.segment.length = 0) +
  labs(x = "Effect size (CLR difference / dispersion)",
       y = expression(-log[10]("FDR q-value"),color = "Enriched in")) + theme_classic()+
  theme(legend.position="none",
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        axis.text.x= element_text(colour="black", size=8, family="sans"), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "bac.z.mulch.volcano.labeled.pdf", device="pdf", units="mm", dpi=300, width=175, height=200, path="plots/")
## fun root ####
# ALDEx2 effect sizes
eff <- aldex.effect(mulch.fun.r.clr)
# Build plotting table
volc <- eff %>% 
  data.frame() %>% 
  rownames_to_column("Feature") %>%
  left_join(ttest.mulch.fun.r %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(wi.eBH),sig = wi.eBH < 0.05, direction = ifelse(rab.win.mulch > rab.win.no, "mulch", "no"),
         color_group = case_when(
           sig & direction == "mulch" ~ "mulch",
           sig & direction == "no"    ~ "no",
           TRUE                       ~ "ns"))
write.table(volc, file="data_tables/taxa_enrichments/genus.mulch.fun.r.eff.csv", sep=",")

ggplot(volc, aes(x = effect, y = neglog10_q)) +
  geom_point(aes(color = color_group), size = 3.5,alpha=0.9) +
  scale_color_manual(values = c(mulch = "#0065a2", no = "#00b0ba", ns = "gray30"),
                     breaks = c("mulch", "no"),   # legend shows only these
                     labels = c("Mulch-enriched", "Mulch–depleted")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color="gray50") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2,color="gray50") +
  ggrepel::geom_text_repel(
    data = subset(volc, sig),
    aes(label = Feature),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.25,
    point.padding = 1,
    min.segment.length = 0) +
  labs(x = "Effect size (CLR difference / dispersion)",
       y = expression(-log[10]("FDR q-value"),color = "Enriched in")) + theme_classic()+
  theme(legend.position="none",
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        axis.text.x= element_text(colour="black", size=8, family="sans"), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "fun.r.mulch.volcano.labeled.pdf", device="pdf", units="mm", dpi=300, width=75, height=100, path="plots/")
## fun rhizo ####
# ALDEx2 effect sizes
eff <- aldex.effect(mulch.fun.z.clr)
# Build plotting table
volc <- eff %>% 
  data.frame() %>% 
  rownames_to_column("Feature") %>%
  left_join(ttest.mulch.fun.z %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(wi.eBH),sig = wi.eBH < 0.05, direction = ifelse(rab.win.mulch > rab.win.no, "mulch", "no"),
         color_group = case_when(
           sig & direction == "mulch" ~ "mulch",
           sig & direction == "no"    ~ "no",
           TRUE                       ~ "ns"))
write.table(volc, file="data_tables/taxa_enrichments/genus.mulch.fun.z.eff.csv", sep=",")

ggplot(volc, aes(x = effect, y = neglog10_q)) +
  geom_point(aes(color = color_group), size = 3.5,alpha=0.9) +
  scale_color_manual(values = c(mulch = "#0065a2", no = "#00b0ba", ns = "gray30"),
                     breaks = c("mulch", "no"),   # legend shows only these
                     labels = c("Mulch-enriched", "Mulch–depleted")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color="gray50") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2,color="gray50") +
  #ggrepel::geom_text_repel(data = subset(volc, sig),aes(label = Feature),size = 3,max.overlaps = 30,box.padding = 0.25,point.padding = 1,min.segment.length = 0) +
  labs(x = "Effect size (CLR difference / dispersion)",
       y = expression(-log[10]("FDR q-value"),color = "Enriched in")) + theme_classic()+
  theme(legend.position="none",
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        axis.text.x= element_text(colour="black", size=8, family="sans"), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "fun.z.mulch.volcano.pdf", device="pdf", units="mm", dpi=300, width=75, height=100, path="plots/")

##### Take a look at the read counts for Pleurostoma ASVs ####
fun.z.ps.raw.treat.mulch<-merge_samples(fun.z.ps.raw.treat, "Mulch")
fun.z.raw.treat.mulch.OTU<-t(data.frame(otu_table(fun.z.ps.raw.treat.mulch)))
Pleurostoma.ASVs <- c("fASV_711","fASV_762","fASV_1326","fASV_6546","fASV_6594")
fun.z.raw.treat.mulch.OTU.Pleuristoma <- fun.z.raw.treat.mulch.OTU[rownames(fun.z.raw.treat.mulch.OTU) %in% Pleurostoma.ASVs, ]
