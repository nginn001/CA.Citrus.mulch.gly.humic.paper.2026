# Script created December 3, 2025 by Nichole Ginnan, Assistant Project Scientist (nginn001@ucr.edu)
# Last updated July 13, 2026 (revised post-peer reviewed)
# Amplicon sequencing microbiome data analysis
# ECDRE Project, Lead PI: Caroline Roper, funded by the USDA
# Lindcove Research and Extension Center, Exeter, California | 91C plot experiment
# Field experiment testing the impacts of mulch, glyphosate, and humic acid citrus tree microbiome
# 16S bacterial and ITS fungal communities
###############################################################################
# Load Libraries ####
# Script was written and executed using R v 4.5.0
library(tidyverse) # v 2.0.0
library(writexl) # v 1.5.4
library(readxl) # v 1.4.5
library(phyloseq) # v 1.53.0
library(lme4) # v 1.1-37
library(lmerTest) # v 3.1-3
library(scales) # v 1.4.0
library(emmeans) # v 1.11.2
library(multcompView) # v 0.1-10
library(multcomp) # v 1.4-28
library(vegan) # v 2.7-1
library(ALDEx2) # v 1.41.0
library(MuMIn) #v 1.48.1

set.seed(4444)
###############################################################################
# Load data objects; centered log-ratio transformed phyloseq objects ####
setwd("/Users/nicholeginnan/Documents/UCR- Current/CA.citrus.paper") # set working directory
bac.clr.ps<-readRDS("Revisions.ISME.Comm/phyloseq_objects/CA.CLR.16S.greengenes.phyloseq.obj.July26.RDS")
fun.clr.ps<-readRDS("Revisions.ISME.Comm/phyloseq_objects/CA.CLR.ITS.UNITE.phyloseq.obj.July26.RDS")
# Load data objects; raw counts phyloseq objects ####
bac.raw.ps<-readRDS("Revisions.ISME.Comm/phyloseq_objects/CA.raw.16S.greengenes.phyloseq.obj.July26.RDS")
fun.raw.ps<-readRDS("Revisions.ISME.Comm/phyloseq_objects/CA.raw.ITS.UNITE.phyloseq.obj.July26.RDS")
###############################################################################
# subset data CLR #####
# tissue type #
bac.r.ps<- subset_samples(bac.clr.ps, Tissue_type == "Roots")
bac.z.ps<- subset_samples(bac.clr.ps, Tissue_type == "Rhizosphere")
fun.r.ps<- subset_samples(fun.clr.ps, Tissue_type == "Roots")
fun.z.ps<- subset_samples(fun.clr.ps, Tissue_type == "Rhizosphere")

# fix rep IDS that have re-sequenced "R" indicators
fix_treat_rep <- function(ps) {
  sd <- data.frame(sample_data(ps))
  sd$Treat_rep_clean <- gsub("R$", "", as.character(sd$Treat_rep))
  sample_data(ps) <- sample_data(sd)
  ps}
bac.r.ps <- fix_treat_rep(bac.r.ps)
bac.z.ps <- fix_treat_rep(bac.z.ps)
fun.r.ps <- fix_treat_rep(fun.r.ps)
fun.z.ps <- fix_treat_rep(fun.z.ps)

# Subset year zero (pre-treatment year)
bac.r.ps.T0 <- subset_samples(bac.r.ps, Timepoint == "T0") 
bac.z.ps.T0 <- subset_samples(bac.z.ps, Timepoint == "T0") 
fun.r.ps.T0 <- subset_samples(fun.r.ps, Timepoint == "T0") 
fun.z.ps.T0 <- subset_samples(fun.z.ps, Timepoint == "T0") 
# keep all years
bac.r.ps.T0_T3 <- subset_samples(bac.r.ps, Timepoint == "T0" | Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
bac.z.ps.T0_T3 <- subset_samples(bac.z.ps, Timepoint == "T0" | Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.r.ps.T0_T3 <- subset_samples(fun.r.ps, Timepoint == "T0" | Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.z.ps.T0_T3 <- subset_samples(fun.z.ps, Timepoint == "T0" | Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3")
# remove year zero (pre-treatment year and only use years 1-3) #
bac.r.ps.treat <- subset_samples(bac.r.ps, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
bac.z.ps.treat <- subset_samples(bac.z.ps, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.r.ps.treat <- subset_samples(fun.r.ps, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.z.ps.treat <- subset_samples(fun.z.ps, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 

# subset data - raw counts #####
# tissue type #
bac.r.ps.raw<- subset_samples(bac.raw.ps, Tissue_type == "Roots")
bac.z.ps.raw<- subset_samples(bac.raw.ps, Tissue_type == "Rhizosphere")
fun.r.ps.raw<- subset_samples(fun.raw.ps, Tissue_type == "Roots")
fun.z.ps.raw<- subset_samples(fun.raw.ps, Tissue_type == "Rhizosphere")

# fix rep IDS that have re-sequenced "R" indicators
fix_treat_rep <- function(ps) {
  sd <- data.frame(sample_data(ps))
  sd$Treat_rep_clean <- gsub("R$", "", as.character(sd$Treat_rep))
  sample_data(ps) <- sample_data(sd)
  ps}
bac.r.ps.raw <- fix_treat_rep(bac.r.ps.raw)
bac.z.ps.raw <- fix_treat_rep(bac.z.ps.raw)
fun.r.ps.raw <- fix_treat_rep(fun.r.ps.raw)
fun.z.ps.raw <- fix_treat_rep(fun.z.ps.raw)
# T0-T3 (pre-treatment and treatment years) #
bac.r.ps.raw.T0_T3 <- subset_samples(bac.r.ps.raw, Timepoint == "T0" | Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
bac.z.ps.raw.T0_T3 <- subset_samples(bac.z.ps.raw, Timepoint == "T0" | Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.r.ps.raw.T0_T3 <- subset_samples(fun.r.ps.raw, Timepoint == "T0" | Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.z.ps.raw.T0_T3 <- subset_samples(fun.z.ps.raw, Timepoint == "T0" | Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3")
# remove year zero (pre-treatment) #
bac.r.ps.raw.treat <- subset_samples(bac.r.ps.raw, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
bac.z.ps.raw.treat <- subset_samples(bac.z.ps.raw, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.r.ps.raw.treat <- subset_samples(fun.r.ps.raw, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
fun.z.ps.raw.treat <- subset_samples(fun.z.ps.raw, Timepoint == "T1" | Timepoint == "T2" | Timepoint == "T3") 
###############################################################################
###### Alpha Diversity - Change from Timepoint Zero (pre-treatment) ###########
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
##
## Bac root - Delta Shannon ####
md.bac.r<-data.frame(sample_data(bac.r.ps.T0_T3)) # extract metadata from the phyloseq object
md.bac.r$z_logObs<-scale(md.bac.r$logObs) # create new col with logObs z-score to control for sequencing depth variability
# calculate Delta Shannon for each tree at each
T0 <- md.bac.r %>%filter(Timepoint == "T0") %>%dplyr::select(Treat_rep_clean, Shannon) %>%rename(Shannon_T0 = Shannon)
md.bac.r.delta <- md.bac.r %>%left_join(T0, by = "Treat_rep_clean")
md.bac.r.delta <- md.bac.r.delta %>%mutate(Delta_Shannon = Shannon - Shannon_T0)
md.bac.r.delta <- md.bac.r.delta %>%filter(Timepoint != "T0")
# Mixed Effects Model
m0<-lmer(Delta_Shannon~Mulch*Glyphosate*Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta)) # model 1 
# No 3-way interaction
m1 <- lmer(Delta_Shannon ~Mulch + Glyphosate + Humic + Mulch*Glyphosate +Mulch*Humic+Glyphosate*Humic+z_logObs+(1|Timepoint)+(1|Treat_rep_clean),data = md.bac.r.delta)
# mulch interactions only
m2<-lmer(Delta_Shannon~Mulch * Glyphosate + Mulch * Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta)) # model 1 
# Mulch:humic interaction only
m3<-lmer(Delta_Shannon~ Mulch * Humic +Glyphosate + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta)) # model 1 
# Gly:humic interaction only
m4<-lmer(Delta_Shannon~Mulch + Glyphosate * Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta)) # model 1 
# Additive
m5<-lmer(Delta_Shannon~Mulch + Glyphosate + Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta)) # model 1 
# Check fit
plot(resid(m0)~fitted(m0)) # looks good
plot(resid(m1)~fitted(m1)) # looks good
plot(resid(m2)~fitted(m2)) # looks good
plot(resid(m3)~fitted(m3)) # looks good
plot(resid(m4)~fitted(m4)) # looks good
plot(resid(m5)~fitted(m5)) # looks good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3, m4,m5)
#.    (Int) Gly Hmc Mlc  z_lgO Gly:Hmc Gly:Mlc Hmc:Mlc Gly:Hmc:Mlc df  logLik  AICc delta weight
#m3 0.7114   +   +   + 0.4031                       +              9 -75.956 172.2  0.00  0.419 <--
#m5 0.4866   +   +   + 0.4025                                      8 -78.063 173.9  1.74  0.176
#m4 0.6455   +   +   + 0.4002       +                              9 -76.980 174.2  2.05  0.150
#m2 0.7107   +   +   + 0.4041               +       +             10 -75.920 174.7  2.47  0.122
#m1 0.8393   +   +   + 0.4019       +       +       +             11 -75.093 175.6  3.42  0.076
#m0 0.9329   +   +   + 0.4017       +       +       +           + 12 -74.032 176.2  3.98  0.057
AIC(m0,m1,m2,m3, m4,m5)
#.   df      AIC
# m0 12 172.0637
# m1 11 172.1863
# m2 10 171.8402
# m3  9 169.9116 <--
# m4  9 171.9604
# m5  8 172.1263
# Test
anova(m3, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Mulch        0.4306  0.4306     1 23.961  2.7156   0.11243    
#Humic        0.5969  0.5969     1 31.450  3.7644   0.06136 .  
#Glyphosate   0.2453  0.2453     1 24.684  1.5469   0.22528    
#z_logObs    10.2379 10.2379     1 49.970 64.5619 1.474e-10 ***
#Mulch:Humic  0.6967  0.6967     1 30.951  4.3936   0.04435 *  

#emm<-emmeans(m3, specs=~Humic:Mulch|Glyphosate, type="response")
emm<-emmeans(m3, specs=~Humic:Mulch:Glyphosate, type="response")

cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>% # plotting emmeans #
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic, group = interaction(Glyphosate, Humic))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.75), color = "black") +
  geom_point(data=md.bac.r.delta,mapping=aes(x=Mulch, y= Delta_Shannon, fill = Glyphosate, color = Glyphosate,shape=Humic), alpha=0.3, position=position_jitterdodge(0.001), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.75)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 3, vjust = -4, hjust = 0.5,
            position = position_dodge(width = 0.75)) +
  ylab("Delta Shannon Index") +
  xlab("Treatment") +
  scale_x_discrete(labels=c("mulch"="mulch", "no"="no mulch"))+
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "bac.delta.shannon.rt.mulch.humic.glyphosate.zlogObs.points.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

## Bac rhizos - Delta Shannon ####
md.bac.z<-data.frame(sample_data(bac.z.ps.T0_T3)) # extract metadata from the phyloseq object
md.bac.z$z_logObs<-scale(md.bac.z$logObs) # create new col with logObs z-score to control for sequencing depth variability
# calculate Delta Shannon for each tree at each
T0 <- md.bac.z %>%filter(Timepoint == "T0") %>%dplyr::select(Treat_rep_clean, Shannon) %>%rename(Shannon_T0 = Shannon)
md.bac.z.delta <- md.bac.z %>%left_join(T0, by = "Treat_rep_clean")
md.bac.z.delta <- md.bac.z.delta %>%mutate(Delta_Shannon = Shannon - Shannon_T0)
md.bac.z.delta <- md.bac.z.delta %>%filter(Timepoint != "T0")
# Mixed Effects Model
m0<-lmer(Delta_Shannon~Mulch*Glyphosate*Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# No 3-way interaction
m1 <- lmer(Delta_Shannon ~Mulch + Glyphosate + Humic + Mulch*Glyphosate +Mulch*Humic+Glyphosate*Humic+z_logObs+(1|Timepoint)+(1|Treat_rep_clean),data = md.bac.z.delta)
# mulch interactions only
m2<-lmer(Delta_Shannon~Mulch * Glyphosate + Mulch * Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# Mulch:humic interaction only
m3<-lmer(Delta_Shannon~ Mulch * Humic +Glyphosate + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# Gly:humic interaction only
m4<-lmer(Delta_Shannon~Mulch + Glyphosate * Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# Additive
m5<-lmer(Delta_Shannon~Mulch + Glyphosate + Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# check fit
plot(resid(m0)~fitted(m0)) # looks good
plot(resid(m1)~fitted(m1)) # looks good
plot(resid(m2)~fitted(m2)) # looks good
plot(resid(m3)~fitted(m3)) # looks good
plot(resid(m4)~fitted(m4)) # looks good
plot(resid(m5)~fitted(m5)) # looks good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3, m4,m5)
#.    (Int) Gly Hmc Mlc  z_lgO Gly:Hmc Gly:Mlc Hmc:Mlc Gly:Hmc:Mlc df  logLik  AICc delta weight
#m3 0.35420   +   +   + 0.7117                       +              9 -65.280 150.8  0.00  0.560 <--
#m2 0.27820   +   +   + 0.7119               +       +             10 -65.001 152.8  1.97  0.209
#m0 0.22880   +   +   + 0.7197       +       +       +           + 12 -62.886 153.8  2.99  0.126
#m1 0.36670   +   +   + 0.7156       +       +       +             11 -64.611 154.6  3.78  0.084
#m5 0.03516   +   +   + 0.7229                                      8 -70.168 158.1  7.31  0.014
#m4 0.12680   +   +   + 0.7256       +                              9 -69.722 159.7  8.88  0.007
AIC(m0,m1,m2,m3, m4,m5)
#.   df      AIC
# 12 149.7712
#m1 11 151.2224
#m2 10 150.0020
#m3  9 148.5597 <--
#m4  9 157.4440
#m5  8 156.3367
# Test
anova(m3, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq  NumDF  DenDF F value   Pr(>F)    
#Mulch        0.0004  0.0004     1 26.915   0.0034 0.95368    
#Humic        0.0229  0.0229     1 26.993   0.2223 0.64108    
#Glyphosate   0.1151  0.1151     1 27.293   1.1160 0.30004    
#z_logObs    14.4094 14.4094     1 65.874 139.6865 < 2e-16 ***
#Mulch:Humic  1.1826  1.1826     1 26.968  11.4647 0.00219 **  

emm<-emmeans(m3, specs=~Humic:Mulch:Glyphosate, type="response")
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>% # plotting emmeans #
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic, group = interaction(Glyphosate, Humic))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.75), color = "black") +
  geom_point(data=md.bac.z.delta,mapping=aes(x=Mulch, y= Delta_Shannon, fill = Glyphosate, color = Glyphosate,shape=Humic), alpha=0.3, position=position_jitterdodge(0.001), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.75)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 3, vjust = -4, hjust = 0.5,
            position = position_dodge(width = 0.75)) +
  ylab("Delta Shannon Index") +
  xlab("Treatment") +
  scale_x_discrete(labels=c("mulch"="mulch", "no"="no mulch"))+
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "bac.delta.shannon.z.mulch.humic.glyphosate.zlogObs.points.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

## Fun root - Delta Shannon ####
md.fun.r<-data.frame(sample_data(fun.r.ps.T0_T3)) # extract metadata from the phyloseq object
md.fun.r$z_logObs<-scale(md.fun.r$logObs) # create new col with logObs z-score to control for sequencing depth variability
# calculate Delta Shannon for each tree at each
T0 <- md.fun.r %>%filter(Timepoint == "T0") %>%dplyr::select(Treat_rep_clean, Shannon) %>%rename(Shannon_T0 = Shannon)
md.fun.r.delta <- md.fun.r %>%left_join(T0, by = "Treat_rep_clean")
md.fun.r.delta <- md.fun.r.delta %>%mutate(Delta_Shannon = Shannon - Shannon_T0)
md.fun.r.delta <- md.fun.r.delta %>%filter(Timepoint != "T0")
# Mixed Effects Model
m0<-lmer(Delta_Shannon~Mulch*Glyphosate*Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta)) # model 1 
# No 3-way interaction
m1 <- lmer(Delta_Shannon ~Mulch + Glyphosate + Humic_acid + Mulch*Glyphosate +Mulch*Humic_acid+Glyphosate*Humic_acid+z_logObs+(1|Timepoint)+(1|Treat_rep_clean),data = md.fun.r.delta)
# mulch interactions only
m2<-lmer(Delta_Shannon~Mulch * Glyphosate + Mulch * Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta)) # model 1 
# Mulch:Humic_acid interaction only
m3<-lmer(Delta_Shannon~ Mulch * Humic_acid +Glyphosate + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta)) # model 1 
# Gly:Humic_acid interaction only
m4<-lmer(Delta_Shannon~Mulch + Glyphosate * Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta)) # model 1 
# Additive
m5<-lmer(Delta_Shannon~Mulch + Glyphosate + Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta)) # model 1 
# Check fit
plot(resid(m0)~fitted(m0)) # looks good
plot(resid(m1)~fitted(m1)) # looks good
plot(resid(m2)~fitted(m2)) # looks good
plot(resid(m3)~fitted(m3)) # looks good
plot(resid(m4)~fitted(m4)) # looks good
plot(resid(m5)~fitted(m5)) # looks good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3, m4,m5)
#.    (Int) Gly Hmc_acd Mlc  z_lgO Gly:Hmc_acd Gly:Mlc Hmc_acd:Mlc Gly:Hmc_acd:Mlc df   logLik  AICc delta weight
#m0 -0.40440   +       +   + 0.4647           +       +           +               + 12 -125.689 279.3  0.00  0.304 <--
#m5  0.13660   +       +   + 0.4628                                                  8 -130.970 279.7  0.35  0.255
#m4 -0.04779   +       +   + 0.4548           +                                      9 -129.940 280.1  0.75  0.209
#m3  0.15900   +       +   + 0.4643                               +                  9 -130.485 281.2  1.84  0.121
#m2  0.09982   +       +   + 0.4607                   +           +                 10 -129.886 282.5  3.16  0.063
#m1 -0.06817   +       +   + 0.4540           +       +           +                 11 -128.854 283.0  3.68  0.048
AIC(m0,m1,m2,m3, m4,m5)
#.   df      AIC
#m0 12 275.3779 <--
#m1 11 279.7085
#m2 10 279.7730
#m3  9 278.9695
#m4  9 277.8798
#m5  8 277.9409
# Test
anova(m0, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Mulch                        0.0846  0.0846     1 23.984  0.1481   0.70377    
#Glyphosate                   0.4621  0.4621     1 23.998  0.8087   0.37744    
#Humic_acid                   0.0336  0.0336     1 28.226  0.0589   0.81006    
#z_logObs                    12.5134 12.5134     1 68.727 21.8993 1.395e-05 ***
#Mulch:Glyphosate             0.0730  0.0730     1 24.623  0.1277   0.72385    
#Mulch:Humic_acid             0.0044  0.0044     1 27.894  0.0078   0.93038    
#Glyphosate:Humic_acid        0.8832  0.8832     1 27.638  1.5456   0.22423    
#Mulch:Glyphosate:Humic_acid  2.4569  2.4569     1 27.679  4.2998   0.04754 *  

emm<-emmeans(m0, specs=~Humic_acid:Mulch:Glyphosate, type="response")
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>% # plotting emmeans #
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic_acid, group = interaction(Glyphosate, Humic_acid))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.75), color = "black") +
  geom_point(data=md.fun.r.delta,mapping=aes(x=Mulch, y= Delta_Shannon, fill = Glyphosate, color = Glyphosate,shape=Humic_acid), alpha=0.3, position=position_jitterdodge(0.001), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.75)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 3, vjust = -4, hjust = 0.5,
            position = position_dodge(width = 0.75)) +
  ylab("Delta Shannon Index") +
  xlab("Treatment") +
  scale_x_discrete(labels=c("mulch"="mulch", "no"="no mulch"))+
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "fun.delta.shannon.rt.mulch.Humic_acid.glyphosate.zlogObs.points.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

## Fun rhizos - Delta Shannon ####
md.fun.z<-data.frame(sample_data(fun.z.ps.T0_T3)) # extract metadata from the phyloseq object
md.fun.z$z_logObs<-scale(md.fun.z$logObs) # create new col with logObs z-score to control for sequencing depth variability
# calculate Delta Shannon for each tree at each
T0 <- md.fun.z %>%filter(Timepoint == "T0") %>%dplyr::select(Treat_rep_clean, Shannon) %>%rename(Shannon_T0 = Shannon)
md.fun.z.delta <- md.fun.z %>%left_join(T0, by = "Treat_rep_clean")
md.fun.z.delta <- md.fun.z.delta %>%mutate(Delta_Shannon = Shannon - Shannon_T0)
md.fun.z.delta <- md.fun.z.delta %>%filter(Timepoint != "T0")
# Mixed Effects Model
m0<-lmer(Delta_Shannon~Mulch*Glyphosate*Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta)) # model 1 
# No 3-way interaction
m1 <- lmer(Delta_Shannon ~Mulch + Glyphosate + Humic_acid + Mulch*Glyphosate +Mulch*Humic_acid+Glyphosate*Humic_acid+z_logObs+(1|Timepoint)+(1|Treat_rep_clean),data = md.fun.z.delta)
# mulch interactions only
m2<-lmer(Delta_Shannon~Mulch * Glyphosate + Mulch * Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta)) # model 1 
# Mulch:Humic_acid interaction only
m3<-lmer(Delta_Shannon~ Mulch * Humic_acid +Glyphosate + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta)) # model 1 
# Gly:Humic_acid interaction only
m4<-lmer(Delta_Shannon~Mulch + Glyphosate * Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta)) # model 1 
# Additive
m5<-lmer(Delta_Shannon~Mulch + Glyphosate + Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta)) # model 1 
# check fit
plot(resid(m0)~fitted(m0)) # looks good
plot(resid(m1)~fitted(m1)) # looks good
plot(resid(m2)~fitted(m2)) # looks good
plot(resid(m3)~fitted(m3)) # looks good
plot(resid(m4)~fitted(m4)) # looks good
plot(resid(m5)~fitted(m5)) # looks good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3, m4,m5)
#.    (Int) Gly Hmc_acd Mlc  z_lgO Gly:Hmc_acd Gly:Mlc Hmc_acd:Mlc Gly:Hmc_acd:Mlc df   logLik  AICc delta weight
#m0 -0.53890   +       +   + 0.5977           +       +           +               + 12 -136.153 300.5  0.00  0.307 <--
#m1 -0.38910   +       +   + 0.6009           +       +           +                 11 -137.836 301.1  0.68  0.219
#m2 -0.14030   +       +   + 0.5871                   +           +                 10 -139.303 301.5  1.00  0.187
#m4  0.08401   +       +   + 0.6153           +                                      9 -141.068 302.4  1.98  0.114
#m5  0.30180   +       +   + 0.6027                                                  8 -142.424 302.7  2.20  0.102
#m3  0.19370   +       +   + 0.6039                               +                  9 -141.551 303.4  2.94  0.071
AIC(m0,m1,m2,m3, m4,m5)
#.   df      AIC
#m0 12 296.3062 <--
#m1 11 297.6714
#m2 10 298.6050
#m3  9 301.1019
#m4  9 300.1354
#m5  8 300.8471
# Test
anova(m0, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#                            Sum Sq  Mean Sq  NumDF  DenDF F value   Pr(>F)    
#Mulch                        0.0980  0.0980     1 23.990  0.1086   0.74459    
#Glyphosate                   1.4685  1.4685     1 22.820  1.6271   0.21494    
#Humic_acid                   0.2912  0.2912     1 22.831  0.3226   0.57558    
#z_logObs                    26.7347 26.7347     1 64.243 29.6220 8.787e-07 ***
#Mulch:Glyphosate             3.0711  3.0711     1 22.925  3.4028   0.07805 .  
#Mulch:Humic_acid             0.3860  0.3860     1 22.736  0.4277   0.51969    
#Glyphosate:Humic_acid        1.3529  1.3529     1 23.184  1.4990   0.23313    
#Mulch:Glyphosate:Humic_acid  0.5707  0.5707     1 22.806  0.6324   0.43469    

emm<-emmeans(m0, specs=~Humic_acid:Mulch:Glyphosate, type="response")
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>% # plotting emmeans #
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic_acid, group = interaction(Glyphosate, Humic_acid))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.75), color = "black") +
  geom_point(data=md.fun.z.delta,mapping=aes(x=Mulch, y= Delta_Shannon, fill = Glyphosate, color = Glyphosate,shape=Humic_acid), alpha=0.3, position=position_jitterdodge(0.001), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.75)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 3, vjust = -4, hjust = 0.5,
            position = position_dodge(width = 0.75)) +
  ylab("Delta Shannon Index") +
  xlab("Treatment") +
  scale_x_discrete(labels=c("mulch"="mulch", "no"="no mulch"))+
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "fun.delta.shannon.z.mulch.humic.glyphosate.zlogObs.points.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")


## Bac root - Delta InvSimpson ####
md.bac.r<-data.frame(sample_data(bac.r.ps.T0_T3)) # extract metadata from the phyloseq object
md.bac.r$z_logObs<-scale(md.bac.r$logObs) # create new col with logObs z-score to control for sequencing depth variability
# calculate Delta InvSimpson for each tree at each
T0 <- md.bac.r %>%filter(Timepoint == "T0") %>%dplyr::select(Treat_rep_clean, InvSimpson) %>%rename(InvSimpson_T0 = InvSimpson)
md.bac.r.delta <- md.bac.r %>%left_join(T0, by = "Treat_rep_clean")
md.bac.r.delta <- md.bac.r.delta %>%mutate(Delta_InvSimpson = InvSimpson - InvSimpson_T0)
md.bac.r.delta <- md.bac.r.delta %>%filter(Timepoint != "T0")
# check for outliers
ggplot(md.bac.r.delta.2, aes(y=Delta_InvSimpson, x=Mulch, color=Glyphosate, shape=Humic, label = rownames(md.bac.r.delta.2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("30","83")) 
md.bac.r.delta.2<-md.bac.r.delta %>% filter(!rownames(md.bac.r.delta) %in% outliers)
outliers<-as.character(c("79")) 
md.bac.r.delta.3<-md.bac.r.delta.2 %>% filter(!rownames(md.bac.r.delta.2) %in% outliers)
# Mixed Effects Model
m0<-lmer(Delta_InvSimpson~Mulch*Glyphosate*Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta.3)) # model 1 
# No 3-way interaction
m1 <- lmer(Delta_InvSimpson ~Mulch + Glyphosate + Humic + Mulch*Glyphosate +Mulch*Humic+Glyphosate*Humic+z_logObs+(1|Timepoint)+(1|Treat_rep_clean),data = md.bac.r.delta.3)
# mulch interactions only
m2<-lmer(Delta_InvSimpson~Mulch * Glyphosate + Mulch * Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta.3)) # model 1 
# Mulch:humic interaction only
m3<-lmer(Delta_InvSimpson~ Mulch * Humic +Glyphosate + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta.3)) # model 1 
# Gly:humic interaction only
m4<-lmer(Delta_InvSimpson~Mulch + Glyphosate * Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta.3)) # model 1 
# Additive
m5<-lmer(Delta_InvSimpson~Mulch + Glyphosate + Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.r.delta.3)) # model 1 
# Check fit
plot(resid(m0)~fitted(m0)) # looks ok
text(fitted(m0),resid(m0),labels = rownames(md.bac.r.delta.2),pos = 4,cex = 0.6)
qqnorm(resid(m0));qqline(resid(m0))
plot(resid(m1)~fitted(m1)) # looks ok
plot(resid(m2)~fitted(m2)) # looks ok
plot(resid(m3)~fitted(m3)) # looks ok
plot(resid(m4)~fitted(m4)) # looks ok
plot(resid(m5)~fitted(m5)) # looks ok
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3, m4,m5)
#.    (Int) Gly Hmc Mlc  z_lgO Gly:Hmc Gly:Mlc Hmc:Mlc Gly:Hmc:Mlc df  logLik  AICc delta weight
#m0 20.37   +   +   + 5.795       +       +       +           + 12 -358.263 744.8  0.00  0.987 <--
#m1 13.79   +   +   + 6.002       +       +       +             11 -364.083 753.7  8.93  0.011
#m2 15.15   +   +   + 5.881               +       +             10 -367.708 758.3 13.55  0.001
#m3 17.93   +   +   + 5.729                       +              9 -371.579 763.5 18.73  0.000
#m4 10.90   +   +   + 5.795       +                              9 -372.717 765.8 21.00  0.000
#m5 12.00   +   +   + 5.697                                      8 -376.308 770.5 25.69  0.000
AIC(m0,m1,m2,m3, m4,m5)
#.   df      AIC
#m0 12 740.5254 <--
#m1 11 750.1658
#m2 10 755.4156
#m3  9 761.1573
#m4  9 763.4348
#m5  8 768.6151
# Test
anova(m0, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Mulch                   493.75  493.75     1 23.745  1.4427 0.24154  
#Glyphosate              265.61  265.61     1 25.066  0.7761 0.38670  
#Humic                  1307.41 1307.41     1 27.570  3.8202 0.06085 .
#z_logObs               2300.83 2300.83     1 27.801  6.7230 0.01501 *
#Mulch:Glyphosate        213.60  213.60     1 23.550  0.6241 0.43740  
#Mulch:Humic            1034.04 1034.04     1 26.581  3.0214 0.09374 .
#Glyphosate:Humic         33.51   33.51     1 26.826  0.0979 0.75678  
#Mulch:Glyphosate:Humic 1172.55 1172.55     1 26.611  3.4262 0.07530 . 

emm<-emmeans(m0, specs=~Humic:Mulch:Glyphosate, type="response")

cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>% # plotting emmeans #
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic, group = interaction(Glyphosate, Humic))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.75), color = "black") +
  geom_point(data=md.bac.r.delta.3,mapping=aes(x=Mulch, y= Delta_InvSimpson, fill = Glyphosate, color = Glyphosate,shape=Humic), alpha=0.3, position=position_jitterdodge(0.001), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.75)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 3, vjust = -4, hjust = 0.5,
            position = position_dodge(width = 0.75)) +
  ylab("Delta InvSimpson Index") +
  xlab("Treatment") +
  scale_x_discrete(labels=c("mulch"="mulch", "no"="no mulch"))+
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "bac.delta.InvSimpson.rt.mulch.humic.glyphosate.zlogObs.points.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

## Bac rhizos - Delta InvSimpson ####
md.bac.z<-data.frame(sample_data(bac.z.ps.T0_T3)) # extract metadata from the phyloseq object
md.bac.z$z_logObs<-scale(md.bac.z$logObs) # create new col with logObs z-score to control for sequencing depth variability
# calculate Delta InvSimpson for each tree at each
T0 <- md.bac.z %>%filter(Timepoint == "T0") %>%dplyr::select(Treat_rep_clean, InvSimpson) %>%rename(InvSimpson_T0 = InvSimpson)
md.bac.z.delta <- md.bac.z %>%left_join(T0, by = "Treat_rep_clean")
md.bac.z.delta <- md.bac.z.delta %>%mutate(Delta_InvSimpson = InvSimpson - InvSimpson_T0)
md.bac.z.delta <- md.bac.z.delta %>%filter(Timepoint != "T0")
# Mixed Effects Model
m0<-lmer(Delta_InvSimpson~Mulch*Glyphosate*Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# No 3-way interaction
m1 <- lmer(Delta_InvSimpson ~Mulch + Glyphosate + Humic + Mulch*Glyphosate +Mulch*Humic+Glyphosate*Humic+z_logObs+(1|Timepoint)+(1|Treat_rep_clean),data = md.bac.z.delta)
# mulch interactions only
m2<-lmer(Delta_InvSimpson~Mulch * Glyphosate + Mulch * Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# Mulch:humic interaction only
m3<-lmer(Delta_InvSimpson~ Mulch * Humic +Glyphosate + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# Gly:humic interaction only
m4<-lmer(Delta_InvSimpson~Mulch + Glyphosate * Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# Additive
m5<-lmer(Delta_InvSimpson~Mulch + Glyphosate + Humic + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.bac.z.delta)) # model 1 
# check fit
plot(resid(m0)~fitted(m0)) # looks good
plot(resid(m1)~fitted(m1)) # looks good
plot(resid(m2)~fitted(m2)) # looks good
plot(resid(m3)~fitted(m3)) # looks good
plot(resid(m4)~fitted(m4)) # looks good
plot(resid(m5)~fitted(m5)) # looks good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3, m4,m5)
#.    (Int) Gly Hmc Mlc  z_lgO Gly:Hmc Gly:Mlc Hmc:Mlc Gly:Hmc:Mlc df  logLik  AICc delta weight
#m0 40.69   +   +   + 55.91       +       +       +           + 12 -482.879  993.8  0.00  0.997 <--
#m1 66.07   +   +   + 54.85       +       +       +             11 -490.048 1005.4 11.68  0.003
#m2 48.56   +   +   + 54.08               +       +             10 -495.738 1014.2 20.47  0.000
#m3 57.86   +   +   + 54.04                       +              9 -500.964 1022.1 28.39  0.000
#m4 50.80   +   +   + 55.40       +                              9 -501.635 1023.5 29.73  0.000
#m5 32.90   +   +   + 54.64                                      8 -507.331 1032.4 38.66  0.000
AIC(m0,m1,m2,m3, m4,m5)
#.   df      AIC
#m0 12  989.7577 <--
#m1 11 1002.0959
#m2 10 1011.4755
#m3  9 1019.9271
#m4  9 1021.2691
#m5  8 1030.6614
# Test
anova(m0, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq  NumDF  DenDF F value   Pr(>F)    
#Mulch                    9709    9709     1 23.647  3.1862    0.0871 .  
#Glyphosate               2072    2072     1 24.067  0.6799    0.4177    
#Humic                      11      11     1 23.729  0.0037    0.9520    
#z_logObs                89728   89728     1 66.147 29.4466 8.777e-07 ***
#Mulch:Glyphosate         1285    1285     1 23.672  0.4218    0.5223    
#Mulch:Humic              8670    8670     1 23.706  2.8452    0.1048    
#Glyphosate:Humic         4334    4334     1 23.826  1.4224    0.2447    
#Mulch:Glyphosate:Humic   9206    9206     1 23.802  3.0213    0.0951 .  

emm<-emmeans(m0, specs=~Humic:Mulch:Glyphosate, type="response")
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>% # plotting emmeans #
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic, group = interaction(Glyphosate, Humic))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.75), color = "black") +
  geom_point(data=md.bac.z.delta,mapping=aes(x=Mulch, y= Delta_InvSimpson, fill = Glyphosate, color = Glyphosate,shape=Humic), alpha=0.3, position=position_jitterdodge(0.001), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.75)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 3, vjust = -4, hjust = 0.5,
            position = position_dodge(width = 0.75)) +
  ylab("Delta InvSimpson Index") +
  xlab("Treatment") +
  scale_x_discrete(labels=c("mulch"="mulch", "no"="no mulch"))+
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "bac.delta.InvSimpson.z.mulch.humic.glyphosate.zlogObs.points.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

## Fun root - Delta InvSimpson ####
md.fun.r<-data.frame(sample_data(fun.r.ps.T0_T3)) # extract metadata from the phyloseq object
md.fun.r$z_logObs<-scale(md.fun.r$logObs) # create new col with logObs z-score to control for sequencing depth variability
# calculate Delta InvSimpson for each tree at each
T0 <- md.fun.r %>%filter(Timepoint == "T0") %>%dplyr::select(Treat_rep_clean, InvSimpson) %>%rename(InvSimpson_T0 = InvSimpson)
md.fun.r.delta <- md.fun.r %>%left_join(T0, by = "Treat_rep_clean")
md.fun.r.delta <- md.fun.r.delta %>%mutate(Delta_InvSimpson = InvSimpson - InvSimpson_T0)
md.fun.r.delta <- md.fun.r.delta %>%filter(Timepoint != "T0")
# check for outliers
ggplot(md.fun.r.delta.2, aes(y=Delta_InvSimpson, x=Mulch, color=Glyphosate, shape=Humic_acid, label = rownames(md.fun.r.delta.2))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("61","64")) 
md.fun.r.delta.2<-md.fun.r.delta %>% filter(!rownames(md.fun.r.delta) %in% outliers)
# Mixed Effects Model
m0<-lmer(Delta_InvSimpson~Mulch*Glyphosate*Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta.2)) # model 1 
# No 3-way interaction
m1 <- lmer(Delta_InvSimpson ~Mulch + Glyphosate + Humic_acid + Mulch*Glyphosate +Mulch*Humic_acid+Glyphosate*Humic_acid+z_logObs+(1|Timepoint)+(1|Treat_rep_clean),data = md.fun.r.delta.2)
# mulch interactions only
m2<-lmer(Delta_InvSimpson~Mulch * Glyphosate + Mulch * Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta.2)) # model 1 
# Mulch:Humic_acid interaction only
m3<-lmer(Delta_InvSimpson~ Mulch * Humic_acid +Glyphosate + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta.2)) # model 1 
# Gly:Humic_acid interaction only
m4<-lmer(Delta_InvSimpson~Mulch + Glyphosate * Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta.2)) # model 1 
# Additive
m5<-lmer(Delta_InvSimpson~Mulch + Glyphosate + Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.r.delta.2)) # model 1 
# Check fit
plot(resid(m0)~fitted(m0)) # looks good
plot(resid(m1)~fitted(m1)) # looks good
plot(resid(m2)~fitted(m2)) # looks good
plot(resid(m3)~fitted(m3)) # looks good
plot(resid(m4)~fitted(m4)) # looks good
plot(resid(m5)~fitted(m5)) # looks good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3, m4,m5)
#.    (Int) Gly Hmc_acd Mlc  z_lgO Gly:Hmc_acd Gly:Mlc Hmc_acd:Mlc Gly:Hmc_acd:Mlc df   logLik  AICc delta weight
#m0 -1.1770   +       +   + 0.9757           +       +           +               + 12 -216.302 460.7  0.00  0.852 <--
#m1 -0.3653   +       +   + 0.9359           +       +           +                 11 -220.029 465.4  4.79  0.078
#m4 -0.4252   +       +   + 0.9257           +                                      9 -222.972 466.2  5.54  0.053
#m2  0.6085   +       +   + 0.9877                   +           +                 10 -223.727 470.2  9.58  0.007
#m3  0.5033   +       +   + 0.9766                               +                  9 -225.292 470.8 10.18  0.005
#m5  0.5968   +       +   + 0.9840                                                  8 -226.809 471.4 10.74  0.004
AIC(m0,m1,m2,m3, m4,m5)
#.   df      AIC
#m0 12 456.6049 <--
#m1 11 462.0576
#m2 10 467.4545
#m3  9 468.5839
#m4  9 463.9445
#m5  8 469.6188
# Test
anova(m0, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Mulch                        9.214   9.214     1 23.839  1.3066 0.264364   
#Glyphosate                   9.143   9.143     1 23.838  1.2966 0.266158   
#Humic_acid                   0.817   0.817     1 27.219  0.1158 0.736237   
#z_logObs                    59.378  59.378     1 75.762  8.4202 0.004856 **
#Mulch:Glyphosate             0.450   0.450     1 25.019  0.0638 0.802631   
#Mulch:Humic_acid             0.309   0.309     1 26.828  0.0439 0.835697   
#Glyphosate:Humic_acid       39.620  39.620     1 26.402  5.6183 0.025350 * 
#Mulch:Glyphosate:Humic_acid 23.946  23.946     1 26.477  3.3957 0.076590 . 

emm<-emmeans(m0, specs=~Humic_acid:Mulch:Glyphosate, type="response")
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>% # plotting emmeans #
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic_acid, group = interaction(Glyphosate, Humic_acid))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.75), color = "black") +
  geom_point(data=md.fun.r.delta.2,mapping=aes(x=Mulch, y= Delta_InvSimpson, fill = Glyphosate, color = Glyphosate,shape=Humic_acid), alpha=0.3, position=position_jitterdodge(0.001), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.75)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 3, vjust = -4, hjust = 0.5,
            position = position_dodge(width = 0.75)) +
  ylab("Delta InvSimpson Index") +
  xlab("Treatment") +
  scale_x_discrete(labels=c("mulch"="mulch", "no"="no mulch"))+
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "fun.delta.InvSimpson.rt.mulch.Humic_acid.glyphosate.zlogObs.points.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

## Fun rhizos - Delta InvSimpson ####
md.fun.z<-data.frame(sample_data(fun.z.ps.T0_T3)) # extract metadata from the phyloseq object
md.fun.z$z_logObs<-scale(md.fun.z$logObs) # create new col with logObs z-score to control for sequencing depth variability
# calculate Delta InvSimpson for each tree at each
T0 <- md.fun.z %>%filter(Timepoint == "T0") %>%dplyr::select(Treat_rep_clean, InvSimpson) %>%rename(InvSimpson_T0 = InvSimpson)
md.fun.z.delta <- md.fun.z %>%left_join(T0, by = "Treat_rep_clean")
md.fun.z.delta <- md.fun.z.delta %>%mutate(Delta_InvSimpson = InvSimpson - InvSimpson_T0)
md.fun.z.delta <- md.fun.z.delta %>%filter(Timepoint != "T0")
# check for outliers
ggplot(md.fun.z.delta, aes(y=Delta_InvSimpson, x=Mulch, color=Glyphosate, shape=Humic_acid, label = rownames(md.fun.z.delta))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("33")) 
md.fun.z.delta.2<-md.fun.z.delta %>% filter(!rownames(md.fun.z.delta) %in% outliers)

# Mixed Effects Model
m0<-lmer(Delta_InvSimpson~Mulch*Glyphosate*Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta.2)) # model 1 
# No 3-way interaction
m1 <- lmer(Delta_InvSimpson ~Mulch + Glyphosate + Humic_acid + Mulch*Glyphosate +Mulch*Humic_acid+Glyphosate*Humic_acid+z_logObs+(1|Timepoint)+(1|Treat_rep_clean),data =data.frame(md.fun.z.delta.2))
# mulch interactions only
m2<-lmer(Delta_InvSimpson~Mulch * Glyphosate + Mulch * Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta.2)) # model 1 
# Mulch:Humic_acid interaction only
m3<-lmer(Delta_InvSimpson~ Mulch * Humic_acid +Glyphosate + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta.2)) # model 1 
# Gly:Humic_acid interaction only
m4<-lmer(Delta_InvSimpson~Mulch + Glyphosate * Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta.2)) # model 1 
# Additive
m5<-lmer(Delta_InvSimpson~Mulch + Glyphosate + Humic_acid + z_logObs + (1|Timepoint)+(1|Treat_rep_clean), data=data.frame(md.fun.z.delta.2)) # model 1 
# check fit
plot(resid(m0)~fitted(m0)) # looks good
plot(resid(m1)~fitted(m1)) # looks good
plot(resid(m2)~fitted(m2)) # looks good
plot(resid(m3)~fitted(m3)) # looks good
plot(resid(m4)~fitted(m4)) # looks good
plot(resid(m5)~fitted(m5)) # looks good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3, m4,m5)
#.    (Int) Gly Hmc_acd Mlc  z_lgO Gly:Hmc_acd Gly:Mlc Hmc_acd:Mlc Gly:Hmc_acd:Mlc df   logLik  AICc delta weight
#m0 -4.8900   +       +   + 3.377           +       +           +               + 12 -300.275 628.8  0.00  0.888 <--
#m1 -4.1670   +       +   + 3.378           +       +           +                 11 -303.859 633.2  4.47  0.095
#m2 -2.7890   +       +   + 3.299                   +           +                 10 -306.995 636.9  8.12  0.015
#m4 -1.5620   +       +   + 3.539           +                                      9 -310.877 642.1 13.32  0.001
#m3 -0.2776   +       +   + 3.470                               +                  9 -311.171 642.7 13.91  0.001
#m5 -0.3800   +       +   + 3.457                                                  8 -313.968 645.8 17.01  0.000
AIC(m0,m1,m2,m3, m4,m5)
#.   df      AIC
#m0 12 624.5504 <--
#m1 11 629.7176
#m2 10 633.9907
#m3  9 640.3419
#m4  9 639.7535
#m5  8 643.9353
# Test
anova(m0, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#                            Sum Sq  Mean Sq  NumDF  DenDF F value   Pr(>F)    
#Mulch                        40.28   40.28     1 24.490  0.5986 0.4465276    
#Glyphosate                  146.88  146.88     1 23.059  2.1828 0.1530890    
#Humic_acid                   32.16   32.16     1 23.003  0.4780 0.4962498    
#z_logObs                    913.87  913.87     1 53.991 13.5811 0.0005306 ***
#Mulch:Glyphosate            190.92  190.92     1 23.305  2.8373 0.1054466    
#Mulch:Humic_acid              0.00    0.00     1 23.013  0.0001 0.9937428    
#Glyphosate:Humic_acid        44.98   44.98     1 23.335  0.6685 0.4218519    
#Mulch:Glyphosate:Humic_acid  14.48   14.48     1 23.050  0.2152 0.6470380    

emm<-emmeans(m0, specs=~Humic_acid:Mulch:Glyphosate, type="response")
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>% # plotting emmeans #
  ggplot(aes(x = Mulch, y = emmean, fill = Glyphosate, color = Glyphosate,shape=Humic_acid, group = interaction(Glyphosate, Humic_acid))) +
  geom_pointrange(aes(ymin = emmean-SE, ymax = emmean+SE),size = 0.8, position = position_dodge(width = 0.75), color = "black") +
  geom_point(data=md.fun.z.delta.2,mapping=aes(x=Mulch, y= Delta_InvSimpson, fill = Glyphosate, color = Glyphosate,shape=Humic_acid), alpha=0.3, position=position_jitterdodge(0.001), size=1.5) +
  geom_point(size = 3.9, color = "black", position = position_dodge(width = 0.75)) +
  geom_point(size = 3.1, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = emmean + SE),
            color = "black", size = 3, vjust = -4, hjust = 0.5,
            position = position_dodge(width = 0.75)) +
  ylab("Delta InvSimpson Index") +
  xlab("Treatment") +
  scale_x_discrete(labels=c("mulch"="mulch", "no"="no mulch"))+
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("Glyphosate","no"),labels = c("Glyphosate","no Glyphosate")) +
  theme_bw()+ plot.theme
ggsave(filename= "fun.delta.InvSimpson.z.mulch.humic.glyphosate.zlogObs.points.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

###############################################################################
###### Beta Diversity - Comparing Timepoint/Year ##############################
###############################################################################
# Bacteria roots - T0-T3 ####
# Constrained ordination, n=121
cap.bac.r <- ordinate(bac.r.ps.T0_T3, method='CAP',distance='euclidean',formula=~Timepoint + Condition(scale(logObs)))
anova.cca(cap.bac.r)
#           Df Variance      F Pr(>F)    
#Model      3   313.63 4.0922  0.001 ***
#Residual 116  2963.46                  
scores.bac.r <- data.frame(scores(cap.bac.r, c(1,2,3,4),display='sites')) # Extract axis scores (points)
cap.bac.r.df<- cbind(sample_data(bac.r.ps.T0_T3),scores.bac.r) # bind scores to metadata
summary(cap.bac.r) # CAP1: 5.72%; CAP2: 2.76% (Proportion Explained)
RsquareAdj(cap.bac.r) # 0.08634053
# plot
ggplot(cap.bac.r.df,aes(x=CAP1,y=CAP2,color=Timepoint,fill = Timepoint))+
  geom_point(size = 3) + stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05) +
  scale_color_manual(values=c("#1b9e77","#f768a1","#ae017e","#49006a"),limits=c("T0","T1","T2","T3")) +
  scale_fill_manual(values=c("#1b9e77","#f768a1","#ae017e","#49006a"),limits=c("T0","T1","T2","T3")) +
  xlab("CAP1 [5.72%]") + ylab("CAP2 [2.76%]") + theme_bw() + theme(panel.grid = element_blank())
ggsave(filename= "T0_T3_CAP_timepoint_bac_r.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")
# Bacteria rhizosphere - T0-T3 ####
# Constrained ordination
cap.bac.z <- ordinate(bac.z.ps.T0_T3, method='CAP',distance='euclidean',formula=~Timepoint + Condition(scale(logObs)))
anova.cca(cap.bac.z)
#           Df Variance      F Pr(>F)    
#Model      3    658.0 4.5641  0.001 ***
#Residual 118   5670.3                  
scores.bac.z <- data.frame(scores(cap.bac.z, c(1,2,3,4),display='sites')) # Extract axis scores (points)
cap.bac.z.df<- cbind(sample_data(bac.z.ps.T0_T3),scores.bac.z) # bind scores to metadata
summary(cap.bac.z) # CAP1: 5.75%; CAP2: 3.05% (Proportion Explained)
RsquareAdj(cap.bac.z) # 0.09100988
# plot
ggplot(cap.bac.z.df,aes(x=CAP1,y=CAP2,color=Timepoint,fill = Timepoint))+
  geom_point(size = 3) + stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05) +
  scale_color_manual(values=c("#1b9e77","#f768a1","#ae017e","#49006a"),limits=c("T0","T1","T2","T3")) +
  scale_fill_manual(values=c("#1b9e77","#f768a1","#ae017e","#49006a"),limits=c("T0","T1","T2","T3")) +
  xlab("CAP1 [5.75%]") + ylab("CAP2 [3.052%]") + theme_bw() + theme(panel.grid = element_blank())
ggsave(filename= "T0_T3_CAP_timepoint.bac.z.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

# Fungi roots - T0-T3 ####
# Constrained ordination
cap.fun.r <- ordinate(fun.r.ps.T0_T3, method='CAP',distance='euclidean',formula=~Timepoint + Condition(scale(logObs)))
anova.cca(cap.fun.r)
#           Df Variance      F Pr(>F)    
#Model      3   72.89 2.9091  0.001 ***
#Residual 119   993.85                  
scores.fun.r <- data.frame(scores(cap.fun.r, c(1,2,3,4),display='sites')) # Extract axis scores (points)
cap.fun.r.df<- cbind(sample_data(fun.r.ps.T0_T3),scores.fun.r) # bind scores to metadata
summary(cap.fun.r) # CAP1: 4.61%; CAP2: 1.16% (Proportion Explained)
RsquareAdj(cap.fun.r) #0.0600901
# plot
ggplot(cap.fun.r.df,aes(x=CAP1,y=CAP2,color=Timepoint,fill = Timepoint))+
  geom_point(size = 3) + stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05) +
  scale_color_manual(values=c("#1b9e77","#f768a1","#ae017e","#49006a"),limits=c("T0","T1","T2","T3")) +
  scale_fill_manual(values=c("#1b9e77","#f768a1","#ae017e","#49006a"),limits=c("T0","T1","T2","T3")) +
  xlab("CAP1 [4.61%]") + ylab("CAP2 [1.16%]") + theme_bw() + theme(panel.grid = element_blank())
ggsave(filename= "T0_T3_CAP_timepoint.fun.r.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")
# Fungi rhizosphere - T0-T3 ####
# Constrained ordination
cap.fun.z <- ordinate(fun.z.ps.T0_T3, method='CAP',distance='euclidean',formula=~Timepoint + Condition(scale(logObs)))
anova.cca(cap.fun.z)
#           Df Variance      F Pr(>F)    
#Model      3   229.34 4.1577  0.001 ***
#Residual 117  2151.23                  
scores.fun.z <- data.frame(scores(cap.fun.z, c(1,2,3,4),display='sites')) # Extract axis scores (points)
cap.fun.z.df<- cbind(sample_data(fun.z.ps.T0_T3),scores.fun.z) # bind scores to metadata
summary(cap.fun.z) # CAP1: 6.17%; CAP2: 1.85% (Proportion Explained)
RsquareAdj(cap.fun.z) # 0.07976562
# plot
ggplot(cap.fun.z.df,aes(x=CAP1,y=CAP2,color=Timepoint,fill = Timepoint))+
  geom_point(size = 3) + stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05) +
  scale_color_manual(values=c("#1b9e77","#f768a1","#ae017e","#49006a"),limits=c("T0","T1","T2","T3")) +
  scale_fill_manual(values=c("#1b9e77","#f768a1","#ae017e","#49006a"),limits=c("T0","T1","T2","T3")) +
  xlab("CAP1 [6.17%]") + ylab("CAP2 [1.85%]") + theme_bw() + theme(panel.grid = element_blank())
ggsave(filename= "T0_T3_CAP_timepoint.fun.z.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

###############################################################################
###### Beta Diversity - Comparing plots at Timepoint Zero #####################
###############################################################################
# Bacteria roots - T0 only ####
cap.bac.r <- ordinate(bac.r.ps.T0, method='CAP',distance='euclidean',formula=~Mulch * Glyphosate * Humic + Condition(scale(logObs)))
anova.cca(cap.bac.r)
#           Df Variance      F Pr(>F)    
#Model       7   835.86 1.1812  0.029 *
#Residual 23  2324.99    
anova.cca(cap.bac.r, by="term")
#                         Df Variance      F Pr(>F) 
# Mulch                   1   146.19 1.4462  0.040 *
# Glyphosate              1    90.93 0.8995  0.685  
# Humic                   1   132.18 1.3076  0.081 .
# Mulch:Glyphosate        1   148.67 1.4707  0.038 *
# Mulch:Humic             1   122.32 1.2100  0.122  
# Glyphosate:Humic        1    86.44 0.8551  0.800  
# Mulch:Glyphosate:Humic  1   109.13 1.0796  0.270  
# Residual               23  2324.99  
RsquareAdj(cap.bac.r) # 0.2447975
scores.bac.r <- data.frame(scores(cap.bac.r, c(1,2,3,4),display='sites')) # Extract axis scores (points)
cap.bac.r.df<- cbind(sample_data(bac.r.ps.T0),scores.bac.r) # bind scores to metadata
summary(cap.bac.r) # Check Proportion Explained for CAP1 and CAP2
# plot
ggplot(cap.bac.r.df,aes(x=CAP1,y=CAP2,color=Treatment_code,fill = Treatment_code))+
  geom_point(size = 3) +
  #stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  xlab("CAP1 [8.25%]") + ylab("CAP2 [4.36%]") + theme_bw() + theme(panel.grid = element_blank())
ggsave(filename= "T0_CAP_treat_bac.r.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

# Bacteria rhizosphere - T0 only ####
cap.bac.z <- ordinate(bac.z.ps.T0, method='CAP',distance='euclidean',formula=~Mulch * Glyphosate * Humic + Condition(scale(logObs)))
anova.cca(cap.bac.z)
#           Df Variance      F Pr(>F)    
#Model       7   1552.9 1.2586  0.004 **
#Residual   23   4053.9    
anova.cca(cap.bac.z, by="term")
#                         Df Variance      F Pr(>F) 
#Mulch                   1    253.2 1.4363  0.043 * 
#Glyphosate              1    217.9 1.2362  0.110   
#Humic                   1    234.5 1.3307  0.064 . 
#Mulch:Glyphosate        1    284.2 1.6124  0.007 **
#Mulch:Humic             1    221.2 1.2552  0.120   
#Glyphosate:Humic        1    160.2 0.9087  0.637   
#Mulch:Glyphosate:Humic  1    181.7 1.0309  0.330   
#Residual               23   4053.9                 
RsquareAdj(cap.bac.z) # 0.2367647
scores.bac.z <- data.frame(scores(cap.bac.z, c(1,2,3,4),display='sites')) # Extract axis scores (points)
cap.bac.z.df<- cbind(sample_data(bac.z.ps.T0),scores.bac.z) # bind scores to metadata
summary(cap.bac.z) # Check Proportion Explained for CAP1 and CAP2
# plot
ggplot(cap.bac.z.df,aes(x=CAP1,y=CAP2,color=Treatment_code,fill = Treatment_code))+
  geom_point(size = 3) +
  #stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  xlab("CAP1 [9.03%]") + ylab("CAP2 [5.25%]") + theme_bw() + theme(panel.grid = element_blank())
ggsave(filename= "T0_CAP_treat_bac.z_updated.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

# fungi roots - T0 only ####
cap.fun.r <- ordinate(fun.r.ps.T0, method='CAP',distance='euclidean',formula=~Mulch * Glyphosate * Humic_acid + Condition(scale(logObs)))
anova.cca(cap.fun.r)
#           Df Variance      F Pr(>F)    
#Model       7   181.86 0.9291   0.68
#Residual   23   643.15    
anova.cca(cap.fun.r, by="term")
#                            Df Variance      F Pr(>F) 
#Mulch                        1    26.34 0.9421  0.496
#Glyphosate                   1    30.17 1.0789  0.265
#Humic_acid                   1    25.04 0.8954  0.597
#Mulch:Glyphosate             1    30.62 1.0949  0.249
#Mulch:Humic_acid             1    19.49 0.6970  0.948
#Glyphosate:Humic_acid        1    24.49 0.8758  0.627
#Mulch:Glyphosate:Humic_acid  1    25.72 0.9198  0.534
#Residual                    23   643.15   
RsquareAdj(cap.fun.r) # 0.1936004
scores.fun.r <- data.frame(scores(cap.fun.r, c(1,2,3,4),display='sites')) # Extract axis scores (points)
cap.fun.r.df<- cbind(sample_data(fun.r.ps.T0),scores.fun.r) # bind scores to metadata
summary(cap.fun.r) # Check Proportion Explained for CAP1 and CAP2
# plot
ggplot(cap.fun.r.df,aes(x=CAP1,y=CAP2,color=Treatment_code,fill = Treatment_code))+
  geom_point(size = 3) +
  #stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  xlab("CAP1 [6.24%]") + ylab("CAP2 [4.26%]") + theme_bw() + theme(panel.grid = element_blank())
ggsave(filename= "T0_CAP_treat_fun.r.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

# fungi rhizosphere - T0 only ####
cap.fun.z <- ordinate(fun.z.ps.T0, method='CAP',distance='euclidean',formula=~Mulch * Glyphosate * Humic_acid + Condition(scale(logObs)))
anova.cca(cap.fun.z)
#           Df Variance      F Pr(>F)    
# Model     7   439.78 1.1534  0.084 .
# Residual 22  1198.34               
anova.cca(cap.fun.z, by="term")
#                         Df Variance      F Pr(>F) 
#Mulch                        1    73.73 1.3535  0.088 .
#Glyphosate                   1    64.43 1.1828  0.194  
#Humic_acid                   1    45.54 0.8361  0.769  
#Mulch:Glyphosate             1    75.07 1.3781  0.070 .
#Mulch:Humic_acid             1    74.53 1.3683  0.084 .
#Glyphosate:Humic_acid        1    53.73 0.9863  0.398  
#Mulch:Glyphosate:Humic_acid  1    52.76 0.9686  0.401  
#Residual                    22  1198.34                
RsquareAdj(cap.fun.z) # 0.2209673
scores.fun.z <- data.frame(scores(cap.fun.z, c(1,2,3,4),display='sites')) # Extract axis scores (points)
cap.fun.z.df<- cbind(sample_data(fun.z.ps.T0),scores.fun.z) # bind scores to metadata
summary(cap.fun.z) # Check Proportion Explained for CAP1 and CAP2
# plot
ggplot(cap.fun.z.df,aes(x=CAP1,y=CAP2,color=Treatment_code,fill = Treatment_code))+
  geom_point(size = 3) +
  #stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  xlab("CAP1 [7.83%]") + ylab("CAP2 [4.89%]") + theme_bw() + theme(panel.grid = element_blank())
ggsave(filename= "T0_CAP_treat_fun.z_updated.pdf", device="pdf", units="mm", dpi=300, width=100, height=75, path="Revisions.ISME.Comm/plots/")

###############################################################################
# Beta Diversity - Comparing treatments and controlling for pre-treat diffs ###
###############################################################################
# Bacterial roots ####
# T0 partial ordination: removes sequencing depth, no treatment constraints
cap.bac.r.T0.partial <- ordinate(bac.r.ps.T0,method = "CAP",distance = "euclidean",formula = ~ Condition(scale(logObs)))
# Extract first 3 baseline axes
T0_scores <- as.data.frame(scores(cap.bac.r.T0.partial, display = "sites", choices = 1:3)) %>%rownames_to_column("SampleID.1") %>%
  rename(T0_MDS1 = MDS1,T0_MDS2 = MDS2,T0_MDS3 = MDS3)
# Get tree IDs from T0 metadata
T0_meta <- data.frame(sample_data(bac.r.ps.T0)) %>%
  tibble::rownames_to_column("SampleID.1") %>%
  dplyr::filter(!is.na(SampleID.1), !is.na(Treat_rep_clean)) %>%
  dplyr::select(SampleID.1, Treat_rep_clean)
# Combine T0 scores with tree IDs
T0_baseline_axes <- T0_scores %>%
  left_join(T0_meta, by = "SampleID.1") %>%
  dplyr::select(Treat_rep_clean, T0_MDS1, T0_MDS2, T0_MDS3)
# Add T0 baseline axes to T1-T3 metadata
treat_meta <- data.frame(sample_data(bac.r.ps.treat)) %>%
  rownames_to_column("SampleID.1") %>%
  left_join(T0_baseline_axes, by = "Treat_rep_clean")
# Check matches
table(is.na(treat_meta$T0_MDS1))
# Replace metadata in T1-T3 phyloseq object
rownames(treat_meta) <- treat_meta$SampleID.1
treat_meta$SampleID.1 <- NULL
sample_data(bac.r.ps.treat) <- sample_data(treat_meta)
# T1-T3 treatment CAP while conditioning on: sequencing depth, timepoint, and baseline T0 community axes
cap.bac.r.baseline <- ordinate(bac.r.ps.treat,method = "CAP",distance = "euclidean",formula =~
                                 Mulch * Glyphosate + Humic + # no significant humic interaction terms, model simplified
                                 Condition(scale(logObs)) +
                                 Condition(Timepoint) +
                                 Condition(T0_MDS1) +
                                 Condition(T0_MDS2) +
                                 Condition(T0_MDS3))
# Stats
set.seed(444)
anova.cca(cap.bac.r.baseline)
#          Df Variance      F Pr(>F)    
#Model     4   189.17 1.4544  0.001 ***
#Residual 78  2536.27          
anova.cca(cap.bac.r.baseline, by = "terms")
##          Df Variance      F Pr(>F)
#Mulch             1    65.09 2.0018  0.001 ***
#Glyphosate        1    39.19 1.2052  0.073 .  
#Humic             1    43.11 1.3259  0.026 *  
#Mulch:Glyphosate  1    41.78 1.2847  0.043 *  
#Residual         78  2536.27  
summary(cap.bac.r.baseline) # 0.05392 R2
RsquareAdj(cap.bac.r.baseline) #0.05392367
# Plotting
RDAscores.bac.r <- data.frame(scores(cap.bac.r.baseline, c(1,2,3,4),display='sites'))
cap.r.df<- cbind(sample_data(bac.r.ps.treat),RDAscores.bac.r)
# Plot
ggplot(cap.r.df,aes(x=CAP1,y=CAP2,color=Mulch,fill = Mulch, shape=Humic))+
  geom_point(size = 3) +
  stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05,aes(fill = Mulch)) +
  scale_shape_manual(values=c(1, 17),limits=c("humic","no"),labels = c("humic acid","no humic acid")) +
  scale_color_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  xlab("CAP 1 [2.49%]") +
  facet_wrap(.~Glyphosate) +
  ylab("CAP 2 [1.80%]") +
  theme_bw() +
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
ggsave(filename= "CAP.bac.r.humic.gly.mulch.baseline.ctrl.treat.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="Revisions.ISME.Comm/plots/")

# Bacterial rhizosphere ####
# T0 partial ordination: removes sequencing depth, no treatment constraints
cap.bac.z.T0.partial <- ordinate(bac.z.ps.T0,method = "CAP",distance = "euclidean",formula = ~ Condition(scale(logObs)))
# Extract first 3 baseline axes
T0_scores <- as.data.frame(scores(cap.bac.z.T0.partial, display = "sites", choices = 1:3)) %>%rownames_to_column("SampleID.1") %>%
  rename(T0_MDS1 = MDS1,T0_MDS2 = MDS2,T0_MDS3 = MDS3)
# Get tree IDs from T0 metadata
T0_meta <- data.frame(sample_data(bac.z.ps.T0)) %>%
  tibble::rownames_to_column("SampleID.1") %>%
  dplyr::filter(!is.na(SampleID.1), !is.na(Treat_rep_clean)) %>%
  dplyr::select(SampleID.1, Treat_rep_clean)
# Combine T0 scores with tree IDs
T0_baseline_axes <- T0_scores %>%
  left_join(T0_meta, by = "SampleID.1") %>%
  dplyr::select(Treat_rep_clean, T0_MDS1, T0_MDS2, T0_MDS3)
# Add T0 baseline axes to T1-T3 metadata
treat_meta <- data.frame(sample_data(bac.z.ps.treat)) %>%
  rownames_to_column("SampleID.1") %>%
  left_join(T0_baseline_axes, by = "Treat_rep_clean")
# Check matches
table(is.na(treat_meta$T0_MDS1))
# Replace metadata in T1-T3 phyloseq object
rownames(treat_meta) <- treat_meta$SampleID.1
treat_meta$SampleID.1 <- NULL
sample_data(bac.z.ps.treat) <- sample_data(treat_meta)
# T1-T3 treatment CAP while conditioning on: sequencing depth, timepoint, and baseline T0 community axes
cap.bac.z.baseline <- ordinate(bac.z.ps.treat,method = "CAP",distance = "euclidean",formula =~
                                 Mulch * Glyphosate + Glyphosate* Humic + # no significant 3-way interaction , model simplified
                                 Condition(scale(logObs)) +
                                 Condition(Timepoint) +
                                 Condition(T0_MDS1) +
                                 Condition(T0_MDS2) +
                                 Condition(T0_MDS3))
# Stats
set.seed(444)
anova.cca(cap.bac.z.baseline)
#          Df Variance      F Pr(>F)    
#Model     5    557.6 1.8443  0.001 ***
#Residual 79   4777.0          
anova.cca(cap.bac.z.baseline, by = "terms")
##          Df Variance      F Pr(>F)
#Mulch             1    239.0 3.9526  0.001 ***
#Glyphosate        1     69.8 1.1547  0.149    
#Humic             1     82.0 1.3554  0.033 *  
#Mulch:Glyphosate  1     87.0 1.4394  0.012 *  
#Glyphosate:Humic  1     79.8 1.3193  0.031 *  
#Residual         79   4777.0     
summary(cap.bac.z.baseline) # 0.07797 R2
# Plotting
RDAscores.bac.z <- data.frame(scores(cap.bac.z.baseline, c(1,2,3,4),display='sites'))
cap.z.df<- cbind(sample_data(bac.z.ps.treat),RDAscores.bac.z)
# Plot
ggplot(cap.z.df,aes(x=CAP1,y=CAP2,color=Mulch,fill = Mulch, shape=Humic))+
  geom_point(size = 3) +
  stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05,aes(fill = Mulch)) +
  scale_shape_manual(values=c(1, 17),limits=c("humic","no"),labels = c("humic acid","no humic acid")) +
  scale_color_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  xlab("CAP 1 [4.58%]") +
  facet_wrap(.~Glyphosate) +
  ylab("CAP 2 [1.74%]") +
  theme_bw() +
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
ggsave(filename= "CAP.bac.z.humic.gly.mulch.baseline.ctrl.treat.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="Revisions.ISME.Comm/plots/")

# Fungi roots ####
# T0 partial ordination: removes sequencing depth, no treatment constraints
cap.fun.r.T0.partial <- ordinate(fun.r.ps.T0,method = "CAP",distance = "euclidean",formula = ~ Condition(scale(logObs)))
# Extract first 3 baseline axes
T0_scores <- as.data.frame(scores(cap.fun.r.T0.partial, display = "sites", choices = 1:3)) %>%rownames_to_column("SampleID.1") %>%
  rename(T0_MDS1 = MDS1,T0_MDS2 = MDS2,T0_MDS3 = MDS3)
# Get tree IDs from T0 metadata
T0_meta <- data.frame(sample_data(fun.r.ps.T0)) %>%
  tibble::rownames_to_column("SampleID.1") %>%
  dplyr::filter(!is.na(SampleID.1), !is.na(Treat_rep_clean)) %>%
  dplyr::select(SampleID.1, Treat_rep_clean)
# Combine T0 scores with tree IDs
T0_baseline_axes <- T0_scores %>%
  left_join(T0_meta, by = "SampleID.1") %>%
  dplyr::select(Treat_rep_clean, T0_MDS1, T0_MDS2, T0_MDS3)
# Add T0 baseline axes to T1-T3 metadata
treat_meta <- data.frame(sample_data(fun.r.ps.treat)) %>%
  rownames_to_column("SampleID.1") %>%
  left_join(T0_baseline_axes, by = "Treat_rep_clean")
# Check matches
table(is.na(treat_meta$T0_MDS1))
# Replace metadata in T1-T3 phyloseq object
rownames(treat_meta) <- treat_meta$SampleID.1
treat_meta$SampleID.1 <- NULL
sample_data(fun.r.ps.treat) <- sample_data(treat_meta)
# T1-T3 treatment CAP while conditioning on: sequencing depth, timepoint, and baseline T0 community axes
cap.fun.r.baseline <- ordinate(fun.r.ps.treat,method = "CAP",distance = "euclidean",formula =~
                                 Mulch + Glyphosate + Humic_acid + # no significant humic interaction terms, model simplified
                                 Condition(scale(logObs)) +
                                 Condition(Timepoint) +
                                 Condition(T0_MDS1) +
                                 Condition(T0_MDS2) +
                                 Condition(T0_MDS3))
# Stats
set.seed(444)
anova.cca(cap.fun.r.baseline)
#          Df Variance      F Pr(>F)    
#Model     3    50.08 1.4387  0.006 **
#Residual 82   951.49          
anova.cca(cap.fun.r.baseline, by = "terms")
##          Df Variance      F Pr(>F)
#Mulch       1    23.87 2.0568  0.004 **
#Glyphosate  1    13.09 1.1284  0.208   
#Humic_acid  1    13.12 1.1310  0.207   
#Residual   82   951.49                 
summary(cap.fun.r.baseline) # 0.03873 R2
# Plotting
RDAscores.fun.r <- data.frame(scores(cap.fun.r.baseline, c(1,2,3,4),display='sites'))
cap.r.df<- cbind(sample_data(fun.r.ps.treat),RDAscores.fun.r)
# Plot
ggplot(cap.r.df,aes(x=CAP1,y=CAP2,color=Mulch,fill = Mulch, shape=Humic_acid))+
  geom_point(size = 3) +
  stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05,aes(fill = Mulch)) +
  scale_shape_manual(values=c(1, 17),limits=c("humic","no"),labels = c("humic acid","no humic acid")) +
  scale_color_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  xlab("CAP 1 [2.56%]") +
  facet_wrap(.~Glyphosate) +
  ylab("CAP 2 [1.39%]") +
  theme_bw() +
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
ggsave(filename= "CAP.fun.r.humic.gly.mulch.baseline.ctrl.treat.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="Revisions.ISME.Comm/plots/")

# Fungi rhizosphere ####
# T0 partial ordination: removes sequencing depth, no treatment constraints
cap.fun.z.T0.partial <- ordinate(fun.z.ps.T0,method = "CAP",distance = "euclidean",formula = ~ Condition(scale(logObs)))
# Extract first 3 baseline axes
T0_scores <- as.data.frame(scores(cap.fun.z.T0.partial, display = "sites", choices = 1:3)) %>%rownames_to_column("SampleID.1") %>%
  rename(T0_MDS1 = MDS1,T0_MDS2 = MDS2,T0_MDS3 = MDS3)
# Get tree IDs from T0 metadata
T0_meta <- data.frame(sample_data(fun.z.ps.T0)) %>%
  tibble::rownames_to_column("SampleID.1") %>%
  dplyr::filter(!is.na(SampleID.1), !is.na(Treat_rep_clean)) %>%
  dplyr::select(SampleID.1, Treat_rep_clean)
# Combine T0 scores with tree IDs
T0_baseline_axes <- T0_scores %>%
  left_join(T0_meta, by = "SampleID.1") %>%
  dplyr::select(Treat_rep_clean, T0_MDS1, T0_MDS2, T0_MDS3)
# Add T0 baseline axes to T1-T3 metadata
treat_meta <- data.frame(sample_data(fun.z.ps.treat)) %>%
  rownames_to_column("SampleID.1") %>%
  left_join(T0_baseline_axes, by = "Treat_rep_clean")
# Check matches
table(is.na(treat_meta$T0_MDS1))
treat_meta %>%filter(is.na(T0_MDS1)) %>%dplyr::select(Treat_rep_clean, Timepoint, Treatment_code) # check the TRUE
#Z.T0.CCC4 missing, can not match to post-treatment CCC4 samples; remove CCC4
# Remove samples without T0 baseline match
treat_meta.filt <- treat_meta %>%filter(!is.na(T0_MDS1))
rownames(treat_meta.filt) <- treat_meta.filt$SampleID.1
treat_meta.filt$SampleID.1 <- NULL
fun.z.ps.treat.baseline <- prune_samples(rownames(treat_meta.filt),fun.z.ps.treat)
treat_meta.filt <- treat_meta.filt[sample_names(fun.z.ps.treat.baseline),]
sample_data(fun.z.ps.treat.baseline) <- sample_data(treat_meta.filt)
# T1-T3 treatment CAP while conditioning on: sequencing depth, timepoint, and baseline T0 community axes
cap.fun.z.baseline <- ordinate(fun.z.ps.treat.baseline,method = "CAP",distance = "euclidean",formula =~
                                 Mulch + Glyphosate+ Humic_acid + # no significant 3-way interaction , model simplified
                                 Condition(scale(logObs)) +
                                 Condition(Timepoint) +
                                 Condition(T0_MDS1) +
                                 Condition(T0_MDS2) +
                                 Condition(T0_MDS3))
# Stats
set.seed(444)
anova.cca(cap.fun.z.baseline)
#          Df Variance      F Pr(>F)    
#Model     3   135.59 1.7624  0.002 **
#Residual 78  2000.31          
anova.cca(cap.fun.z.baseline, by = "terms")
##          Df Variance      F Pr(>F)
#Mulch       1    73.25 2.8562  0.002 **
#Glyphosate  1    39.13 1.5258  0.067 . 
#Humic_acid  1    23.21 0.9051  0.515   
#Residual   78  2000.31  
summary(cap.fun.z.baseline) # 0.04438 R2
# Plotting
RDAscores.fun.z <- data.frame(scores(cap.fun.z.baseline, c(1,2,3,4),display='sites'))
cap.z.df<- cbind(sample_data(fun.z.ps.treat.baseline),RDAscores.fun.z)
# Plot
ggplot(cap.z.df,aes(x=CAP1,y=CAP2,color=Mulch,fill = Mulch, shape=Humic_acid))+
  geom_point(size = 3) +
  stat_ellipse(geom="polygon",linewidth=0.4,alpha = 0.05,aes(fill = Mulch)) +
  scale_shape_manual(values=c(1, 17),limits=c("humic","no"),labels = c("humic acid","no humic acid")) +
  scale_color_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  xlab("CAP 1 [3.49%]") +
  facet_wrap(.~Glyphosate) +
  ylab("CAP 2 [1.83%]") +
  theme_bw() +
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
ggsave(filename= "CAP.fun.z.humic.gly.mulch.baseline.ctrl.treat.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="Revisions.ISME.Comm/plots/")
###############################################################################
###### ALDEX2 Genera Differential Abundance Pre vs. Post Treatment ############
###############################################################################
# Add pre vs. post treatment column ####
# bac.r
meta.bac.r <- data.frame(sample_data(bac.r.ps.raw.T0_T3))
meta.bac.r$PrePostTreat <- ifelse(meta.bac.r$Timepoint == "T0","1.pre-treatment","2.post-treatment")
meta.bac.r$PrePostTreat <- factor(meta.bac.r$PrePostTreat,levels = c("1.pre-treatment", "2.post-treatment"))
sample_data(bac.r.ps.raw.T0_T3) <- sample_data(meta.bac.r)
# bac.z
meta.bac.z <- data.frame(sample_data(bac.z.ps.raw.T0_T3))
meta.bac.z$PrePostTreat <- ifelse(meta.bac.z$Timepoint == "T0","1.pre-treatment","2.post-treatment")
meta.bac.z$PrePostTreat <- factor(meta.bac.z$PrePostTreat,levels = c("1.pre-treatment", "2.post-treatment"))
sample_data(bac.z.ps.raw.T0_T3) <- sample_data(meta.bac.z)
# fun.r
meta.fun.r <- data.frame(sample_data(fun.r.ps.raw.T0_T3))
meta.fun.r$PrePostTreat <- ifelse(meta.fun.r$Timepoint == "T0","1.pre-treatment","2.post-treatment")
meta.fun.r$PrePostTreat <- factor(meta.fun.r$PrePostTreat,levels = c("1.pre-treatment", "2.post-treatment"))
sample_data(fun.r.ps.raw.T0_T3) <- sample_data(meta.fun.r)
# fun.z
meta.fun.z <- data.frame(sample_data(fun.z.ps.raw.T0_T3))
meta.fun.z$PrePostTreat <- ifelse(meta.fun.z$Timepoint == "T0","1.pre-treatment","2.post-treatment")
meta.fun.z$PrePostTreat <- factor(meta.fun.z$PrePostTreat,levels = c("1.pre-treatment", "2.post-treatment"))
sample_data(fun.z.ps.raw.T0_T3) <- sample_data(meta.fun.z)
#### agglomerate to Genus level ####
set.seed(4444)
# bac root #
bac.r.ps.raw.T0_T3.genus <- tax_glom(bac.r.ps.raw.T0_T3, taxrank = "Genus", NArm = FALSE)
tax <- tax_table(bac.r.ps.raw.T0_T3.genus)
taxa_names(bac.r.ps.raw.T0_T3.genus) <- make.unique(as.character(tax[, "Genus"])) # Some genus do not have names ("NA"). Number the NAs so they are unique (NA.1, NA.2,...)
taxa_names(bac.r.ps.raw.T0_T3.genus)[is.na(taxa_names(bac.r.ps.raw.T0_T3.genus))] <- "NA.0" # add a .0 to the first NA
# bac rhzo #
bac.z.ps.raw.T0_T3.genus <- tax_glom(bac.z.ps.raw.T0_T3, taxrank = "Genus", NArm = FALSE)
tax <- tax_table(bac.z.ps.raw.T0_T3.genus)
taxa_names(bac.z.ps.raw.T0_T3.genus) <- make.unique(as.character(tax[, "Genus"]))
taxa_names(bac.z.ps.raw.T0_T3.genus)[is.na(taxa_names(bac.z.ps.raw.T0_T3.genus))] <- "NA.0"
# fun root #
fun.r.ps.raw.T0_T3.genus <- tax_glom(fun.r.ps.raw.T0_T3, taxrank = "Genus", NArm = FALSE)
tax <- tax_table(fun.r.ps.raw.T0_T3.genus)
taxa_names(fun.r.ps.raw.T0_T3.genus) <- make.unique(as.character(tax[, "Genus"]))
taxa_names(fun.r.ps.raw.T0_T3.genus)[is.na(taxa_names(fun.r.ps.raw.T0_T3.genus))] <- "NA.0"
# fun rhzo #
fun.z.ps.raw.T0_T3.genus <- tax_glom(fun.z.ps.raw.T0_T3, taxrank = "Genus", NArm = FALSE)
tax <- tax_table(fun.z.ps.raw.T0_T3.genus)
taxa_names(fun.z.ps.raw.T0_T3.genus) <- make.unique(as.character(tax[, "Genus"]))
taxa_names(fun.z.ps.raw.T0_T3.genus)[is.na(taxa_names(fun.z.ps.raw.T0_T3.genus))] <- "NA.0"
## Bac root - testing and plotting - All 8 treatment groups ####
# subseting treatments
bac.r.ps.raw.T0_T3.genus.CCC <- subset_samples(bac.r.ps.raw.T0_T3.genus, Treatment_code == "CCC")
bac.r.ps.raw.T0_T3.genus.CCH <- subset_samples(bac.r.ps.raw.T0_T3.genus, Treatment_code == "CCH") 
bac.r.ps.raw.T0_T3.genus.CGH <- subset_samples(bac.r.ps.raw.T0_T3.genus, Treatment_code == "CGH") 
bac.r.ps.raw.T0_T3.genus.CGC <- subset_samples(bac.r.ps.raw.T0_T3.genus, Treatment_code == "CGC") 
bac.r.ps.raw.T0_T3.genus.MCC <- subset_samples(bac.r.ps.raw.T0_T3.genus, Treatment_code == "MCC") 
bac.r.ps.raw.T0_T3.genus.MCH <- subset_samples(bac.r.ps.raw.T0_T3.genus, Treatment_code == "MCH") 
bac.r.ps.raw.T0_T3.genus.MGH <- subset_samples(bac.r.ps.raw.T0_T3.genus, Treatment_code == "MGH") 
bac.r.ps.raw.T0_T3.genus.MGC <- subset_samples(bac.r.ps.raw.T0_T3.genus, Treatment_code == "MGC") 
# Set DF for test variable
bac.r.CCC.pp <- sample_data(bac.r.ps.raw.T0_T3.genus.CCC)$PrePostTreat %>% as.character() # 2 levels
bac.r.CCH.pp <- sample_data(bac.r.ps.raw.T0_T3.genus.CCH)$PrePostTreat %>% as.character() # 2 levels
bac.r.CGH.pp <- sample_data(bac.r.ps.raw.T0_T3.genus.CGH)$PrePostTreat %>% as.character() # 2 levels
bac.r.CGC.pp <- sample_data(bac.r.ps.raw.T0_T3.genus.CGC)$PrePostTreat %>% as.character() # 2 levels
bac.r.MCC.pp <- sample_data(bac.r.ps.raw.T0_T3.genus.MCC)$PrePostTreat %>% as.character() # 2 levels
bac.r.MCH.pp <- sample_data(bac.r.ps.raw.T0_T3.genus.MCH)$PrePostTreat %>% as.character() # 2 levels
bac.r.MGH.pp <- sample_data(bac.r.ps.raw.T0_T3.genus.MGH)$PrePostTreat %>% as.character() # 2 levels
bac.r.MGC.pp <- sample_data(bac.r.ps.raw.T0_T3.genus.MGC)$PrePostTreat %>% as.character() # 2 levels
# Extract genus table
genus.bac.r.CCC <- as(otu_table(bac.r.ps.raw.T0_T3.genus.CCC, taxa_are_rows = TRUE), 'matrix') 
genus.bac.r.CCC <- t(genus.bac.r.CCC)
genus.bac.r.CCH <- as(otu_table(bac.r.ps.raw.T0_T3.genus.CCH, taxa_are_rows = TRUE), 'matrix') 
genus.bac.r.CCH <- t(genus.bac.r.CCH)
genus.bac.r.CGH <- as(otu_table(bac.r.ps.raw.T0_T3.genus.CGH, taxa_are_rows = TRUE), 'matrix') 
genus.bac.r.CGH <- t(genus.bac.r.CGH)
genus.bac.r.CGC <- as(otu_table(bac.r.ps.raw.T0_T3.genus.CGC, taxa_are_rows = TRUE), 'matrix') 
genus.bac.r.CGC <- t(genus.bac.r.CGC)
genus.bac.r.MCC <- as(otu_table(bac.r.ps.raw.T0_T3.genus.MCC, taxa_are_rows = TRUE), 'matrix') 
genus.bac.r.MCC <- t(genus.bac.r.MCC)
genus.bac.r.MCH <- as(otu_table(bac.r.ps.raw.T0_T3.genus.MCH, taxa_are_rows = TRUE), 'matrix') 
genus.bac.r.MCH <- t(genus.bac.r.MCH)
genus.bac.r.MGH <- as(otu_table(bac.r.ps.raw.T0_T3.genus.MGH, taxa_are_rows = TRUE), 'matrix') 
genus.bac.r.MGH <- t(genus.bac.r.MGH)
genus.bac.r.MGC <- as(otu_table(bac.r.ps.raw.T0_T3.genus.MGC, taxa_are_rows = TRUE), 'matrix') 
genus.bac.r.MGC <- t(genus.bac.r.MGC)
#### 2-level factors aldex CLR/MC sampling ##
bac.r.CCC.pp.clr <- aldex.clr(genus.bac.r.CCC, conds=bac.r.CCC.pp, mc.samples=200) #
bac.r.CCH.pp.clr <- aldex.clr(genus.bac.r.CCH, conds=bac.r.CCH.pp, mc.samples=200) # 
bac.r.CGH.pp.clr <- aldex.clr(genus.bac.r.CGH, conds=bac.r.CGH.pp, mc.samples=200) # 
bac.r.CGC.pp.clr <- aldex.clr(genus.bac.r.CGC, conds=bac.r.CGC.pp, mc.samples=200) # 
bac.r.MCC.pp.clr <- aldex.clr(genus.bac.r.MCC, conds=bac.r.MCC.pp, mc.samples=200) #
bac.r.MCH.pp.clr <- aldex.clr(genus.bac.r.MCH, conds=bac.r.MCH.pp, mc.samples=200) # 
bac.r.MGH.pp.clr <- aldex.clr(genus.bac.r.MGH, conds=bac.r.MGH.pp, mc.samples=200) # 
bac.r.MGC.pp.clr <- aldex.clr(genus.bac.r.MGC, conds=bac.r.MGC.pp, mc.samples=200) # 
#### 2 level factor statistical test ##
ttest.bac.r.CCC.pp <- aldex.ttest(bac.r.CCC.pp.clr)
ttest.bac.r.CCH.pp <- aldex.ttest(bac.r.CCH.pp.clr)
ttest.bac.r.CGH.pp <- aldex.ttest(bac.r.CGH.pp.clr)
ttest.bac.r.CGC.pp <- aldex.ttest(bac.r.CGC.pp.clr)
ttest.bac.r.MCC.pp <- aldex.ttest(bac.r.MCC.pp.clr)
ttest.bac.r.MCH.pp <- aldex.ttest(bac.r.MCH.pp.clr)
ttest.bac.r.MGH.pp <- aldex.ttest(bac.r.MGH.pp.clr)
ttest.bac.r.MGC.pp <- aldex.ttest(bac.r.MGC.pp.clr)
#### check for significant genera; BH = adjusted p-value ##
sum(ttest.bac.r.CCC.pp$wi.eBH < 0.05) # Wilcoxon test # 0
sum(ttest.bac.r.CCC.pp$we.eBH < 0.05) # Welch’s t test # 2
sum(ttest.bac.r.CCH.pp$wi.eBH < 0.05) # 0
sum(ttest.bac.r.CCH.pp$we.eBH < 0.05) # 0
sum(ttest.bac.r.CGH.pp$wi.eBH < 0.05) # 0
sum(ttest.bac.r.CGH.pp$we.eBH < 0.05) # 0
sum(ttest.bac.r.CGC.pp$wi.eBH < 0.05) # 0
sum(ttest.bac.r.CGC.pp$we.eBH < 0.05) # 0
sum(ttest.bac.r.MCC.pp$wi.eBH < 0.05) # 0
sum(ttest.bac.r.MCC.pp$we.eBH < 0.05) # 10
sum(ttest.bac.r.MCH.pp$wi.eBH < 0.05) # 0
sum(ttest.bac.r.MCH.pp$we.eBH < 0.05) # 2
sum(ttest.bac.r.MGH.pp$wi.eBH < 0.05) # 0
sum(ttest.bac.r.MGH.pp$we.eBH < 0.05) # 10
sum(ttest.bac.r.MGC.pp$wi.eBH < 0.05) # 34
sum(ttest.bac.r.MGC.pp$we.eBH < 0.05) # 47
# merge results with sample data and taxa table and write tables
tax.CCC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CCC)) %>% rownames_to_column("Taxa")
ttest.bac.r.CCC.pp.df <- ttest.bac.r.CCC.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.r.CCC.pp <- dplyr::left_join(tax.CCC,ttest.bac.r.CCC.pp.df,by = "Taxa")
sig.bac.r.CCC.pp<- results.bac.r.CCC.pp %>% filter(we.eBH < 0.05)

tax.CCH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CCH)) %>% rownames_to_column("Taxa")
ttest.bac.r.CCH.pp.df <- ttest.bac.r.CCH.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.r.CCH.pp <- dplyr::left_join(tax.CCH,ttest.bac.r.CCH.pp.df,by = "Taxa")
sig.bac.r.CCH.pp<- results.bac.r.CCH.pp %>% filter(we.eBH < 0.05)

tax.CGH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CGH)) %>% rownames_to_column("Taxa")
ttest.bac.r.CGH.pp.df <- ttest.bac.r.CGH.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.r.CGH.pp <- dplyr::left_join(tax.CGH,ttest.bac.r.CGH.pp.df,by = "Taxa")
sig.bac.r.CGH.pp<- results.bac.r.CGH.pp %>% filter(we.eBH < 0.05)

tax.CGC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CGC)) %>% rownames_to_column("Taxa")
ttest.bac.r.CGC.pp.df <- ttest.bac.r.CGC.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.r.CGC.pp <- dplyr::left_join(tax.CGC,ttest.bac.r.CGC.pp.df,by = "Taxa")
sig.bac.r.CGC.pp<- results.bac.r.CGC.pp %>% filter(we.eBH < 0.05)

tax.MCC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MCC)) %>% rownames_to_column("Taxa")
ttest.bac.r.MCC.pp.df <- ttest.bac.r.MCC.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.r.MCC.pp <- dplyr::left_join(tax.MCC,ttest.bac.r.MCC.pp.df,by = "Taxa")
sig.bac.r.MCC.pp<- results.bac.r.MCC.pp %>% filter(we.eBH < 0.05)

tax.MCH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MCH)) %>% rownames_to_column("Taxa")
ttest.bac.r.MCH.pp.df <- ttest.bac.r.MCH.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.r.MCH.pp <- dplyr::left_join(tax.MCH,ttest.bac.r.MCH.pp.df,by = "Taxa")
sig.bac.r.MCH.pp<- results.bac.r.MCH.pp %>% filter(we.eBH < 0.05)

tax.MGH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MGH)) %>% rownames_to_column("Taxa")
ttest.bac.r.MGH.pp.df <- ttest.bac.r.MGH.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.r.MGH.pp <- dplyr::left_join(tax.MGH,ttest.bac.r.MGH.pp.df,by = "Taxa")
sig.bac.r.MGH.pp<- results.bac.r.MGH.pp %>% filter(we.eBH < 0.05)

tax.MGC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MGC)) %>% rownames_to_column("Taxa")
ttest.bac.r.MGC.pp.df <- ttest.bac.r.MGC.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.r.MGC.pp <- dplyr::left_join(tax.MGC,ttest.bac.r.MGC.pp.df,by = "Taxa")
sig.bac.r.MGC.pp<- results.bac.r.MGC.pp %>% filter(we.eBH < 0.05)
# combind all of the results tables
sig.bac.r.pp.all <- dplyr::bind_rows(
  CCC = sig.bac.r.CCC.pp,
  CCH = sig.bac.r.CCH.pp,
  CGC = sig.bac.r.CGC.pp,
  CGH = sig.bac.r.CGH.pp,
  MCC = sig.bac.r.MCC.pp,
  MCH = sig.bac.r.MCH.pp,
  MGC = sig.bac.r.MGC.pp,
  MGH = sig.bac.r.MGH.pp,
  .id = "Treatment_code")
write.table(sig.bac.r.pp.all, file="Revisions.ISME.Comm/taxa_enrichments/Genus.sig.bac.r.PrePostTreat.csv", sep=",")
# summarize
sig.bac.r.pp.summary <- sig.bac.r.pp.all %>%group_by(Taxa) %>%summarise(
  n_treatments = n_distinct(Treatment_code),
  treatments = paste(sort(unique(Treatment_code)), collapse = ", "),
  Phylum = first(Phylum),Class = first(Class),Order = first(Order),Family = first(Family),Genus = first(Genus),.groups = "drop") %>%
  arrange(desc(n_treatments))
write.table(sig.bac.r.pp.summary, file="Revisions.ISME.Comm/taxa_enrichments/summary.Genus.sig.bac.r.PrePostTreat.csv", sep=",")

##### Calculate effect size ###
eff.bac.r.MCC <- aldex.effect(bac.r.MCC.pp.clr)
eff.bac.r.MGC <- aldex.effect(bac.r.MGC.pp.clr)
eff.bac.r.MGH <- aldex.effect(bac.r.MGH.pp.clr)
eff.bac.r.MCH <- aldex.effect(bac.r.MCH.pp.clr)
eff.bac.r.CCC <- aldex.effect(bac.r.CCC.pp.clr)
eff.bac.r.CGC <- aldex.effect(bac.r.CGC.pp.clr)
eff.bac.r.CGH <- aldex.effect(bac.r.CGH.pp.clr)
eff.bac.r.CCH <- aldex.effect(bac.r.CCH.pp.clr)
# get tax tables
tax.bac.r.MCC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MCC)) %>% rownames_to_column("Feature")
tax.bac.r.MGC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MGC)) %>% rownames_to_column("Feature")
tax.bac.r.MGH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MGH)) %>% rownames_to_column("Feature")
tax.bac.r.MCH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MCH)) %>% rownames_to_column("Feature")
tax.bac.r.CCC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CCC)) %>% rownames_to_column("Feature")
tax.bac.r.CGC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CGC)) %>% rownames_to_column("Feature")
tax.bac.r.CGH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CGH)) %>% rownames_to_column("Feature")
tax.bac.r.CCH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CCH)) %>% rownames_to_column("Feature")
# Build plotting table
volc.bac.r.MCC <- eff.bac.r.MCC %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.r.MCC.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.r.MCC, by="Feature")
volc.bac.r.MGC <- eff.bac.r.MGC %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.r.MGC.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.r.MGC, by="Feature")
volc.bac.r.MGH <- eff.bac.r.MGH %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.r.MGH.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.r.MGH, by="Feature")
volc.bac.r.MCH <- eff.bac.r.MCH %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.r.MCH.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.r.MCH, by="Feature")
volc.bac.r.CCC <- eff.bac.r.CCC %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.r.CCC.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.r.CCC, by="Feature")
volc.bac.r.CGC <- eff.bac.r.CGC %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.r.CGC.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.r.CGC, by="Feature")
volc.bac.r.CGH <- eff.bac.r.CGH %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.r.CGH.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.r.CGH, by="Feature")
volc.bac.r.CCH <- eff.bac.r.CCH %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.r.CCH.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.r.CCH, by="Feature")
##### Combine tables ###
res.list.bac.r <- list(
  CCC = volc.bac.r.CCC,
  CCH = volc.bac.r.CCH,
  CGC = volc.bac.r.CGC,
  CGH = volc.bac.r.CGH,
  MCC = volc.bac.r.MCC,
  MCH = volc.bac.r.MCH,
  MGC = volc.bac.r.MGC,
  MGH = volc.bac.r.MGH)
# Row-bind all tables and save
all.bac.r <- bind_rows(res.list.bac.r,.id = "Treatment_code")
write_xlsx(all.bac.r,path = "Revisions.ISME.Comm/taxa_enrichments/bacterial_root_PrePost_ALDEx2_full_results.xlsx")
########### dot plot #####################
sig.dot.df.bac.r <- imap_dfr(res.list.bac.r,~ .x %>%mutate(Treatment_code = .y) %>%filter(we.eBH < 0.05) %>%
                               dplyr::select(Feature, Treatment_code, effect, we.eBH,direction, Phylum, Class, Order, Family, Genus)) %>%
  group_by(Feature) %>%mutate(max_abs_effect = max(abs(effect), na.rm = TRUE)) %>%ungroup()
feature_order <- sig.dot.df.bac.r %>%mutate(treatment_priority = case_when(
  Treatment_code == "MGC" ~ 1L,
  Treatment_code == "MCC" ~ 2L,
  Treatment_code == "MCH" ~ 3L,
  Treatment_code == "CCC" ~ 4L,
  Treatment_code == "MGH" ~ 5L,
  TRUE                    ~ 6L)) %>%
  group_by(Feature) %>%arrange(treatment_priority,desc(effect), .by_group = TRUE) %>%slice(1) %>%ungroup() %>%arrange(treatment_priority, desc(effect))
feature_levels <- feature_order$Feature
sig.dot.df.bac.r <- sig.dot.df.bac.r %>%
  mutate(Feature = factor(Feature,levels = rev(feature_levels)))

ggplot(sig.dot.df.bac.r, aes(x = effect, y = Feature,fill = Treatment_code)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(shape = 21, color = "black",stroke=0.4,size = 3.5, alpha = 0.75, aes(fill=Treatment_code)) +
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  labs(x = "Effect size (CLR difference/dispersion)\nPre-treatment enriched \u2190    \u2192 Post-treatment enriched",y = "Genera",color = "Treatment",) +
  theme_bw() +
  xlim(-3.2,2.1)+
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = "gray95"),
        panel.grid.major.x = element_line(color = "gray95"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8,face = "italic"),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "right")
ggsave(filename= "genus.bac.r.treat.PrePost.enriched.pdf", device="pdf", units="mm", dpi=300, width=130, height=175, path="Revisions.ISME.Comm/plots/")
# "NA" genus names were manually updated to family or order name using Adobe Illustrator 

## Bac rhizo - testing and plotting - All 8 treatment groups ####
set.seed(4444)
# subseting treatments
bac.z.ps.raw.T0_T3.genus.CCC <- subset_samples(bac.z.ps.raw.T0_T3.genus, Treatment_code == "CCC") 
bac.z.ps.raw.T0_T3.genus.CCH <- subset_samples(bac.z.ps.raw.T0_T3.genus, Treatment_code == "CCH") 
bac.z.ps.raw.T0_T3.genus.CGH <- subset_samples(bac.z.ps.raw.T0_T3.genus, Treatment_code == "CGH") 
bac.z.ps.raw.T0_T3.genus.CGC <- subset_samples(bac.z.ps.raw.T0_T3.genus, Treatment_code == "CGC") 
bac.z.ps.raw.T0_T3.genus.MCC <- subset_samples(bac.z.ps.raw.T0_T3.genus, Treatment_code == "MCC") 
bac.z.ps.raw.T0_T3.genus.MCH <- subset_samples(bac.z.ps.raw.T0_T3.genus, Treatment_code == "MCH") 
bac.z.ps.raw.T0_T3.genus.MGH <- subset_samples(bac.z.ps.raw.T0_T3.genus, Treatment_code == "MGH") 
bac.z.ps.raw.T0_T3.genus.MGC <- subset_samples(bac.z.ps.raw.T0_T3.genus, Treatment_code == "MGC") 
# Set DF for test variable
bac.z.CCC.pp <- sample_data(bac.z.ps.raw.T0_T3.genus.CCC)$PrePostTreat %>% as.character() # 2 levels
bac.z.CCH.pp <- sample_data(bac.z.ps.raw.T0_T3.genus.CCH)$PrePostTreat %>% as.character() # 2 levels
bac.z.CGH.pp <- sample_data(bac.z.ps.raw.T0_T3.genus.CGH)$PrePostTreat %>% as.character() # 2 levels
bac.z.CGC.pp <- sample_data(bac.z.ps.raw.T0_T3.genus.CGC)$PrePostTreat %>% as.character() # 2 levels
bac.z.MCC.pp <- sample_data(bac.z.ps.raw.T0_T3.genus.MCC)$PrePostTreat %>% as.character() # 2 levels
bac.z.MCH.pp <- sample_data(bac.z.ps.raw.T0_T3.genus.MCH)$PrePostTreat %>% as.character() # 2 levels
bac.z.MGH.pp <- sample_data(bac.z.ps.raw.T0_T3.genus.MGH)$PrePostTreat %>% as.character() # 2 levels
bac.z.MGC.pp <- sample_data(bac.z.ps.raw.T0_T3.genus.MGC)$PrePostTreat %>% as.character() # 2 levels
# Extract Genus table
genus.bac.z.CCC <- as(otu_table(bac.z.ps.raw.T0_T3.genus.CCC, taxa_are_rows = TRUE), 'matrix') 
genus.bac.z.CCC <- t(genus.bac.z.CCC)
genus.bac.z.CCH <- as(otu_table(bac.z.ps.raw.T0_T3.genus.CCH, taxa_are_rows = TRUE), 'matrix') 
genus.bac.z.CCH <- t(genus.bac.z.CCH)
genus.bac.z.CGH <- as(otu_table(bac.z.ps.raw.T0_T3.genus.CGH, taxa_are_rows = TRUE), 'matrix') 
genus.bac.z.CGH <- t(genus.bac.z.CGH)
genus.bac.z.CGC <- as(otu_table(bac.z.ps.raw.T0_T3.genus.CGC, taxa_are_rows = TRUE), 'matrix') 
genus.bac.z.CGC <- t(genus.bac.z.CGC)
genus.bac.z.MCC <- as(otu_table(bac.z.ps.raw.T0_T3.genus.MCC, taxa_are_rows = TRUE), 'matrix') 
genus.bac.z.MCC <- t(genus.bac.z.MCC)
genus.bac.z.MCH <- as(otu_table(bac.z.ps.raw.T0_T3.genus.MCH, taxa_are_rows = TRUE), 'matrix') 
genus.bac.z.MCH <- t(genus.bac.z.MCH)
genus.bac.z.MGH <- as(otu_table(bac.z.ps.raw.T0_T3.genus.MGH, taxa_are_rows = TRUE), 'matrix') 
genus.bac.z.MGH <- t(genus.bac.z.MGH)
genus.bac.z.MGC <- as(otu_table(bac.z.ps.raw.T0_T3.genus.MGC, taxa_are_rows = TRUE), 'matrix') 
genus.bac.z.MGC <- t(genus.bac.z.MGC)
#### 2-level factors aldex CLR/MC sampling ##
bac.z.CCC.pp.clr <- aldex.clr(genus.bac.z.CCC, conds=bac.z.CCC.pp, mc.samples=200) #
bac.z.CCH.pp.clr <- aldex.clr(genus.bac.z.CCH, conds=bac.z.CCH.pp, mc.samples=200) # 
bac.z.CGH.pp.clr <- aldex.clr(genus.bac.z.CGH, conds=bac.z.CGH.pp, mc.samples=200) # 
bac.z.CGC.pp.clr <- aldex.clr(genus.bac.z.CGC, conds=bac.z.CGC.pp, mc.samples=200) # 
bac.z.MCC.pp.clr <- aldex.clr(genus.bac.z.MCC, conds=bac.z.MCC.pp, mc.samples=200) #
bac.z.MCH.pp.clr <- aldex.clr(genus.bac.z.MCH, conds=bac.z.MCH.pp, mc.samples=200) # 
bac.z.MGH.pp.clr <- aldex.clr(genus.bac.z.MGH, conds=bac.z.MGH.pp, mc.samples=200) # 
bac.z.MGC.pp.clr <- aldex.clr(genus.bac.z.MGC, conds=bac.z.MGC.pp, mc.samples=200) # 
#### 2 level factor statistical test ##
ttest.bac.z.CCC.pp <- aldex.ttest(bac.z.CCC.pp.clr)
ttest.bac.z.CCH.pp <- aldex.ttest(bac.z.CCH.pp.clr)
ttest.bac.z.CGH.pp <- aldex.ttest(bac.z.CGH.pp.clr)
ttest.bac.z.CGC.pp <- aldex.ttest(bac.z.CGC.pp.clr)
ttest.bac.z.MCC.pp <- aldex.ttest(bac.z.MCC.pp.clr)
ttest.bac.z.MCH.pp <- aldex.ttest(bac.z.MCH.pp.clr)
ttest.bac.z.MGH.pp <- aldex.ttest(bac.z.MGH.pp.clr)
ttest.bac.z.MGC.pp <- aldex.ttest(bac.z.MGC.pp.clr)
#### check for significant genera; BH = adjusted p-value ##
sum(ttest.bac.z.CCC.pp$wi.eBH < 0.05) # Wilcoxon test # 1
sum(ttest.bac.z.CCC.pp$we.eBH < 0.05) # Welch’s t test # 0
sum(ttest.bac.z.CCH.pp$wi.eBH < 0.05) # 0
sum(ttest.bac.z.CCH.pp$we.eBH < 0.05) # 0
sum(ttest.bac.z.CGH.pp$wi.eBH < 0.05) # 0
sum(ttest.bac.z.CGH.pp$we.eBH < 0.05) # 0
sum(ttest.bac.z.CGC.pp$wi.eBH < 0.05) # 2
sum(ttest.bac.z.CGC.pp$we.eBH < 0.05) # 0
sum(ttest.bac.z.MCC.pp$wi.eBH < 0.05) # 3
sum(ttest.bac.z.MCC.pp$we.eBH < 0.05) # 1
sum(ttest.bac.z.MCH.pp$wi.eBH < 0.05) # 1
sum(ttest.bac.z.MCH.pp$we.eBH < 0.05) # 0
sum(ttest.bac.z.MGH.pp$wi.eBH < 0.05) # 5
sum(ttest.bac.z.MGH.pp$we.eBH < 0.05) # 10
sum(ttest.bac.z.MGC.pp$wi.eBH < 0.05) # 54
sum(ttest.bac.z.MGC.pp$we.eBH < 0.05) # 54
# merge results with sample data and taxa table and write tables
tax.CCC <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.CCC)) %>% rownames_to_column("Taxa")
ttest.bac.z.CCC.pp.df <- ttest.bac.z.CCC.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.z.CCC.pp <- dplyr::left_join(tax.CCC,ttest.bac.z.CCC.pp.df,by = "Taxa")
sig.bac.z.CCC.pp<- results.bac.z.CCC.pp %>% filter(we.eBH < 0.05)

tax.CCH <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.CCH)) %>% rownames_to_column("Taxa")
ttest.bac.z.CCH.pp.df <- ttest.bac.z.CCH.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.z.CCH.pp <- dplyr::left_join(tax.CCH,ttest.bac.z.CCH.pp.df,by = "Taxa")
sig.bac.z.CCH.pp<- results.bac.z.CCH.pp %>% filter(we.eBH < 0.05)

tax.CGH <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.CGH)) %>% rownames_to_column("Taxa")
ttest.bac.z.CGH.pp.df <- ttest.bac.z.CGH.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.z.CGH.pp <- dplyr::left_join(tax.CGH,ttest.bac.z.CGH.pp.df,by = "Taxa")
sig.bac.z.CGH.pp<- results.bac.z.CGH.pp %>% filter(we.eBH < 0.05)

tax.CGC <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.CGC)) %>% rownames_to_column("Taxa")
ttest.bac.z.CGC.pp.df <- ttest.bac.z.CGC.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.z.CGC.pp <- dplyr::left_join(tax.CGC,ttest.bac.z.CGC.pp.df,by = "Taxa")
sig.bac.z.CGC.pp<- results.bac.z.CGC.pp %>% filter(we.eBH < 0.05)

tax.MCC <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.MCC)) %>% rownames_to_column("Taxa")
ttest.bac.z.MCC.pp.df <- ttest.bac.z.MCC.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.z.MCC.pp <- dplyr::left_join(tax.MCC,ttest.bac.z.MCC.pp.df,by = "Taxa")
sig.bac.z.MCC.pp<- results.bac.z.MCC.pp %>% filter(we.eBH < 0.05)

tax.MCH <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.MCH)) %>% rownames_to_column("Taxa")
ttest.bac.z.MCH.pp.df <- ttest.bac.z.MCH.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.z.MCH.pp <- dplyr::left_join(tax.MCH,ttest.bac.z.MCH.pp.df,by = "Taxa")
sig.bac.z.MCH.pp<- results.bac.z.MCH.pp %>% filter(we.eBH < 0.05)

tax.MGH <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.MGH)) %>% rownames_to_column("Taxa")
ttest.bac.z.MGH.pp.df <- ttest.bac.z.MGH.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.z.MGH.pp <- dplyr::left_join(tax.MGH,ttest.bac.z.MGH.pp.df,by = "Taxa")
sig.bac.z.MGH.pp<- results.bac.z.MGH.pp %>% filter(we.eBH < 0.05)

tax.MGC <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.MGC)) %>% rownames_to_column("Taxa")
ttest.bac.z.MGC.pp.df <- ttest.bac.z.MGC.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.bac.z.MGC.pp <- dplyr::left_join(tax.MGC,ttest.bac.z.MGC.pp.df,by = "Taxa")
sig.bac.z.MGC.pp<- results.bac.z.MGC.pp %>% filter(we.eBH < 0.05)
# combind all of the results tables
sig.bac.z.pp.all <- dplyr::bind_rows(
  CCC = sig.bac.z.CCC.pp,
  CCH = sig.bac.z.CCH.pp,
  CGC = sig.bac.z.CGC.pp,
  CGH = sig.bac.z.CGH.pp,
  MCC = sig.bac.z.MCC.pp,
  MCH = sig.bac.z.MCH.pp,
  MGC = sig.bac.z.MGC.pp,
  MGH = sig.bac.z.MGH.pp,
  .id = "Treatment_code")
write.table(sig.bac.z.pp.all, file="Revisions.ISME.Comm/taxa_enrichments/Genus.sig.bac.z.PrePostTreat.csv", sep=",")
# summarize
sig.bac.z.pp.summary <- sig.bac.z.pp.all %>%group_by(Taxa) %>%summarise(
  n_treatments = n_distinct(Treatment_code),
  treatments = paste(sort(unique(Treatment_code)), collapse = ", "),
  Phylum = first(Phylum),Class = first(Class),Order = first(Order),Family = first(Family),Genus = first(Genus),.groups = "drop") %>%
  arrange(desc(n_treatments))
write.table(sig.bac.z.pp.summary, file="Revisions.ISME.Comm/taxa_enrichments/summary.Genus.sig.bac.z.PrePostTreat.csv", sep=",")

# Calculate effect size #
eff.bac.z.MCC <- aldex.effect(bac.z.MCC.pp.clr)
eff.bac.z.MGC <- aldex.effect(bac.z.MGC.pp.clr)
eff.bac.z.MGH <- aldex.effect(bac.z.MGH.pp.clr)
eff.bac.z.MCH <- aldex.effect(bac.z.MCH.pp.clr)
eff.bac.z.CCC <- aldex.effect(bac.z.CCC.pp.clr)
eff.bac.z.CGC <- aldex.effect(bac.z.CGC.pp.clr)
eff.bac.z.CGH <- aldex.effect(bac.z.CGH.pp.clr)
eff.bac.z.CCH <- aldex.effect(bac.z.CCH.pp.clr)
# get tax tables
tax.bac.z.MCC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MCC)) %>% rownames_to_column("Feature")
tax.bac.z.MGC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MGC)) %>% rownames_to_column("Feature")
tax.bac.z.MGH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.MGH)) %>% rownames_to_column("Feature")
tax.bac.z.MCH <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.MCH)) %>% rownames_to_column("Feature")
tax.bac.z.CCC <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.CCC)) %>% rownames_to_column("Feature")
tax.bac.z.CGC <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CGC)) %>% rownames_to_column("Feature")
tax.bac.z.CGH <- data.frame(tax_table(bac.r.ps.raw.T0_T3.genus.CGH)) %>% rownames_to_column("Feature")
tax.bac.z.CCH <- data.frame(tax_table(bac.z.ps.raw.T0_T3.genus.CCH)) %>% rownames_to_column("Feature")
# Build plotting table
volc.bac.z.MCC <- eff.bac.z.MCC %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.z.MCC.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.z.MCC, by="Feature")
volc.bac.z.MGC <- eff.bac.z.MGC %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.z.MGC.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.z.MGC, by="Feature")
volc.bac.z.MGH <- eff.bac.z.MGH %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.z.MGH.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>%left_join(tax.bac.z.MGH, by="Feature")
volc.bac.z.MCH <- eff.bac.z.MCH %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.z.MCH.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.z.MCH, by="Feature")
volc.bac.z.CCC <- eff.bac.z.CCC %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.z.CCC.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.z.CCC, by="Feature")
volc.bac.z.CGC <- eff.bac.z.CGC %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.z.CGC.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.z.CGC, by="Feature")
volc.bac.z.CGH <- eff.bac.z.CGH %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.z.CGH.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>%left_join(tax.bac.z.CGH, by="Feature")
volc.bac.z.CCH <- eff.bac.z.CCH %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.bac.z.CCH.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.bac.z.CCH, by="Feature")
##### Combine tables ###
res.list.bac.z <- list(
  CCC = volc.bac.z.CCC,
  CCH = volc.bac.z.CCH,
  CGC = volc.bac.z.CGC,
  CGH = volc.bac.z.CGH,
  MCC = volc.bac.z.MCC,
  MCH = volc.bac.z.MCH,
  MGC = volc.bac.z.MGC,
  MGH = volc.bac.z.MGH)
# Row-bind all tables and save
all.bac.z <- bind_rows(res.list.bac.z,.id = "Treatment_code")
write_xlsx(all.bac.z,path = "Revisions.ISME.Comm/taxa_enrichments/bacterial_rhizosphere_PrePost_ALDEx2_full_results.xlsx")
########### dot plot #####################
sig.dot.df.bac.z <- imap_dfr(res.list.bac.z,~ .x %>%mutate(Treatment_code = .y) %>%filter(we.eBH < 0.05) %>%
                               dplyr::select(Feature, Treatment_code, effect, we.eBH,direction, Phylum, Class, Order, Family, Genus)) %>%
  group_by(Feature) %>%mutate(max_abs_effect = max(abs(effect), na.rm = TRUE)) %>% ungroup()
feature_order <- sig.dot.df.bac.z %>%mutate(treatment_priority = case_when(
  Treatment_code == "MGC" ~ 1L,
  Treatment_code == "MGH" ~ 2L,
  Treatment_code == "MCC" ~ 3L,
  TRUE                    ~ 4L)) %>%
  group_by(Feature) %>%arrange(treatment_priority, desc(effect), .by_group = TRUE) %>%slice(1) %>%ungroup() %>%arrange(treatment_priority, desc(effect))
feature_levels <- feature_order$Feature
sig.dot.df.bac.z <- sig.dot.df.bac.z %>%
  mutate(Feature = factor(Feature,levels = rev(feature_levels)))

ggplot(sig.dot.df.bac.z, aes(x = effect, y = Feature,fill = Treatment_code)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(shape = 21, color = "black",stroke=0.4,size = 3.5, alpha = 0.75, aes(fill=Treatment_code)) +
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  labs(x = "Effect size (CLR difference/dispersion)\nPre-treatment enriched \u2190    \u2192 Post-treatment enriched",y = "Genera",color = "Treatment",) +
  theme_bw() +
  xlim(-3.2,2.1)+
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = "gray95"),
        panel.grid.major.x = element_line(color = "gray95"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8,face = "italic"),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "right")
ggsave(filename= "genus.bac.z.treat.PrePost.enriched.pdf", device="pdf", units="mm", dpi=300, width=133, height=175, path="Revisions.ISME.Comm/plots/")
# "NA" genus names were manually updated to family or order name using Adobe Illustrator 

## fun root - testing and plotting - Mulch/no mulch ####
# The fungal communities had no significant genera enriched pre vs. post treatment when each of the 8 treatment groups were test
# The beta diversity analysis showed that only mulch treatment significantly explained community composition shifts
# Therefore, we looked specifically at genera enriched pre vs. post mulch/no mulch treatments
set.seed(4444)
# subset treatments
fun.r.ps.raw.T0_T3.genus.c <- subset_samples(fun.r.ps.raw.T0_T3.genus, Mulch == "no")
fun.r.ps.raw.T0_T3.genus.m <- subset_samples(fun.r.ps.raw.T0_T3.genus, Mulch == "mulch")
# Set DF for test variable
fun.r.c.pp <- sample_data(fun.r.ps.raw.T0_T3.genus.c)$PrePostTreat %>% as.character() # 2 levels
fun.r.m.pp <- sample_data(fun.r.ps.raw.T0_T3.genus.m)$PrePostTreat %>% as.character() # 2 levels
# Extract genus table
genus.fun.r.c <- as(otu_table(fun.r.ps.raw.T0_T3.genus.c, taxa_are_rows = TRUE), 'matrix') 
genus.fun.r.c <- t(genus.fun.r.c)
genus.fun.r.m <- as(otu_table(fun.r.ps.raw.T0_T3.genus.m, taxa_are_rows = TRUE), 'matrix') 
genus.fun.r.m <- t(genus.fun.r.m)
#### 2-level factors aldex CLR/MC sampling ##
fun.r.c.pp.clr <- aldex.clr(genus.fun.r.c, conds=fun.r.c.pp, mc.samples=200) #
fun.r.m.pp.clr <- aldex.clr(genus.fun.r.m, conds=fun.r.m.pp, mc.samples=200) #
#### 2 level factor statistical test ##
ttest.fun.r.c.pp <- aldex.ttest(fun.r.c.pp.clr)
ttest.fun.r.m.pp <- aldex.ttest(fun.r.m.pp.clr)
#### check for significant genera; BH = adjusted p-value ##
sum(ttest.fun.r.c.pp$wi.eBH < 0.1) # Wilcoxon test # 0
sum(ttest.fun.r.c.pp$we.eBH < 0.1) # Welch’s t test # 0
sum(ttest.fun.r.m.pp$wi.eBH < 0.1) # Wilcoxon test # 1
sum(ttest.fun.r.m.pp$we.eBH < 0.1) # Welch’s t test # 1
# merge results with sample data and taxa table and write tables
tax.c <- data.frame(tax_table(fun.r.ps.raw.T0_T3.genus.c)) %>% rownames_to_column("Taxa")
ttest.fun.r.c.pp.df <- ttest.fun.r.c.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.fun.r.c.pp <- dplyr::left_join(tax.c,ttest.fun.r.c.pp.df,by = "Taxa")
sig.fun.r.c.pp<- results.fun.r.c.pp %>% filter(we.eBH < 0.1)

tax.m <- data.frame(tax_table(fun.r.ps.raw.T0_T3.genus.m)) %>% rownames_to_column("Taxa")
ttest.fun.r.m.pp.df <- ttest.fun.r.m.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.fun.r.m.pp <- dplyr::left_join(tax.m,ttest.fun.r.m.pp.df,by = "Taxa")
sig.fun.r.m.pp<- results.fun.r.m.pp %>% filter(we.eBH < 0.1)
# combind all of the results tables
sig.fun.r.pp.all <- dplyr::bind_rows(
  no_mulch = sig.fun.r.c.pp,
  mulch = sig.fun.r.m.pp,
  .id = "Mulch")
write.table(sig.fun.r.pp.all, file="Revisions.ISME.Comm/taxa_enrichments/Genus.sig.fun.r.PrePostTreat.mulch.csv", sep=",")
# summarize
sig.fun.r.pp.summary <- sig.fun.r.pp.all %>%group_by(Taxa) %>%summarise(
  n_treatments = n_distinct(Mulch),
  treatments = paste(sort(unique(Mulch)), collapse = ", "),
  Phylum = first(Phylum),Class = first(Class),Order = first(Order),Family = first(Family),Genus = first(Genus),.groups = "drop") %>%
  arrange(desc(n_treatments))
write.table(sig.fun.r.pp.summary, file="Revisions.ISME.Comm/taxa_enrichments/summary.Genus.sig.fun.r.PrePostTreat.mulch.csv", sep=",")

# Calculate effect size ###
eff.fun.r.c <- aldex.effect(fun.r.c.pp.clr)
eff.fun.r.m <- aldex.effect(fun.r.m.pp.clr)
# get taxa tables
tax.fun.r.c <- data.frame(tax_table(fun.r.ps.raw.T0_T3.genus.c)) %>% rownames_to_column("Feature")
tax.fun.r.m <- data.frame(tax_table(fun.r.ps.raw.T0_T3.genus.m)) %>% rownames_to_column("Feature")
# build plotting table
volc.fun.r.c <- eff.fun.r.c %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.fun.r.c.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.fun.r.c, by="Feature")
volc.fun.r.m <- eff.fun.r.m %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.fun.r.m.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.fun.r.m, by="Feature")
# Combine tables #
res.list.fun.r <- list(
  no_mulch = volc.fun.r.c,
  mulch = volc.fun.r.m)
# Row-bind all tables and save
all.fun.r <- bind_rows(res.list.fun.r,.id = "Treatment_code")
write_xlsx(all.fun.r,path = "Revisions.ISME.Comm/taxa_enrichments/fungal_root_PrePost_ALDEx2_full_results.xlsx")
################ dot plot ##############
sig.dot.df.fun.r <- imap_dfr(res.list.fun.r,~ .x %>%mutate(Treatment_code = .y) %>%filter(we.eBH < 0.1) %>%
                               dplyr::select(Feature, Treatment_code, effect, we.eBH,direction, Phylum, Class, Order, Family, Genus)) %>%
  group_by(Feature) %>%mutate(max_abs_effect = max(abs(effect), na.rm = TRUE)) %>%
  ungroup() 

ggplot(sig.dot.df.fun.r, aes(x = effect, y = Feature,fill = Treatment_code)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(shape = 21, color = "black",stroke=0.4,size = 3.5, alpha = 0.75, aes(fill=Treatment_code)) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no_mulch"))+
  labs(x = "Effect size (CLR difference/dispersion)\nPre-treatment enriched \u2190    \u2192 Post-treatment enriched",y = "Genera",fill = "Treatment",) +
  theme_bw() +
  xlim(-3.2,2.1)+
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = "gray95"),
        panel.grid.major.x = element_line(color = "gray95"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8,face = "italic"),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "right")
ggsave(filename= "genus.fun.r.mulch.PrePost.enriched.pdf", device="pdf", units="mm", dpi=300, width=130, height=25, path="Revisions.ISME.Comm/plots/")
# "NA" genus names were manually updated to family or order name using Adobe Illustrator 

## fun rhizo - testing and plotting - Mulch/no mulch ####
# The fungal communities had no significant genera enriched pre vs. post treatment when each of the 8 treatment groups were test
# The beta diversity analysis showed that only mulch treatment significantly explained community composition shifts
# Therefore, we looked specifically at genera enriched pre vs. post mulch/no mulch treatments
set.seed(4444)
# subset treatments
fun.z.ps.raw.T0_T3.genus.c <- subset_samples(fun.z.ps.raw.T0_T3.genus, Mulch == "no")
fun.z.ps.raw.T0_T3.genus.m <- subset_samples(fun.z.ps.raw.T0_T3.genus, Mulch == "mulch")
# Set DF for test variable
fun.z.c.pp <- sample_data(fun.z.ps.raw.T0_T3.genus.c)$PrePostTreat %>% as.character() # 2 levels
fun.z.m.pp <- sample_data(fun.z.ps.raw.T0_T3.genus.m)$PrePostTreat %>% as.character() # 2 levels
# Extract genus table
genus.fun.z.c <- as(otu_table(fun.z.ps.raw.T0_T3.genus.c, taxa_are_rows = TRUE), 'matrix') 
genus.fun.z.c <- t(genus.fun.z.c)
genus.fun.z.m <- as(otu_table(fun.z.ps.raw.T0_T3.genus.m, taxa_are_rows = TRUE), 'matrix') 
genus.fun.z.m <- t(genus.fun.z.m)
#### 2-level factors aldex CLR/MC sampling ##
fun.z.c.pp.clr <- aldex.clr(genus.fun.z.c, conds=fun.z.c.pp, mc.samples=200) #
fun.z.m.pp.clr <- aldex.clr(genus.fun.z.m, conds=fun.z.m.pp, mc.samples=200) #
#### 2 level factor statistical test ##
ttest.fun.z.c.pp <- aldex.ttest(fun.z.c.pp.clr)
ttest.fun.z.m.pp <- aldex.ttest(fun.z.m.pp.clr)
#### check for significant genera; BH = adjusted p-value ##
sum(ttest.fun.z.c.pp$wi.eBH < 0.1) # Wilcoxon test # 0
sum(ttest.fun.z.c.pp$we.eBH < 0.1) # Welch’s t test # 1
sum(ttest.fun.z.m.pp$wi.eBH < 0.1) # Wilcoxon test # 17
sum(ttest.fun.z.m.pp$we.eBH < 0.1) # Welch’s t test # 15
# merge results with sample data and taxa table and write tables
tax.c <- data.frame(tax_table(fun.z.ps.raw.T0_T3.genus.c)) %>% rownames_to_column("Taxa")
ttest.fun.z.c.pp.df <- ttest.fun.z.c.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.fun.z.c.pp <- dplyr::left_join(tax.c,ttest.fun.z.c.pp.df,by = "Taxa")
sig.fun.z.c.pp<- results.fun.z.c.pp %>% filter(we.eBH < 0.1)

tax.m <- data.frame(tax_table(fun.z.ps.raw.T0_T3.genus.m)) %>% rownames_to_column("Taxa")
ttest.fun.z.m.pp.df <- ttest.fun.z.m.pp %>% data.frame() %>% rownames_to_column("Taxa")
results.fun.z.m.pp <- dplyr::left_join(tax.m,ttest.fun.z.m.pp.df,by = "Taxa")
sig.fun.z.m.pp<- results.fun.z.m.pp %>% filter(we.eBH < 0.1)
# combind all of the results tables
sig.fun.z.pp.all <- dplyr::bind_rows(
  no_mulch = sig.fun.z.c.pp,
  mulch = sig.fun.z.m.pp,
  .id = "Mulch")
write.table(sig.fun.z.pp.all, file="Revisions.ISME.Comm/taxa_enrichments/Genus.sig.fun.z.PrePostTreat.mulch.csv", sep=",")
# summarize
sig.fun.z.pp.summary <- sig.fun.z.pp.all %>%group_by(Taxa) %>%summarise(
  n_treatments = n_distinct(Mulch),
  treatments = paste(sort(unique(Mulch)), collapse = ", "),
  Phylum = first(Phylum),Class = first(Class),Order = first(Order),Family = first(Family),Genus = first(Genus),.groups = "drop") %>%
  arrange(desc(n_treatments))
write.table(sig.fun.z.pp.summary, file="Revisions.ISME.Comm/taxa_enrichments/summary.Genus.sig.fun.z.PrePostTreat.mulch.csv", sep=",")

# Calculate effect size #
eff.fun.z.c <- aldex.effect(fun.z.c.pp.clr)
eff.fun.z.m <- aldex.effect(fun.z.m.pp.clr)
# get taxa tables #
tax.fun.z.c <- data.frame(tax_table(fun.z.ps.raw.T0_T3.genus.c)) %>% rownames_to_column("Feature")
tax.fun.z.m <- data.frame(tax_table(fun.z.ps.raw.T0_T3.genus.m)) %>% rownames_to_column("Feature")
# build plotting table #
volc.fun.z.c <- eff.fun.z.c %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.fun.z.c.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.fun.z.c, by="Feature")
volc.fun.z.m <- eff.fun.z.m %>% data.frame() %>% rownames_to_column("Feature") %>% left_join(ttest.fun.z.m.pp %>% data.frame() %>% rownames_to_column("Feature"),by = "Feature") %>%
  mutate(neglog10_q = -log10(we.eBH),sig = we.eBH < 0.05, direction = case_when(effect > 0 ~ "post-treatment",effect < 0 ~ "pre-treatment",TRUE ~ "no change"),
         color_group = case_when(sig & effect > 0 ~ "post-treatment",sig & effect < 0 ~ "pre-treatment",TRUE ~ "ns")) %>% left_join(tax.fun.z.m, by="Feature")
# Combine tables #
res.list.fun.z <- list(
  no_mulch = volc.fun.z.c,
  mulch = volc.fun.z.m)
# Row-bind all tables and save
all.fun.z <- bind_rows(res.list.fun.z,.id = "Treatment_code")
write_xlsx(all.fun.z,path = "Revisions.ISME.Comm/taxa_enrichments/fungal_rhizosphere_PrePost_ALDEx2_full_results.xlsx")
########### dot plot #####################
sig.dot.df.fun.z <- imap_dfr(res.list.fun.z,~ .x %>%mutate(Treatment_code = .y) %>%filter(we.eBH < 0.1) %>%
                               dplyr::select(Feature, Treatment_code, effect, we.eBH,direction, Phylum, Class, Order, Family, Genus)) %>%
  group_by(Feature) %>%mutate(max_abs_effect = max(abs(effect), na.rm = TRUE)) %>%ungroup() %>% mutate(Feature = forcats::fct_reorder(Feature,effect,.fun = max,.desc = FALSE))

ggplot(sig.dot.df.fun.z, aes(x = effect, y = Feature,fill = Treatment_code)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(shape = 21, color = "black",stroke=0.4,size = 3.5, alpha = 0.75, aes(fill=Treatment_code)) +
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no_mulch"))+
  labs(x = "Effect size (CLR difference/dispersion)\nPre-treatment enriched \u2190    \u2192 Post-treatment enriched",y = "Genera",fill = "Treatment",) +
  theme_bw() +
  xlim(-3.2,2.1)+
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = "gray95"),
        panel.grid.major.x = element_line(color = "gray95"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8,face = "italic"),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "right")
ggsave(filename= "genus.fun.z.mulch.PrePost.enriched.pdf", device="pdf", units="mm", dpi=300, width=130, height=65, path="Revisions.ISME.Comm/plots/")
# "NA" genus names were manually updated to family or order name using Adobe Illustrator 


