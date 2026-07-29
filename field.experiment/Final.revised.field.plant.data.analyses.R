# Script created November 20, 2025 by Riley Jones, PhD Student
# Last updated July 15, 2026 by Nichole Ginnan, Assistant Project Scientist
# Analysis to understand the impacts of mulch, glyphosate, and humic acid, on field tree traits and weed growth
# ECDRE Project, Lead PI: Caroline Roper, funded by the USDA
# Lindcove Research and Extension Center, Exeter, California | 91C plot experiment
# Includes Leaf nutrients LICOR, weed biomass, and fruit packline data
###############################################################################
# Load Libraries ####
# analyses completed in R v 4.5.0
library(multcomp) # v 1.4-28
library(ggplot2) # v 3.5.2
library(emmeans) # v 1.11.2
library(lme4) # v 1.1-37
library(lmerTest) # v 3.1-3
library(tidyverse) # v 2.0.0
library(readxl) # v 1.4.5
library(MuMIn) # v 1.48.19

#### Load data ############
setwd("/Users/nicholeginnan/Documents/UCR- Current/CA.citrus.paper/") # set working directory
licor_cleaned<- read_excel("Revisions.ISME.Comm/Sup.Table.S3.FieldPlantData.xlsx", sheet = "Assimilation_and_transpiration")
final_packline<- read_excel("Revisions.ISME.Comm/Sup.Table.S3.FieldPlantData.xlsx", sheet = "Fruit_weight")
weed_biomass<- read_excel("Revisions.ISME.Comm/Sup.Table.S3.FieldPlantData.xlsx", sheet = "Weed_biomass")
nutrients<- read_excel("Revisions.ISME.Comm/Sup.Table.S3.FieldPlantData.xlsx", sheet = "Leaf_nutrients")
set.seed(4444)

######################################################
######################################################
#### Leaf nutrients T0.5-T3 - Temporal analysis ########
######################################################
######################################################
# Plot raw data T0.5-T3, switch out y= "nutrient"
ggplot(nutrients, aes(x = Timepoint, y = Lf_B_ppm, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) 

#### Leaf nutrients of interest with visual differences and similarities across treatments ###
# diff.         # similar
# `Lf_K_%`        `Lf_Cu_ppm`
# `Lf_Ca_%`.      `Lf_Zn_ppm`
# `Lf_Fe_ppm`
# `Lf_Na_%`
# Lf_N_%
# `Lf_P_%`
# Lf_B_ppm

# Calculate the relative nutrients for T1, T2, and T3 or percent change from T0.5
T0.5 <- nutrients %>%filter(Timepoint == "T0.5") %>%dplyr::select(Treat_rep,`Lf_K_%`,`Lf_Ca_%`,`Lf_Na_%`,`Lf_N_%`,`Lf_P_%`,Lf_Fe_ppm,Lf_B_ppm,Lf_Cu_ppm,Lf_Zn_ppm) %>% 
  rename(`Lf_K_%_T0.5`=`Lf_K_%`,`Lf_Ca_%_T0.5`=`Lf_Ca_%`,`Lf_Na_%_T0.5`=`Lf_Na_%`,`Lf_N_%_T0.5`=`Lf_N_%`,
         `Lf_P_%_T0.5`=`Lf_P_%`,`Lf_Fe_ppm_T0.5`=Lf_Fe_ppm,`Lf_B_ppm_T0.5`=Lf_B_ppm,`Lf_Cu_ppm_T0.5`=Lf_Cu_ppm,`Lf_Zn_ppm_T0.5`=Lf_Zn_ppm)
nutrients.relative <- nutrients %>%left_join(T0.5, by = "Treat_rep")
nutrients.relative <- nutrients.relative %>%mutate(`rel.Lf_K_%`=(`Lf_K_%`-`Lf_K_%_T0.5`)/`Lf_K_%_T0.5`,
                                                   `rel.Lf_Ca_%`=(`Lf_Ca_%`-`Lf_Ca_%_T0.5`)/`Lf_Ca_%_T0.5`,
                                                   `rel.Lf_Na_%`=(`Lf_Na_%`-`Lf_Na_%_T0.5`)/`Lf_Na_%_T0.5`,
                                                   `rel.Lf_N_%`=(`Lf_N_%`-`Lf_N_%_T0.5`)/`Lf_N_%_T0.5`,
                                                   `rel.Lf_P_%`=(`Lf_P_%`-`Lf_P_%_T0.5`)/`Lf_P_%_T0.5`,
                                                   `rel.Lf_Fe_ppm`=(`Lf_Fe_ppm`-`Lf_Fe_ppm_T0.5`)/`Lf_Fe_ppm_T0.5`,
                                                   `rel.Lf_B_ppm`=(`Lf_B_ppm`-`Lf_B_ppm_T0.5`)/`Lf_B_ppm_T0.5`,
                                                   `rel.Lf_Cu_ppm`=(`Lf_Cu_ppm`-`Lf_Cu_ppm_T0.5`)/`Lf_Cu_ppm_T0.5`,
                                                   `rel.Lf_Zn_ppm`=(`Lf_Zn_ppm`-`Lf_Zn_ppm_T0.5`)/`Lf_Zn_ppm_T0.5`)
nutrients.relative <- nutrients.relative %>%filter(Timepoint != "T0.5") # remove T0.5
#### leaf nitrogen ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = `rel.Lf_N_%`, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)
# Mixed Effects Model
m0<-lmer(`rel.Lf_N_%`~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(`rel.Lf_N_%`~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(`rel.Lf_N_%`~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(`rel.Lf_N_%`~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(`rel.Lf_N_%`~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(`rel.Lf_N_%`~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(`rel.Lf_N_%`~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(`rel.Lf_N_%`~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(`rel.Lf_N_%`~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(`rel.Lf_N_%`~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(`rel.Lf_N_%`~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(`rel.Lf_N_%`~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(`rel.Lf_N_%`~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(`rel.Lf_N_%`~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # good
plot(resid(m1)~fitted(m1)) # good
plot(resid(m2)~fitted(m2)) # good
plot(resid(m3)~fitted(m3)) # not good
plot(resid(m4)~fitted(m4)) # not good
plot(resid(m5)~fitted(m5)) # ok
plot(resid(m6)~fitted(m6)) # ok
plot(resid(m7)~fitted(m7)) # not good
plot(resid(m8)~fitted(m8)) # not good
plot(resid(m9)~fitted(m9)) # not good
plot(resid(m10)~fitted(m10)) # not good
plot(resid(m11)~fitted(m11)) # not good
plot(resid(m12)~fitted(m12)) # not good
plot(resid(m13)~fitted(m13)) # not good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m5,m6) # including only the models with good to ok fit

#.      (Int) Gly Hmc Mlc Tmp Gly:Hmc Gly:Mlc Gly:Tmp Hmc:Mlc Hmc:Tmp Mlc:Tmp Gly:Hmc:Mlc Gly:Hmc:Tmp Gly:Mlc:Tmp Hmc:Mlc:Tmp Gly:Hmc:Mlc:Tmp df logLik   AICc delta weight
#m1 -0.13190   +   +   +   +               +       +                       +                                   +                             15 98.590 -158.6  0.00  0.935
#m6 -0.05624   +   +   +   +                                               +                                                                 10 88.445 -153.3  5.33  0.065
#m5 -0.07078   +   +   +   +                       +               +       +                                                                 14 78.165 -121.0 37.65  0.000
#m2 -0.05010   +   +   +   +                               +       +       +                                               +                 15 76.698 -114.8 43.78  0.000
#m0 -0.13250   +   +   +   +       +       +       +       +       +       +           +           +           +           +               + 26 80.935  -78.7 79.94  0.000
AIC(m0,m1,m2,m5,m6)
#.   df      AIC
#m0 26 -109.8700
#m1 15 -167.1802
#m2 15 -123.3960
#m5 14 -128.3305
#m6 10 -156.8890
# Test
anova(m1, type = "III") 
#.                            Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
#Mulch                      0.01822 0.01822     1    19  34.6687 1.140e-05 ***
#Glyphosate                 0.00347 0.00347     1    19   6.6079  0.018722 *  
#Timepoint                  0.80104 0.40052     2    40 761.9612 < 2.2e-16 ***
#Humic                      0.00687 0.00687     1    19  13.0614  0.001848 ** 
#Mulch:Glyphosate           0.00487 0.00487     1    19   9.2598  0.006692 ** 
#Mulch:Timepoint            0.17615 0.08807     2    40 167.5563 < 2.2e-16 ***
#Glyphosate:Timepoint       0.00427 0.00214     2    40   4.0623  0.024766 *  
#Mulch:Glyphosate:Timepoint 0.02915 0.01457     2    40  27.7266 2.789e-08 ***
emm<-emmeans(m1, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response")

pairs(emm)
emmip(m1,Mulch ~ Glyphosate)
emmip(m1,Glyphosate ~ Mulch | Timepoint|Humic)
emm_3way <- emmeans(m1, ~ Mulch * Glyphosate * Timepoint)
joint_tests(emm_3way, by = "Timepoint")
pairs_analysis <- pairs(emm_3way, by = c("Timepoint", "Mulch"))
summary(pairs_analysis, adjust = "tukey")
pairs(emm_3way, by = c("Timepoint", "Glyphosate"))

cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
      Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
      Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
      Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
      Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
      Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
      Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
      Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
      Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.N <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=`rel.Lf_N_%`),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Total Nitrogen (%)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_N_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")
#ggsave(filename= "EMM.rel.Lf_N_Timepoint.line_no_fct.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

#### leaf Potassium ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = `rel.Lf_K_%`, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)

# Mixed Effects Model
m0<-lmer(`rel.Lf_K_%`~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(`rel.Lf_K_%`~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(`rel.Lf_K_%`~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(`rel.Lf_K_%`~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(`rel.Lf_K_%`~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(`rel.Lf_K_%`~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(`rel.Lf_K_%`~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(`rel.Lf_K_%`~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(`rel.Lf_K_%`~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(`rel.Lf_K_%`~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(`rel.Lf_K_%`~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(`rel.Lf_K_%`~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(`rel.Lf_K_%`~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(`rel.Lf_K_%`~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # ok
plot(resid(m1)~fitted(m1)) # ok
plot(resid(m2)~fitted(m2)) # ok
plot(resid(m3)~fitted(m3)) # ok 
plot(resid(m4)~fitted(m4)) # not good
plot(resid(m5)~fitted(m5)) # ok
plot(resid(m6)~fitted(m6)) # ok
plot(resid(m7)~fitted(m7)) # ok
plot(resid(m8)~fitted(m8)) # ok
plot(resid(m9)~fitted(m9)) # not good
plot(resid(m10)~fitted(m10)) # not good
plot(resid(m11)~fitted(m11)) # not good
plot(resid(m12)~fitted(m12)) # not good
plot(resid(m13)~fitted(m13)) # not good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3,m5,m6,m7,m8) # including only the models with good to ok fit

#.      Int) Gly Hmc Mlc Tmp Gly:Hmc Gly:Mlc Gly:Tmp Hmc:Mlc Hmc:Tmp Mlc:Tmp Gly:Hmc:Mlc Gly:Hmc:Tmp Gly:Mlc:Tmp Hmc:Mlc:Tmp Gly:Hmc:Mlc:Tmp df logLik   AICc delta weight
#m7 0.8045   +   +   +   +                       +                                                                                         10 -15.015 53.6  0.00  0.698
#m1 1.0150   +   +   +   +               +       +                       +                                   +                             15  -8.641 55.9  2.22  0.230
#m5 0.7818   +   +   +   +                       +               +       +                                                                 14 -12.330 60.0  6.39  0.029
#m8 0.5469   +   +   +   +                                       +                                                                         10 -18.494 60.6  6.96  0.022
#m3 0.7141   +   +   +   +       +               +               +                               +                                         15 -11.333 61.2  7.60  0.016
#m6 0.7408   +   +   +   +                                               +                                                                 10 -19.792 63.2  9.55  0.006
#m2 0.6357   +   +   +   +                               +       +       +                                               +                 15 -18.127 74.8 21.19  0.000
#m0 0.9600   +   +   +   +       +       +       +       +       +       +           +           +           +           +               + 26   0.223 82.8 29.12  0.000
AIC(m0,m1,m2,m3,m5,m6,m7,m8)
#.   df      AIC
#m0 26 51.55500
#m1 15 47.28188 <-- going with m1 
#m2 15 66.25331
#m3 15 52.66586
#m5 14 52.65985
#m6 10 59.58350
#m7 10 50.03090
#m8 10 56.98733
# Test
anova(m1, type = "III") 
#.                            Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
#Mulch                      0.2952 0.29520     1    19  8.9708 0.0074414 ** 
#Glyphosate                 0.0587 0.05869     1    19  1.7836 0.1974894    
#Timepoint                  4.5112 2.25558     2    40 68.5459 1.195e-13 ***
#Humic                      0.0014 0.00145     1    19  0.0439 0.8361857    
#Mulch:Glyphosate           0.0014 0.00136     1    19  0.0415 0.8408210    
#Mulch:Timepoint            0.2637 0.13183     2    40  4.0064 0.0259462 *  
#Glyphosate:Timepoint       0.8031 0.40155     2    40 12.2028 7.290e-05 ***
#Mulch:Glyphosate:Timepoint 0.6448 0.32238     2    40  9.7968 0.0003445 ***
emm<-emmeans(m1, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response")
pairs(emm)
#emmip(m1,Mulch ~ Timepoint | Glyphosate)
#emmip(m1,Mulch:Glyphosate:Humic ~ Timepoint)
emmip(m1,Mulch:Glyphosate ~ Timepoint|Humic)
#emmip(m1,Glyphosate ~ Timepoint | Mulch)
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
  Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
  Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.K <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=`rel.Lf_K_%`),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Potassium (%)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_K_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

#### leaf Calcium ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = `rel.Lf_Ca_%`, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)

# Mixed Effects Model
m0<-lmer(`rel.Lf_Ca_%`~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(`rel.Lf_Ca_%`~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(`rel.Lf_Ca_%`~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(`rel.Lf_Ca_%`~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(`rel.Lf_Ca_%`~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(`rel.Lf_Ca_%`~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(`rel.Lf_Ca_%`~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(`rel.Lf_Ca_%`~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(`rel.Lf_Ca_%`~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(`rel.Lf_Ca_%`~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(`rel.Lf_Ca_%`~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(`rel.Lf_Ca_%`~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(`rel.Lf_Ca_%`~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(`rel.Lf_Ca_%`~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # ok
plot(resid(m1)~fitted(m1)) # good
plot(resid(m2)~fitted(m2)) # good
plot(resid(m3)~fitted(m3)) #good
plot(resid(m4)~fitted(m4)) # ok
plot(resid(m5)~fitted(m5)) # good
plot(resid(m6)~fitted(m6)) # good
plot(resid(m7)~fitted(m7)) # good
plot(resid(m8)~fitted(m8)) # good
plot(resid(m9)~fitted(m9)) #good
plot(resid(m10)~fitted(m10)) # good
plot(resid(m11)~fitted(m11)) # ok
plot(resid(m12)~fitted(m12)) # ok
plot(resid(m13)~fitted(m13)) # ok
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13) # including only the models with good to ok fit
# #(Int) Gly Hmc Mlc Tmp Gly:Hmc Gly:Mlc Gly:Tmp Hmc:Mlc Hmc:Tmp Mlc:Tmp Gly:Hmc:Mlc Gly:Hmc:Tmp Gly:Mlc:Tmp Hmc:Mlc:Tmp Gly:Hmc:Mlc:Tmp df logLik  AICc delta weight
# m13 0.12650   +   +   +   +                                                                                                                  8 38.201 -58.1  0.00  0.771
# m11 0.10410   +   +   +   +                               +                                                                                  9 37.299 -53.7  4.42  0.084
# m10 0.12340   +   +   +   +               +                                                                                                  9 37.052 -53.2  4.91  0.066
# m12 0.12580   +   +   +   +       +                                                                                                          9 37.048 -53.2  4.92  0.066
# m7  0.09803   +   +   +   +                       +                                                                                         10 35.962 -48.3  9.80  0.006
# m6  0.15240   +   +   +   +                                               +                                                                 10 35.577 -47.5 10.57  0.004
# m8  0.14900   +   +   +   +                                       +                                                                         10 35.138 -46.7 11.45  0.003
# m9  0.10020   +   +   +   +       +       +               +                                                                                 11 35.051 -43.7 14.41  0.001
# m4  0.10530   +   +   +   +       +       +               +                           +                                                     12 34.672 -40.1 18.06  0.000
# m1  0.07017   +   +   +   +               +       +                       +                                   +                             15 37.701 -36.8 21.29  0.000
# m3  0.08341   +   +   +   +       +               +               +                               +                                         15 32.565 -26.6 31.56  0.000
# m5  0.14640   +   +   +   +                       +               +       +                                                                 14 30.356 -25.3 32.77  0.000
# m2  0.15340   +   +   +   +                               +       +       +                                               +                 15 28.611 -18.7 39.46  0.000
# m0  0.01350   +   +   +   +       +       +       +       +       +       +           +           +           +           +               + 26 33.457  16.3 74.40  0.000
AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13)
#.    df      AIC
# m0  26 -14.91317
# m1  15 -45.40177
# m2  15 -27.22231
# m3  15 -35.13048
# m4  12 -45.34335
# m5  14 -32.71193
# m6  10 -51.15324
# m7  10 -51.92354
# m8  10 -50.27513
# m9  11 -48.10129
# m10  9 -56.10417
# m11  9 -56.59710
# m12  9 -56.09530
# m13  8 -60.40123 <--

# Test
anova(m13, type = "III") 
#.                            Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
# Mulch      0.009023 0.0090231     1    20  1.2988 0.26791  
# Glyphosate 0.003817 0.0038166     1    20  0.5493 0.46720  
# Humic      0.001948 0.0019478     1    20  0.2804 0.60229  
# Timepoint  0.048928 0.0244641     2    46  3.5213 0.03776 *
emm<-emmeans(m13, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response")
pairs(emm)
#emmip(m1,Mulch ~ Timepoint | Glyphosate)
#emmip(m1,Mulch:Glyphosate:Humic ~ Timepoint)
emmip(m13,Mulch:Glyphosate ~ Timepoint|Humic)
#emmip(m1,Glyphosate ~ Timepoint | Mulch)
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
  Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
  Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.Ca <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=`rel.Lf_Ca_%`),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Calcium (%)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_Ca_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

#### leaf Sodium ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = `rel.Lf_Na_%`, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)

# Mixed Effects Model
m0<-lmer(`rel.Lf_Na_%`~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(`rel.Lf_Na_%`~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(`rel.Lf_Na_%`~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(`rel.Lf_Na_%`~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(`rel.Lf_Na_%`~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(`rel.Lf_Na_%`~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(`rel.Lf_Na_%`~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(`rel.Lf_Na_%`~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(`rel.Lf_Na_%`~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(`rel.Lf_Na_%`~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(`rel.Lf_Na_%`~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(`rel.Lf_Na_%`~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(`rel.Lf_Na_%`~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(`rel.Lf_Na_%`~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # good
plot(resid(m1)~fitted(m1)) # good
plot(resid(m2)~fitted(m2)) # good
plot(resid(m3)~fitted(m3)) #  good
plot(resid(m4)~fitted(m4)) # ok
plot(resid(m5)~fitted(m5)) # ok
plot(resid(m6)~fitted(m6)) # good
plot(resid(m7)~fitted(m7)) #  good
plot(resid(m8)~fitted(m8)) # ok
plot(resid(m9)~fitted(m9)) #  good
plot(resid(m10)~fitted(m10)) #  good
plot(resid(m11)~fitted(m11)) #  good
plot(resid(m12)~fitted(m12)) #  good
plot(resid(m13)~fitted(m13)) #  good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13) # including only the models with good to ok fit
# 
#         (Int) Gly Hmc Mlc Tmp Gly:Hmc Gly:Mlc Gly:Tmp Hmc:Mlc Hmc:Tmp Mlc:Tmp Gly:Hmc:Mlc Gly:Hmc:Tmp Gly:Mlc:Tmp Hmc:Mlc:Tmp Gly:Hmc:Mlc:Tmp df logLik  AICc delta weight
# m8  -0.21830   +   +   +   +                                       +                                                                         10 35.391 -47.2  0.00  0.743
# m10 -0.08550   +   +   +   +               +                                                                                                  9 32.357 -43.8  3.36  0.138
# m13 -0.15440   +   +   +   +                                                                                                                  8 30.721 -43.2  4.02  0.100
# m11 -0.13970   +   +   +   +                               +                                                                                  9 29.660 -38.4  8.76  0.009
# m12 -0.15890   +   +   +   +       +                                                                                                          9 29.558 -38.2  8.96  0.008
# m9  -0.07522   +   +   +   +       +       +               +                                                                                 11 29.942 -33.5 13.69  0.001
# m7  -0.12700   +   +   +   +                       +                                                                                         10 28.076 -32.5 14.63  0.000
# m6  -0.16650   +   +   +   +                                               +                                                                 10 26.968 -30.3 16.85  0.000
# m4  -0.07373   +   +   +   +       +       +               +                           +                                                     12 29.390 -29.5 17.68  0.000
# m5  -0.20290   +   +   +   +                       +               +       +                                                                 14 28.970 -22.6 24.60  0.000
# m3  -0.19010   +   +   +   +       +               +               +                               +                                         15 28.883 -19.2 27.98  0.000
# m2  -0.20840   +   +   +   +                               +       +       +                                               +                 15 27.517 -16.5 30.71  0.000
# m1  -0.06592   +   +   +   +               +       +                       +                                   +                             15 23.371  -8.2 39.00  0.000
# m0  -0.12530   +   +   +   +       +       +       +       +       +       +           +           +           +           +               + 26 18.470  46.3 93.43  0.000
AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13)
# #.   df      AIC
# m0  26  15.05965
# m1  15 -16.74105
# m2  15 -25.03396
# m3  15 -27.76628
# m4  12 -34.78047
# m5  14 -29.93914
# m6  10 -33.93526
# m7  10 -36.15146
# m8  10 -50.78112 <--
# m9  11 -37.88349
# m10  9 -46.71319
# m11  9 -41.31930
# m12  9 -41.11662
# m13  8 -45.44270
# Test
anova(m8, type = "III") 
#.                            Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
# Humic           0.000043 0.000043     1    20  0.0063 0.9375560    
# Timepoint       0.177006 0.088503     2    44 12.8362 4.061e-05 ***
# Mulch           0.020106 0.020106     1    20  2.9162 0.1031737    
# Glyphosate      0.013782 0.013782     1    20  1.9989 0.1727943    
# Humic:Timepoint 0.146717 0.073358     2    44 10.6397 0.0001702 ***
emm<-emmeans(m8, specs=~Humic | Timepoint, type="response") 
pairs(emm, adjust="holm") # Humic | Timepoint not significant 

emm<-emmeans(m8, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response") # using this contrast to match the other plots
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
  Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
  Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.Na <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=`rel.Lf_Na_%`),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Sodium (%)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_Na_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

#### leaf Phosphorus ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = `rel.Lf_P_%`, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)
# Mixed Effects Model
m0<-lmer(`rel.Lf_P_%`~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(`rel.Lf_P_%`~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(`rel.Lf_P_%`~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(`rel.Lf_P_%`~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(`rel.Lf_P_%`~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(`rel.Lf_P_%`~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(`rel.Lf_P_%`~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(`rel.Lf_P_%`~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(`rel.Lf_P_%`~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(`rel.Lf_P_%`~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(`rel.Lf_P_%`~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(`rel.Lf_P_%`~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(`rel.Lf_P_%`~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(`rel.Lf_P_%`~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # not good
plot(resid(m1)~fitted(m1)) # not good
plot(resid(m2)~fitted(m2)) # ok
plot(resid(m3)~fitted(m3)) # ok
plot(resid(m4)~fitted(m4)) # ok
plot(resid(m5)~fitted(m5)) # not good
plot(resid(m6)~fitted(m6)) # not good
plot(resid(m7)~fitted(m7)) # not good
plot(resid(m8)~fitted(m8)) # not good
plot(resid(m9)~fitted(m9)) # not good
plot(resid(m10)~fitted(m10)) # not good
plot(resid(m11)~fitted(m11)) # not good
plot(resid(m12)~fitted(m12)) # not good
plot(resid(m13)~fitted(m13)) # not good
# Select most parsimonious model
#MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13) 
MuMIn::model.sel(m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13) # including only the models with good to ok fit

#        (Int) Gly Hmc Mlc Tmp Hmc:Mlc Hmc:Tmp Mlc:Tmp Hmc:Mlc:Tmp Gly:Hmc Gly:Tmp Gly:Hmc:Tmp Gly:Mlc Gly:Hmc:Mlc df logLik AICc delta weight
# m7  0.4069   +   +   +   +                                                   +                                 10 14.063 -4.5  0.00  0.944
# m13 0.2972   +   +   +   +                                                                                      8  8.043  2.2  6.72  0.033
# m10 0.3467   +   +   +   +                                                                       +              9  8.384  4.1  8.65  0.012
# m11 0.3323   +   +   +   +       +                                                                              9  7.559  5.8 10.30  0.005
# m12 0.3076   +   +   +   +                                           +                                          9  6.856  7.2 11.71  0.003
# m6  0.3630   +   +   +   +                       +                                                             10  7.916  7.8 12.29  0.002
# m5  0.4736   +   +   +   +               +       +                           +                                 14 11.480 12.4 16.93  0.000
# m9  0.3922   +   +   +   +       +                                   +                           +             11  6.687 13.0 17.55  0.000
# m8  0.2981   +   +   +   +               +                                                                     10  5.089 13.4 17.95  0.000
# m4  0.4201   +   +   +   +       +                                   +                           +           + 12  6.612 16.1 20.58  0.000
# m3  0.4566   +   +   +   +               +                           +       +           +                     15  8.952 20.7 25.19  0.000
# m2  0.4143   +   +   +   +       +       +       +           +                                                 15  2.940 32.7 37.21  0.000
AIC(m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13)
#.   df      AIC
# m2  15 24.1205261
# m3  15 12.0960075
# m4  12 10.7760501
# m5  14  5.0404529
# m6  10  4.1676267
# m7  10 -8.1253562 <--
# m8  10  9.8223727
# m9  11  8.6263478
# m10  9  1.2322038
# m11  9  2.8829727
# m12  9  4.2884271
# m13  8 -0.0857875
# Test
anova(m7, type = "III") 
#.                    Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
#Glyphosate           0.08872 0.08872     1    20  4.5415 0.0456784 *  
#Timepoint            0.72953 0.36477     2    44 18.6730 1.345e-06 ***
#Mulch                0.18795 0.18795     1    20  9.6214 0.0056213 ** 
#Humic                0.00628 0.00628     1    20  0.3216 0.5769725    
#Glyphosate:Timepoint 0.43294 0.21647     2    44 11.0814 0.0001266 ***
emm<-emmeans(m7, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response")
pairs(emm)
#emmip(m1,Mulch ~ Timepoint | Glyphosate)
#emmip(m1,Mulch:Glyphosate:Humic ~ Timepoint)
emmip(m7,Mulch:Glyphosate ~ Timepoint|Humic)
#emmip(m1,Glyphosate ~ Timepoint | Mulch)
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
  Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
  Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.P <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=`rel.Lf_P_%`),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Phosphorus (%)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_P_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

#### leaf Iron ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = rel.Lf_Fe_ppm, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)
# Mixed Effects Model
m0<-lmer(rel.Lf_Fe_ppm~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(rel.Lf_Fe_ppm~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(rel.Lf_Fe_ppm~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(rel.Lf_Fe_ppm~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(rel.Lf_Fe_ppm~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(rel.Lf_Fe_ppm~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(rel.Lf_Fe_ppm~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(rel.Lf_Fe_ppm~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(rel.Lf_Fe_ppm~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(rel.Lf_Fe_ppm~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(rel.Lf_Fe_ppm~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(rel.Lf_Fe_ppm~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(rel.Lf_Fe_ppm~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(rel.Lf_Fe_ppm~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # good
plot(resid(m1)~fitted(m1)) # good
plot(resid(m2)~fitted(m2)) # good
plot(resid(m3)~fitted(m3)) # good
plot(resid(m4)~fitted(m4)) # good
plot(resid(m5)~fitted(m5)) # good
plot(resid(m6)~fitted(m6)) # good
plot(resid(m7)~fitted(m7)) # good
plot(resid(m8)~fitted(m8)) # ok
plot(resid(m9)~fitted(m9)) # ok
plot(resid(m10)~fitted(m10)) # ok
plot(resid(m11)~fitted(m11)) # ok
plot(resid(m12)~fitted(m12)) # ok
plot(resid(m13)~fitted(m13)) # ok
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13) # including only the models with good to ok fit

#        (Int) Gly Hmc Mlc Tmp Gly:Hmc Gly:Mlc Gly:Tmp Hmc:Mlc Hmc:Tmp Mlc:Tmp Gly:Hmc:Mlc Gly:Hmc:Tmp Gly:Mlc:Tmp Hmc:Mlc:Tmp Gly:Hmc:Mlc:Tmp df logLik  AICc delta weight
# m6  -0.2629   +   +   +   +                                               +                                                                 10 33.934 -44.3  0.00  0.990 <--
# m5  -0.3543   +   +   +   +                       +               +       +                                                                 14 35.127 -34.9  9.38  0.009
# m13 -0.1612   +   +   +   +                                                                                                                  8 23.542 -28.8 15.46  0.000
# m1  -0.3314   +   +   +   +               +       +                       +                                   +                             15 33.571 -28.6 15.69  0.000
# m2  -0.2805   +   +   +   +                               +       +       +                                               +                 15 32.877 -27.2 17.08  0.000
# m12 -0.1339   +   +   +   +       +                                                                                                          9 22.793 -24.7 19.58  0.000
# m11 -0.1863   +   +   +   +                               +                                                                                  9 22.674 -24.4 19.82  0.000
# m10 -0.1484   +   +   +   +               +                                                                                                  9 22.225 -23.5 20.71  0.000
# m7  -0.2163   +   +   +   +                       +                                                                                         10 23.341 -23.1 21.19  0.000
# m8  -0.1976   +   +   +   +                                       +                                                                         10 21.457 -19.3 24.95  0.000
# m9  -0.1461   +   +   +   +       +       +               +                                                                                 11 20.599 -14.8 29.46  0.000
# m4  -0.1568   +   +   +   +       +       +               +                           +                                                     12 19.928 -10.6 33.69  0.000
# m3  -0.2187   +   +   +   +       +               +               +                               +                                         15 18.429   1.7 45.98  0.000
# m0  -0.3427   +   +   +   +       +       +       +       +       +       +           +           +           +           +               + 26 29.270  24.7 68.92  0.000
AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13)
#.   df      AIC
# m0  26  -6.539566
# m1  15 -37.142926
# m2  15 -35.753518
# m3  15  -6.857261
# m4  12 -15.855766
# m5  14 -42.253489
# m6  10 -47.867777 <--
# m7  10 -26.681251
# m8  10 -22.914489
# m9  11 -19.197069
# m10  9 -26.449462
# m11  9 -27.347576
# m12  9 -27.586944
# m13  8 -31.083951
# Test
anova(m6, type = "III") 
#.                Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
#Mulch           0.11377 0.113774     1    20 11.7691  0.002648 ** 
#Timepoint       0.57265 0.286325     2    44 29.6183 7.106e-09 ***
#Glyphosate      0.02123 0.021231     1    20  2.1962  0.153939    
#Humic           0.01052 0.010517     1    20  1.0879  0.309376    
#Mulch:Timepoint 0.37218 0.186091     2    44 19.2498 9.863e-07 ***
emm<-emmeans(m6, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response")
pairs(emm)
emmip(m6,Mulch:Glyphosate ~ Timepoint|Humic)
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
  Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
  Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.Fe <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=rel.Lf_Fe_ppm),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Iron (ppm)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_Fe_ppm_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

#### leaf Boron ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = rel.Lf_B_ppm, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)
# Mixed Effects Model
m0<-lmer(rel.Lf_B_ppm~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(rel.Lf_B_ppm~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(rel.Lf_B_ppm~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(rel.Lf_B_ppm~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(rel.Lf_B_ppm~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(rel.Lf_B_ppm~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(rel.Lf_B_ppm~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(rel.Lf_B_ppm~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(rel.Lf_B_ppm~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(rel.Lf_B_ppm~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(rel.Lf_B_ppm~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(rel.Lf_B_ppm~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(rel.Lf_B_ppm~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(rel.Lf_B_ppm~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # good
plot(resid(m1)~fitted(m1)) # good
plot(resid(m2)~fitted(m2)) # good
plot(resid(m3)~fitted(m3)) # good
plot(resid(m4)~fitted(m4)) # good
plot(resid(m5)~fitted(m5)) # good
plot(resid(m6)~fitted(m6)) # good
plot(resid(m7)~fitted(m7)) # good
plot(resid(m8)~fitted(m8)) # good
plot(resid(m9)~fitted(m9)) # good
plot(resid(m10)~fitted(m10)) # good
plot(resid(m11)~fitted(m11)) # good
plot(resid(m12)~fitted(m12)) # good
plot(resid(m13)~fitted(m13)) # good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13) # including only the models with good to ok fit

#       (Int) Gly Hmc Mlc Tmp Gly:Hmc Gly:Mlc Gly:Tmp Hmc:Mlc Hmc:Tmp Mlc:Tmp Gly:Hmc:Mlc Gly:Hmc:Tmp Gly:Mlc:Tmp Hmc:Mlc:Tmp Gly:Hmc:Mlc:Tmp df logLik  AICc delta weight
# m7  0.4603   +   +   +   +                       +                                                                                         10 33.068 -42.5  0.00  0.691
# m13 0.3950   +   +   +   +                                                                                                                  8 28.682 -39.1  3.45  0.123
# m5  0.4529   +   +   +   +                       +               +       +                                                                 14 37.111 -38.9  3.68  0.110
# m12 0.3551   +   +   +   +       +                                                                                                          9 28.353 -35.8  6.73  0.024
# m10 0.4311   +   +   +   +               +                                                                                                  9 28.195 -35.5  7.04  0.020
# m8  0.3450   +   +   +   +                                       +                                                                         10 29.245 -34.9  7.65  0.015
# m11 0.3775   +   +   +   +                               +                                                                                  9 27.673 -34.4  8.09  0.012
# m6  0.4376   +   +   +   +                                               +                                                                 10 27.920 -32.2 10.30  0.004
# m9  0.3737   +   +   +   +       +       +               +                                                                                 11 26.846 -27.3 15.24  0.000
# m3  0.3576   +   +   +   +       +               +               +                               +                                         15 32.422 -26.3 16.26  0.000
# m4  0.3729   +   +   +   +       +       +               +                           +                                                     12 26.367 -23.4 19.08  0.000
# m1  0.5383   +   +   +   +               +       +                       +                                   +                             15 29.848 -21.1 21.41  0.000
# m2  0.3518   +   +   +   +                               +       +       +                                               +                 15 25.827 -13.1 29.45  0.000
# m0  0.3749   +   +   +   +       +       +       +       +       +       +           +           +           +           +               + 26 26.835  29.5 72.06  0.000
AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13)
#.   df      AIC
# m0  26  -1.669414
# m1  15 -29.695532
# m2  15 -21.654007
# m3  15 -34.843490
# m4  12 -28.733763
# m5  14 -46.222654
# m6  10 -35.840013
# m7  10 -46.136398 <--
# m8  10 -38.490564
# m9  11 -31.691466
# m10  9 -38.389185
# m11  9 -37.346512
# m12  9 -38.706741
# m13  8 -41.363604
# Test
anova(m7, type = "III") 
#.                            Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
# Glyphosate           0.019326 0.019326     1    20  2.5447 0.1263433    
# Timepoint            0.204636 0.102318     2    44 13.4726 2.727e-05 ***
# Mulch                0.125352 0.125352     1    20 16.5056 0.0006076 ***
# Humic                0.010994 0.010994     1    20  1.4477 0.2429504    
# Glyphosate:Timepoint 0.153494 0.076747     2    44 10.1056 0.0002446 ***
emm<-emmeans(m7, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response")
pairs(emm)
#emmip(m1,Mulch ~ Timepoint | Glyphosate)
#emmip(m1,Mulch:Glyphosate:Humic ~ Timepoint)
emmip(m7,Mulch:Glyphosate ~ Timepoint|Humic)
#emmip(m1,Glyphosate ~ Timepoint | Mulch)
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
  Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
  Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.B <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=rel.Lf_B_ppm),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Boron (ppm)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_B_ppm_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

#### leaf Copper ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = rel.Lf_Cu_ppm, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)
# Mixed Effects Model
m0<-lmer(rel.Lf_Cu_ppm~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(rel.Lf_Cu_ppm~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(rel.Lf_Cu_ppm~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(rel.Lf_Cu_ppm~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(rel.Lf_Cu_ppm~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(rel.Lf_Cu_ppm~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(rel.Lf_Cu_ppm~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(rel.Lf_Cu_ppm~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(rel.Lf_Cu_ppm~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(rel.Lf_Cu_ppm~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(rel.Lf_Cu_ppm~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(rel.Lf_Cu_ppm~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(rel.Lf_Cu_ppm~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(rel.Lf_Cu_ppm~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # good
plot(resid(m1)~fitted(m1)) # good
plot(resid(m2)~fitted(m2)) # not good
plot(resid(m3)~fitted(m3)) # ok
plot(resid(m4)~fitted(m4)) # not good
plot(resid(m5)~fitted(m5)) # not good
plot(resid(m6)~fitted(m6)) # not good
plot(resid(m7)~fitted(m7)) # not good
plot(resid(m8)~fitted(m8)) # not good
plot(resid(m9)~fitted(m9)) # not good
plot(resid(m10)~fitted(m10)) # not good
plot(resid(m11)~fitted(m11)) # not good
plot(resid(m12)~fitted(m12)) # not good
plot(resid(m13)~fitted(m13)) # not good
# Select most parsimonious model
MuMIn::model.sel(m0,m1,m3) # including only the models with good to ok fit
# 
#     (Int) Gly Hmc Mlc Tmp Gly:Hmc Gly:Mlc Gly:Tmp Hmc:Mlc Hmc:Tmp Mlc:Tmp Gly:Hmc:Mlc Gly:Hmc:Tmp Gly:Mlc:Tmp Hmc:Mlc:Tmp Gly:Hmc:Mlc:Tmp df   logLik  AICc delta weight
# m0 5.356   +   +   +   +       +       +       +       +       +       +           +           +           +           +               + 26 -123.667 330.5  0.00      1
# m1 7.646   +   +   +   +               +       +                       +                                   +                             15 -163.936 366.4 35.91      0
# m3 5.927   +   +   +   +       +               +               +                               +                                         15 -164.957 368.5 37.95      0
AIC(m0,m1,m3)
#.   df      AIC
# m0 26 299.3339
# m1 15 357.8719
# m3 15 359.9141
# Test
anova(m0, type = "III") 
#                                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# Mulch                             6.951   6.951     1    16  4.0021  0.062709 .  
# Glyphosate                        1.961   1.961     1    16  1.1289  0.303791    
# Humic                             0.566   0.566     1    16  0.3258  0.576101    
# Timepoint                        97.618  48.809     2    32 28.1020 9.007e-08 ***
# Mulch:Glyphosate                  4.132   4.132     1    16  2.3788  0.142533    
# Mulch:Humic                       1.190   1.190     1    16  0.6853  0.419931    
# Glyphosate:Humic                  2.739   2.739     1    16  1.5769  0.227230    
# Mulch:Timepoint                   0.060   0.030     2    32  0.0172  0.982942    
# Glyphosate:Timepoint             11.529   5.765     2    32  3.3189  0.048999 *  
# Humic:Timepoint                  24.257  12.129     2    32  6.9830  0.003044 ** 
# Mulch:Glyphosate:Humic            0.896   0.896     1    16  0.5159  0.482940    
# Mulch:Glyphosate:Timepoint       44.044  22.022     2    32 12.6792 8.808e-05 ***
# Mulch:Humic:Timepoint            54.836  27.418     2    32 15.7860 1.699e-05 ***
# Glyphosate:Humic:Timepoint       14.628   7.314     2    32  4.2112  0.023793 *  
# Mulch:Glyphosate:Humic:Timepoint 15.392   7.696     2    32  4.4310  0.020012 *  
emm<-emmeans(m0, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response")
#emm<-emmeans(m0, specs=~Mulch:Glyphosate:Humic|Timepoint, type="response")

pairs(emm)
#emmip(m1,Mulch ~ Timepoint | Glyphosate)
#emmip(m1,Mulch:Glyphosate:Humic ~ Timepoint)
emmip(m1,Mulch:Glyphosate ~ Timepoint|Humic)
#emmip(m1,Glyphosate ~ Timepoint | Mulch)
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
  Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
  Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.Cu <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=rel.Lf_Cu_ppm),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Copper (ppm)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_Cu_ppm_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

#### leaf Zinc ####
# Plot raw points 
ggplot(nutrients.relative, aes(x = Timepoint, y = rel.Lf_Zn_ppm, color = Treatment_code, fill = Treatment_code)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = Treatment_code),fun = mean,geom = "point",size = 3) +
  stat_summary(aes(group = Treatment_code),fun.data = mean_se,geom = "errorbar",width = 0.1) +
  geom_hline(linetype="dashed",yintercept = 0) +
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank()) +
  facet_wrap(Mulch~.)
# Mixed Effects Model
m0<-lmer(rel.Lf_Zn_ppm~ Mulch*Glyphosate*Humic*Timepoint + (1|Treat_rep), data=nutrients.relative) # 4-way interactions
# single 3-way interaction
m1<-lmer(rel.Lf_Zn_ppm~ Mulch * Glyphosate * Timepoint + Humic + (1|Treat_rep), data=nutrients.relative) # 
m2<-lmer(rel.Lf_Zn_ppm~ Mulch * Humic * Timepoint + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
m3<-lmer(rel.Lf_Zn_ppm~ Glyphosate * Humic * Timepoint + Mulch + (1|Treat_rep), data=nutrients.relative) # 
m4<-lmer(rel.Lf_Zn_ppm~ Mulch * Glyphosate * Humic +Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# No treatment interactions, just treat X timepoints
m5<-lmer(rel.Lf_Zn_ppm~ Mulch*Timepoint + Glyphosate*Timepoint + Humic*Timepoint + (1|Treat_rep), data=nutrients.relative)
m6<-lmer(rel.Lf_Zn_ppm~ Mulch * Timepoint + Glyphosate + Humic + (1|Treat_rep), data=nutrients.relative) # 
m7<-lmer(rel.Lf_Zn_ppm~ Glyphosate*Timepoint +Mulch + Humic + (1|Treat_rep), data=nutrients.relative) # 
m8<-lmer(rel.Lf_Zn_ppm~ Humic*Timepoint +Mulch + Glyphosate + (1|Treat_rep), data=nutrients.relative) # 
# interactions consistent across timepoints
m9<-lmer(rel.Lf_Zn_ppm~ Mulch*Glyphosate + Mulch*Humic + Glyphosate*Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
m10<-lmer(rel.Lf_Zn_ppm~ Mulch * Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m11<-lmer(rel.Lf_Zn_ppm~ Mulch * Humic + Glyphosate + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
m12<-lmer(rel.Lf_Zn_ppm~ Glyphosate * Humic + Mulch + Timepoint + (1|Treat_rep), data=nutrients.relative) # # single 2-way interactions
# additive 
m13<-lmer(rel.Lf_Zn_ppm~ Mulch + Glyphosate + Humic + Timepoint + (1|Treat_rep), data=nutrients.relative) # 
# Check fit
plot(resid(m0)~fitted(m0)) # ok
plot(resid(m1)~fitted(m1)) # not good
plot(resid(m2)~fitted(m2)) # not good
plot(resid(m3)~fitted(m3)) # not good
plot(resid(m4)~fitted(m4)) # not good
plot(resid(m5)~fitted(m5)) # not good
plot(resid(m6)~fitted(m6)) # not good
plot(resid(m7)~fitted(m7)) # not good
plot(resid(m8)~fitted(m8)) # not good
plot(resid(m9)~fitted(m9)) # not good
plot(resid(m10)~fitted(m10)) # not good
plot(resid(m11)~fitted(m11)) # not good
plot(resid(m12)~fitted(m12)) # not good
plot(resid(m13)~fitted(m13)) # not good
# Select most parsimonious model
MuMIn::model.sel(m0) # including only the models with good to ok fit
AIC(m0)
# Test
anova(m0, type = "III") 
#                                  Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
# Mulch                              0.44    0.44     1    16   0.3085 0.5862916    
# Glyphosate                         0.48    0.48     1    16   0.3347 0.5709846    
# Humic                              0.04    0.04     1    16   0.0304 0.8637139    
# Timepoint                        874.66  437.33     2    32 306.9310 < 2.2e-16 ***
# Mulch:Glyphosate                   3.47    3.47     1    16   2.4336 0.1383191    
# Mulch:Humic                        4.74    4.74     1    16   3.3260 0.0869270 .  
# Glyphosate:Humic                  16.70   16.70     1    16  11.7234 0.0034798 ** 
# Mulch:Timepoint                    2.48    1.24     2    32   0.8687 0.4291507    
# Glyphosate:Timepoint               4.50    2.25     2    32   1.5788 0.2218567    
# Humic:Timepoint                    3.42    1.71     2    32   1.2001 0.3143629    
# Mulch:Glyphosate:Humic             1.39    1.39     1    16   0.9729 0.3386365    
# Mulch:Glyphosate:Timepoint        30.87   15.44     2    32  10.8335 0.0002553 ***
# Mulch:Humic:Timepoint             24.98   12.49     2    32   8.7661 0.0009209 ***
# Glyphosate:Humic:Timepoint        38.15   19.07     2    32  13.3873 5.962e-05 ***
# Mulch:Glyphosate:Humic:Timepoint   3.78    1.89     2    32   1.3258 0.2797931    
emm<-emmeans(m0, specs=~Mulch:Glyphosate|Humic|Timepoint, type="response")
pairs(emm)
#emmip(m1,Mulch ~ Timepoint | Glyphosate)
#emmip(m1,Mulch:Glyphosate:Humic ~ Timepoint)
emmip(m0,Mulch:Glyphosate ~ Timepoint|Humic)
#emmip(m1,Glyphosate ~ Timepoint | Mulch)
cld<-cld(emm, alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# add treatment code to table
cld.df <- cld %>%as.data.frame() %>%mutate(Treatment_code = case_when(
  Mulch == "no"    & Glyphosate == "no"          & Humic == "no"    ~ "CCC",
  Mulch == "no"    & Glyphosate == "no"          & Humic == "humic" ~ "CCH",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "no"    ~ "CGC",
  Mulch == "no"    & Glyphosate == "Glyphosate" & Humic == "humic" ~ "CGH",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "no"    ~ "MCC",
  Mulch == "mulch" & Glyphosate == "no"          & Humic == "humic" ~ "MCH",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "no"    ~ "MGC",
  Mulch == "mulch" & Glyphosate == "Glyphosate" & Humic == "humic" ~ "MGH"))
cld.summary.Zn <- cld.df %>%dplyr::select(Timepoint, Treatment_code, .group) %>%mutate(Treatment_code = factor(Treatment_code,
                                                                                                              levels = c("CCC","CGC","MCC","MGC","CCH","CGH","MCH","MGH"))) %>%arrange(Timepoint, Treatment_code)
# plot EMMs and raw data points
ggplot(cld.df,aes(x = Timepoint, y = emmean, color = Treatment_code, fill = Treatment_code, group = Treatment_code)) +
  geom_point(data=nutrients.relative,mapping=aes(x=Timepoint,y=rel.Lf_Zn_ppm),position = position_jitter(width = 0.07), alpha = 0.5, size=2) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width=0.1,size = 1) +
  geom_line(size=0.8) + geom_point(shape=18,size=4) +
  geom_hline(linetype="dashed",yintercept = 0) +
  ylab("Rel. Leaf Zinc (ppm)")+
  scale_color_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH"))+
  scale_fill_manual(values=c("#277da1","#284b63","#43aa8b","#90be6d","#f94144","#9d0208","#f9844a","#f9c74f"),limits=c("CCC","CCH","CGC","CGH","MCC","MCH","MGC","MGH")) +
  theme_bw()+theme(panel.grid = element_blank(),strip.background=element_blank()) +
  facet_wrap(Humic~., labeller = labeller(Humic=c("humic"="humic","no"="no humic")))
ggsave(filename= "EMM.rel.Lf_Zn_ppm_Timepoint.line.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")
#ggsave(filename= "EMM.rel.Lf_Zn_ppm_Timepoint.line_no_fct.pdf", device="pdf", units="mm", dpi=300, width=125, height=75, path="Revisions.ISME.Comm/plots/")

######### Combine CLD table, manually add significance to plots in Illustrator #####
cld.all <- bind_rows(list(
    Zn = cld.summary.Zn,
    N  = cld.summary.N,
    Na = cld.summary.Na,
    K  = cld.summary.K,
    P  = cld.summary.P,
    Cu = cld.summary.Cu,
    Ca = cld.summary.Ca,
    B  = cld.summary.B,
    Fe = cld.summary.Fe),
  .id = "Nutrient")

writexl::write_xlsx(cld.all,"Revisions.ISME.Comm/nutrient_EMM_summary.xlsx")

######################################################
######################################################
#### LI-COR Analysis #################################
######################################################
######################################################
#########################
####### Licor E #########
#########################
# check for outliers
ggplot(licor_cleaned, aes(y=licor_E, x=mulch, color=glyphosate, shape=humic_acid, label = ID)) +
  geom_boxplot(outliers = FALSE) +
  #facet_grid(.~year) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("2024_CCC1","2024_CGH1","2024_MCC6")) 
licor_cleaned.2<-licor_cleaned %>% filter(!ID %in% outliers)

m0 <- lm(log10(licor_E) ~ mulch*glyphosate*humic_acid + year, data = licor_cleaned.2)
m1 <- lm(log10(licor_E) ~ mulch*glyphosate + mulch*humic_acid + glyphosate*humic_acid + year, data = licor_cleaned.2)
m2 <- lm(log10(licor_E) ~ mulch*glyphosate + humic_acid + year, data = licor_cleaned.2)
m3 <- lm(log10(licor_E) ~ mulch*humic_acid + glyphosate + year, data = licor_cleaned.2)
m4 <- lm(log10(licor_E) ~ glyphosate*humic_acid + mulch + year, data = licor_cleaned.2)
m5 <- lm(log10(licor_E) ~ glyphosate*humic_acid + glyphosate*mulch + year, data = licor_cleaned.2)
m6 <- lm(log10(licor_E) ~ mulch*humic_acid + glyphosate*mulch + year, data = licor_cleaned.2)
m7 <- lm(log10(licor_E) ~ mulch*humic_acid + glyphosate*humic_acid + year, data = licor_cleaned.2)
m8 <- lm(log10(licor_E) ~ mulch + humic_acid + glyphosate + year, data = licor_cleaned.2)
m9 <- lm(log10(licor_E) ~ mulch*humic_acid + year, data = licor_cleaned.2)
m10 <- lm(log10(licor_E) ~ mulch*glyphosate + year, data = licor_cleaned.2)
m11 <- lm(log10(licor_E) ~ humic_acid*glyphosate + year, data = licor_cleaned.2)
m12 <- lm(log10(licor_E) ~ mulch + glyphosate + year, data = licor_cleaned.2)
m13 <- lm(log10(licor_E) ~ mulch + humic_acid + year, data = licor_cleaned.2)
m14 <- lm(log10(licor_E) ~ humic_acid + glyphosate + year, data = licor_cleaned.2)
m15 <- lm(log10(licor_E) ~ mulch + year, data = licor_cleaned.2)
m16 <- lm(log10(licor_E) ~ humic_acid + year, data = licor_cleaned.2)
m17 <- lm(log10(licor_E) ~ glyphosate + year, data = licor_cleaned.2)

plot(resid(m0)~fitted(m0))
plot(resid(m1)~fitted(m1))
plot(resid(m2)~fitted(m2))
plot(resid(m3)~fitted(m3))
plot(resid(m4)~fitted(m4))
plot(resid(m5)~fitted(m5))
plot(resid(m6)~fitted(m6))
plot(resid(m7)~fitted(m7))
plot(resid(m8)~fitted(m8))
plot(resid(m9)~fitted(m9))
plot(resid(m10)~fitted(m10))
plot(resid(m11)~fitted(m11))
plot(resid(m12)~fitted(m12))
plot(resid(m13)~fitted(m13))
plot(resid(m14)~fitted(m14))
plot(resid(m15)~fitted(m15))
plot(resid(m16)~fitted(m16))
plot(resid(m17)~fitted(m17))
# model selection
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17) # including only the models with good to ok fit
#       (Int) gly hmc_acd mlc     yer gly:hmc_acd gly:mlc hmc_acd:mlc gly:hmc_acd:mlc df logLik   AICc delta weight
# m12 -75.72   +           + 0.03601                                                  5 71.848 -133.4  0.00  0.298 <--
# m8  -75.53   +       +   + 0.03591                                                  6 72.415 -132.4  1.00  0.181
# m10 -75.89   +           + 0.03609                   +                              6 72.219 -132.0  1.39  0.149
# m2  -75.70   +       +   + 0.03599                   +                              7 72.792 -131.0  2.40  0.090
# m4  -75.91   +       +   + 0.03610           +                                      7 72.649 -130.7  2.69  0.078
# m3  -75.55   +       +   + 0.03592                               +                  7 72.418 -130.2  3.15  0.062
# m5  -76.09   +       +   + 0.03618           +       +                              8 73.031 -129.3  4.11  0.038
# m6  -75.72   +       +   + 0.03600                   +           +                  8 72.797 -128.8  4.57  0.030
# m7  -75.93   +       +   + 0.03611           +                   +                  8 72.653 -128.5  4.86  0.026
# m17 -75.21   +             0.03574                                                  4 67.751 -127.3  6.08  0.014
# m1  -76.11   +       +   + 0.03619           +       +           +                  9 73.036 -127.1  6.30  0.013
# m14 -75.03   +       +     0.03564                                                  5 68.278 -126.2  7.14  0.008
# m0  -76.21   +       +   + 0.03625           +       +           +               + 10 73.208 -125.2  8.19  0.005
# m11 -75.40   +       +     0.03583           +                                      6 68.492 -124.5  8.84  0.004
# m15 -75.11               + 0.03569                                                  4 66.046 -123.9  9.49  0.003
# m13 -74.92           +   + 0.03558                                                  5 66.636 -122.9 10.42  0.002
# m9  -74.94           +   + 0.03559                               +                  6 66.640 -120.8 12.55  0.001
# m16 -74.41           +     0.03531                                                  4 62.691 -117.2 16.20  0.000

anova(m12)
#             Df Sum Sq Mean Sq F value    Pr(>F)    
# mulch        1 0.2299 0.22989  8.2206 0.0046224 ** 
# glyphosate   1 0.3250 0.32497 11.6204 0.0008006 ***
# year         1 0.1487 0.14868  5.3165 0.0222339 *  
# Residuals  185 5.1736 0.02797                      

## VISUALISATION OF LICOR E: ####
emm_licor_e<-emmeans(m12, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#cld 
cld.licor_e<-cld(emm_licor_e, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#Visualize the results
ggplot(data.frame(cld.licor_e), aes(x=mulch, y=response, fill=glyphosate, color=glyphosate))+
  geom_point(data=licor_cleaned.2,aes(x=mulch,y=licor_E, fill=glyphosate, color=glyphosate), , position=position_jitterdodge(dodge.width = 1.0, jitter.width = 0.2), alpha = 0.2, size=2) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), linewidth = 1.5,size=1.2, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_x_discrete(labels=c("mulch"="mulch", "control"="no mulch")) + 
  geom_text(aes(label=.group, y = upper.CL, group=glyphosate),color="black",size=3, vjust = -14,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Transpiration rate (mol m-2 s-1)") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.licor.e.pdf", device="pdf", units="mm", dpi=300, width=97, height=100, path="Revisions.ISME.Comm/plots/")

#########################
####### Licor A #########
#########################
# check for outliers
ggplot(licor_cleaned, aes(y=licor_A, x=mulch, color=glyphosate, shape=humic_acid, label = ID)) +
  geom_boxplot(outliers = FALSE) +
  #facet_grid(.~year) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("2024_CCH8","2024_CGH1")) 
licor_cleaned.2<-licor_cleaned %>% filter(!ID %in% outliers)

m0 <- lm(sqrt(licor_A) ~ mulch*glyphosate*humic_acid + year, data = licor_cleaned.2)
m1 <- lm(sqrt(licor_A) ~ mulch*glyphosate + mulch*humic_acid + glyphosate*humic_acid + year, data = licor_cleaned.2)
m2 <- lm(sqrt(licor_A) ~ mulch*glyphosate + humic_acid + year, data = licor_cleaned.2)
m3 <- lm(sqrt(licor_A) ~ mulch*humic_acid + glyphosate + year, data = licor_cleaned.2)
m4 <- lm(sqrt(licor_A) ~ glyphosate*humic_acid + mulch + year, data = licor_cleaned.2)
m5 <- lm(sqrt(licor_A) ~ glyphosate*humic_acid + glyphosate*mulch + year, data = licor_cleaned.2)
m6 <- lm(sqrt(licor_A) ~ mulch*humic_acid + glyphosate*mulch + year, data = licor_cleaned.2)
m7 <- lm(sqrt(licor_A) ~ mulch*humic_acid + glyphosate*humic_acid + year, data = licor_cleaned.2)
m8 <- lm(sqrt(licor_A) ~ mulch + humic_acid + glyphosate + year, data = licor_cleaned.2)
m9 <- lm(sqrt(licor_A) ~ mulch*humic_acid + year, data = licor_cleaned.2)
m10 <- lm(sqrt(licor_A) ~ mulch*glyphosate + year, data = licor_cleaned.2)
m11 <- lm(sqrt(licor_A) ~ humic_acid*glyphosate + year, data = licor_cleaned.2)
m12 <- lm(sqrt(licor_A) ~ mulch + glyphosate + year, data = licor_cleaned.2)
m13 <- lm(sqrt(licor_A) ~ mulch + humic_acid + year, data = licor_cleaned.2)
m14 <- lm(sqrt(licor_A) ~ humic_acid + glyphosate + year, data = licor_cleaned.2)
m15 <- lm(sqrt(licor_A) ~ mulch + year, data = licor_cleaned.2)
m16 <- lm(sqrt(licor_A) ~ humic_acid + year, data = licor_cleaned.2)
m17 <- lm(sqrt(licor_A) ~ glyphosate + year, data = licor_cleaned.2)

plot(resid(m0)~fitted(m0))
plot(resid(m1)~fitted(m1))
plot(resid(m2)~fitted(m2))
plot(resid(m3)~fitted(m3))
plot(resid(m4)~fitted(m4))
plot(resid(m5)~fitted(m5))
plot(resid(m6)~fitted(m6))
plot(resid(m7)~fitted(m7))
plot(resid(m8)~fitted(m8))
plot(resid(m9)~fitted(m9))
plot(resid(m10)~fitted(m10))
plot(resid(m11)~fitted(m11))
plot(resid(m12)~fitted(m12))
plot(resid(m13)~fitted(m13))
plot(resid(m14)~fitted(m14))
plot(resid(m15)~fitted(m15))
plot(resid(m16)~fitted(m16))
plot(resid(m17)~fitted(m17))
# model selection
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17) # including only the models with good to ok fit
#       (Int) gly hmc_acd mlc     yer gly:hmc_acd gly:mlc hmc_acd:mlc gly:hmc_acd:mlc df  logLik  AICc delta weight
# m12 -177.5   +           + 0.08873                                                  5 -69.961 150.2  0.00  0.220 <--
# m10 -177.6   +           + 0.08873                   +                              6 -68.914 150.3  0.04  0.216 
# m4  -177.2   +       +   + 0.08858           +                                      7 -68.602 151.8  1.57  0.100
# m5  -177.2   +       +   + 0.08858           +       +                              8 -67.514 151.8  1.57  0.100
# m8  -177.2   +       +   + 0.08858                                                  6 -69.882 152.2  1.97  0.082
# m2  -177.2   +       +   + 0.08858                   +                              7 -68.834 152.3  2.03  0.080
# m7  -176.6   +       +   + 0.08828           +                   +                  8 -68.329 153.5  3.21  0.044
# m1  -176.6   +       +   + 0.08828           +       +           +                  9 -67.238 153.5  3.23  0.044
# m3  -176.6   +       +   + 0.08828                               +                  7 -69.613 153.8  3.59  0.037
# m6  -176.6   +       +   + 0.08828                   +           +                  8 -68.561 153.9  3.67  0.035
# m0  -176.6   +       +   + 0.08828           +       +           +               + 10 -67.094 155.4  5.17  0.017
# m15 -177.6               + 0.08873                                                  4 -73.668 155.6  5.30  0.016
# m13 -177.3           +   + 0.08858                                                  5 -73.592 157.5  7.26  0.006
# m9  -176.7           +   + 0.08828                               +                  6 -73.333 159.1  8.88  0.003
# m17 -174.2   +             0.08701                                                  4 -79.274 166.8 16.51  0.000
# m11 -173.8   +       +     0.08684           +                                      6 -78.024 168.5 18.26  0.000
# m14 -173.8   +       +     0.08684                                                  5 -79.184 168.7 18.44  0.000
# m16 -173.9           +     0.08684                                                  4 -82.553 173.3 23.07  0.000

anova(m12)
#             Df Sum Sq Mean Sq F value    Pr(>F)    
# mulch        1  2.3580 2.35803 18.8779 2.289e-05 ***
# glyphosate   1  0.9244 0.92437  7.4003   0.00714 ** 
# year         1  0.9089 0.90887  7.2762   0.00763 ** 
# Residuals  186 23.2332 0.12491                      

## VISUALISATION OF LICOR A: ####
emm_licor_A<-emmeans(m12, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#cld 
cld.licor_A<-cld(emm_licor_A, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#Visualize the results
ggplot(data.frame(cld.licor_A), aes(x=mulch, y=response, fill=glyphosate, color=glyphosate))+
  geom_point(data=licor_cleaned.2,aes(x=mulch,y=licor_A, fill=glyphosate, color=glyphosate), , position=position_jitterdodge(dodge.width = 1.0, jitter.width = 0.2), alpha = 0.2, size=2) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), linewidth = 1.5,size=1.2, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_x_discrete(labels=c("mulch"="mulch", "control"="no mulch")) + 
  geom_text(aes(label=.group, y = upper.CL, group=glyphosate),color="black",size=3, vjust = -14,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Assimilation rate (µmol m-2 s-1)") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.licor.a.pdf", device="pdf", units="mm", dpi=300, width=90, height=100, path="Revisions.ISME.Comm/plots/")

######################################################
######################################################
#### Weed Biomass Analysis ###########################
######################################################
######################################################
weed_biomass # data frame
# check for outliers
ggplot(weed_biomass, aes(y=weed_dry_weight, x=mulch, color=glyphosate, shape=humic_acid, label = ID)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c()) 
weed_biomass.2<-weed_biomass %>% filter(!ID %in% outliers)

m0 <- lm(sqrt(weed_dry_weight) ~ mulch*glyphosate*humic_acid + year, data = weed_biomass.2)
m1 <- lm(sqrt(weed_dry_weight) ~ mulch*glyphosate + mulch*humic_acid + glyphosate*humic_acid + year, data = weed_biomass.2)
m2 <- lm(sqrt(weed_dry_weight) ~ mulch*glyphosate + humic_acid + year, data = weed_biomass.2)
m3 <- lm(sqrt(weed_dry_weight) ~ mulch*humic_acid + glyphosate + year, data = weed_biomass.2)
m4 <- lm(sqrt(weed_dry_weight) ~ glyphosate*humic_acid + mulch + year, data = weed_biomass.2)
m5 <- lm(sqrt(weed_dry_weight) ~ glyphosate*humic_acid + glyphosate*mulch + year, data = weed_biomass.2)
m6 <- lm(sqrt(weed_dry_weight) ~ mulch*humic_acid + glyphosate*mulch + year, data = weed_biomass.2)
m7 <- lm(sqrt(weed_dry_weight) ~ mulch*humic_acid + glyphosate*humic_acid + year, data = weed_biomass.2)
m8 <- lm(sqrt(weed_dry_weight) ~ mulch + humic_acid + glyphosate + year, data = weed_biomass.2)
m9 <- lm(sqrt(weed_dry_weight) ~ mulch*humic_acid + year, data = weed_biomass.2)
m10 <- lm(sqrt(weed_dry_weight) ~ mulch*glyphosate + year, data = weed_biomass.2)
m11 <- lm(sqrt(weed_dry_weight) ~ humic_acid*glyphosate + year, data = weed_biomass.2)
m12 <- lm(sqrt(weed_dry_weight) ~ mulch + glyphosate + year, data = weed_biomass.2)
m13 <- lm(sqrt(weed_dry_weight) ~ mulch + humic_acid + year, data = weed_biomass.2)
m14 <- lm(sqrt(weed_dry_weight) ~ humic_acid + glyphosate + year, data = weed_biomass.2)
m15 <- lm(sqrt(weed_dry_weight) ~ mulch + year, data = weed_biomass.2)
m16 <- lm(sqrt(weed_dry_weight) ~ humic_acid + year, data = weed_biomass.2)
m17 <- lm(sqrt(weed_dry_weight) ~ glyphosate + year, data = weed_biomass.2)

plot(resid(m0)~fitted(m0))
plot(resid(m1)~fitted(m1))
plot(resid(m2)~fitted(m2))
plot(resid(m3)~fitted(m3))
plot(resid(m4)~fitted(m4))
plot(resid(m5)~fitted(m5))
plot(resid(m6)~fitted(m6))
plot(resid(m7)~fitted(m7))
plot(resid(m8)~fitted(m8))
plot(resid(m9)~fitted(m9))
plot(resid(m10)~fitted(m10))
plot(resid(m11)~fitted(m11))
plot(resid(m12)~fitted(m12))
plot(resid(m13)~fitted(m13))
plot(resid(m14)~fitted(m14))
plot(resid(m15)~fitted(m15))
plot(resid(m16)~fitted(m16))
plot(resid(m17)~fitted(m17))
# model selection
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17) # including only the models with good to ok fit
#       (Int) gly hmc_acd mlc     yer gly:hmc_acd gly:mlc hmc_acd:mlc gly:hmc_acd:mlc df  logLik  AICc delta weight
# m10 -80.32   +           + 0.04016                   +                              6   6.601  0.8  0.00  0.337 <--
# m12 -80.37   +           + 0.04016                                                  5   4.998  1.4  0.59  0.252
# m2  -80.34   +       +   + 0.04016                   +                              7   6.763  3.3  2.43  0.100
# m8  -80.39   +       +   + 0.04016                                                  6   5.149  3.7  2.90  0.079
# m6  -80.38   +       +   + 0.04016                   +           +                  8   7.615  4.5  3.62  0.055
# m3  -80.43   +       +   + 0.04016                               +                  7   5.945  4.9  4.06  0.044
# m15 -80.46               + 0.04016                                                  4   1.656  5.6  4.77  0.031
# m5  -80.32   +       +   + 0.04016           +       +                              8   6.990  5.7  4.87  0.030
# m4  -80.37   +       +   + 0.04016           +                                      7   5.361  6.1  5.23  0.025
# m1  -80.36   +       +   + 0.04016           +       +           +                  9   7.850  7.0  6.19  0.015
# m7  -80.41   +       +   + 0.04016           +                   +                  8   6.164  7.4  6.52  0.013
# m13 -80.48           +   + 0.04016                                                  5   1.788  7.9  7.01  0.010
# m9  -80.52           +   + 0.04016                               +                  6   2.478  9.1  8.25  0.005
# m0  -80.37   +       +   + 0.04016           +       +           +               + 10   7.970 10.0  9.16  0.003
# m17 -80.68   +             0.04016                                                  4 -21.589 52.1 51.26  0.000
# m16 -80.79           +     0.04016                                                  4 -22.698 54.3 53.48  0.000
# m14 -80.70   +       +     0.04016                                                  5 -21.539 54.5 53.66  0.000
# m11 -80.68   +       +     0.04016           +                                      6 -21.470 57.0 56.14  0.000

anova(m10)
#             Df Sum Sq Mean Sq F value    Pr(>F)    
# mulch             1 4.6273  4.6273 93.2143 2.48e-12 ***
# glyphosate        1 0.3409  0.3409  6.8680  0.01208 *  
# year              1 0.0194  0.0194  0.3899  0.53567    
# mulch:glyphosate  1 0.1475  0.1475  2.9709  0.09196 .  
# Residuals        43 2.1346  0.0496                    

## VISUALISATION OF WEED DRY WT: ####
emm_weed_dry_weight<-emmeans(m10, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#cld 
cld.weed_dry_weight<-cld(emm_weed_dry_weight, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#Visualize the results
ggplot(data.frame(cld.weed_dry_weight), aes(x=mulch, y=response, fill=glyphosate, color=glyphosate))+
  geom_point(data=weed_biomass.2,aes(x=mulch,y=weed_dry_weight, fill=glyphosate, color=glyphosate), , position=position_jitterdodge(dodge.width = 1.0, jitter.width = 0.2), alpha = 0.2, size=2) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), linewidth = 1.5,size=1.2, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_x_discrete(labels=c("mulch"="mulch", "control"="no mulch")) + 
  geom_text(aes(label=.group, y = upper.CL, group=glyphosate),color="black",size=3, vjust = -8,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Weed dry weight") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.weed.dry.wt.pdf", device="pdf", units="mm", dpi=300, width=95, height=100, path="Revisions.ISME.Comm/plots/")

######################################################
######################################################
#### Pack-line fruit Analysis ########################
######################################################
######################################################
################################################
######### Average fruit weight #################
################################################
# check for outliers
ggplot(final_packline.2, aes(y=mean_fruit_weight, x=mulch, color=glyphosate, shape=humic, label = ID)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("2023_CCH8","2025_CCH6")) 
final_packline.2<-final_packline %>% filter(!ID %in% outliers)

m0 <- lm(mean_fruit_weight ~ mulch*glyphosate*humic + Year, data = final_packline.2)
m1 <- lm(mean_fruit_weight ~ mulch*glyphosate + mulch*humic + glyphosate*humic + Year, data = final_packline.2)
m2 <- lm(mean_fruit_weight ~ mulch*glyphosate + humic + Year, data = final_packline.2)
m3 <- lm(mean_fruit_weight ~ mulch*humic + glyphosate + Year, data = final_packline.2)
m4 <- lm(mean_fruit_weight ~ glyphosate*humic + mulch + Year, data = final_packline.2)
m5 <- lm(mean_fruit_weight ~ glyphosate*humic + glyphosate*mulch + Year, data = final_packline.2)
m6 <- lm(mean_fruit_weight ~ mulch*humic + glyphosate*mulch + Year, data = final_packline.2)
m7 <- lm(mean_fruit_weight ~ mulch*humic + glyphosate*humic + Year, data = final_packline.2)
m8 <- lm(mean_fruit_weight ~ mulch + humic + glyphosate + Year, data = final_packline.2)
m9 <- lm(mean_fruit_weight ~ mulch*humic + Year, data = final_packline.2)
m10 <- lm(mean_fruit_weight ~ mulch*glyphosate + Year, data = final_packline.2) # <--
m11 <- lm(mean_fruit_weight ~ humic*glyphosate + Year, data = final_packline.2)
m12 <- lm(mean_fruit_weight ~ mulch + glyphosate + Year, data = final_packline.2)
m13 <- lm(mean_fruit_weight ~ mulch + humic + Year, data = final_packline.2)
m14 <- lm(mean_fruit_weight ~ humic + glyphosate + Year, data = final_packline.2)
m15 <- lm(mean_fruit_weight ~ mulch + Year, data = final_packline.2)
m16 <- lm(mean_fruit_weight ~ humic + Year, data = final_packline.2)
m17 <- lm(mean_fruit_weight ~ glyphosate + Year, data = final_packline.2)

plot(resid(m0)~fitted(m0))
plot(resid(m1)~fitted(m1))
plot(resid(m2)~fitted(m2))
plot(resid(m3)~fitted(m3))
plot(resid(m4)~fitted(m4))
plot(resid(m5)~fitted(m5))
plot(resid(m6)~fitted(m6))
plot(resid(m7)~fitted(m7))
plot(resid(m8)~fitted(m8))
plot(resid(m9)~fitted(m9))
plot(resid(m10)~fitted(m10))
plot(resid(m11)~fitted(m11))
plot(resid(m12)~fitted(m12))
plot(resid(m13)~fitted(m13))
plot(resid(m14)~fitted(m14))
plot(resid(m15)~fitted(m15))
plot(resid(m16)~fitted(m16))
plot(resid(m17)~fitted(m17))
# model selection
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17) # including only the models with good to ok fit
#       (Int) gly hmc mlc   Yer gly:hmc gly:mlc hmc:mlc gly:hmc:mlc df   logLik   AICc delta weight
# m10 -18170   +       + 9.017               +                      6 -701.971 1416.4  0.00  0.275 <--
# m15 -18170           + 9.017                                      4 -704.369 1417.0  0.55  0.208
# m2  -18170   +   +   + 9.017               +                      7 -701.906 1418.4  2.02  0.100
# m5  -18170   +   +   + 9.017       +       +                      8 -700.865 1418.5  2.12  0.095
# m13 -18170       +   + 9.017                                      5 -704.296 1418.9  2.52  0.078
# m12 -18170   +       + 9.017                                      5 -704.368 1419.1  2.66  0.073
# m6  -18170   +   +   + 9.017               +       +              8 -701.903 1420.6  4.20  0.034
# m1  -18170   +   +   + 9.017       +       +       +              9 -700.864 1420.7  4.33  0.032
# m9  -18170       +   + 9.017                       +              6 -704.291 1421.0  4.64  0.027
# m8  -18170   +   +   + 9.017                                      6 -704.295 1421.0  4.65  0.027
# m4  -18170   +   +   + 9.017       +                              7 -703.315 1421.2  4.84  0.024
# m0  -18170   +   +   + 9.017       +       +       +           + 10 -700.857 1422.9  6.54  0.010
# m3  -18170   +   +   + 9.017                       +              7 -704.290 1423.2  6.79  0.009
# m7  -18170   +   +   + 9.017       +               +              8 -703.311 1423.4  7.02  0.008
# m16 -18160       +     9.017                                      4 -712.111 1432.4 16.04  0.000
# m17 -18160   +         9.017                                      4 -712.159 1432.5 16.13  0.000
# m14 -18160   +   +     9.017                                      5 -712.107 1434.5 18.14  0.000
# m11 -18170   +   +     9.017       +                              6 -711.146 1434.8 18.35  0.000

anova(m10)
#                         Df Sum Sq Mean Sq F value    Pr(>F)    
#mulch              1  1578.8  1578.8  16.2234  8.21e-05 ***
#glyphosate         1     0.3     0.3   0.0029   0.95736    
#Year               1 10245.4 10245.4 105.2797 < 2.2e-16 ***
#mulch:glyphosate   1   459.9   459.9   4.7259   0.03098 *  
#Residuals        185 18003.4    97.3                       

## VISUALISATION OF AVG FRUIT WT: ####
emm_mean_fruit_weight<-emmeans(m10, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#cld 
cld.mean_fruit_weight<-cld(emm_mean_fruit_weight, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#Visualize the results
#Visualize the results
ggplot(data.frame(cld.mean_fruit_weight), aes(x=mulch, y=emmean, fill=glyphosate, color=glyphosate))+
  geom_point(data=final_packline.2,aes(x=mulch,y=mean_fruit_weight, fill=glyphosate, color=glyphosate), , position=position_jitterdodge(dodge.width = 1.0, jitter.width = 0.2), alpha = 0.2, size=2) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), linewidth = 1.5,size=1.2, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_x_discrete(labels=c("mulch"="mulch", "control"="no mulch")) + 
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  geom_text(aes(label=.group, y = upper.CL, group=glyphosate),color="black",size=3, vjust = -8,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Mean fruit weight (g)") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.fruit.wt.2023_2025.pdf", device="pdf", units="mm", dpi=300, width=98, height=100, path="Revisions.ISME.Comm/plots/")

################################################
######### Yield (number of fruit) ##############
################################################
final_packline # data frame 
# check for outliers
ggplot(final_packline.2, aes(y=N_of_Fruits, x=mulch, color=glyphosate, shape=humic, label = ID)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(0.1)) + theme_classic() + geom_text(position = position_jitterdodge(0.1),check_overlap = FALSE)
# remove outliers
outliers<-as.character(c("2025_CCH6","2025_CCH4","2025_CGC6", "2025_CGC4")) 
final_packline.2<-final_packline %>% filter(!ID %in% outliers)

m0 <- lm(log(N_of_Fruits) ~ mulch*glyphosate*humic + Year, data = final_packline.2)
m1 <- lm(log(N_of_Fruits) ~ mulch*glyphosate + mulch*humic + glyphosate*humic + Year, data = final_packline.2)
m2 <- lm(log(N_of_Fruits) ~ mulch*glyphosate + humic + Year, data = final_packline.2)
m3 <- lm(log(N_of_Fruits) ~ mulch*humic + glyphosate + Year, data = final_packline.2)
m4 <- lm(log(N_of_Fruits) ~ glyphosate*humic + mulch + Year, data = final_packline.2)
m5 <- lm(log(N_of_Fruits) ~ glyphosate*humic + glyphosate*mulch + Year, data = final_packline.2)
m6 <- lm(log(N_of_Fruits) ~ mulch*humic + glyphosate*mulch + Year, data = final_packline.2)
m7 <- lm(log(N_of_Fruits) ~ mulch*humic + glyphosate*humic + Year, data = final_packline.2)
m8 <- lm(log(N_of_Fruits) ~ mulch + humic + glyphosate + Year, data = final_packline.2)
m9 <- lm(log(N_of_Fruits) ~ mulch*humic + Year, data = final_packline.2)
m10 <- lm(log(N_of_Fruits) ~ mulch*glyphosate + Year, data = final_packline.2) # <--
m11 <- lm(log(N_of_Fruits) ~ humic*glyphosate + Year, data = final_packline.2)
m12 <- lm(log(N_of_Fruits) ~ mulch + glyphosate + Year, data = final_packline.2)
m13 <- lm(log(N_of_Fruits) ~ mulch + humic + Year, data = final_packline.2)
m14 <- lm(log(N_of_Fruits) ~ humic + glyphosate + Year, data = final_packline.2)
m15 <- lm(log(N_of_Fruits) ~ mulch + Year, data = final_packline.2)
m16 <- lm(log(N_of_Fruits) ~ humic + Year, data = final_packline.2)
m17 <- lm(log(N_of_Fruits) ~ glyphosate + Year, data = final_packline.2)

plot(resid(m0)~fitted(m0))
text(x = fitted(m0), y = resid(m0), labels = final_packline.2$ID)
plot(resid(m1)~fitted(m1))
plot(resid(m2)~fitted(m2))
plot(resid(m3)~fitted(m3))
plot(resid(m4)~fitted(m4))
plot(resid(m5)~fitted(m5))
plot(resid(m6)~fitted(m6))
plot(resid(m7)~fitted(m7))
plot(resid(m8)~fitted(m8))
plot(resid(m9)~fitted(m9))
plot(resid(m10)~fitted(m10))
plot(resid(m11)~fitted(m11))
plot(resid(m12)~fitted(m12))
plot(resid(m13)~fitted(m13))
plot(resid(m14)~fitted(m14))
plot(resid(m15)~fitted(m15))
plot(resid(m16)~fitted(m16))
plot(resid(m17)~fitted(m17))
# model selection
MuMIn::model.sel(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17) # including only the models with good to ok fit
#      (Int) gly hmc mlc     Yer gly:hmc gly:mlc hmc:mlc gly:hmc:mlc df   logLik  AICc delta weight
# m12  1212   +       + -0.5954                                      5 -143.473 297.3  0.00  0.195 <--
# m8   1212   +   +   + -0.5954                                      6 -142.629 297.7  0.45  0.156
# m10  1212   +       + -0.5954               +                      6 -142.851 298.2  0.89  0.125
# m15  1212           + -0.5954                                      4 -145.190 298.6  1.32  0.101
# m2   1212   +   +   + -0.5954               +                      7 -142.050 298.7  1.45  0.095
# m13  1212       +   + -0.5954                                      5 -144.438 299.2  1.93  0.074
# m3   1212   +   +   + -0.5954                       +              7 -142.442 299.5  2.23  0.064
# m4   1212   +   +   + -0.5953       +                              7 -142.626 299.9  2.60  0.053
# m6   1212   +   +   + -0.5954               +       +              8 -141.883 300.6  3.30  0.037
# m5   1212   +   +   + -0.5953       +       +                      8 -142.046 300.9  3.63  0.032
# m9   1212       +   + -0.5954                       +              6 -144.291 301.1  3.77  0.030
# m7   1212   +   +   + -0.5953       +               +              8 -142.438 301.7  4.41  0.022
# m1   1212   +   +   + -0.5953       +       +       +              9 -141.880 302.8  5.50  0.012
# m0   1213   +   +   + -0.5958       +       +       +           + 10 -141.844 304.9  7.67  0.004
# m17  1235   +         -0.6067                                      4 -158.702 325.6 28.35  0.000
# m14  1235   +   +     -0.6067                                      5 -157.986 326.3 29.03  0.000
# m16  1235       +     -0.6067                                      4 -159.522 327.3 29.99  0.000
# m11  1235   +   +     -0.6069       +                              6 -157.980 328.4 31.15  0.000

anova(m12)
#                         Df Sum Sq Mean Sq F value    Pr(>F)    
# mulch        1 10.782  10.782  38.5237 3.562e-09 ***
# glyphosate   1  0.949   0.949   3.3918   0.06715 .  
# Year         1 43.111  43.111 154.0376 < 2.2e-16 ***
# Residuals  182 50.937   0.280                                           

## VISUALISATION OF # OF FRUIT: ####
emm_N_of_Fruits<-emmeans(m12, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
emm_N_of_Fruits.m<-emmeans(m12, specs=~mulch, type="response") # calculate estimated marginal means
#cld 
cld.N_of_Fruits<-cld(emm_N_of_Fruits, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#Visualize the results
ggplot(data.frame(cld.N_of_Fruits), aes(x=mulch, y=response, fill=glyphosate, color=glyphosate))+
  geom_point(data=final_packline.2,aes(x=mulch,y=N_of_Fruits, fill=glyphosate, color=glyphosate), , position=position_jitterdodge(dodge.width = 1.0, jitter.width = 0.2), alpha = 0.2, size=2) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), linewidth = 1.5,size=1.2, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_fill_manual(values=c("#25533f", "#aacc96"),limits=c("glyphosate","control"), labels = c("glyphosate","no glyphosate")) +
  scale_x_discrete(labels=c("mulch"="mulch", "control"="no mulch")) + 
  #scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  geom_text(aes(label=.group, y = upper.CL, group=glyphosate),color="black",size=3, vjust = -8,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Yield (# of Fruit)") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.fruit.number_yield.2023_2025.pdf", device="pdf", units="mm", dpi=300, width=98, height=100, path="Revisions.ISME.Comm/plots/")

