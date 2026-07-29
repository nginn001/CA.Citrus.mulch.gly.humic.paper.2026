# Script created February 20, 2026 by Nichole Ginnan, Assistant Project Scientist (nginn001@ucr.edu)
# Last updated July 18, 2026
# Greenhouse soil microbiome functional effects on citrus seedlings data analysis
# ECDRE Project, Lead PI: Caroline Roper, funded by the USDA
# University of California, Riverside, California | Plant Research Facility 1 Greenhouse
# Include seeding establishment and root and shoot growth data
###############################################################################
# Load Libraries ####
library(tidyverse)
library(ggplot2)
library(scales)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(brglm2)
library(readxl)
library(grid)
library(MuMIn)
###################################################
##### Greenhouse experiment final measurements #####
###################################################
#### Load data ####
plant<-read_excel("/Users/nicholeginnan/Desktop/Sup.Table.S4. Greenhouse Plant Data.xlsx",col_names = TRUE)
plant.2 <- plant %>% filter(end_status != "none")
## established #### binary response model ####
plant$established <- ifelse(plant$end_status == "none", 1, 0)
plant_A <- subset(plant, Active.or.killed == "A")

m_establish_firth.0 <- glm(established ~ mulch * gly *humic,data = plant_A, family = binomial("logit"),method = "brglmFit")
m_establish_firth.1 <- glm(established ~ mulch*gly + gly*humic + mulch*humic ,data = plant_A, family = binomial("logit"),method = "brglmFit")
m_establish_firth.2 <- glm(established ~ mulch*gly + gly*humic,data = plant_A, family = binomial("logit"),method = "brglmFit")
m_establish_firth.3 <- glm(established ~ mulch*gly + humic,data = plant_A, family = binomial("logit"),method = "brglmFit")
m_establish_firth.4 <- glm(established ~ mulch + gly + humic,data = plant_A, family = binomial("logit"),method = "brglmFit")

MuMIn::model.sel(m_establish_firth.0,m_establish_firth.1,m_establish_firth.2,m_establish_firth.3,m_establish_firth.4)
#                      Int) gly hmc mlc gly:hmc gly:mlc hmc:mlc gly:hmc:mlc df  logLik  AICc delta weight
#m_establish_firth.3 2.358   +   +   +               +                      5 -47.292 105.1  0.00  0.373 
#m_establish_firth.2 2.667   +   +   +       +       +                      6 -46.336 105.4  0.29  0.323
#m_establish_firth.0 1.758   +   +   +       +       +       +           +  8 -44.709 106.6  1.55  0.172 <--
#m_establish_firth.1 2.175   +   +   +       +       +       +              7 -46.113 107.2  2.08  0.132
#m_establish_firth.4 3.482   +   +   +                                      4 -54.242 116.8 11.73  0.001

# the full model is <2 different than the top. ranked model (meaning they are practically the same). I am moving forward with the full model

summary(m_establish_firth.0)
anova(m_establish_firth.0, test = "Chisq")
#Analysis of Deviance Table
#.                Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
#NULL                              127    174.308              
#mulch            1   50.010       126    124.298 1.529e-12 ***
#gly              1   10.242       125    114.055 0.0013725 ** 
#humic            1    5.572       124    108.484 0.0182537 *  
#mulch:gly        1   13.899       123     94.585 0.0001929 ***
#mulch:humic      1    1.489       122     93.096 0.2223886    
#gly:humic        1    0.869       121     92.227 0.3513010    
#mulch:gly:humic  1    2.810       120     89.417 0.0936870 . 
emm <- emmeans(m_establish_firth.0, ~ mulch:gly:humic, type = "response")
pairs(emm, adjust = "fdr")
cld_soil <- cld(emm, alpha = 0.05, adjust = "fdr", Letters = LETTERS)
#### plot % plants established ####
established.cnt<-plant.2%>%count(mulch,gly,humic,Active.or.killed)
established.cnt <- established.cnt %>%mutate(Treatment = paste(mulch, gly, humic, sep = "_"))

ggplot(established.cnt, aes(x = Treatment,y = (n/16), fill=mulch)) +
  geom_col(position = position_dodge(width = 0.8),color="black") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(.~Active.or.killed,labeller = labeller(Active.or.killed = c("A" = "Active", "K" = "Autoclaved"))) +
  scale_x_discrete(labels = c("no_no_no"="CCC","no_no_humic"="CCH","no_gly_no"="CGC","no_gly_humic"="CGH",
                              "mulch_no_no"="MCC","mulch_no_humic"="MCH","mulch_gly_no"="MGC","mulch_gly_humic"="MGH"),
                   limits = c("no_no_no","no_no_humic","no_gly_no","no_gly_humic",
                              "mulch_no_no","mulch_no_humic","mulch_gly_no","mulch_gly_humic"))+
  scale_fill_manual(values=c("#0065a2", "#00b0ba"),limits=c("mulch","no"),labels = c("mulch","no mulch")) +
  labs(x = "Soil Microbiota", y = "Plants established (%)") + 
  theme_classic() +
  theme(legend.title = element_text(size=10,family="sans"), 
        legend.text=element_text(size=10,family="sans"),
        legend.position="right",
        legend.box.spacing = unit(0.0, "pt"),
        legend.spacing.x = unit(5.0, 'pt'),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        panel.background = element_rect(colour = "black",linewidth=1),
        #strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.3, "lines"),
        axis.text.x= element_text(colour="black", size=8, family="sans", angle=25, hjust=1), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "GH.bar.established.legend.pdf", device="pdf", units="mm", dpi=300, width=150, height=80, path="plots/")

################################################
##### Root biomass #############################
################################################
# Plot raw data to identify major outliers
ggplot(plant.2, aes(x=Soil, y= rt_mass_g,color=Active.or.killed)) + 
  geom_boxplot(outliers = FALSE) +
  geom_point(alpha=0.2, position=position_jitterdodge(0.2), size=2) +
  geom_text(aes(label=SampleID), position=position_jitterdodge(0.2), 
            vjust=-0.5, size=3, alpha=0.7, check_overlap=TRUE) +
  theme_classic()
# remove major outliers from the dataset to improve model fit
outliers<-as.character(c("D1_27","C1_23")) 
plant.3<-plant.2 %>% filter(!SampleID %in% outliers)# model
# Root mass model
#m_root <- lmer((rt_mass_g) ~ mulch*gly*Active.or.killed*humic + (1|Block),data = plant.3) #over fitting
#ranova(m_root)
m_root.0 <- lm(rt_mass_g ~ mulch*gly*Active.or.killed*humic,data = plant.3)
m_root.1 <- lm(rt_mass_g ~ mulch*gly*Active.or.killed+humic,data = plant.3) # reduced model
MuMIn::model.sel(m_root.0,m_root.1)
#           (Int) Act.or.kll gly hmc mlc Act.or.kll:gly Act.or.kll:hmc Act.or.kll:mlc gly:hmc gly:mlc hmc:mlc Act.or.kll:gly:hmc Act.or.kll:gly:mlc Act.or.kll:hmc:mlc gly:hmc:mlc Act.or.kll:gly:hmc:mlc df  logLik   AICc delta weight
#m_root.1 0.07822          +   +   +   +              +                             +               +                                             +                                                       10 243.024 -464.7  0.00  0.987
#m_root.0 0.14000          +   +   +   +              +              +              +       +       +       +                  +                  +                  +           +                      + 17 246.928 -456.1  8.67  0.013
plot(m_root.1)
qqnorm(resid(m_root.1));qqline(resid(m_root.1))
anova(m_root.1)
#.                             Df  Sum Sq Mean Sq F value    Pr(>F)    
# mulch                        1 0.00121 0.00121  0.2922  0.589529    
# gly                          1 0.00039 0.00039  0.0935  0.760162    
# Active.or.killed             1 0.34150 0.34150 82.4715 2.563e-16 ***
# humic                        1 0.05486 0.05486 13.2494  0.000361 ***
# mulch:gly                    1 0.00124 0.00124  0.2996  0.584876    
# mulch:Active.or.killed       1 0.03178 0.03178  7.6746  0.006220 ** 
# gly:Active.or.killed         1 0.00395 0.00395  0.9535  0.330222    
# mulch:gly:Active.or.killed   1 0.01669 0.01669  4.0294  0.046290 *  
# Residuals                  171 0.70809 0.00414                      
## Plot Root mass ####
emm_main<-emmeans(m_root.1, specs=~humic, type="response") # significant
emm_interact<-emmeans(m_root.1, specs=~mulch:gly:Active.or.killed, type="response") # significant
cld_main<-cld(emm_main, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld_interact<-cld(emm_interact, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# Plot the data
# main effects #
ggplot(cld_main, aes(x = humic, y = emmean, color=humic)) +
  geom_point(data=plant.3,mapping=aes(x=humic, y=rt_mass_g), alpha=0.2, position=position_jitter(width = 0.2,height=0.03), size=3) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1, linewidth = 1.1) +
  theme_classic(base_size = 10) +
  #facet_grid(.~Active.or.killed,labeller = labeller(Active.or.killed = c("A" = "Active", "K" = "Autoclaved"))) +
  geom_text(aes(label=.group, y = upper.CL+0.1),color="black",size=4, vjust=-3,hjust=0.5) +
  scale_x_discrete(labels=c("humic"="HA","no"="no HA")) +
  scale_color_manual(values=c("#704776", "#f0be39"),limits=c("humic","no"),labels = c("HA","no HA")) +
  xlab("Treatment")+
  ylab("Dry Root Mass (g)")+
    theme(legend.title = element_text(size=10,family="sans"), 
      legend.text=element_text(size=10,family="sans"),
      legend.position="none",
      legend.box.spacing = unit(0.0, "pt"),
      legend.spacing.x = unit(5.0, 'pt'),
      axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
      strip.text.x = element_text(size = 10),
      axis.title.y=element_text(size=10, family="sans",vjust = 1), 
      panel.background = element_rect(colour = "black",linewidth=1),
      #strip.background=element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.spacing=unit(0.5, "lines"),
      axis.text.x= element_text(colour="black", size=8, family="sans", angle=0, hjust=0.5), 
      axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "GH.root_mass.humic.marginal.effects.emm.pdf", device="pdf", units="mm", dpi=300, width=40, height=100, path="/Users/nicholeginnan/Documents/UCR- Current/CA.citrus.paper/Revisions.ISME.Comm/plots/")

# Interactions #
ggplot(cld_interact, aes(x = mulch, y = emmean,color=gly)) +
  geom_point(data=plant.3,mapping=aes(x=mulch, y=rt_mass_g,color=gly), alpha=0.2,
             position=position_jitterdodge(dodge.width = 1,jitter.width = 0.2,jitter.height = 0.03), size=3) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1, linewidth = 1.1,position = position_dodge2(1.0)) +
  theme_classic(base_size = 10) +
  facet_grid(.~Active.or.killed,labeller = labeller(Active.or.killed = c("A" = "Active", "K" = "Autoclaved"))) +
  geom_text(aes(label=.group, y = upper.CL+0.1),color="black",size=4, vjust=-0.5,hjust=0.5,position = position_dodge2(1.0)) +
  scale_color_manual(values=c("#25533f", "#aacc96"),limits=c("gly","no"),labels = c("glyphosate","no glyphosate")) +
  xlab("Treatment")+
  ylab("Dry Root Mass (g)")+
  theme(legend.title = element_text(size=10,family="sans"), 
        legend.text=element_text(size=10,family="sans"),
        legend.position="none",
        legend.box.spacing = unit(0.0, "pt"),
        legend.spacing.x = unit(5.0, 'pt'),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        panel.background = element_rect(colour = "black",linewidth=1),
        #strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"),
        axis.text.x= element_text(colour="black", size=8, family="sans", angle=0, hjust=0.5), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "GH.root_mass.mulch.gly.active.emm.pdf", device="pdf", units="mm", dpi=300, width=80, height=100, path="plots/")

########################
### Shoot biomass #####
#######################
# Plot raw data to identify major outliers
ggplot(plant.3, aes(x=Soil, y= sht_mass_g,color=Active.or.killed)) + 
  geom_boxplot(outliers = FALSE) +
  geom_point(alpha=0.2, position=position_jitterdodge(0.2), size=2) +
  geom_text(aes(label=SampleID), position=position_jitterdodge(0.2), 
            vjust=-0.5, size=3, alpha=0.7, check_overlap=TRUE) +
  theme_classic()
# remove major outliers from the dataset to improve model fit
outliers<-as.character(c("D1_27","D3_25","C2_3","D3_23")) 
plant.3<-plant.2 %>% filter(!SampleID %in% outliers)# model
# model
m_sht.0 <- lmer(sht_mass_g ~ mulch*gly*Active.or.killed*humic+ (1|Block),data = plant.3, REML = FALSE) 
# Remove only the four-way interaction; retain all three-way interactions
m_sht.1 <- lmer(sht_mass_g ~ mulch*gly*Active.or.killed +mulch*gly*humic+ gly*Active.or.killed*humic+ mulch*humic*Active.or.killed+ (1|Block),data = plant.3, REML = FALSE)
# Retain the significant three-way interactions
m_sht.2 <- lmer(sht_mass_g ~mulch*Active.or.killed*humic+gly +(1 | Block),data = plant.3, REML = FALSE)
# Retain the significant three-way interactions and drop no-sig variable
m_sht.3 <- lmer(sht_mass_g ~mulch*Active.or.killed*humic+(1 | Block),data = plant.3, REML = FALSE)
# Retain all two-way interactions, but no three- or four-way interactions
m_sht.4 <- lmer(sht_mass_g ~mulch*gly + mulch*Active.or.killed + mulch*humic + gly*Active.or.killed + humic*Active.or.killed+ gly*humic+(1 | Block),data = plant.3, REML = FALSE)
# Additive model
m_sht.5 <- lmer(sht_mass_g ~ mulch + gly + Active.or.killed + humic +(1 | Block),data = plant.3, REML = FALSE)

MuMIn::model.sel(m_sht.0,m_sht.1,m_sht.2,m_sht.3,m_sht.4,m_sht.5)
#.         (Int) Act.or.kll gly hmc mlc Act.or.kll:gly Act.or.kll:hmc Act.or.kll:mlc gly:hmc gly:mlc hmc:mlc Act.or.kll:gly:hmc Act.or.kll:gly:mlc Act.or.kll:hmc:mlc gly:hmc:mlc Act.or.kll:gly:hmc:mlc df  logLik   AICc delta weight
#m_sht.3 0.11640          +       +   +                             +              +                       +                                                        +                                    10 290.883 -560.4  0.00  0.676 <--
#m_sht.2 0.11720          +   +   +   +                             +              +                       +                                                        +                                    11 290.934 -558.3  2.17  0.228
#m_sht.1 0.13200          +   +   +   +              +              +              +       +       +       +                  +                  +                  +           +                        17 296.063 -554.3  6.15  0.031
#m_sht.5 0.11120          +   +   +   +                                                                                                                                                                   7 284.333 -554.0  6.44  0.027
#m_sht.4 0.09126          +   +   +   +              +              +              +       +       +       +                                                                                             13 290.768 -553.3  7.13  0.019
#m_sht.0 0.14620          +   +   +   +              +              +              +       +       +       +                  +                  +                  +           +                      + 18 296.748 -553.2  7.25  0.018
ranova(m_sht.3) # block is significant
plot(m_sht.3)
qqnorm(resid(m_sht.3));qqline(resid(m_sht.3))
anova(m_sht.3)
#.                                   Df   Sum Sq Mean Sq F value    Pr(>F)    
# mulch                        0.014034 0.014034     1 165.30  6.4381   0.01210 *  
# Active.or.killed             0.051825 0.051825     1 165.38 23.7741 2.526e-06 ***
# humic                        0.038529 0.038529     1 166.13 17.6750 4.272e-05 ***
# mulch:Active.or.killed       0.008686 0.008686     1 165.67  3.9844   0.04756 *  
# mulch:humic                  0.002827 0.002827     1 165.84  1.2970   0.25640    
# Active.or.killed:humic       0.004746 0.004746     1 166.57  2.1773   0.14195    
# mulch:Active.or.killed:humic 0.013430 0.013430     1 166.20  6.1609   0.01405 * 
## Plot Sht mass ####
emm_main<-emmeans(m_sht.3, specs=~mulch:humic:Active.or.killed, type="response") # significant
cld_main<-cld(emm_main, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
#cld_interact<-cld(emm_interact, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# Plot the data
# main effects #
ggplot(cld_main, aes(x = mulch, y = response,color=humic)) +
  geom_point(data=plant.3,mapping=aes(x=mulch, y=sht_mass_g,color=humic), alpha=0.2,
             position=position_jitterdodge(dodge.width = 1,jitter.width = 0.2,jitter.height = 0.01), size=3) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1, linewidth = 1.1,position = position_dodge2(1.0)) +
  theme_classic(base_size = 10) +
  facet_grid(.~Active.or.killed,labeller = labeller(Active.or.killed = c("A" = "Active", "K" = "Autoclaved"))) +
  geom_text(aes(label=.group, y = upper.CL+0.08),color="black",size=4, vjust=-0.5,hjust=0.5,position = position_dodge2(1.0)) +
  scale_color_manual(values=c("#704776", "#f0be39"),limits=c("humic","no"),labels = c("humic","no humic")) +
  ylab("Dry Shoot Mass (g)")+
  xlab("Treatment")+
  theme(legend.title = element_text(size=10,family="sans"), 
        legend.text=element_text(size=10,family="sans"),
        legend.position="none",
        legend.box.spacing = unit(0.0, "pt"),
        legend.spacing.x = unit(5.0, 'pt'),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        panel.background = element_rect(colour = "black",linewidth=1),
        #strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"),
        axis.text.x= element_text(colour="black", size=8, family="sans", angle=0, hjust=0.5), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "GH.sht_mass.mulch.humic.active.emm.pdf", device="pdf", units="mm", dpi=300, width=80, height=100, path="plots/")

########################
### Shoot height #####
#######################
# Plot raw data to identify major outliers
ggplot(plant.3, aes(x=Soil, y= ht_cm,color=Active.or.killed)) + 
  geom_boxplot(outliers = FALSE) +
  geom_point(alpha=0.2, position=position_jitterdodge(0.2), size=2) +
  geom_text(aes(label=SampleID), position=position_jitterdodge(0.2), 
            vjust=-0.5, size=3, alpha=0.7, check_overlap=TRUE) +
  theme_classic()
# remove major outliers from the dataset to improve model fit
outliers<-as.character(c("D3_25","C2_3","C4_9","C2_3","D2_7")) 
plant.3<-plant.2 %>% filter(!SampleID %in% outliers)# model
# model
m_sht.0 <- lmer(ht_cm ~ mulch*Active.or.killed*humic+gly + (1|Block),data = plant.3, REML=FALSE)
# Remove only the four-way interaction; retain all three-way interactions
m_sht.1 <- lmer(ht_cm ~ mulch*gly*Active.or.killed +mulch*gly*humic+ gly*Active.or.killed*humic+ mulch*humic*Active.or.killed+ (1|Block),data = plant.3, REML = FALSE)
# Retain the significant three-way interactions
m_sht.2 <- lmer(ht_cm ~mulch*Active.or.killed*humic+gly +(1 | Block),data = plant.3, REML = FALSE)
# Retain the significant three-way interactions and drop no-sig variable
m_sht.3 <- lmer(ht_cm ~mulch*Active.or.killed*humic+(1 | Block),data = plant.3, REML = FALSE)
# Retain all two-way interactions, but no three- or four-way interactions
m_sht.4 <- lmer(ht_cm ~mulch*gly + mulch*Active.or.killed + mulch*humic + gly*Active.or.killed + humic*Active.or.killed+ gly*humic+(1 | Block),data = plant.3, REML = FALSE)
# Additive model
m_sht.5 <- lmer(ht_cm ~ mulch + gly + Active.or.killed + humic +(1 | Block),data = plant.3, REML = FALSE)

MuMIn::model.sel(m_sht.0,m_sht.1,m_sht.2,m_sht.3,m_sht.4,m_sht.5)
#         (Int) Act.or.kll gly hmc mlc Act.or.kll:hmc Act.or.kll:mlc hmc:mlc Act.or.kll:hmc:mlc Act.or.kll:gly gly:hmc gly:mlc Act.or.kll:gly:hmc Act.or.kll:gly:mlc gly:hmc:mlc df   logLik  AICc delta weight
#m_sht.3 7.141          +       +   +              +              +       +                  +                                                                                  10 -296.871 615.1  0.00  0.554 <--
#m_sht.0 7.181          +   +   +   +              +              +       +                  +                                                                                  11 -296.654 616.9  1.84  0.221
#m_sht.2 7.181          +   +   +   +              +              +       +                  +                                                                                  11 -296.654 616.9  1.84  0.221
#m_sht.1 7.188          +   +   +   +              +              +       +                  +              +       +       +                  +                  +           + 17 -293.715 625.3 10.19  0.003
#m_sht.4 5.473          +   +   +   +              +              +       +                                 +       +       +                                                   13 -303.967 636.2 21.09  0.000
#m_sht.5 6.741          +   +   +   +                                                                                                                                            7 -313.256 641.2 26.11  0.000
ranova(m_sht.3) # block is significant
plot(m_sht.3)
qqnorm(resid(m_sht.3));qqline(resid(m_sht.3))
anova(m_sht.3)
#.                             Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
# mulch                        18.897  18.897     1 165.24 11.7604  0.000764 ***
# Active.or.killed             34.333  34.333     1 165.17 21.3670 7.610e-06 ***
# humic                        47.157  47.157     1 166.03 29.3481 2.094e-07 ***
# mulch:Active.or.killed       14.585  14.585     1 165.72  9.0768  0.002995 ** 
# mulch:humic                  14.294  14.294     1 165.94  8.8956  0.003290 ** 
# Active.or.killed:humic        5.543   5.543     1 166.37  3.4499  0.065022 .  
# mulch:Active.or.killed:humic 29.380  29.380     1 166.23 18.2844 3.200e-05 ***
## Plot Sht height ####
emm_main<-emmeans(m_sht.3, specs=~mulch:humic:Active.or.killed, type="response") # significant
cld_main<-cld(emm_main, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
#cld_interact<-cld(emm_interact, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
# Plot the data
# main effects #
ggplot(cld_main, aes(x = mulch, y = response,color=humic)) +
  geom_point(data=plant.3,mapping=aes(x=mulch, y=ht_cm,color=humic), alpha=0.2,
             position=position_jitterdodge(dodge.width = 1,jitter.width = 0.2,jitter.height = 0.01), size=3) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1, linewidth = 1.1,position = position_dodge2(1.0)) +
  theme_classic(base_size = 10) +
  facet_grid(.~Active.or.killed,labeller = labeller(Active.or.killed = c("A" = "Active", "K" = "Autoclaved"))) +
  geom_text(aes(label=.group, y = upper.CL+0.08),color="black",size=4, vjust=-0.5,hjust=0.5,position = position_dodge2(1.0)) +
  scale_color_manual(values=c("#704776", "#f0be39"),limits=c("humic","no"),labels = c("humic","no humic")) +
  ylab("Shoot Height (cm)")+
  xlab("Treatment")+
  theme(legend.title = element_text(size=10,family="sans"), 
        legend.text=element_text(size=10,family="sans"),
        legend.position="none",
        legend.box.spacing = unit(0.0, "pt"),
        legend.spacing.x = unit(5.0, 'pt'),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x = element_text(size = 10),
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        panel.background = element_rect(colour = "black",linewidth=1),
        #strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"),
        axis.text.x= element_text(colour="black", size=8, family="sans", angle=0, hjust=0.5), 
        axis.text.y= element_text(colour="black", size=8, family="sans"))
ggsave(filename= "GH.sht_height.mulch.humic.active.emm.pdf", device="pdf", units="mm", dpi=300, width=80, height=100, path="plots/")

