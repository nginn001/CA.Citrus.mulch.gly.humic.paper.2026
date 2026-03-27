
# Script created November 20, 2025 by Riley Jones, PhD Student
# Last updated March 26, 2026 by Nichole Ginnan
# Analysis to understand the impacts of mulch, glyphosate, and humic acid, or field tree traits and weed growth
# ECDRE Project, Lead PI: Caroline Roper, funded by the USDA
# Lindcove Research and Extension Center, Exeter, California | 91C plot experiment
# Includes LICOR, weed biomass, and fruit packline data
###############################################################################
# Load Libraries ####
library(bestNormalize)
library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)
library(janitor)

#### Load data ############
licor_cleaned<- read_excel("Sup.Table.S3. Field Plant Data.xlsx", sheet = "Assimilation_and_transpiration")
final_packline<- read_excel("Sup.Table.S3. Field Plant Data.xlsx", sheet = "Fruit_weight")
weed_biomass<- read_excel("Sup.Table.S3. Field Plant Data.xlsx", sheet = "Weed_biomass")
set.seed(4444)

######################################################
######################################################
#### LI-COR Analysis #################################
######################################################
######################################################

#########################
####### Licor E #########
#########################
#use bestNormalize:
normalized.licor.E <- bestNormalize(licor_cleaned$licor_e)
#plot normalized residuals:
# Step 1: Run bestNormalize on your response
normalized.licor.e <- bestNormalize(licor_cleaned$licor_e)
# Step 2: Add the transformed variable to your dataset
licor_cleaned$e_trans <- predict(normalized.licor.e)
# Step 3: Fit a linear mixed model
#model_trans_licor_e <- lmer(e_trans ~ mulch*humic_acid*glyphosate + (1 | year), data = licor_cleaned)
#remove humic acid bc not significant
model_trans_licor_e <- lmer(e_trans ~ mulch*glyphosate + (1 | year), data = licor_cleaned)
# Step 4: Plot residuals vs fitted values
plot(fitted(model_licor_e), resid(model_licor_e),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")
#step 5: stats
anova(model_trans_licor_e) # significant
ranova(model_trans_licor_e) # significant
#####
#VISUALISATION OF LICOR E:
emm_licor_e<-emmeans(model_trans_licor_e, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#because I used best normalize to create the model, have to do the back transformation by hand using this code:
# Back-transform the EMMs and their CIs
emm_licor_e_back_transformed <- as.data.frame(emm_licor_e) %>%
  mutate(response = predict(normalized.licor.e, emmean, inverse = TRUE),
         bt_lowerCL= predict(normalized.licor.e, lower.CL, inverse = TRUE),
         bt_upperCL = predict(normalized.licor.e, upper.CL, inverse = TRUE))
#cld 
cld.licor_e<-cld(emm_licor_e, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#merge the tables
emm_licor_e_back_transformed$treatment <- paste(emm_licor_e_back_transformed$mulch, emm_licor_e_back_transformed$glyphosate)
cld.licor_e$treatment <- paste(cld.licor_e$mulch, cld.licor_e$glyphosate)
#use this for visualization: plot response column as y variable, mulch.x for mulch
cld_licor_e_for_plot <- left_join(emm_licor_e_back_transformed, data.frame(cld.licor_e), by = "treatment")
#Visualize the results
ggplot(data.frame(cld_licor_e_for_plot), aes(x=mulch.x, y=response, fill=glyphosate.x, color=glyphosate.x))+
  facet_wrap(.~glyphosate.x)+
  geom_pointrange(aes(ymin=bt_lowerCL, ymax=bt_upperCL), size=0.9, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_fill_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_x_discrete(labels=c("mulch.x"="Mulch", "control"="Control")) + 
  geom_text(aes(label=.group, y = bt_upperCL, group=glyphosate.x),color="black",size=3, vjust = -0.5,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Transpiration rate") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.licor.e.pdf", device="pdf", units="mm", dpi=300, width=100, height=100, path="/Users/rileyjones/Library/Mobile Documents/com~apple~CloudDocs/91 C data")

#########################
####### Licor A #########
#########################
# Step 1: Run bestNormalize on your response
normalized.licor.a <- bestNormalize(licor_cleaned$licor_a)
# Step 2: Add the transformed variable to your dataset
licor_cleaned$licor_a_trans <- predict(normalized.licor.a)
# Step 3: Fit a linear mixed model
#mulch*humic_acid*glyphosate:
#model_trans_licor_a<- lmer(licor_a_trans ~ mulch*humic_acid*glyphosate + (1 | year), data = licor_cleaned)
#remove humic b/c not significant:
model_trans_licor_a<- lmer(licor_a_trans ~ mulch*glyphosate + (1 | year), data = licor_cleaned)
# Step 4: Plot residuals vs fitted values
plot(fitted(model_trans_licor_a), resid(model_trans_licor_a),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")

anova(model_trans_licor_a)
ranova(model_trans_licor_a) #year is significant
#VISUALISATION OF LICOR A:
emm_licor_a<-emmeans(model_trans_licor_a, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#because I used best normalize to create the model, have to do the back transformation by hand using this code:
# Back-transform the EMMs and their CIs
emm_licor_a_back_transformed <- as.data.frame(emm_licor_a) %>%
  mutate(response = predict(normalized.licor.a, emmean, inverse = TRUE),
         bt_lowerCL= predict(normalized.licor.a, lower.CL, inverse = TRUE),
         bt_upperCL = predict(normalized.licor.a, upper.CL, inverse = TRUE))
#cld 
cld.licor_a<-cld(emm_licor_a, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#merge the tables
emm_licor_a_back_transformed$treatment <- paste(emm_licor_a_back_transformed$mulch, emm_licor_a_back_transformed$glyphosate)
cld.licor_a$treatment <- paste(cld.licor_a$mulch, cld.licor_a$glyphosate)
#use this for visualization: plot response column as y variable, mulch.x for mulch
cld_licor_a_for_plot <- left_join(emm_licor_a_back_transformed, data.frame(cld.licor_a), by = "treatment")
#Visualize the results
ggplot(data.frame(cld_licor_a_for_plot), aes(x=mulch.x, y=response, fill=glyphosate.x, color=glyphosate.x))+
  facet_wrap(.~glyphosate.x)+
  geom_pointrange(aes(ymin=bt_lowerCL, ymax=bt_upperCL), size=0.9, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_fill_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_x_discrete(labels=c("mulch.x"="Mulch", "control"="Control")) + 
  geom_text(aes(label=.group, y = bt_upperCL, group=glyphosate.x),color="black",size=3, vjust = -0.5,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Carbon assimilation rate") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.licor.a.pdf", device="pdf", units="mm", dpi=300, width=100, height=100, path="/Users/rileyjones/Library/Mobile Documents/com~apple~CloudDocs/91 C data")

######################
#### Licor gsw #######
######################
# Step 1: Run bestNormalize on your response
normalized.licor.gsw <- bestNormalize(licor_cleaned$licor_gsw)
# Step 2: Add the transformed variable to your dataset
licor_cleaned$licor_gsw_trans <- predict(normalized.licor.gsw)
# Step 3: Fit a linear mixed model
#model_trans_licor_gsw<- lmer(licor_gsw_trans ~ mulch*humic_acid*glyphosate + (1 | year), data = licor_cleaned)
#remove humic_acid b/c not significant:
model_trans_licor_gsw<- lmer(licor_gsw_trans ~ mulch*glyphosate + (1 | year), data = licor_cleaned)
# Step 4: Plot residuals vs fitted values
plot(fitted(model_trans_licor_gsw), resid(model_trans_licor_gsw),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")
anova(model_trans_licor_gsw)
ranova(model_trans_licor_gsw) #year is significant
#VISUALISATION OF LICOR GSW:
emm_licor_gsw<-emmeans(model_trans_licor_gsw, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#because I used best normalize to create the model, have to do the back transformation by hand using this code:
# Back-transform the EMMs and their CIs
emm_licor_gsw_back_transformed <- as.data.frame(emm_licor_gsw) %>%
  mutate(response = predict(normalized.licor.gsw, emmean, inverse = TRUE),
         bt_lowerCL= predict(normalized.licor.gsw, lower.CL, inverse = TRUE),
         bt_upperCL = predict(normalized.licor.gsw, upper.CL, inverse = TRUE))
#cld 
cld.licor_gsw<-cld(emm_licor_gsw, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#merge the tables
emm_licor_gsw_back_transformed$treatment <- paste(emm_licor_gsw_back_transformed$mulch, emm_licor_gsw_back_transformed$glyphosate)
cld.licor_gsw$treatment <- paste(cld.licor_gsw$mulch, cld.licor_gsw$glyphosate)
#use this for visualization: plot response column as y variable, mulch.x for mulch
cld_licor_gsw_for_plot <- left_join(emm_licor_gsw_back_transformed, data.frame(cld.licor_gsw), by = "treatment")
#Visualize the results
ggplot(data.frame(cld_licor_gsw_for_plot), aes(x=mulch.x, y=response, fill=glyphosate.x, color=glyphosate.x))+
  facet_wrap(.~glyphosate.x)+
  geom_pointrange(aes(ymin=bt_lowerCL, ymax=bt_upperCL), size=0.9, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_fill_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_x_discrete(labels=c("mulch.x"="Mulch", "control"="Control")) + 
  geom_text(aes(label=.group, y = bt_upperCL, group=glyphosate.x),color="black",size=3, vjust = -0.5,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Stomatal Conductance") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.licor.gsw.pdf", device="pdf", units="mm", dpi=300, width=100, height=100, path="/Users/rileyjones/Library/Mobile Documents/com~apple~CloudDocs/91 C data")

######################################################
######################################################
#### Weed Biomass Analysis ###########################
######################################################
######################################################
weed_biomass # data frame
# clean-up spaces between column names
weed_biomass_cleaned <- weed_biomass %>% clean_names()
# Step 1: Run bestNormalize on your response
normalized.weed_biomass_dry_weight <- bestNormalize(weed_biomass_cleaned$weed_biomass_dry_weight)
# Step 2: Add the transformed variable to your dataset
weed_biomass_cleaned$weed_biomass_dry_weight_trans <- predict(normalized.weed_biomass_dry_weight)
# Step 3: Fit a linear mixed model
#model_trans_weed_biomass_dry_weight<- lm(weed_biomass_dry_weight_trans ~ mulch*humic_acid*glyphosate, data = weed_biomass_cleaned)
#remove humic_acid bc not sig
model_trans_weed_biomass_dry_weight<- lm(weed_biomass_dry_weight_trans ~ mulch*glyphosate, data = weed_biomass_cleaned)
# Step 4: Plot residuals vs fitted values
plot(fitted(model_trans_weed_biomass_dry_weight), resid(model_trans_weed_biomass_dry_weight),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")

anova(model_trans_weed_biomass_dry_weight)
ranova(model_trans_weed_biomass_dry_weight) #year is not significant

#Visualization of weed biomass dry weight
emm_weed_biomass_dry_weight<-emmeans(model_trans_weed_biomass_dry_weight, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#because I used best normalize to create the model, have to do the back transformation by hand using this code:
# Back-transform the EMMs and their CIs
emm_weed_biomass_dry_weight_back_transformed <- as.data.frame(emm_weed_biomass_dry_weight) %>%
  mutate(response = predict(normalized.weed_biomass_dry_weight, emmean, inverse = TRUE),
         bt_lowerCL= predict(normalized.weed_biomass_dry_weight, lower.CL, inverse = TRUE),
         bt_upperCL = predict(normalized.weed_biomass_dry_weight, upper.CL, inverse = TRUE))
#cld 
cld.weed_biomass_dry_weight<-cld(emm_weed_biomass_dry_weight, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#merge the tables
emm_weed_biomass_dry_weight_back_transformed$treatment <- paste(emm_weed_biomass_dry_weight_back_transformed$mulch, emm_weed_biomass_dry_weight_back_transformed$glyphosate)
cld.weed_biomass_dry_weight$treatment <- paste(cld.weed_biomass_dry_weight$mulch, cld.weed_biomass_dry_weight$glyphosate)
#use this for visualization: plot response column as y variable, mulch.x for mulch
cld_dry_weight_for_plot <- left_join(emm_weed_biomass_dry_weight_back_transformed, data.frame(cld.weed_biomass_dry_weight), by = "treatment")
#Visualize the results
ggplot(data.frame(cld_dry_weight_for_plot), aes(x=mulch.x, y=response, fill=glyphosate.x, color=glyphosate.x))+
  facet_wrap(.~glyphosate.x)+
  geom_pointrange(aes(ymin=bt_lowerCL, ymax=bt_upperCL), size=0.9, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_fill_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_x_discrete(labels=c("mulch"="Mulch", "control"="Control")) + 
  geom_text(aes(label=.group, y = bt_upperCL, group=glyphosate.x),color="black",size=3, vjust = -0.5,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Weed biomass dry weight") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.weed.biomass.dry.weight.pdf", device="pdf", units="mm", dpi=300, width=100, height=100, path="/Users/rileyjones/Library/Mobile Documents/com~apple~CloudDocs/91 C data")

######################################################
######################################################
#### Pack-line fruit Analysis ########################
######################################################
######################################################
final_packline # data frame
colnames(final_packline)[colnames(final_packline) == "Fruit_weight"] <- "mean_weight"

# Model mean fruit weight
model_mean_weight <- lmer(mean_weight ~ mulch*humic*glyphosate + (1 | year), data = final_packline)
#simplified because only mulch and glyphosate had an effect:
model_mean_weight<- lmer(mean_weight ~ mulch*glyphosate + (1 | year), data = final_packline)
# Plot residuals vs fitted values
plot(fitted(model_mean_weight), resid(model_mean_weight),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")

anova(model_mean_weight)
ranova(model_mean_weight)

#emm
emm_mean_weight<-emmeans(model_mean_weight, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means
#cld 
cld.mean_weight<-cld(emm_mean_weight, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display
#Visualize the results
ggplot(data.frame(cld.mean_weight), aes(x=mulch, y=emmean, fill=glyphosate, color=glyphosate))+
  facet_wrap(.~glyphosate)+
  geom_pointrange(aes(ymin=emmean - SE, ymax=emmean + SE), size=0.9, position=position_dodge2(1.0)) +
  scale_color_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_fill_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("glyphosate","control"), labels = c("Glyphosate","Control")) +
  scale_x_discrete(labels=c("mulch"="Mulch", "control"="Control")) + 
  geom_text(aes(label=.group, y = emmean + SE + 5, group=glyphosate),color="black",size=3, vjust = -0.5,hjust=0.5, position=position_dodge2(1.0))+
  ylab("Fruit Weight") +
  xlab("Treatment") +
  theme_bw()
# save plot
ggsave(filename= "emm.final.fruit.weight.pdf", device="pdf", units="mm", dpi=300, width=100, height=100, path="/Users/rileyjones/Library/Mobile Documents/com~apple~CloudDocs/91 C data")

