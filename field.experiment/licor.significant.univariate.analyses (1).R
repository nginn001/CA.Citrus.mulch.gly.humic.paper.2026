library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)

#Licor E:

#####


library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)

set.seed(4444)
#use bestNormalize:
normalized.licor.E <- bestNormalize(licor_cleaned$licor_e)

#plot normalized residuals:
library(bestNormalize)
library(lme4)

# Step 1: Run bestNormalize on your response
normalized.licor.e <- bestNormalize(licor_cleaned$licor_e)

# Step 2: Add the transformed variable to your dataset
licor_cleaned$e_trans <- predict(normalized.licor.e)

# Step 3: Fit a linear mixed model
#mulch*humic_acid*glyphosate
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

#VISUALISATION OF LICOR E: getting a weird output
#####

#Load libraries 
library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)

emm_licor_e<-emmeans(model_trans_licor_e, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means

#because i used best normalize to create the model, have to do the back transformation by hand using this code:

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



#####

#Licor A:

#####

# de-bugged code:

library(multcomp)
library(ggplot2)
library(emmeans)
library(bestNormalize)
library(lme4)

set.seed(4444)

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

#Significant
#####

#VISUALISATION OF LICOR A:

#####

#Load libraries 
library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)

emm_licor_a<-emmeans(model_trans_licor_a, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means

#because i used best normalize to create the model, have to do the back transformation by hand using this code:

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




#####

#Licor gsw:

#####
#de-bugged code:

library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)

set.seed(4444)
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
#significant
#####

#VISUALISATION OF LICOR GSW:
#####

#Load libraries 
library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)

emm_licor_gsw<-emmeans(model_trans_licor_gsw, specs=~mulch:glyphosate, type="response") # calculate estimated marginal means

#because i used best normalize to create the model, have to do the back transformation by hand using this code:

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


#####

#Licor Ca: Ignore

#####
library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)

set.seed(4444)
# Step 1: Run bestNormalize on your response
normalized.licor.ca <- bestNormalize(licor_cleaned$licor_ca)

# Step 2: Add the transformed variable to your dataset
licor_cleaned$licor_ca_trans <- predict(normalized.licor.ca)

# Step 3: Fit a linear mixed model
model_trans_licor_ca<- lmer(licor_ca_trans ~ mulch*humic_acid*glyphosate + (1 | year), data = licor_cleaned)
model_trans_licor_ca<- lmer(licor_ca_trans ~ mulch*glyphosate + (1 | year), data = licor_cleaned)

# Step 4: Plot residuals vs fitted values: bimodal and very difficult to fix - leave out of final analysis?
plot(fitted(model_trans_licor_ca), resid(model_trans_licor_ca),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")

anova(model_trans_licor_ca)
ranova(model_trans_licor_ca) #year is significant

#####

#Licor Ci: NOT significant

#####

library(multcomp)
library(ggplot2)
library(emmeans)
library(bestNormalize)
library(lme4)

set.seed(4444)
# Step 1: Run bestNormalize on your response
normalized.licor.ci <- bestNormalize(licor_cleaned$licor_ci)

# Step 2: Add the transformed variable to your dataset
licor_cleaned$licor_ci_trans <- predict(normalized.licor.ci)

# Step 3: Fit a linear mixed model
model_trans_licor_ci<- lmer(licor_ci_trans ~ treatment + (1 | year), data = licor_cleaned)
#mulch*humic_acid*glyphosate:
model_trans_licor_ci<- lmer(licor_ci_trans ~ mulch*humic_acid*glyphosate + (1 | year), data = licor_cleaned)
#remove humic b/c not significant:
model_trans_licor_ci<- lmer(licor_ci_trans ~ mulch*glyphosate + (1 | year), data = licor_cleaned)


# Step 4: Plot residuals vs fitted values
plot(fitted(model_trans_licor_ci), resid(model_trans_licor_ci),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")

anova(model_trans_licor_ci)
ranova(model_trans_licor_ci) #year is significant

#####

#Licor iWUE: NOT significant

#####
library(multcomp)
library(ggplot2)
library(emmeans)
library(bestNormalize)
library(lme4)

set.seed(4444)

# Step 1: Run bestNormalize on your response
normalized.licor.iWUE <- bestNormalize(licor_cleaned$licor_i_wue)

# Step 2: Add the transformed variable to your dataset
licor_cleaned$licor_i_wue_trans <- predict(normalized.licor.iWUE)

# Step 3: Fit a linear mixed model
model_trans_licor_iWUE<- lmer(licor_i_wue_trans ~ mulch*humic_acid*glyphosate + (1 | year), data = licor_cleaned)


# Step 4: Plot residuals vs fitted values
plot(fitted(model_trans_licor_iWUE), resid(model_trans_licor_iWUE),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")

anova(model_trans_licor_iWUE)
ranova(model_trans_licor_iWUE)
####

#visualization of: E, A, gsw, Ca
#trying to make a loop: 

##### 

licor.list <- c(model_trans_licor_ca, model_trans_licor_gsw, model_trans_licor_a, model_trans_licor_e) 
for(i in licor.list
){emm_licor<-emmeans(i, specs=~mulch:glyphosate) # calculate estimated marginal means

cld.licor<-cld(emm_licor, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # Tukey posthoc and compact letter display

#Visualize the results
ggplot(data.frame(cld.licor), aes(x=mulch, y=emmean, fill=glyphosate, color=glyphosate))+
  #facet_wrap(.~glyphosate)+
  geom_pointrange(aes(ymin=emmean-SE, ymax=emmean+SE), size=0.9, position=position_dodge2(1.0)) +
  #scale_color_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("humic_acid","control"), labels = c("Humic Acid","Control")) +
 # scale_fill_manual(values=c("#8d1c5c", "#f7cf49"),limits=c("humic_acid","control"), labels = c("Humic Acid","Control")) +
  scale_x_discrete(labels=c("mulch"="Mulch", "control"="Control")) + 
  geom_text(aes(label=.group, y = emmean+SE, group=glyphosate),color="black",size=3, vjust = -0.5,hjust=0.5, position=position_dodge2(1.0))+
  ylab(paste0(i)) +
  xlab("treatment") +
  theme_bw()
# save plot
ggsave(filename= paste0("emm.",i,".pdf"), device="pdf", units="mm", dpi=300, width=30, height=100, path="/Users/rileyjones/Library/Mobile Documents/com~apple~CloudDocs/91 C data")
print(i)
}

 