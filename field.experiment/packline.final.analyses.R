

final_packline <- read_excel("final_packline_data_for_analysis.xlsx")

model_mean_weight <- lmer(mean_weight ~ mulch*humic*glyphosate + (1 | year), data = final_packline)
#simplified because only mulch and glyphosate had an effect:
model_mean_weight<- lmer(mean_weight ~ mulch*glyphosate + (1 | year), data = final_packline)

# Step 4: Plot residuals vs fitted values
plot(fitted(model_mean_weight), resid(model_mean_weight),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")

anova(model_mean_weight)
ranova(model_mean_weight)


#Visualization of fruit number:

#####
#Load libraries 
library(multcomp)
library(ggplot2)
library(emmeans)
library(lme4)

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
