#!/usr/local/bin/Rscript

#Install packages
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  tidyverse,
  data.table)

#Association analysis of Protsig with DM,BP markers
trait <- fread("${lipid_trait}")
prs_pcs <- fread("${ug_prs_pcs}")
traits <- trait %>% 
  dplyr::filter(!(is.na(sangerid)|sangerid == "")) %>% 
  rename(FID = sangerid) %>% 
  setnames(., names(.), c(gsub("-","_",colnames(.))))
traits_cov <- merge(prs_pcs, traits, by="FID")

#Transform to same scale
traits_cov_z <- traits_cov %>%
  mutate(
    HbA1c = ifelse(HbA1c <= 0, NA, HbA1c),
    TC = ifelse(TC <= 0, NA, TC),
    LDL_C = ifelse(LDL_C <= 0, NA, LDL_C),
    HDL_C = ifelse(HDL_C <= 0, NA, HDL_C),
    systolic = ifelse(systolic <= 0, NA,systolic),
    TG = ifelse(TG <= 0, NA,TG),
  ) %>%
  mutate(
    HbA1c_z = scale(log(HbA1c)),
    TC_z = scale(log(TC)),
    LDL_C_z = scale(log(LDL_C)),
    HDL_C_z = scale(log(HDL_C)),
    systolic_z = scale(log( systolic)),
    TG_z = scale(log(TG)),
  )

#Prediction - Lipid traits
predictor <- "s1"
outcome_variables <- c("TC_z", "HbA1c_z", "LDL_C_z" , "HDL_C_z", "systolic_z", "TG_z")
models <- list()
for (outcome_variable in outcome_variables) {
  formula <- as.formula(paste0(outcome_variable, "~", predictor, "+ sex + age + PC1 + PC2 + PC3"))
  model <- lm(formula, data=traits_cov_z)
  models[[outcome_variable]] <- model
  cat("Linear model for", outcome_variable, ":\\n")
  print(summary(model))
  conf_intervals <- confint(model)
  cat("Confidence Intervals:\\n")
  print(conf_intervals)
  cat("\\n")
}

# Generate a forest plot
results <- data.frame(
  Outcome = character(),
  Estimate = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each model to extract coefficients for 's1' and store in the dataframe
for (outcome_variable in outcome_variables) {
  model <- models[[outcome_variable]]
  
  # Check if the 's1' term exists in the model summary
  term <- "s1"  # This is the main effect of interest
  if (term %in% rownames(summary(model)\$coefficients)) {
    est <- coef(summary(model))[term, "Estimate"]
    ci <- confint(model, term, level = 0.95)
    
    # Append results to dataframe
    results <- rbind(results, data.frame(
      Outcome = outcome_variable,
      Estimate = est,
      LowerCI = ci[1],
      UpperCI = ci[2]
    ))
  } else {
    # Handle the case where the main effect is not present
    results <- rbind(results, data.frame(
      Outcome = outcome_variable,
      Estimate = NA,
      LowerCI = NA,
      UpperCI = NA
    ))
  }
}
write.table(results, file ="${blood_trait}_results.csv", sep = "\\t", col.names =F, quote=F, row.names=F)

# Adjust the plot with modified 'fatten' parameter
vitc_forest_plot <- ggplot(results, aes(x=reorder(Outcome,Estimate), 
                                        y=Estimate, ymin=LowerCI, ymax=UpperCI)) + 
  geom_linerange(size=5,position=position_dodge(width = 0.5),color="blue4") + 
  geom_hline(yintercept=0, lty=2) + 
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,
             position=position_dodge(width = 0.5)) +
  coord_flip() + 
  theme_classic() +
  theme(axis.title.y = element_blank()) + 
  theme( panel.background = element_blank(), panel.grid = element_blank()) +
  theme(text=element_text(size = 14)) +ylab("Beta(95%CI)")


pdf("${blood_trait}_forest_plot.pdf", width = 6, height = 4)
print(vitc_forest_plot)
dev.off()

