#!/usr/local/bin/Rscript

#Install packages
#if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  readr,
  dplyr,
  data.table,
  glmnet,
  purrr,
  tidyverse,
  caret,
  ggplot2)

#Step 0 Read in data
#Data frame with FID and PRS scores
prs = read.table("${score_files}", header = TRUE, sep = "\\t")
#Dataframe with FID and pheno
phenotype = read.table("${pheno}", header = TRUE, sep = "\\t")
pheno = data.frame(FID = phenotype\$FID, ${blood_trait} = phenotype\$${blood_trait})

#Step 1: Build Elastic Net Model
#1.1 Transform prs data
prstransform <- prs %>% mutate(across(-1, ~ if(is.numeric(.)) { scale(log1p(.)) } else { . }))

#1.2 merge transformed data to pheno data
pheno_prs_transformed <- merge(pheno, prstransform,  by ="FID")
#remove FID and IID columns
pheno_prs_transformednoFID <- pheno_prs_transformed[,-1]

#1.3 Create training and test data
set.seed(123)
# Check for missing values or NaNs
any_missing <- anyNA(pheno_prs_transformednoFID\$${blood_trait})
if (any_missing) {
  pheno_prs_transformednoFID <- pheno_prs_transformednoFID[complete.cases(pheno_prs_transformednoFID\$${blood_trait}), ]
}
#trainingx <- createDataPartition(pheno_prs_transformednoFID\$mcv, p = 0.8, list = FALSE)
trainingx <- pheno_prs_transformednoFID\$${blood_trait} %>% createDataPartition(p=0.8, list=FALSE)
traindta  <- pheno_prs_transformednoFID[trainingx, ]
testdta  <- pheno_prs_transformednoFID[-trainingx, ]

#1.4 Run Elasticnet
set.seed(123)
modelenet <- train(${blood_trait} ~., data = traindta, method ="glmnet", trControl = trainControl("cv", number =10), tuneLenght=10)

x <- coef(modelenet\$finalModel, modelenet\$bestTune\$lambda)
dense_matrix <- as.matrix(x)
predicted <- as.data.frame(dense_matrix)
write.table(predicted, file ="${blood_trait}_predicted.txt", sep = "\\t", col.names =T, quote=F, row.names=T)

#1.5 make predictions
predictions <- modelenet %>% predict(testdta) 

#1.6 Model prediction performance
data.frame(RMSE = RMSE(predictions,
 testdta\$${blood_trait}), Rsquare = R2(predictions, testdta\$${blood_trait}))


#Step 2: Generate intergrated scores

# 2.1 Make sure the dimensions of the matrix and the coefficients vector match
coefficients <- as.matrix(coef(modelenet\$finalModel, modelenet\$bestTune\$lambda)[-1, 1, drop = FALSE])
predictor_variables <- as.matrix(prstransform[, -1])

#2.2 Check if dimensions match before matrix multiplication
if (ncol(predictor_variables) == length(coefficients)) {
 integrated_scores <- predictor_variables %*% coefficients
  #Combine the scores with the participant IDs
  final_scores <- data.frame(FID = prstransform[, 1], IntegratedScore = integrated_scores)
} else {
  stop("The number of predictor variables does not match the number of coefficients.")
}

#2.3 Write out the intergrated scores which is the prot signature
write.table(final_scores, file ="${blood_trait}_protsignature.txt", col.names =T, quote=F, row.names=F)

#Step 3
# 3.1 Correlation of intergrated PRS and phenotype
protsigdta <-  merge(final_scores, pheno, by = "FID")
protsigdtax <- na.omit(protsigdta)
cor.test.result <- suppressWarnings(cor.test(protsigdtax\$s1, protsigdtax\$${blood_trait}, method = "pearson"))
pearson_cor <- cor.test.result\$estimate
p_value <- cor.test.result\$p.value# Create the base plot
p <- ggplot(protsigdtax, aes(x = s1, y = ${blood_trait})) +
  geom_point() +  # This adds the scatter plot points
  geom_smooth(method = "lm", col = "red") +  # This adds the regression line
  labs(x = "intergrated_score", y = "${blood_trait}") +
  annotate("text", x = Inf, y = Inf, label = sprintf("Pearson: %.2f\np-value: %.3f", pearson_cor, p_value),
             hjust = 1.1, vjust = 1.1, size = 4, color = "blue")

pdf("${blood_trait}_correlation.pdf", width = 6, height = 4)
print(p)
dev.off()

#Step 4 
as_input = read.table("${cov}", header = TRUE)

ug_prs_pcs <- merge(protsigdta, as_input, by = "FID")
 
alm <- lm(${blood_trait} ~ s1+age+sex+PC1+PC2+PC3, data = ug_prs_pcs)
summary_alm <- summary(alm)
almx <- lm(${blood_trait}~ age+sex+PC1+PC2+PC3, data = ug_prs_pcs)
summary_almx = summary(almx)

pred_score <- summary_alm\$adj.r.squared-summary_almx\$adj.r.squared
write.table(pred_score, file ="${blood_trait}_pred_score.txt", col.names =F, quote=F, row.names=F)
write.table(ug_prs_pcs, file ="${blood_trait}_ug_prs_pcs.txt", col.names =T, quote=F, row.names=F)
