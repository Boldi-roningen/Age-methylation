#Set up libraries
library(readr)
library(tidyverse)
library(caret)
library(tidymodels)
library(mice)
#Import Data
data <- read_rds("Results/Cov_ZF_SampleInfo_Feb2024_Subset.rds")
#Creating the data frame for modelling, one column with Age at Collection, the remaining Columns with cpg sites with their methylation percentage.
model_data <- data_frame(data[c("cpg_site", "MethylationPercentage", "AgeCollectionDays")])%>%
  pivot_wider(names_from = cpg_site, values_from = MethylationPercentage)
print("Data pivoted")
#Imputing missing values using the package mice
#init = mice(model_data, maxit=0)
#meth = init$method
#predM = init$predictorMatrix
#colnames(model_data)
#predM[, c("AgeCollectionDays")]=0
#meth[c("AgeCollectionDays")]=""
#set.seed(123)
#imputed = mice(model_data, method = meth, predictorMatrix=predM, m=5)
#print("missing valeus imputed")

#Data transformation: removing cpgs with little variation in TRAINING data
nearzerovar_cpg_list <- nearZeroVar(model_data, freqCut = 85/15, uniqueCut = 50, saveMetrics = FALSE) 
if (length(nearzerovar_cpg_list) == 0) {
  nearzerovar_data <- model_data
} else {
  nearzerovar_data <- model_data[,-nearzerovar_cpg_list]
}
print("Cpg sites with little to no variation removed")

#Data transformation: identifying cpg sites that are highly correlated
highlyCorDescr <- findCorrelation(nearzerovar_data, cutoff = 0.8)
nearzerovar_nocor_data <- nearzerovar_data[,-highlyCorDescr]
print("highly correlated cpg sites removed")

#Data transformation: standardizing and scaling data
preProcValues <- preProcess(nearzerovar_nocor_data, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, nearzerovar_nocor_data)
print("Data centered on mean 0, std 1")

#Splitting data
set.seed(123)
splits <- createDataPartition(trainTransformed$AgeCollectionDays, p = 0.8, list = FALSE, times = 1)
training_data <- trainTransformed[splits,]
testing_data <- trainTransformed[-splits,]
print("Data split")
#Model tuning

fitControl <- trainControl(method = 'repeatedcv' ,
                           number=10, repeats=10)

lasso_model <- train(AgeCollectionDays ~ ., 
                     data = training_data, 
                     method = "glmnet", 
                     trControl = fitControl, 
                     tuneGrid = data.frame(
                       alpha = 1,
                       lambda = 10^seq(-3, 3, length = 100)), 
                     tuneLength = 10)
print("lasso model created")

elastic_model <- train(AgeCollectionDays ~ ., 
                       data = training_data, 
                       method ="glmnet", 
                       trControl = fitControl, 
                       tuneLength = 10)
print("elastic model created")

elastic_model_05 <- train(AgeCollectionDays ~ ., 
                          data = training_data,
                          method = "glmnet", 
                          trControl = fitControl, 
                          tuneGrid = data.frame(
                            alpha = 0.5, 
                            lambda = 10^seq(-3, 3, length = 100)), 
                          tuneLength = 10)
print("elasti model with alpha 0.5 created")

#Comparing Models
set.seed(123)
models_compare <- resamples(list(LM=lasso_model, EM=elastic_model, EM05=elastic_model_05))
sink(file = "model_comparison.txt")
summary(models_compare)
sink(file = NULL)
print("model comparison created")

#Model evaluations
set.seed(123)
sink(file = "lasso_model.txt")
print("Training data")
predicted.age <- predict.train(lasso_model)
postResample(pred = predicted.age, training_data$AgeCollectionDays)
cor.test(predicted.age, training_data$AgeCollectionDays)
print("Testing data")
predict.lasso.test <- predict(lasso_model, testTransformed)
postResample(pred = predict.lasso.test, testing_data$AgeCollectionDays)
cor.test(predict.lasso.test, testing_data$AgeCollectionDays)
sink(file = NULL)
print("lasso model_Evaluated")

set.seed(123)
sink(file = "elastic_model.txt")
print("Training data")
predicted.age <- predict.train(elastic_model)
postResample(pred = predicted.age, training_data$AgeCollectionDays)
cor.test(predicted.age, training_data$AgeCollectionDays)
print("Testing data")
predict.elastic.test <- predict(elastic_model, testTransformed)
postResample(pred = predict.elastic.test, testing_data$AgeCollectionDays)
cor.test(predict.elastic.test, testing_data$AgeCollectionDays)
sink(file = NULL)
print("elastic model evaluated")

set.seed(123)
sink(file = "elastic_model_05.txt")
print("Training data")
predicted.age <- predict.train(elastic_model_05)
postResample(pred = predicted.age, training_data$AgeCollectionDays)
cor.test(predicted.age, training_data$AgeCollectionDays)
print("Testing data")
predict.elastic.test <- predict(elastic_model_05, testTransformed)
postResample(pred = predict.elastic.test, testing_data$AgeCollectionDays)
cor.test(predict.elastic.test, testing_data$AgeCollectionDays)
sink(file = NULL)
print("elastic model with alpha 0.5 evaluated")

#Testing models on TESTING data
set.seed(123)
sink(file = "lasso_testing.txt")
predict_ridge_test <- predict(lasso_model, testTransformed)
postResample(pred = predict.ridge.test, testing_data$AgeCollectionDays)
cor.test(predict.ridge.test, testing_data$AgeCollectionDays)
sink(file = NULL)

set.seed(123)
sink(file = "elastic_testing.txt")
predict_ridge_test <- predict(elastic_model, testTransformed)
postResample(pred = predict.ridge.test, testing_data$AgeCollectionDays)
cor.test(predict.ridge.test, testing_data$AgeCollectionDays)
sink(file = NULL)

set.seed(123)
sink(file = "elastic05_testing.txt")
predict_ridge_test <- predict(elastic_model_05, testTransformed)
postResample(pred = predict.ridge.test, testing_data$AgeCollectionDays)
cor.test(predict.ridge.test, testing_data$AgeCollectionDays)
sink(file = NULL)
print("models tested")
