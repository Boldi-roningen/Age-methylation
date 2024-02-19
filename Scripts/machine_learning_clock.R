#Set up libraries
library(tidymodels)
library(readr)
library(caret)
library(tidyverse)

#Import Data
data <- read_rds("")

#Create a data frame for modelling
head(merge_Covn)
model_data <- data[c("Chromosome", "MethylationPercentage", "BirdID", "AgeCollectionDays")]
#random seed for model reproductability
set.seed(2000)

#Split the data into training (80%) and testing (20%)
splits <- initial_split(model_data, strata = FinalAgeDays)
training_data <- training(splits)
testing_data <- testing(splits)

#Data transformation: removing cpgs with little variation in TRAINING data

nearzerovar_cpg_list <- nearZeroVar(training_data, freqCut = 85/15, uniqueCut = 50) 
nearzerovar_training_data <- training_data[, -nearzerovar_cpg_list]

#Data transformation: identifying cpg sites that are highly correlated

highlyCorDescr <- findCorrelation(nearzerovar_training_data, cutoff = 0.8)
training_data_model <- nearzerovar_training_data[,-highlyCorDescr]

#Data transformation: standardizing and scaling data
preProcValues <- preProcess(training_data_model, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, training_data_model)

#Data transformation: removing cpgs with little variation in TESTING data
nearzerovar_cpg_list2 <- nearZeroVar(testing_data, freqCut = 85/15, uniqueCut = 50) 
nearzerovar_testing_data <- testing_data[, -nearzerovar_cpg_list2]

#Data transformation: identifying cpg sites that are highly correlated

highlyCorDescr2 <- findCorrelation(nearzerovar_testing_data, cutoff = 0.8)
testing_data_model <- nearzerovar_testing_data[,-highlyCorDescr2]

#Data transformation: standardizing and scaling data
preProcValues2 <- preProcess(testing_data_model, method = c("center", "scale"))
testTransformed <- predict(preProcValues2, testing_data_model)

#Model tuning
fitControl <- trainControl(method = 'repeatedcv' ,
                           number=10, repeats=10)

ridge_model <- train(AgeCollectionDays ~ ., 
                     data = trainTransformed, 
                     method = "glmnet", 
                     trControl = fitControl, 
                     tuneGrid = data.frame(
                     alpha = 0,
                     lambda = 10^seq(-3, 3, length = 100)), 
                     tuneLength = 10)

lasso_model <- train(AgeCollectionDays ~ ., 
                     data = trainTransformed, 
                     method = "glmnet", 
                     trControl = fitControl, 
                     tuneGrid = data.frame(
                       alpha = 1,
                       lambda = 10^seq(-3, 3, length = 100)), 
                     tuneLength = 10)
elastic_model <- train(AgeCollectionDays ~ ., 
                       data = trainTransformed, 
                       method ="glmnet", 
                       trControl = fitControl, 
                       tuneLength = 10)

elastic_model_05 <- train(AgeCollectionDays ~ ., 
                          data = trainTransformed,
                          method = "glmnet", 
                          trControl = fitControl, 
                          tuneGrid = data.frame(
                            alpha = 0.5, 
                            lambda = 10^seq(-3, 3, length = 100)), 
                          tuneLength = 10)
#Comparing Models
models_compare <- resamples(list(R=ridge_model,
                                 LM=lasso_model, EM=elastic_model, EM05=elastic_model.05))
sink(file = "model_comparison.txt")
summary(models_compare)
sink(file = NULL)
#Model evaluations
sink(file = "ridge_model.txt")
predicted_age <- predict.train(ridge_model)
postResample(pred = predicted_age, trainTransformed$AgeCollectionDays)
cor.test(predicted_age, trainTransformed$AgeCollectionDays)
sink(file = NULL)

sink(file = "lasso_model.txt")
predicted.age <- predict.train(lasso_model)
postResample(pred = predicted.age, trainTransformed$AgeCollectionDays)
cor.test(predicted.age, trainTransformed$AgeCollectionDays)
sink(file = NULL)

sink(file = "elastic_model.txt")
predicted.age <- predict.train(elastic_model)
postResample(pred = predicted.age, trainTransformed$AgeCollectionDays)
cor.test(predicted.age, trainTransformed$AgeCollectionDays)
sink(file = NULL)

sink(file = "elastic_model_05.txt")
predicted.age <- predict.train(elastic_model_05)
postResample(pred = predicted.age, trainTransformed$AgeCollectionDays)
cor.test(predicted.age, trainTransformed$AgeCollectionDays)
sink(file = NULL)

#Testing models on transformed TESTING data
sink(file = "ridge_testing.txt")
predict_ridge_test <- predict(ridge_model, testTransformed)
postResample(pred = predict.ridge.test, testTransformed$AgeCollectionDays)
cor.test(predict.ridge.test, testTransformed$AgeCollectionDays)
sink(file = NULL)

sink(file = "lasso_testing.txt")
predict_ridge_test <- predict(lasso_model, testTransformed)
postResample(pred = predict.ridge.test, testTransformed$AgeCollectionDays)
cor.test(predict.ridge.test, testTransformed$AgeCollectionDays)
sink(file = NULL)


sink(file = "elastic_testing.txt")
predict_ridge_test <- predict(elastic_model, testTransformed)
postResample(pred = predict.ridge.test, testTransformed$AgeCollectionDays)
cor.test(predict.ridge.test, testTransformed$AgeCollectionDays)
sink(file = NULL)

sink(file = "elastic05_testing.txt")
predict_ridge_test <- predict(elastic_model_05, testTransformed)
postResample(pred = predict.ridge.test, testTransformed$AgeCollectionDays)
cor.test(predict.ridge.test, testTransformed$AgeCollectionDays)
sink(file = NULL)



#Looking at different predictors