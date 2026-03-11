##**************************
## BINF*6970 - Assignment 2: Problem 2
##
## Student Names: Cyrus Akbarally & William Feinman 
##
## Due Date: 2026-03-13
##
##**************************

## _Packages used----
library(glmnet)
library(ggplot2)
library(caret)
library(pROC)
library(readxl)
library(dplyr)
## setwd("~/BINF6970 Statistical Bioinformatics/Assignment 2")



## Read and prepare data----
immunology <- read_excel("Resources/Immunologic profiles of patients with COVID-19.xlsx")

X <- model.matrix(~ . - 1, data = immunology[, c(2, 6:32)])
Y <- as.factor(immunology$Severirty)


# Split data into training and test data (75/25)

set.seed(1717)

train_index <- createDataPartition(Y, p = 0.75, list = FALSE)
X_train <- X[train_index, ]
X_test <- X[-train_index, ]
Y_train <- Y[train_index]
Y_test <- Y[-train_index]


# Cross-Validation----

# NOTE: I used "deviance" instead of "auc" because the R warning message said too few observations per fold when using "auc".  Try out both methods, though, if you wish.

# Set seed for cross-validation (for sanity purposes)
set.seed(1717)

# Define alpha search grid
alpha_grid <- seq(0, 1, .01)

# Create data frames for collecting cross-validation results
result_10_fold <- data.frame()
result_20_fold <- data.frame()


# 10-Fold
for (a in alpha_grid) {
  cv_10fold <- cv.glmnet(X_train, Y_train, 
                         alpha = a, 
                         nfolds = 10, 
                         family = "binomial", 
                         type.measure = "deviance", 
                         keep = TRUE)
  
  result_10_fold <- rbind(result_10_fold, data.frame(
    alpha = a, 
    lambda_min = cv_10fold$lambda.min, 
    lambda_1se = cv_10fold$lambda.1se, 
    cvm_min = min(cv_10fold$cvm), 
    cvm_1se = cv_10fold$cvm[which(cv_10fold$lambda == cv_10fold$lambda.1se)]
  ))
}

# Extract best alpha and lambda values (10-Fold)
best_alpha_min_10 <- result_10_fold$alpha[which.min(result_10_fold$cvm_min)]
best_alpha_1se_10 <- result_10_fold$alpha[which.min(result_10_fold$cvm_1se)]

best_lambda_min_10 <- result_10_fold$lambda_min[result_10_fold$alpha == best_alpha_min_10]
best_lambda_1se_10 <- result_10_fold$lambda_1se[result_10_fold$alpha == best_alpha_1se_10]


plot(cv_10fold, main = "Elastic Net (10-Fold) - CV Error")

# 20-Fold
for (a in alpha_grid) {
  cv_20fold <- cv.glmnet(X_train, Y_train, 
                         alpha = a, 
                         nfolds = 20, 
                         family = "binomial", 
                         type.measure = "deviance", 
                         keep = TRUE)
  
  result_20_fold <- rbind(result_20_fold, data.frame(
    alpha = a, 
    lambda_min = cv_20fold$lambda.min, 
    lambda_1se = cv_20fold$lambda.1se, 
    cvm_min = min(cv_20fold$cvm), 
    cvm_1se = cv_20fold$cvm[which(cv_20fold$lambda == cv_20fold$lambda.1se)]
  ))
}

# Extract best alpha and lambda values (20-Fold)
best_alpha_min_20 <- result_20_fold$alpha[which.min(result_20_fold$cvm_min)]
best_alpha_1se_20 <- result_20_fold$alpha[which.min(result_20_fold$cvm_1se)]

best_lambda_min_20 <- result_20_fold$lambda_min[result_20_fold$alpha == best_alpha_min_20]
best_lambda_1se_20 <- result_20_fold$lambda_1se[result_20_fold$alpha == best_alpha_1se_20]

plot(cv_20fold, main = "Elastic Net (20-Fold) - CV Error")


# Fitting the final Elastic-Net models based on training data

# 10-Fold
model_10fold_min <- glmnet(X_train, Y_train, 
                         alpha = best_alpha_min_10, 
                         lambda = best_lambda_min_10,
                         family = "binomial")


model_10fold_1se <- glmnet(X_train, Y_train, 
                         alpha = best_alpha_1se_10, 
                         lambda = best_lambda_1se_10,
                         family = "binomial")

# 20-Fold
model_20fold_min <- glmnet(X_train, Y_train, 
                         alpha = best_alpha_min_20, 
                         lambda = best_lambda_min_20,
                         family = "binomial")


model_20fold_1se <- glmnet(X_train, Y_train, 
                         alpha = best_alpha_1se_20, 
                         lambda = best_lambda_1se_20,
                         family = "binomial")
