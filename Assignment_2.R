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
immunology <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")

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

set.seed(1717)
alpha_grid <- seq(0, 1, .01)


# 10-Fold
cv_10fold <- cv.glmnet_(X_train, Y_train, 
                        alpha = 0.5, 
                        nfolds = 10, 
                        family = "binomial", 
                        type.measure = "auc")

# 20-Fold
cv_20fold <- cv.glmnet(X_train, Y_train, 
                       alpha = 0.5, 
                       nfolds = 20, 
                       family = "binomial", 
                       type.measure = "auc")