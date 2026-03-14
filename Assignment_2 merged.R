##**************************
## BINF*6970 - Assignment 2
##
## Student Names: Cyrus Akbarally, Iroayo Toki, Sodiq Dada & William Feinman 
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
library("stringr")
library(tidyverse)
library(ggfortify)
## setwd("~/BINF6970 Statistical Bioinformatics/Assignment 2")

##**************************


#Problem 1
#1. Data acquisition and preprocessing----
#Load Data frame
load("Resources/geneexpression2.rda")
Df <- dat
Df
#Rownames are metadata to be converted before plotting PCA
row.names(Df)

#Preprocessing

#Converting rownames into cell type and status for visualization
df_clean <- Df %>%
  tibble::rownames_to_column("ID") %>%
  separate(ID, into = c("Status", "Cell_type", "Replicate"), sep = "_") %>%
  mutate(
    Status = str_remove_all(Status, "\\d+")
  ) %>% select(-Replicate)
#2. Computing and plotting PCA----
#Compute PCA
PCA <-  prcomp(df_clean[,-(1:2)], center = T)

# Auto plot PCA
#PLot shows distribution of dataset by cell type across PC1
autoplot( PCA , data = df_clean, colour = "Cell_type", shape = "Status", main = "PCA: PC1 vs PC2" )
#Visualizing eigenvectors
autoplot( PCA , data = df_clean, colour = "Cell_type", shape = "Status", loadings= T,main = "PCA: PC1 vs PC2 with Eigenvectors" )
#Adding confidence interval ellipses
autoplot( PCA , data = df_clean, colour = "Cell_type", shape = "Status", main = "PCA: PC1 vs PC2") + stat_ellipse(aes(group = Cell_type, colour = Cell_type), level = 0.95)




#Dataframe for PC1 vs PC3
pca_data <- data.frame(PC1 = PCA$x[,1], PC3 = PCA$x[,3])
pca_data <- pca_data %>% mutate(Status = df_clean$Status, 
                                Cell_type = df_clean$Cell_type)
#Plotting PC1 vs PC3
#PC3 shows what seems to be an even distribution among the classes
PC1_PC3 <- ggplot(pca_data, aes(x = PC1, y = PC3, colour = Cell_type, shape = Status)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC3",
       x = "PC1",
       y = "PC3") 
PC1_PC3
#Plot with ellipses 
PC1_PC3 + stat_ellipse(aes(group = Cell_type, colour = Cell_type), level = 0.95)


# PCA Justification using screeplot and result matrix----
## The sample covariance matrix
#Original data frame used as it has the original data structure without string variables
CVmatrix <- var(Df)
round( CVmatrix,2)
#Correlation Matrix
Crmatrix <-  cor(Df) 
## The eigenvalue-eigenvector pairs
eigen(CVmatrix)

#Sum of eigenvalues
sum(diag(CVmatrix)) 
sum(eigen(CVmatrix)$values)


## First, using eigen() directly:

eg.vals <- eigen(CVmatrix)$values
eg.vcs <- eigen(CVmatrix)$vectors 	## columns are eigenvectors

#Creating result matrix with eigenvectors, eigenvalues and cumulative variance for each PC
resultmat <- matrix(NA, ncol = 156, nrow = 158) # 156 columns because we have 156 variables

#Adding eigenvalues
resultmat[157,] <- eg.vals

## cumulative sum of PC variances 
resultmat[158,] <- cumsum(eg.vals/sum(eg.vals)*100) # % of variability explained by the PC

rownames(resultmat) <- c( colnames(Df)  , "lambda.hat"   , "cumulative %" )

colnames(resultmat) <- paste0("PC", 1:156)

resultmat[1:156,1:156] <- eg.vcs 
#PC1 and #PC2 account for 77% of the variation
print(resultmat[158,1:5 ]  )


### The Scree PLOT
#Plotting only a few PCs to better visualize elbow
#Elbow at PC2 showing that most of the variation is represented by PC1 and PC2
plot(as.ts(eg.vals[1:5]) , ylab = "lambda-hat" , xlab = "PC" , main = " The Scree Plot")
#Values for important PCs
resultmat[1:156,1:2]

##**************************


#Problem 2
## Read and prepare data----
immunology <- read_excel("Resources/Immunologic profiles of patients with COVID-19.xlsx")

X <- model.matrix(~ . - 1, data = immunology[, c(2, 6:32)])
Y <- as.factor(immunology$Severirty)



# Set seed for cross-validation (for sanity purposes)
set.seed(1717)

# Split data into training and test data (75/25)

train_index <- createDataPartition(Y, p = 0.75, list = FALSE)
X_train_prescaled <- X[train_index, ]
X_test_prescaled <- X[-train_index, ]
Y_train <- Y[train_index]
Y_test <- Y[-train_index]


# Standardize predictors (glmnet does this by default with standardize=TRUE)
# But we'll scale manually for consistency across methods
X_train <- scale(X_train_prescaled)
X_test <- scale(X_test_prescaled, 
                       center = attr(X_train, "scaled:center"),
                       scale = attr(X_train, "scaled:scale"))

# Remove any columns with NA values (if scaling creates them due to zero variance)
na_cols <- which(colSums(is.na(X_train)) > 0)
if (length(na_cols) > 0) {
  X_train_scaled <- X_train[, -na_cols]
  X_test_scaled <- X_test[, -na_cols]
  cat("Removed", length(na_cols), "columns with zero variance\n")
}

# Cross-Validation----

# NOTE: I used "deviance" instead of "auc" because the R warning message said too few observations per fold when using "auc".  Try out both methods, though, if you wish.

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
        cvm_1se = cv_10fold$cvm[which(cv_10fold$lambda == cv_10fold$lambda.1se)],
    #adding in nzero, as was included in original for loop
    nzero_min = cv_10fold$nzero[which(cv_10fold$lambda == cv_10fold$lambda.min)],
    nzero_1se = cv_10fold$nzero[which(cv_10fold$lambda == cv_10fold$lambda.1se)]
  ))
}

print(result_10_fold)

# Extract best alpha and lambda values (10-Fold) by min CV error
best_alpha_min_10 <- result_10_fold$alpha[which.min(result_10_fold$cvm_min)]
best_alpha_1se_10 <- result_10_fold$alpha[which.min(result_10_fold$cvm_1se)]

best_lambda_min_10 <- result_10_fold$lambda_min[result_10_fold$alpha == best_alpha_min_10]
best_lambda_1se_10 <- result_10_fold$lambda_1se[result_10_fold$alpha == best_alpha_1se_10]


#Plots for visualization and alpha grid search
plot(cv_10fold, main = "Elastic Net (10-Fold) - Deviance")


ggplot(result_10_fold, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha, 10-Fold",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "Deviance",
       color = "Lambda Type") +
  theme_minimal()

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
    cvm_1se = cv_20fold$cvm[which(cv_20fold$lambda == cv_20fold$lambda.1se)],
    #adding in nzero, as was included in original for loop
    nzero_min = cv_20fold$nzero[which(cv_20fold$lambda == cv_20fold$lambda.min)],
    nzero_1se = cv_20fold$nzero[which(cv_20fold$lambda == cv_20fold$lambda.1se)]
  ))
}


print(result_20_fold)

# Extract best alpha and lambda values (20-Fold) by min CV error
best_alpha_min_20 <- result_20_fold$alpha[which.min(result_20_fold$cvm_min)]
best_alpha_1se_20 <- result_20_fold$alpha[which.min(result_20_fold$cvm_1se)]

best_lambda_min_20 <- result_20_fold$lambda_min[result_20_fold$alpha == best_alpha_min_20]
best_lambda_1se_20 <- result_20_fold$lambda_1se[result_20_fold$alpha == best_alpha_1se_20]

#Plots for visualization and alpha grid search
plot(cv_20fold, main = "Elastic Net (20-Fold) - Deviance")

ggplot(result_20_fold, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "Deviance vs Alpha, 20-Fold",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "Deviance",
       color = "Lambda Type") +
  theme_minimal()



# Fitting the final Elastic-Net models based on training data

# 10-Fold
model_10fold_min <- glmnet(X_train, Y_train, 
                         alpha = best_alpha_min_10, 
                         lambda = best_lambda_min_10,
                         lambda_type = "min",
                         family = "binomial")
                        
                        

model_10fold_1se <- glmnet(X_train, Y_train, 
                         alpha = best_alpha_1se_10, 
                         lambda = best_lambda_1se_10,
                         lambda_type = "1se",
                         family = "binomial")

# 20-Fold
model_20fold_min <- glmnet(X_train, Y_train, 
                         alpha = best_alpha_min_20, 
                         lambda = best_lambda_min_20,
                         lambda_type = "min",
                         family = "binomial")


model_20fold_1se <- glmnet(X_train, Y_train, 
                         alpha = best_alpha_1se_20, 
                         lambda = best_lambda_1se_20,
                         lambda_type = "1se",
                         family = "binomial")


# Function to extract top N predictors from a model
extract_top_predictors <- function(model, n = 10, lambda_type = "min") {
  # Get coefficients
  if (lambda_type == "min") {
    coefs <- as.matrix(coef(model))[, 1]
  } else {
    # For models stored differently
    coefs <- as.matrix(coef(model, s = paste0("lambda.", lambda_type)))[, 1]
  }
  
  # Remove intercept
  coefs_no_intercept <- coefs[-1]
  
  # Get absolute values for ranking
  abs_coefs <- abs(coefs_no_intercept)
  
  # Get top N predictors
  top_indices <- order(abs_coefs, decreasing = TRUE)[1:min(n, length(abs_coefs))]
  
  # Create result data frame
  result <- data.frame(
    Predictor = names(coefs_no_intercept)[top_indices],
    Coefficient = coefs_no_intercept[top_indices],
    Absolute_Effect = abs_coefs[top_indices],
    Rank = 1:length(top_indices)
  )
  
  return(result)
}

# Extract top 11 predictors for each model. 11 selected to try and identify impact of both age and the top 10 biomarkers (if that many)
m10_top_min <- extract_top_predictors(model_10fold_min, n = 11, lambda_type = "min")
print(m10_top_min)

m10_top_lse <- extract_top_predictors(model_10fold_1se, n = 11, lambda_type = "1se")
print(m10_top_lse)

m20_top_min <- extract_top_predictors(model_20fold_min, n = 11, lambda_type = "min")
print(m20_top_min)

m20_top_lse <- extract_top_predictors(model_20fold_1se, n = 11, lambda_type = "1se")
print(m20_top_lse)


# Function to calculate performance metrics
calculate_metrics <- function(model, X_test, Y_test, lambda_type = "min") {
  # Make predictions
  if (lambda_type == "min") {
    predictions <- predict(model, newx = X_test)
  } else {
    predictions <- predict(model, newx = X_test, s = paste0("lambda.", lambda_type))
  }
  
  # Calculate metrics
  mse <- mean((Y_test - predictions)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(Y_test - predictions))
  r_squared <- 1 - sum((Y_test - predictions)^2) / sum((Y_test - mean(Y_test))^2)
  
  # Count non-zero coefficients (excluding intercept)
  if (lambda_type == "min") {
    n_coefs <- sum(coef(model)[-1] != 0)
  } else {
    n_coefs <- sum(coef(model, s = paste0("lambda.", lambda_type))[-1] != 0)
  }
  
  return(c(MSE = mse, RMSE = rmse, MAE = mae, R2 = r_squared, NonZero = n_coefs))
}


# Calculate metrics for all models. We attempted to get a comparison table working for further cross-model metrics, but we could only get Non-Zero coefficient comparison working before the deadline.

comparison <- data.frame()

# 10-Fold
comparison <- rbind(comparison, 
                    data.frame(Model = "10-Fold (lambda.min)", 
                               as.list(calculate_metrics(model_10fold_min, X_test, Y_test, "min"))))

comparison <- rbind(comparison, 
                    data.frame(Model = "10-Fold (lambda.1se)", 
                               as.list(calculate_metrics(model_10fold_1se, X_test, Y_test, "1se"))))

# 20-Fold
comparison <- rbind(comparison, 
                    data.frame(Model = "20-Fold (lambda.min)", 
                               as.list(calculate_metrics(model_20fold_min, X_test, Y_test, "min"))))

comparison <- rbind(comparison, 
                    data.frame(Model = "20-Fold (lambda.1se)", 
                               as.list(calculate_metrics(model_20fold_1se, X_test, Y_test, "1se"))))

print(comparison)

