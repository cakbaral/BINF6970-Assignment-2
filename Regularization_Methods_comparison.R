# Regularization Methods:

# Load required libraries
library(lars)       # For diabetes dataset
library(glmnet)     # For regularization methods
library(caret)      # For data splitting and preprocessing
library(ggplot2)    # For visualization
library(dplyr)      # For data manipulation

# Set seed for reproducibility
set.seed(123)

# Load the diabetes dataset
data(diabetes)
# The diabetes dataset in lars has: x (matrix), y (response), x2 (expanded matrix)
# We'll use the expanded version (x2) which includes squared terms and interactions

# Prepare data
X <- as.matrix(diabetes$x2)  # Predictors (expanded features)
y <- diabetes$y              # Response variable

# Check dimensions
cat("Dataset dimensions:", dim(X), "\n")
cat("Number of predictors:", ncol(X), "\n")
cat("Number of observations:", length(y), "\n")

# Split data into training (70%) and testing (30%)
train_index <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_index, ]
X_test <- X[-train_index, ]
y_train <- y[train_index]
y_test <- y[-train_index]

# Standardize predictors (glmnet does this by default with standardize=TRUE)
# But we'll scale manually for consistency across methods
X_train_scaled <- scale(X_train)
X_test_scaled <- scale(X_test, 
                       center = attr(X_train_scaled, "scaled:center"),
                       scale = attr(X_train_scaled, "scaled:scale"))

# Remove any columns with NA values (if scaling creates them due to zero variance)
na_cols <- which(colSums(is.na(X_train_scaled)) > 0)
if (length(na_cols) > 0) {
  X_train_scaled <- X_train_scaled[, -na_cols]
  X_test_scaled <- X_test_scaled[, -na_cols]
  cat("Removed", length(na_cols), "columns with zero variance\n")
}

# ----------------------------------------------------------------
# 1. RIDGE REGRESSION (alpha = 0)
# ----------------------------------------------------------------
cat("\n=== RIDGE REGRESSION (alpha=0) ===\n")

# Cross-validation for ridge
cv_ridge <- cv.glmnet(X_train_scaled, y_train, 
                      alpha = 0, 
                      nfolds = 10,
                      type.measure = "mse")

# Plot cross-validation results
plot(cv_ridge, main = "Ridge Regression - CV Error")

# Extract lambda values
ridge_lambda_min <- cv_ridge$lambda.min
ridge_lambda_1se <- cv_ridge$lambda.1se

cat("Ridge - Lambda min:", ridge_lambda_min, "\n")
cat("Ridge - Lambda 1se:", ridge_lambda_1se, "\n")

# Fit final model with lambda.min
ridge_model_min <- glmnet(X_train_scaled, y_train, 
                          alpha = 0, 
                          lambda = ridge_lambda_min)

# Fit final model with lambda.1se
ridge_model_1se <- glmnet(X_train_scaled, y_train, 
                          alpha = 0, 
                          lambda = ridge_lambda_1se)

# ----------------------------------------------------------------
# 2. LASSO REGRESSION (alpha = 1)
# ----------------------------------------------------------------
cat("\n=== LASSO REGRESSION (alpha=1) ===\n")

# Cross-validation for lasso
cv_lasso <- cv.glmnet(X_train_scaled, y_train, 
                      alpha = 1, 
                      nfolds = 10,
                      type.measure = "mse")

plot(cv_lasso, main = "Lasso Regression - CV Error")

# Extract lambda values
lasso_lambda_min <- cv_lasso$lambda.min
lasso_lambda_1se <- cv_lasso$lambda.1se

cat("Lasso - Lambda min:", lasso_lambda_min, "\n")
cat("Lasso - Lambda 1se:", lasso_lambda_1se, "\n")

# Fit final models
lasso_model_min <- glmnet(X_train_scaled, y_train, 
                          alpha = 1, 
                          lambda = lasso_lambda_min)

lasso_model_1se <- glmnet(X_train_scaled, y_train, 
                          alpha = 1, 
                          lambda = lasso_lambda_1se)

# ----------------------------------------------------------------
# 3. ELASTIC NET (alpha = 0.5 as example)
# ----------------------------------------------------------------
cat("\n=== ELASTIC NET (alpha=0.5) ===\n")

# Cross-validation for elastic net with alpha=0.5
cv_enet <- cv.glmnet(X_train_scaled, y_train, 
                     alpha = 0.5, 
                     nfolds = 10,
                     type.measure = "mse")

plot(cv_enet, main = "Elastic Net (alpha=0.5) - CV Error")
# how many predictor left?
# Extract lambda values
enet_lambda_min <- cv_enet$lambda.min
enet_lambda_1se <- cv_enet$lambda.1se

cat("Elastic Net - Lambda min:", enet_lambda_min, "\n")
cat("Elastic Net - Lambda 1se:", enet_lambda_1se, "\n")

# Fit final models
enet_model_min <- glmnet(X_train_scaled, y_train, 
                         alpha = 0.5, 
                         lambda = enet_lambda_min)

enet_model_1se <- glmnet(X_train_scaled, y_train, 
                         alpha = 0.5, 
                         lambda = enet_lambda_1se)

# ----------------------------------------------------------------
# 4. COMPARISON OF DIFFERENT ALPHA VALUES (Grid Search)
# ----------------------------------------------------------------
cat("\n=== GRID SEARCH FOR OPTIMAL ALPHA ===\n")

# Define alpha grid
alpha_grid <- seq(0, 1, by = 0.1)

# Create container for results
results <- data.frame()

# Perform grid search
for (alpha_val in alpha_grid) {
  # Cross-validation for current alpha
  cv_model <- cv.glmnet(X_train_scaled, y_train, 
                        alpha = alpha_val, 
                        nfolds = 10,
                        type.measure = "mse",
                        keep = TRUE)
  
  # Store results
  results <- rbind(results, data.frame(
    alpha = alpha_val,
    lambda_min = cv_model$lambda.min,
    lambda_1se = cv_model$lambda.1se,
    cvm_min = min(cv_model$cvm),
    cvm_1se = cv_model$cvm[which(cv_model$lambda == cv_model$lambda.1se)],
    nzero_min = cv_model$nzero[which(cv_model$lambda == cv_model$lambda.min)],
    nzero_1se = cv_model$nzero[which(cv_model$lambda == cv_model$lambda.1se)]
  ))
}

# Print results
print(results)

# Find best alpha based on minimum CV error
best_alpha_min <- results$alpha[which.min(results$cvm_min)]
best_alpha_1se <- results$alpha[which.min(results$cvm_1se)]

cat("\nBest alpha (lambda.min):", best_alpha_min, "\n")
cat("Best alpha (lambda.1se):", best_alpha_1se, "\n")

# Visualize alpha grid search
ggplot(results, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "CV Error (MSE)",
       color = "Lambda Type") +
  theme_minimal()

# ----------------------------------------------------------------
# 5. EXTRACT COEFFICIENTS FOR TOP PREDICTORS
# ----------------------------------------------------------------
cat("\n=== TOP PREDICTORS FOR EACH MODEL ===\n")

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

# Extract top predictors for each model
cat("\n--- Top 10 Predictors (Lambda.min) ---\n")

cat("\nRidge (lambda.min):\n")
ridge_top_min <- extract_top_predictors(ridge_model_min, n = 10, lambda_type = "min")
print(ridge_top_min)

cat("\nLasso (lambda.min):\n")
lasso_top_min <- extract_top_predictors(lasso_model_min, n = 10, lambda_type = "min")
print(lasso_top_min)

cat("\nElastic Net (lambda.min):\n")
enet_top_min <- extract_top_predictors(enet_model_min, n = 10, lambda_type = "min")
print(enet_top_min)

cat("\n--- Top 10 Predictors (Lambda.1se) ---\n")

cat("\nRidge (lambda.1se):\n")
ridge_top_1se <- extract_top_predictors(ridge_model_1se, n = 10, lambda_type = "min")
print(ridge_top_1se)

cat("\nLasso (lambda.1se):\n")
lasso_top_1se <- extract_top_predictors(lasso_model_1se, n = 10, lambda_type = "min")
print(lasso_top_1se)

cat("\nElastic Net (lambda.1se):\n")
enet_top_1se <- extract_top_predictors(enet_model_1se, n = 10, lambda_type = "min")
print(enet_top_1se)

# ----------------------------------------------------------------
# 6. MODEL COMPARISON ON TEST SET
# ----------------------------------------------------------------
cat("\n=== MODEL COMPARISON (TEST SET PERFORMANCE) ===\n")

# Function to calculate performance metrics
calculate_metrics <- function(model, X_test, y_test, lambda_type = "min") {
  # Make predictions
  if (lambda_type == "min") {
    predictions <- predict(model, newx = X_test)
  } else {
    predictions <- predict(model, newx = X_test, s = paste0("lambda.", lambda_type))
  }
  
  # Calculate metrics
  mse <- mean((y_test - predictions)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(y_test - predictions))
  r_squared <- 1 - sum((y_test - predictions)^2) / sum((y_test - mean(y_test))^2)
  
  # Count non-zero coefficients (excluding intercept)
  if (lambda_type == "min") {
    n_coefs <- sum(coef(model)[-1] != 0)
  } else {
    n_coefs <- sum(coef(model, s = paste0("lambda.", lambda_type))[-1] != 0)
  }
  
  return(c(MSE = mse, RMSE = rmse, MAE = mae, R2 = r_squared, NonZero = n_coefs))
}

# Calculate metrics for all models
comparison <- data.frame()

# Ridge
comparison <- rbind(comparison, 
                    data.frame(Model = "Ridge (lambda.min)", 
                               as.list(calculate_metrics(ridge_model_min, X_test_scaled, y_test, "min"))))

comparison <- rbind(comparison, 
                    data.frame(Model = "Ridge (lambda.1se)", 
                               as.list(calculate_metrics(ridge_model_1se, X_test_scaled, y_test, "min"))))

# Lasso
comparison <- rbind(comparison, 
                    data.frame(Model = "Lasso (lambda.min)", 
                               as.list(calculate_metrics(lasso_model_min, X_test_scaled, y_test, "min"))))

comparison <- rbind(comparison, 
                    data.frame(Model = "Lasso (lambda.1se)", 
                               as.list(calculate_metrics(lasso_model_1se, X_test_scaled, y_test, "min"))))

# Elastic Net
comparison <- rbind(comparison, 
                    data.frame(Model = "ElasticNet (lambda.min)", 
                               as.list(calculate_metrics(enet_model_min, X_test_scaled, y_test, "min"))))

comparison <- rbind(comparison, 
                    data.frame(Model = "ElasticNet (lambda.1se)", 
                               as.list(calculate_metrics(enet_model_1se, X_test_scaled, y_test, "min"))))

print(comparison)

# ----------------------------------------------------------------
# 7. VISUALIZATION OF COEFFICIENT PATHS
# ----------------------------------------------------------------
cat("\n=== COEFFICIENT PATHS ===\n")

# Create a plot showing coefficient paths for different alphas
par(mfrow = c(2, 2))

# Ridge coefficient path
ridge_full <- glmnet(X_train_scaled, y_train, alpha = 0)
plot(ridge_full, xvar = "lambda", main = "Ridge (alpha=0) - Coefficient Paths")

# Lasso coefficient path
lasso_full <- glmnet(X_train_scaled, y_train, alpha = 1)
plot(lasso_full, xvar = "lambda", main = "Lasso (alpha=1) - Coefficient Paths")

# Elastic Net coefficient path
enet_full <- glmnet(X_train_scaled, y_train, alpha = 0.5)
plot(enet_full, xvar = "lambda", main = "Elastic Net (alpha=0.5) - Coefficient Paths")

# Compare coefficient magnitudes across models
coef_comparison <- data.frame(
  Predictor = colnames(X_train_scaled),
  Ridge = as.numeric(coef(ridge_model_min)[-1]),  # Remove intercept
  Lasso = as.numeric(coef(lasso_model_min)[-1]),
  ElasticNet = as.numeric(coef(enet_model_min)[-1])
)

# Select top 20 predictors by variance across models
coef_comparison$Variance <- apply(coef_comparison[, -1], 1, var)
top_predictors <- coef_comparison[order(-coef_comparison$Variance), ][1:20, ]

# Melt for plotting
library(reshape2)
top_melted <- melt(top_predictors[, 1:4], id.vars = "Predictor")

# Plot coefficient comparison
ggplot(top_melted, aes(x = Predictor, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Coefficient Comparison Across Models (Top 20)",
       y = "Coefficient Value",
       x = "Predictor",
       fill = "Model")

# ----------------------------------------------------------------
# 8. LAMBDA.MIN VS LAMBDA.1SE EXPLANATION
# ----------------------------------------------------------------
cat("\n=== EXPLANATION: lambda.min vs lambda.1se ===\n")
cat("\nlambda.min: The value of lambda that gives minimum cross-validation error.\n")
cat("  - Typically results in more complex model with more predictors\n")
cat("  - May overfit to the training/validation data\n")
cat("  - Lower bias but higher variance\n\n")

cat("lambda.1se: The largest value of lambda within 1 standard error of the minimum.\n")
cat("  - Results in a simpler, more regularized model\n")
cat("  - More robust and generalizable\n")
cat("  - Higher bias but lower variance\n")
cat("  - Often preferred for model interpretation and deployment\n")

# Show the difference in number of predictors
cat("\nNumber of non-zero coefficients (excluding intercept):\n")
cat("Ridge lambda.min:", sum(coef(ridge_model_min)[-1] != 0), "\n")
cat("Ridge lambda.1se:", sum(coef(ridge_model_1se)[-1] != 0), "\n")
cat("Lasso lambda.min:", sum(coef(lasso_model_min)[-1] != 0), "\n")
cat("Lasso lambda.1se:", sum(coef(lasso_model_1se)[-1] != 0), "\n")
cat("ElasticNet lambda.min:", sum(coef(enet_model_min)[-1] != 0), "\n")
cat("ElasticNet lambda.1se:", sum(coef(enet_model_1se)[-1] != 0), "\n")

# Reset par
par(mfrow = c(1, 1))

# ----------------------------------------------------------------
# 9. SUMMARY
# ----------------------------------------------------------------
cat("\n=== SUMMARY ===\n")
cat("1. Ridge (alpha=0) shrinks coefficients but keeps all predictors.\n")
cat("2. Lasso (alpha=1) performs variable selection by zeroing some coefficients.\n")
cat("3. Elastic Net (alpha between 0 and 1) balances ridge and lasso properties.\n")
cat("4. Grid search across alpha helps find optimal regularization type.\n")
cat("5. lambda.1se typically yields simpler models than lambda.min.\n")
cat("6. Top predictors can be extracted by sorting absolute coefficient values.\n")

