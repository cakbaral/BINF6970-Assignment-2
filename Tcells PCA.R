library("dplyr")
library("stringr")
library(tidyverse)
library(ggplot2)
library(GGally)
library(ggfortify)
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
?autoplot()
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

