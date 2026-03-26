# download GEMSS package
# remotes::install_github("szhua-stat/GEMSS")
library(tidyverse)
library(hetGP)
library(twinning)
library(Metrics)
library(GEMSS)
library(ContourFunctions)
contour_gp_pred <- function(X, Y, index, g = NULL, r = NULL, X_test, Y_test, Cov_Fun, ngrid = 50, zlim = NULL){
  gr1 <- seq(min(X[,1])-0.01,max(X[,1])+0.01, len = ngrid)
  gr2 <- seq(min(X[,2])-0.01,max(X[,2])+0.01, len = ngrid)
  grid <- expand.grid(x1 = gr1, x2 = gr2) %>% as.matrix()

  GP <- mleHomGP(X = X[index,], Z = Y[index], covtype = Cov_Fun)
  pred_test <- predict(GP, X_test)$mean
  nrmspe <-  rmse(pred_test, Y_test) / rmse(Y_test, mean(Y))

  R <- cov_gen(X[index,], theta = GP$theta, type = Cov_Fun)
  cond <- format(kappa(R, exact = T), digits=4)

  prep_grid <- predict(GP, grid)$mean
  cf_grid(gr1, gr2, matrix(prep_grid,ngrid,ngrid),
          bar=TRUE, main=paste0('nRMSPE: ', round(nrmspe, 4), ', cond_num: ', cond), las=0, cex.main = 1.5,zlim = zlim,
          afterplotfunc=function(){
            points(X[index,, drop = F], col=1, pch='.',cex=6)
            points(X[g,, drop = F], col='darkgrey', pch=20,cex=1.3)
            points(X[r,, drop = F], col='red', pch=4, cex=1.3, lwd = 2)

          })
}

setwd("real_estate_valuation_dataset")
Data <- read.csv('dataset.csv', header = T)
colnames(Data) <- c('No', 'transaction.date', 'house.age', 'nearest.MRT',
                    'n.convenience.stores','latitude','longitude','unit.price')

data <- Data[6:8] ; data$unit.price <- log(data$unit.price)

# ### original data-splitting code
# set.seed(1)
# # split data to training set and testing set
# index_test <- twin(data, r = 6)
# index_train <- setdiff(1:nrow(data), index_test)
# data_train <- data[index_train,]
# data_test <- data[index_test,]
#
# X <- data_train[,c(2,1)] %>% as.matrix() ; Y <- data_train[,3]
# X_test <- data_test[,c(2,1)] %>% as.matrix() ; Y_test <- data_test[,3]
#
# #Use twinning to select 20% of the training data to form the validation set.
# i_val <- twin(cbind(X,Y), r = 5)
# i_tr <- setdiff(1:nrow(X), i_val)

### For duplicate same figure in the paper
index_dup <- read.csv('index_for_reproduce.csv')[,-1]

# testing - training - validation
index_test <- index_dup$index_test

index_train <- setdiff(1:nrow(data), index_test)
data_train <- data[index_train,]
data_test <- data[index_test,]

X <- as.matrix(data_train[,c(2,1)]) ; Y <- data_train[,3]
X_test <- as.matrix(data_test[,c(2,1)]) ; Y_test <- data_test[,3]

i_val <- index_dup$i_val
i_tr <- setdiff(1:nrow(X), i_val)

Cov_Fun <- "Gaussian"

# define the grid
ngrid <- 80
gr1 <- seq(min(data$latitude)-0.01,max(data$latitude)+0.01, len = ngrid)
gr2 <- seq(min(data$longitude)-0.01,max(data$longitude)+0.01, len = ngrid)
grid <- expand.grid(x1 = gr1, x2 = gr2) %>% as.matrix()

# prediction performance of training + validation (fig.4a)
contour_gp_pred(X, Y, index = 1:nrow(X), X_test = X_test, Y_test = Y_test, Cov_Fun = Cov_Fun, zlim = c(2.5,4))

# prediction performance of training data (Fig.4b)
contour_gp_pred(X, Y, index = i_tr, g = i_val, X_test = X_test, Y_test = Y_test, Cov_Fun = Cov_Fun, zlim = c(2.5,4))

# Removing process
result <- gemss_remove(X[i_tr, ], Y[i_tr], X[i_val, ], Y[i_val], n_remove = 150, Cov_Fun)

# Predictive R-square plot (Fig.4d)
eval <- result$eval_matrix
par(mfrow = c(1,1), mar = c(4,4,2,2))
plot(eval[,1], eval[,3], type = 'b', pch = 19, yaxt = 'n', xaxt = 'n', xlab = '', ylab ='', cex = 0.8)
axis(1,cex.axis=1.3)
axis(2,cex.axis=1.3)
mtext('R2_pred', side=2, line=2.3, cex=1.4)
mtext('n removals', side=1, line=2, cex=1.4)
abline(v = which.max(eval[,3]), lty = 2, col = 2)
text(x=which.max(eval[,3])-0, y=par("usr")[3], which.max(eval[,3]), pos = 1, xpd = TRUE, col=2, cex = 1.45)

# GP prediction after the elimination of the red cross (Fig.4c)
final_index <- i_tr[result$index]
removal <- i_tr[result$remove]
contour_gp_pred(X, Y, index = final_index, r = removal, X_test = X_test, Y_test = Y_test, Cov_Fun = Cov_Fun, zlim = c(2.5,4))

