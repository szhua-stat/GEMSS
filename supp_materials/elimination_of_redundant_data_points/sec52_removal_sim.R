# download GEMSS package
# remotes::install_github("szhua-stat/GEMSS")
library(tidyverse)
library(hetGP)
library(twinning)
library(Metrics)
library(GEMSS)
library(ContourFunctions) # packages for contour plot
dropwave = function(xx){
  #  DROP-WAVE FUNCTION   http://www.sfu.ca/~ssurjano/drop.html
  mini = c(-2,-2)
  maxi = c(2,2)
  x1 = mini[1] + xx[1]*(maxi[1] - mini[1])
  x2 = mini[2] + xx[2]*(maxi[2] - mini[2])
  frac1 <- 1 + cos(12*sqrt(x1^2+x2^2))
  frac2 <- 0.5*(x1^2+x2^2) + 2
  y <- -frac1/frac2
  return(y)
}
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
##### Dropwave dataset simulation #####
# original generated code
# set.seed(1)
#
# # generate dropwave dataset of size-400
# dat <- generate_data(SNR = 20, X_type = 'Uniform', n_train = 400, n_test = 10000, fun_type = 'DROPWAVE', rho_range = c(0,0))
# X <- dat$X ; Y <- dat$Y
# X_test <- dat$X_test ; Y_test <- dat$Y_test
#
# # 1/10 validation and 9/10 training index
# i_val <- twin(cbind(X,Y), r = 10)
# i_tr <- setdiff(1:nrow(X), i_val)

### For duplicate same figures in the paper
data_train <- as.matrix(read.csv('data_train.csv')[,-1])
data_test <- as.matrix(read.csv('data_test.csv')[,-1])
X <- data_train[,1:2] ; Y <- data_train[,3]
X_test <- data_test[,1:2] ; Y_test <- data_test[,3]

# validation index
i_val <- read.csv('validation_index.csv')[,-1]
i_tr <- setdiff(1:nrow(X), i_val)

# kernel
Cov_Fun <- 'Matern5_2'

ngrid <- 100
grid.x <- seq(0,1, len = ngrid)
grid <- expand.grid(x1 = grid.x, x2 = grid.x) %>% as.matrix()

# contour plot for dropwave function (Fig.3a)
f_value <- grid %>% apply(1, dropwave)
cf_grid(grid.x,grid.x, matrix(f_value,ngrid, ngrid), las = 0, zlim = c(-1,0.2),
        bar=TRUE, main='')

# prediction performance of training + validation (Fig.3b)
contour_gp_pred(X, Y, index = 1:nrow(X), X_test = X_test, Y_test = Y_test, Cov_Fun = Cov_Fun, zlim = c(-1,0.2))

# prediction performance of training data (Fig.4a)
contour_gp_pred(X, Y, index = i_tr, g = i_val, X_test = X_test, Y_test = Y_test, Cov_Fun = Cov_Fun, zlim = c(-1,0.2))

# Removing process
result <- gemss_remove(X[i_tr, ], Y[i_tr], X[i_val, ], Y[i_val], n_remove = 100, Cov_Fun)

# Predictive R-square plot (Fig.4c)
eval <- result$eval_matrix
par(mfrow = c(1,1), mar = c(4,4,2,2))
plot(eval[,1], eval[,3], type = 'b', pch = 19, yaxt = 'n', xaxt = 'n', xlab = '', ylab ='', cex = 0.8)
axis(1,cex.axis=1.3)
axis(2,cex.axis=1.3)
mtext('R2_pred', side=2, line=2.3, cex=1.4)
mtext('n removals', side=1, line=2, cex=1.4)
abline(v = which.max(eval[,3]), lty = 2, col = 2)
text(x=which.max(eval[,3])-0, y=par("usr")[3], which.max(eval[,3]), pos = 1, xpd = TRUE, col=2, cex = 1.45)

# GP prediction after the elimination of the red cross (Fig.4b)
final_index <- i_tr[result$index]
removal <- i_tr[result$remove]
contour_gp_pred(X, Y, index = final_index, r = removal, X_test = X_test, Y_test = Y_test, Cov_Fun = Cov_Fun, zlim = c(-1,0.2))
