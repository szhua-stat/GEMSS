library(MultiRNG)
library(mixAK)
library(tmvtnorm)
library(Metrics)

# Chang(2023)
Index_EI_MSEPV_Sum = function(X, Y, index_now, index_left, num, Cov_Fun, sub_GP = NULL){
  X_sub = X[index_now,] ; X_left = X[index_left,]
  Y_sub = Y[index_now]  ; Y_left = Y[index_left]
  if(is.null(sub_GP)){
    sub_GP = mleHomGP(X_sub, Y_sub, covtype = Cov_Fun)
  }
  Obj_pred = predict(sub_GP, X_left)
  sub_GP_pred = Obj_pred$mean
  sub_GP_predVar = Obj_pred$sd2 + Obj_pred$nugs
  MSE = (sub_GP_pred - Y_left)^2
  Sum = MSE + sub_GP_predVar

  index_add = sort(Sum, decreasing = TRUE, index.return=TRUE)$ix[1:num]
  index = union(index_now,index_left[index_add])
  return(list(index = index))
}
EI_MSEPV_Sum_val = function(X, Y, ns, init_ind, Cov_Fun, print_result = T){
  iter <- 0 ; d <- ncol(X)

  n = nrow(X)
  index_now = init_ind
  index_left = setdiff(1:n,index_now)

  sub_GP = mleHomGP(X[index_now,], Y[index_now], covtype = Cov_Fun)
  sub_GPP_pred = predict(sub_GP, X[index_left,])$mean
  test_error_0 = mean((Y[index_left] - sub_GPP_pred)^2)

  test_error_old = test_error_0
  k = i = 1
  R_sq = size = c()
  while( length(index_now) < ns ){
    iter <- iter + 1

    res_select = Index_EI_MSEPV_Sum(X, Y, index_now, index_left, k, Cov_Fun, sub_GP)
    index_now = res_select$index
    index_left = setdiff(1:n, index_now)

    sub_GP = mleHomGP(X[index_now,], Y[index_now], covtype = Cov_Fun)
    sub_GPP_pred = predict(sub_GP, X[index_left,])$mean
    test_error_now = mean((Y[index_left] - sub_GPP_pred)^2)
    # Decide the update size a_s
    if( test_error_now < test_error_old ){
      k = min(2*k, ns - length(index_now)) # k = min(3*k, ns - length(index_now)) floor(ns/50)
    }else if(test_error_now >= test_error_old){
      k = min(max(ceiling(k/2), 1), ns - length(index_now)) # k = min(max(ceiling(k/3), 1), ns - length(index_now))
    }
    R_sq[i] = 1-test_error_now/test_error_0
    size[i] = length(index_now)
    if(print_result){
      print(paste0(round(size[i]),' ; ', round(R_sq[i],4)))
    }
    i = i + 1
    test_error_old = test_error_now
  }
  return(list(index = index_now, R2 = data.frame(size = size, R_sq = R_sq)))
}

# ASMECr
Est_grad=function(X,y){     # calculate the bandwidth by assuming additive model
  if(is.null(nrow(X))){
    X <- as.matrix(X)}
  alpha=mean(y)
  pp=ncol(X)
  B=rep(0,pp)
  if(pp >= 3){
    Ff=matrix(0,nrow(X),pp)
    for(i in 1:3){
      for(j in 1:pp){
        # backfitting
        yy=y-rep(alpha,nrow(X))-apply(Ff[,-j],1,sum)
        bw=npregbw(xdat=X[,j],ydat=yy,regtype="ll",ckertype="epanechnikov",bwmethod="cv.aic")$bw
        regres=npreg(bws=bw,txdat=X[,j],tydat=yy,ckertype="epanechnikov",gradients=TRUE)
        B[j]=bw
        Ff[,j]=regres$eval[[1]]
        Ff[,j]=Ff[,j]-mean(Ff[,j])
      }
    }
  }else{
    bw=npregbw(xdat=X,ydat=y,regtype="ll",ckertype="epanechnikov",bwmethod="cv.aic")$bw
    regres=npreg(bws=bw,txdat=X,tydat=y,ckertype="epanechnikov",gradients=TRUE)
    B=bw
  }
  return(B)
}
MEDh_selectr=function(X,y,n){
  if(nrow(X)>1e3){
    ib=sample(1:nrow(X),1e3)
  }
  else{
    ib=1:nrow(X)
  }
  bw=Est_grad(X[ib,],y[ib])
  if(nrow(X)>1e4){

    id=sample(1:nrow(X),1e4,replace = FALSE)
  }
  else{
    id=1:nrow(X)
  }
  regres=npreg(bws=bw,txdat=X[id,],tydat=y[id],ckertype="epanechnikov",gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(apply(grad^2,1,sum))

  MEDres=mined::SelectMinED(as.matrix(X)[id,],log(cgrad),n,1,2)
  MEDloc=MEDres$points
  iMED=subsample(as.matrix(X),MEDloc)

  return(iMED)
}

# assessing function
hms_span <- function(start, end) {
  # calculating time consuming
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}
RMSE_GP= function(X_sub, Y_sub, X_test, Y_test, Cov_Fun){
  # calculating root mean square error
  tryCatch({
    sub_GP = mleHomGP(X_sub, Y_sub, covtype = Cov_Fun)
    sub_GP_pred = predict(sub_GP, X_test)$mean
    sub_rmse = rmse(Y_test, sub_GP_pred)
    return(sub_rmse)
  },error = function(err) {
    return(NA)
  })
}


# fun_type = c('michalewicz', 'shubert', 'dropwave', 'borehole', 'welch',
#              'piston', 'wingweight', 'g_fun', 'currin', 'brat',
#              'roos', 'levy_13', 'bukin', 'otl', 'robot')

generate_data <- function(SNR, X_type, n_train, n_test, fun_type, rho_range = c(0,0)){
  if(length(rho_range) == 2){
    rho_min <- rho_range[1] ; rho_max <- rho_range[2]
  }else if(length(rho_range) == 1){
    rho_min <- rho_max <- rho_range
  }

  d <- switch(fun_type, michalewicz = 2, shubert = 2, dropwave = 2, borehole = 8, welch = 20, piston = 7, wingweight = 10,
              g_fun = dim_Gfun, currin = 2, brat = dim_BRAT, roos = dim_ROOS, levy_13 = 2, bukin = 2, otl = 6, robot = 8)

  f <- switch(fun_type, michalewicz = michalewicz, shubert = shubert, dropwave = dropwave, borehole = borehole, welch = welch, piston = piston, wingweight = wingweight,
              g_fun = gfunc, currin = curretal88exp, brat = bratleyetal92, roos = roosarn63, levy_13 = levy13, bukin = bukin6, otl = otlcircuit, robot = robot)

  if(X_type == 'Uniform'){
    ## Uniform ##
    pho = runif(d*(d-1)/2, min=rho_min, max=rho_max)
    cov_matrix <- diag(1/2,d)
    cov_matrix[upper.tri(cov_matrix)] <- pho
    cov_matrix <- cov_matrix + t(cov_matrix)

    X = draw.d.variate.uniform(no.row=n_train, d=d, cov.mat=cov_matrix)
    X_test = draw.d.variate.uniform(no.row=n_test, d=d, cov.mat=cov_matrix)

  }else if(X_type == 'MVNormal'){
    ## Mixture Normal ##
    r_1 = runif(d^2, min=rho_min, max=rho_max)
    cov_matrix_1 = (matrix(r_1, d, d) + t(matrix(r_1, d, d)))/2
    r_2 = runif(d^2, min=rho_min, max=rho_max)
    cov_matrix_2 = (matrix(r_2, d, d) + t(matrix(r_2, d, d)))/2
    diag(cov_matrix_1) = 1 ; diag(cov_matrix_2) = 1
    SSigma <- list(cov_matrix_1, cov_matrix_2)

    X = rMVNmixture(n_train, weight=c(0.5,0.5), mean=matrix(c(rep(-3,d),rep(3,d)),nrow=2,ncol=d,byrow = T), Sigma=SSigma)
    X = apply(X, 2, FUN = function(t){(t-min(t))/(max(t)-min(t))})

    X_test = rMVNmixture(n_test, weight=c(0.5,0.5), mean=matrix(c(rep(-3,d),rep(3,d)),nrow=2,ncol=d,byrow = T), Sigma=SSigma)
    X_test = apply(X_test, 2, FUN = function(t){(t-min(t))/(max(t)-min(t))} )

  }else if(X_type == 'Mvt'){
    ## Truncated correlated T ##
    pho = runif(d^2, min=rho_min, max=rho_max)
    cov_matrix = (matrix(pho, d, d) + t(matrix(pho, d, d)))/2
    diag(cov_matrix) = 1
    # truncate multivariate t distribution
    X = rtmvt(n=n_train, mean=rep(0,d), sigma = cov_matrix, df=10, lower = rep(-2, d), upper = rep(2, d))
    X = apply(X, 2, FUN = function(t){(t-min(t))/(max(t)-min(t))} )

    X_test = rtmvt(n=n_test, mean=rep(0,d), sigma =cov_matrix, df=10)
    X_test = apply(X_test, 2, FUN = function(t){(t-min(t))/(max(t)-min(t))} )
    ####
  }
  Y_signal = apply(as.matrix(X), 1, f)
  sigma = sqrt(var(Y_signal)/SNR)
  error = rnorm(n_train, mean=0, sd=sigma)
  Y = Y_signal + error
  Y_test = apply(as.matrix(X_test), 1, f)
  return(list(X = X,Y = Y, X_test = X_test, Y_test = Y_test))
}

{
  michalewicz = function(xx){
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2)
  # m = constant (optional), with default value 10
  #
  ##########################################################################
  mini = c(0,0)
  maxi = c(pi,pi)
  x1 = mini[1] + xx[1]*(maxi[1] - mini[1])
  x2 = mini[2] + xx[2]*(maxi[2] - mini[2])
  m = 10
  ii <- c(1:2)
  sum <- sin(x1) * (sin(1*x1^2/pi))^(2*m) + sin(x2) * (sin(2*x2^2/pi))^(2*m)
  y <- -sum
  return(y)
}
  shubert = function(xx){
    mini = c(-5.12,-5.12)
    maxi = c(5.12,5.12)
    #mini = c(-10,-10)
    #maxi = c(10,10)
    x1 <- mini[1] + xx[1]*(maxi[1] - mini[1])
    x2 <- mini[2] + xx[2]*(maxi[2] - mini[2])
    ii <- c(1:5)
    sum1 <- sum(ii * cos((ii+1)*x1+ii))
    sum2 <- sum(ii * cos((ii+1)*x2+ii))
    y <- sum1 * sum2
    return(y)
  }
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
  borehole = function(xx){ #  BOREHOLE FUNCTION d=8
    mini = c(0.05,100,63070,990,63.1,700,1120,9855)
    maxi = c(0.15,50000,115600,1110,116,820,1680,12045)
    rw <- mini[1] + xx[1]*(maxi[1] - mini[1])
    r  <- mini[2] + xx[2]*(maxi[2] - mini[2])
    Tu <- mini[3] + xx[3]*(maxi[3] - mini[3])
    Hu <- mini[4] + xx[4]*(maxi[4] - mini[4])
    Tl <- mini[5] + xx[5]*(maxi[5] - mini[5])
    Hl <- mini[6] + xx[6]*(maxi[6] - mini[6])
    L  <- mini[7] + xx[7]*(maxi[7] - mini[7])
    Kw <- mini[8] + xx[8]*(maxi[8] - mini[8])
    frac1 <- 2 * pi * Tu * (Hu-Hl)
    frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
    frac2b <- Tu / Tl
    frac2 <- log(r/rw) * (1+frac2a+frac2b)
    y <- frac1 / frac2
    return(y)
  }
  welch = function(xx){ #  WELCH d=20
    mini = rep(-0.5,20)
    maxi = rep(0.5,20)
    x1 <- mini[1] + xx[1]*(maxi[1] - mini[1])
    x2  <- mini[2] + xx[2]*(maxi[2] - mini[2])
    x3 <- mini[3] + xx[3]*(maxi[3] - mini[3])
    x4 <- mini[4] + xx[4]*(maxi[4] - mini[4])
    x5 <- mini[5] + xx[5]*(maxi[5] - mini[5])
    x6 <- mini[6] + xx[6]*(maxi[6] - mini[6])
    x7  <- mini[7] + xx[7]*(maxi[7] - mini[7])
    x8 <- mini[8] + xx[8]*(maxi[8] - mini[8])
    x9 <- mini[9] + xx[9]*(maxi[9] - mini[9])
    x10 <- mini[10] + xx[10]*(maxi[10] - mini[10])
    x11 <- mini[11] + xx[11]*(maxi[11] - mini[11])
    x12 <- mini[12] + xx[12]*(maxi[12] - mini[12])
    x13 <- mini[13] + xx[13]*(maxi[13] - mini[13])
    x14 <- mini[14] + xx[14]*(maxi[14] - mini[14])
    x15 <- mini[15] + xx[15]*(maxi[15] - mini[15])
    x16 <- mini[16] + xx[16]*(maxi[16] - mini[16])
    x17 <- mini[17] + xx[17]*(maxi[17] - mini[17])
    x18 <- mini[18] + xx[18]*(maxi[18] - mini[18])
    x19 <- mini[19] + xx[19]*(maxi[19] - mini[19])
    x20 <- mini[20] + xx[20]*(maxi[20] - mini[20])
    term1 <- 5*x12 / (1+x1)
    term2 <- 5 * (x4-x20)^2
    term3 <- x5 + 40*x19^3 - 5*x19
    term4 <- 0.05*x2 + 0.08*x3 - 0.03*x6
    term5 <- 0.03*x7 - 0.09*x9 - 0.01*x10
    term6 <- -0.07*x11 + 0.25*x13^2 - 0.04*x14
    term7 <- 0.06*x15 - 0.01*x17 - 0.03*x18
    y <- term1 + term2 + term3 + term4 + term5 + term6 + term7
    return(y)
  }
  piston = function(xx){ #  PISTON SIMULATION d=7
    mini = c(30,0.005,0.002,1000,90000,290,340)
    maxi = c(60,0.02,0.01,5000,110000,296,360)
    M <- mini[1] + xx[1]*(maxi[1] - mini[1])
    S <- mini[2] + xx[2]*(maxi[2] - mini[2])
    V0  <- mini[3] + xx[3]*(maxi[3] - mini[3])
    k <- mini[4] + xx[4]*(maxi[4] - mini[4])
    P0 <- mini[5] + xx[5]*(maxi[5] - mini[5])
    Ta <- mini[6] + xx[6]*(maxi[6] - mini[6])
    T0 <- mini[7] + xx[7]*(maxi[7] - mini[7])
    Aterm1 <- P0 * S
    Aterm2 <- 19.62 * M
    Aterm3 <- -k*V0 / S
    A <- Aterm1 + Aterm2 + Aterm3
    Vfact1 <- S / (2*k)
    Vfact2 <- sqrt(A^2 + 4*k*(P0*V0/T0)*Ta)
    V <- Vfact1 * (Vfact2 - A)
    fact1 <- M
    fact2 <- k + (S^2)*(P0*V0/T0)*(Ta/(V^2))
    C <- 2 * pi * sqrt(fact1/fact2)
    return(C)
  }
  wingweight = function(xx){  # WING WEIGHT FUNCTION  d=10
    mini = c(150,220,6,-10,16,0.5,0.08,2.5,1700,0.025)
    maxi = c(200,300,10,10,45,1,0.18,6,2500,0.08)
    Sw <- mini[1] + xx[1]*(maxi[1] - mini[1])
    Wfw  <- mini[2] + xx[2]*(maxi[2] - mini[2])
    A <- mini[3] + xx[3]*(maxi[3] - mini[3])
    LamCaps <- (mini[4] + xx[4]*(maxi[4] - mini[4]))* (pi/180)
    q <- mini[5] + xx[5]*(maxi[5] - mini[5])
    lam <- mini[6] + xx[6]*(maxi[6] - mini[6])
    tc  <- mini[7] + xx[7]*(maxi[7] - mini[7])
    Nz <- mini[8] + xx[8]*(maxi[8] - mini[8])
    Wdg <- mini[9] + xx[9]*(maxi[9] - mini[9])
    Wp <- mini[10] + xx[10]*(maxi[10] - mini[10])
    fact1 <- 0.036 * Sw^0.758 * Wfw^0.0035
    fact2 <- (A / ((cos(LamCaps))^2))^0.6
    fact3 <- q^0.006 * lam^0.04
    fact4 <- (100*tc / cos(LamCaps))^(-0.3)
    fact5 <- (Nz*Wdg)^0.49
    term1 <- Sw * Wp
    y <- fact1*fact2*fact3*fact4*fact5 + term1
    return(y)
  }

  dim_Gfun = 10
  gfunc = function(xx, a=(c(1:length(xx))-1)/2){
    ##########################################################################
    #
    # G-FUNCTION
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # INPUTS:
    #
    # xx = c(x1, x2, ..., xd)
    # a = c(a1, a2, ..., ad) (optional), with default values given by Crestaux
    #     et al. (2007)
    #
    ##########################################################################

    new1 <- abs(4*xx-2) + a
    new2 <- 1 + a
    prod <- prod(new1/new2)

    y <- prod
    return(y)
  }

  curretal88exp <- function(xx){
    ##########################################################################
    #
    # CURRIN ET AL. (1988) EXPONENTIAL FUNCTION
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # INPUT:
    # xx = c(x1, x2)
    #
    #########################################################################

    x1 <- xx[1]
    x2 <- xx[2]

    fact1 <- 1 - exp(-1/(2*x2))
    fact2 <- 2300*x1^3 + 1900*x1^2 + 2092*x1 + 60
    fact3 <- 100*x1^3 + 500*x1^2 + 4*x1 + 20

    y <- fact1 * fact2/fact3
    return(y)
  }

  dim_BRAT = 5
  bratleyetal92 <- function(xx){
    ##########################################################################
    #
    # BRATLEY ET AL. (1992) FUNCTION
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # INPUT:
    #
    # xx = c(x1, x2, ..., xd)
    #
    #########################################################################

    d <- length(xx)
    ii <- c(1:d)

    xxmat <- matrix(rep(xx,times=d), d, d, byrow=TRUE)
    xxmatlow <- xxmat
    xxmatlow[upper.tri(xxmatlow)] <- 1

    prod <- apply(xxmatlow, 1, prod)
    xxmatlow[upper.tri(xxmatlow)] <- 0
    sum <- sum(prod*(-1)^ii)

    y <- sum
    return(y)
  }

  dim_ROOS = 15
  roosarn63 <- function(xx){
    ##########################################################################
    #
    # ROOS & ARNOLD (1963) FUNCTION
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # INPUT:
    #
    # xx = c(x1, x2, ..., xd)
    #
    ##########################################################################

    prod <- prod(abs(4*xx - 2))

    y <- prod
    return(y)
  }

  levy13 <- function(xx){
    ##########################################################################
    #
    # LEVY FUNCTION N. 13
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # INPUT:
    #
    # xx = c(x1, x2)
    #
    ##########################################################################
    mini = c(-10,-10)
    maxi = c(10,10)
    x1 <- mini[1] + xx[1]*(maxi[1] - mini[1])
    x2  <- mini[2] + xx[2]*(maxi[2] - mini[2])

    term1 <- (sin(3*pi*x1))^2
    term2 <- (x1-1)^2 * (1+(sin(3*pi*x2))^2)
    term3 <- (x2-1)^2 * (1+(sin(2*pi*x2))^2)

    y <- term1 + term2 + term3
    return(y)
  }

  bukin6 <- function(xx){
    ##########################################################################
    #
    # BUKIN FUNCTION N. 6
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # INPUT:
    #
    # xx = c(x1, x2)
    #
    ##########################################################################
    mini = c(-15,-3)
    maxi = c(-5,3)
    x1 <- mini[1] + xx[1]*(maxi[1] - mini[1])
    x2  <- mini[2] + xx[2]*(maxi[2] - mini[2])

    term1 <- 100 * sqrt(abs(x2 - 0.01*x1^2))
    term2 <- 0.01 * abs(x1+10)

    y <- term1 + term2
    return(y)
  }

  otlcircuit <- function(xx){
    ##########################################################################
    #
    # OTL CIRCUIT FUNCTION
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # OUTPUT AND INPUT:
    #
    # Vm = midpoint voltage
    # xx = c(Rb1, Rb2, Rf, Rc1, Rc2, beta)
    #
    ##########################################################################
    mini = c(50,25,0.5,1.2,0.25,50)
    maxi = c(150,70,3,2.5,1.2,300)

    Rb1  <- mini[1] + xx[1]*(maxi[1] - mini[1])
    Rb2  <- mini[2] + xx[2]*(maxi[2] - mini[2])
    Rf   <- mini[3] + xx[3]*(maxi[3] - mini[3])
    Rc1  <- mini[4] + xx[4]*(maxi[4] - mini[4])
    Rc2  <- mini[5] + xx[5]*(maxi[5] - mini[5])
    beta <- mini[6] + xx[6]*(maxi[6] - mini[6])

    Vb1 <- 12*Rb2 / (Rb1+Rb2)
    term1a <- (Vb1+0.74) * beta * (Rc2+9)
    term1b <- beta*(Rc2+9) + Rf
    term1 <- term1a / term1b

    term2a <- 11.35 * Rf
    term2b <- beta*(Rc2+9) + Rf
    term2 <- term2a / term2b

    term3a <- 0.74 * Rf * beta * (Rc2+9)
    term3b <- (beta*(Rc2+9)+Rf) * Rc1
    term3 <- term3a / term3b

    Vm <- term1 + term2 + term3
    return(Vm)
  }

  robot <- function(xx){
    ##########################################################################
    #
    # ROBOT ARM FUNCTION
    #
    # Authors: Sonja Surjanovic, Simon Fraser University
    #          Derek Bingham, Simon Fraser University
    # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
    #
    # Copyright 2013. Derek Bingham, Simon Fraser University.
    #
    # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
    # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
    # derivative works, such modified software should be clearly marked.
    # Additionally, this program is free software; you can redistribute it
    # and/or modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; version 2.0 of the License.
    # Accordingly, this program is distributed in the hope that it will be
    # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # For function details and reference information, see:
    # http://www.sfu.ca/~ssurjano/
    #
    ##########################################################################
    #
    # OUTPUT AND INPUTS:
    #
    # y = distance from the end of the arm to the origin
    # xx = c(theta1, theta2, theta3, theta4, L1, L2, L3, L4)
    #
    #########################################################################

    mini = c(0,0,0,0,0,0,0,0)
    maxi = c(2*pi,2*pi,2*pi,2*pi,1,1,1,1)


    theta <- mini[1] + xx[1:4]*(maxi[1] - mini[1])
    L     <- mini[5] + xx[5:8]*(maxi[5] - mini[5])

    thetamat <- matrix(rep(theta,times=4), 4, 4, byrow=TRUE)
    thetamatlow <- thetamat
    thetamatlow[upper.tri(thetamatlow)] <- 0
    sumtheta <- rowSums(thetamatlow)

    u <- sum(L*cos(sumtheta))
    v <- sum(L*sin(sumtheta))

    y <- (u^2 + v^2)^(0.5)
    return(y)
  }
}
