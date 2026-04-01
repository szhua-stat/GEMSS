# download all packages
# packages = c('pacman', 'Metrics', 'MultiRNG', 'hetGP', 'mvtnorm', 'SPlit', 'class', 'mixAK',
#              'stringr', 'mined', 'twinning', 'supercompress', 'svMisc', 'tidyverse', 'gss',
#              'FNN', 'tmvtnorm', 'np', 'twingp', 'GpGp', 'fields', 'doParallel', 'foreac')
# install.packages(packages)
#
# download GEMSS package
remotes::install_github("szhua-stat/GEMSS")

# load packages
pacman::p_load(Metrics, MultiRNG, hetGP, mvtnorm, SPlit, class, mixAK, stringr, mined, twinning,
               supercompress, svMisc, tidyverse, gss, FNN, tmvtnorm, np, twingp,
               GpGp, fields, doParallel, foreach, GEMSS)

# download mined package (removed from CRAN since 2026)
remotes::install_version("mined", version = "1.0.3")

#### function ####
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

#### Analysis ####
packages <- c('Metrics', 'hetGP', 'SPlit', 'class', 'stringr',
              'mined', 'twinning', 'supercompress', 'svMisc', 'tidyverse',
              'gss', 'FNN', 'np', 'twingp', 'GpGp', 'fields', 'GEMSS')

setwd("sacros_dataset")
#read data
data_train <- read.csv('sarcos_inv.csv')[,-1]
data_test <- read.csv('sarcos_inv_test.csv')[,-1]

X <- data_train[,1:21] %>% as.matrix()
Y <- data_train[,22]
X_test <- data_test[,1:21] %>% as.matrix()
Y_test <- data_test[,22]
n_train <- nrow(X)
n_test <- nrow(X_test)

# ns_set <- c(250,500,750,1000)  # subdata size
ns_set <- c(100,200,300,400)  # subdata size
dup = 5# 50

index1 = 1
Cov_Fun <-  c('Matern3_2', 'Matern5_2')[index1] # GP kernel

time_start <- Sys.time()

# parallel computing with cl cores
cl = 5
registerDoParallel(cl)
RES <- foreach(i = 1:dup, .errorhandling = "stop", .packages=packages, .combine=rbind) %dopar%
  {
    Result <- matrix(NA,length(ns_set),ncol = 23)
    set.seed(i + 10000)

    # estimate parameters
    d <- ncol(X)
    ind_for_est <- twin(cbind(X,Y), floor(n_train / last(ns_set))) # random subsample of size-1000 for the parameter estimation
    est_GP <- mleHomGP(X[ind_for_est, ], Y[ind_for_est], covtype = Cov_Fun, maxit = 10000, noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps), 10000)))
    Pars <- est_GP[c(1,2,3,6)]
    print(round(unlist(Pars),4))

    test_rmse_GEMSS <- rmse(Y_test, predict(est_GP, X_test)$mean)
    RES_full <- test_rmse_GEMSS / rmse(Y_test, mean(Y))

    for (l in 1:length(ns_set)) {
      ns <- ns_set[l]

      #GEMSS
      RES_GEMSS <-tryCatch({
        t1 <- Sys.time()
        Result_GEMSS <- gemss_select(X, Y, ns, Cov_Fun, Pars, verbose  = F)
        Index_GEMSS <- Result_GEMSS$index
        t2 = Sys.time()
        Time <-  difftime(t2, t1, units='mins')

        # Making prediction on testing data using the same parameters
        sub_GP = gp_predict(X_test, X[Index_GEMSS, ], Y[Index_GEMSS], Cov_Fun, Result_GEMSS$parameters)
        test_rmse_GEMSS <- rmse(Y_test, sub_GP$mean)
        c(rmspe_GEMSS = test_rmse_GEMSS, nrmspe_GEMSS = test_rmse_GEMSS / rmse(Y_test, mean(Y)), Time_GEMSS = Time)

      },error = function(err) {
        print(paste("MY_ERROR:  ",err))
        c(rmspe_GEMSS = NA, nrmspe_GEMSS = NA, Time_GEMSS = NA)
      })

      #Chang (2023)
      RES_Chang <-tryCatch({
        t1 <- Sys.time()
        init_ind <- sample(1:nrow(X), floor(ns*0.3))
        Result_Chang <- EI_MSEPV_Sum_val(X, Y, ns, init_ind, Cov_Fun, print_result = F)
        Index_Chang <- Result_Chang$index
        t2 = Sys.time()
        Time <-  difftime(t2, t1, units='mins')
        test_rmse_Chang <- RMSE_GP(X[Index_Chang, ], Y[Index_Chang], X_test, Y_test, Cov_Fun)
        c(rmspe_Chang = test_rmse_Chang, nrmspe_Chang = test_rmse_Chang / rmse(Y_test, mean(Y)), Time_Chang = Time)

      },error = function(err) {
        print(paste("MY_ERROR:  ",err))
        c(rmspe_Chang = NA, nrmspe_Chang = NA, Time_Chang = NA)
      })

      # ASMEC
      RES_ASMEC <-tryCatch({
        t1 <- Sys.time()
        Index_ASMECr = MEDh_selectr(X, Y, ns)
        t2 = Sys.time()
        Time <-  difftime(t2, t1, units='mins')
        test_rmse_ASMECr <- RMSE_GP(X[Index_ASMECr, ], Y[Index_ASMECr], X_test, Y_test, Cov_Fun)
        c(rmspe_ASMECr = test_rmse_ASMECr, nrmspe_ASMECr = test_rmse_ASMECr / rmse(Y_test, mean(Y)), Time_ASMECr = Time)

      },error = function(err) {
        print(paste("MY_ERROR:  ",err))
        c(rmspe_ASMEC = NA, nrmspe_ASMEC = NA, Time_ASMEC = NA)
      })

      # Supercompress
      RES_supercom <-tryCatch({
        t1 <- Sys.time()
        supercom <- supercompress(ns, X, Y,lam=1/(1+d))
        t2 = Sys.time()

        Time <-  difftime(t2, t1, units='mins')
        test_rmse_supercom <- RMSE_GP(supercom$D, supercom$ybar, X_test, Y_test, Cov_Fun)
        c(rmspe_supercom = test_rmse_supercom, nrmspe_supercom = test_rmse_supercom / rmse(Y_test, mean(Y)), Time_supercom = Time)
      },error = function(err) {
        print(paste("MY_ERROR:  ",err))
        c(rmspe_supercom = NA, nrmspe_supercom = NA, Time_supercom = NA)
      })

      # Simple random sampling
      RES_SRS <-tryCatch({
        t1 <- Sys.time()
        Index_SRS <- sample(1:n_train, ns)
        t2 = Sys.time()
        Time <-  difftime(t2, t1, units='mins')
        test_rmse_SRS <- RMSE_GP(X[Index_SRS, ], Y[Index_SRS], X_test, Y_test, Cov_Fun)
        c(rmspe_SRS = test_rmse_SRS, nrmspe_SRS = test_rmse_SRS / rmse(Y_test, mean(Y)), Time_SRS = Time)

      },error = function(err) {
        print(paste("MY_ERROR:  ",err))
        c(rmspe_SRS = NA, nrmspe_SRS = NA, Time_SRS = NA)
      })

      Result[l, 1:17] <- c(ns, RES_full, RES_GEMSS, RES_Chang, RES_ASMEC, RES_supercom, RES_SRS)
    }
    # twinGP (not subsampling approach)
    RES_twinGP <-tryCatch({
      t1 <- Sys.time()
      Result_twinGP <- twingp(X, Y, X_test)
      t2 = Sys.time()
      Time <-  difftime(t2, t1, units='mins')

      test_rmse_twinGP <- rmse(Y_test, Result_twinGP$mu)
      c(rmspe_twinGP = test_rmse_twinGP, nrmspe_twinGP = test_rmse_twinGP / rmse(Y_test, mean(Y)), Time_twinGP = Time)

    },error = function(err) {
      print(paste("MY_ERROR:  ",err))
      c(rmspe_twinGP = NA, nrmspe_twinGP = NA, Time_twinGP = NA)
    })

    # GpGp (not subsampling approach)
    RES_GpGp <-tryCatch({
      t1 <- Sys.time()
      Result_GpGp <- fit_model(y = Y, locs = X, covfun_name = ifelse(Cov_Fun == "Matern3_2", 'matern15_isotropic', 'matern25_isotropic'), silent = T)
      prediction_GpGp <- predictions(fit = Result_GpGp, locs_pred = X_test,X_pred = rep(1, n_test))
      t2 = Sys.time()
      Time <-  difftime(t2, t1, units='mins')

      test_rmse_GpGp <- rmse(Y_test, prediction_GpGp)
      c(rmspe_GpGp = test_rmse_GpGp, nrmspe_GpGp = test_rmse_GpGp / rmse(Y_test, mean(Y)), Time_GpGp = Time)

    },error = function(err) {
      print(paste("MY_ERROR:  ",err))
      c(rmspe_GpGp52 = NA, nrmspe_GpGp= NA, Time_GpGp = NA)
    })

    Result[1, 18:23] <- c(RES_twinGP, RES_GpGp)
    if(i == 1){colnames(Result) <- names(c(ns=ns, nrmspe_full = RES_full, RES_GEMSS, RES_Chang, RES_ASMEC, RES_supercom, RES_SRS, RES_twinGP, RES_GpGp))}
    Result
}
# stopCluster(cl)
print(RES)

time_end <- Sys.time()
print(paste0('time comsuming: ',hms_span(time_start, time_end)))


#### boxplot ####
#nrmspe
RES <- as.data.frame(RES)
res1 <- RES[,c(1,(1:5)*3+1)] %>% group_by(ns) %>% gather(key = 'method', value = 'nrmspe',c(2:6)) %>% as.data.frame()
res2 <- data.frame(ns = 'Full(left)\nGPGP(middle)\nTwinGP(right)', RES[(1:dup)*4-3,c(2,19,22)] %>% gather(key = 'method', value = 'nrmspe',1:3))
Res <- rbind(res1, res2) %>% mutate(ns = factor(ns, levels = c(ns_set, 'Full(left)\nGPGP(middle)\nTwinGP(right)')))

ylim <- quantile(Res$nrmspe, 0.99, na.rm = T)
Res %>% ggplot(aes(x = ns, y = nrmspe, col = method)) +
  geom_boxplot(width=0.35, position = position_dodge(width = 0.5), outlier.shape = T) +
  stat_summary(mapping = aes(group = method, linetype = method),fun = "median", geom = "line",
               position = position_dodge(width = 0.5)) +
  ylim(c(NA,ylim)) +
  scale_color_manual(name = "Method ",
                     labels = c("ASMECr","Chang(2023)","Full", "GEMSS", 'GpGp', "SRS", "Supercom",'twingp'),
                     values = c('blue4',2:8))+
  scale_linetype_manual(name = "Method",
                        labels = c("ASMECr","Chang(2023)","Full", "GEMSS", 'GpGp', "SRS", "Supercom",'twingp'),
                        values = c(1,2,0,3,0,4,5,0)) +
  theme(strip.text = element_text(size = 14), legend.position = 'bottom',
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(color  = "Method ", linetype = "Method")

# log time (in minute)
res1 <- RES[,c(1,(1:4)*3+2)] %>% group_by(ns) %>% gather(key = 'method', value = 'time',c(2:5)) %>% as.data.frame()
res2 <- data.frame(ns = 'GPGP(left)\nTwinGP(right)', RES[(1:dup)*4-3,c(20,23)] %>%  gather(key = 'method', value = 'time',1:2))
Res_time <- rbind(res1, res2) %>% mutate(ns = factor(ns, levels = c(ns_set, 'GPGP(left)\nTwinGP(right)')))


ylim <- quantile(log10(Res_time$time), 0.995, na.rm = T)
Res_time %>% ggplot(aes(x = ns, y = log10(time), col = method)) +
  geom_boxplot(width=0.35, position = position_dodge(width = 0.5), outlier.shape = NA) +
  stat_summary(mapping = aes(group = method, linetype = method),fun = "median", geom = "line",
               position = position_dodge(width = 0.5)) +
  ylim(c(NA,ylim)) +
  scale_color_manual(name = "Method ",
                     labels = c("ASMECr","Chang(2023)","GEMSS",'GpGp', "Supercom",'twingp'),
                     values = c('blue4',2,4,5,7,8))+
  scale_linetype_manual(name = "Method",
                        labels = c("ASMECr","Chang(2023)","GEMSS",'GpGp', "Supercom",'twingp'),
                        values = c(1,2,3,0,5,0)) +
  theme(strip.text = element_text(size = 14), legend.position = 'bottom',
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(color  = "Method ", linetype = "Method", y = expression("log"[10]*"(Time in minutes)"))

