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

#### simulation function ####
source('all_functions.R')

# packages for doParallel
packages <- c('Metrics', 'MultiRNG', 'hetGP', 'mvtnorm', 'SPlit', 'class', 'mixAK',
              'stringr', 'mined', 'twinning', 'supercompress', 'svMisc', 'tidyverse',
              'gss', 'FNN', 'tmvtnorm', 'np', 'twingp', 'GpGp', 'fields', 'GEMSS')

#### Simulation ####
# settings
# fun_type: determine the simulation function
# X_type: control the X distribution, uniform, multivariate truncated t. and multivariate normal
# Cov_Fun: choose the GP Kernel, Matern3_2 or Matern5_2
# SNR: signal-to-noise ratio. Either 2 or 10 in our setting

index1 = 1
index2 = 1
index3 = 1
index4 = 1

fun_type = c("brat", "robot", "g_fun", "welch")[index1]
X_type = c("Uniform", 'Mvt', 'MVNormal')[index2]
Cov_Fun = c('Matern3_2', 'Matern5_2')[index3]
SNR = c(2,10)[index4] # Signal-to-Noise Ratio

# paper setting, which could be time-consuming
# n_train = 100000 # training data size
# n_test = 10000 # testing data size
# ns_set = c(250, 500, 750, 1000) # subdata size

# the smaller setting
n_train = 10000
n_test = 10000
ns_set = c(100, 200, 300, 400)


dup = 10 # 50 Replication

print(paste0("---",n_train, '-', n_test,"---"))
print(paste0('fun_type = ', fun_type,', X_type = ', X_type, ', SNR = ', SNR))

time_start <- Sys.time()

# parallel computing with cl cores
cl = 5 #detectCores()
registerDoParallel(cl)
RES <- foreach(i = 1:dup, .errorhandling = "stop", .packages=packages, .combine=rbind) %dopar%
  {
    Result <- matrix(NA, length(ns_set), ncol = 23)
    set.seed(i + 10000)
    # Draw X
    rho_range <- switch(X_type, Uniform = c(0,0), Mvt = c(0.5,0.5), MVNormal = c(0.3,0.3))
    dat <- generate_data(SNR, X_type, n_train, n_test, fun_type, rho_range = rho_range)
    X <- dat$X ; Y <- dat$Y
    X_test <- dat$X_test ; Y_test <- dat$Y_test

    # estimate parameters via a twinning sample
    d <- ncol(X)
    ind_for_est <- twin(cbind(X,Y), floor(n_train / last(ns_set))) # random subsample of size-1000 for the parameter estimation
    est_GP <- mleHomGP(X[ind_for_est, ], Y[ind_for_est], covtype = Cov_Fun, maxit = 10000, noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps), 10000)))
    Pars <- est_GP[c(1,2,3,6)]

    # calculate the nrmspe for full method
    test_rmse_GEMSS <- rmse(Y_test, predict(est_GP, X_test)$mean)
    RES_full <- test_rmse_GEMSS / rmse(Y_test, mean(Y))

    for (l in 1:length(ns_set)) {
      ns <- ns_set[l]
      print(paste0(fun_type,': ns = ',ns, ', Replicate = ', i))

      # GEMSS
      RES_GEMSS <-tryCatch({
        t1 <- Sys.time()
        # select subdata
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

      # Chang (2023)
      p_init <- 0.3 # initial subdata size for Chang (2023)
      RES_Chang <-tryCatch({
        t1 <- Sys.time()
        init_ind <- sample(1:nrow(X), floor(ns*p_init))
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
# nrmspe
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
