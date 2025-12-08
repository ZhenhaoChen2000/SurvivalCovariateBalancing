##To solve the optimization problem in sbw, the default solver requires 'osqp' packages, which may require to install cmake


rm(list = ls())
#load function========
if (!require("WeightIt")){
  install.packages("WeightIt")
}
if (!require("mlogit")){
  install.packages("mlogit")
}
library(survminer)
library(cobalt)
library(WeightIt)

library(ATE.ncb) 
library(sbw)
if (!require("devtools")){
  install.packages("devtools")
}
devtools::install_github("yimindai0521/MBalance")
library(MBalance)

get_diff<-function(time, status, group, weight=NULL, tau=NULL, alpha=.05, 
                   xaxismin=0, xaxismax=max(time)){
  if(sum(time<0)>0){print("Error: times must be positive.")
  }else{
    if(sum(weight<=0)>0){print("Error: weights must be greater than 0.")
    }else{
      if(sum(status!=0 & status!=1)>0){print("Error: status must be a vector of 0s and/or 1s.")
      }else{
        
        if(is.null(weight)){weight <- rep(1, length(time))}	
        data <- data.frame(time, status, group, weight)
        data <- data[!is.na(data$group) & !is.na(data$time),]
        data <- data[order(group),] 
        
        #--- If tau not specified, use minimum tau from all groups ---
        j=length(unique(data$group))
        
        if(is.null(tau)){
          taui = rep(999, j)
          for (i in (1:j)){
            groupval <- (levels(data$group)[i])
            dat_group <- data[which(data$group==(groupval)),]
            taui[i] <- max(dat_group$time[dat_group$status==1])
          }
          tau <- min(taui)
        }
        
        #--- Calculate AKM RMST in each group ---
        rmst <- rep(999, length(1:j))
        groupval <- rep(999, length(1:j))
        rmst_var <- rep(999, length(1:j))
        rmst_se <- rep(999, length(1:j))
        
        for (i in 1:j){
          groupval[i] <- (levels(data$group)[i])
          dat_group <- data[which(data$group==(groupval[i])),]
          
          #--- AKM ---
          # Based on 'adjusted.KM' function from {IPWsurvival} package
          # Author: F. Le Borgne and Y. Foucher
          tj <- c(0,sort(unique(dat_group$time[dat_group$status==1])))
          dj <- sapply(tj, function(x){sum(dat_group$weight[dat_group$time==x & dat_group$status==1])})
          yj <- sapply(tj, function(x){sum(dat_group$weight[dat_group$time>=x])})
          st <- cumprod(1-(dj/yj))
          m <- sapply(tj, function(x){sum((dat_group$weight[dat_group$time>=x])^2)})
          mj <- ((yj^2)/m)
          #ft <- data.frame(time=tj, n_risk=yj, n_event=dj, survival=st, variable=i, m=mj)
          ft <- data.frame(tj, yj, dj, st, i, mj)
          
          #--- RMST ---
          # Based on 'rmst1 function' from {survRM2} package
          # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
          rtime <- ft$tj<=tau
          tj_r <- sort(c(ft$tj[rtime],tau))
          st_r <- ft$st[rtime]
          yj_r <- ft$yj[rtime]
          dj_r <- ft$dj[rtime]
          time_diff <- diff(c(0, tj_r))
          areas <- time_diff * c(1, st_r)
          rmst[i] <- sum(areas)
          
          #--- Variance ---
          mj_r <- ft$mj[rtime]
          var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
          #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
          var_r <- c(var_r,0)
          rmst_var[i] <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
          rmst_se[i] <- sqrt(rmst_var[i])
          
        }
      }
    }
  }
  
  #--- Compare RMST between groups and compile output---
  results <- data.frame(groupval,rmst,rmst_var,rmst_se,tau)
  pwc <- ((j^2)-j)/2   #number of pairwise comparisons
  
  label_diff <- rep(999,(pwc))
  rmst_diff <- rep(999,(pwc))
  rmst_diff_se <- rep(999,(pwc))
  rmst_diff_low <- rep(999,(pwc))
  rmst_diff_upp <- rep(999,(pwc))
  rmst_diff_pval <- rep(999,(pwc))
  
  output_diff <- data.frame(label_diff,rmst_diff,rmst_diff_se,rmst_diff_low,rmst_diff_upp,rmst_diff_pval)
  #output_ratio <- data.frame(label_ratio,rmst_logratio,rmst_logratio_se,rmst_ratio,rmst_ratio_low,rmst_ratio_upp,rmst_logratio_pval)
  l <- 1
  
  for (i in 1:(j-1)){
    for (j in (i+1):j){
      # Based on 'rmst2 function' from {survRM2} package
      # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
      
      #--- RMST Difference ---
      output_diff[l,]$label_diff <- paste('Groups',results[j,]$groupval,'vs.',results[i,]$groupval,' ')
      output_diff[l,]$rmst_diff <- (results[j,]$rmst - results[i,]$rmst)
      output_diff[l,]$rmst_diff_se <- sqrt(results[j,]$rmst_var + results[i,]$rmst_var)
      output_diff[l,]$rmst_diff_low <- output_diff[l,]$rmst_diff - qnorm(1-alpha/2)*output_diff[l,]$rmst_diff_se
      output_diff[l,]$rmst_diff_upp <- output_diff[l,]$rmst_diff + qnorm(1-alpha/2)*output_diff[l,]$rmst_diff_se
      output_diff[l,]$rmst_diff_pval <- 2*(1-pnorm(abs(output_diff[l,]$rmst_diff)/output_diff[l,]$rmst_diff_se))
      
      l <- l+1 #move to next row
    }
  }
  return(output_diff)
}

################MDABW method
sbwbal <- function(X, Tr){
  data.matrix <- data.frame(X , factor(Tr))
  dimension   <- dim(X)[2]
  bal = list()
  bal$bal_gri = c(1e-05,1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  
  
  character   <- names(data.matrix)
  bal$bal_cov <- character[1:dimension]
  sbw.result <- FALSE
  ######dimension of the output of sbw function is 15, if sbw does not output correct result, search grid minus one
  while(sum(dim(as.matrix(sbw.result))) != 15){
    sbw.result  <- tryCatch(sbw.result <- sbw(dat = data.matrix, ind = character[1 + dimension], bal = bal, out = character[2 + dimension], par = list(par_est = "ate")), error = function(e) { skip_to_next <<- FALSE})
    bal$bal_gri <- bal$bal_gri[-1]
  }
  sbw.weight <- sbw.result$dat_weights$sbw_weights
  return(list(weight = sbw.weight))
}

sbwbal.high.dim <- function(X, Tr){
  data.matrix <- data.frame(X , factor(Tr))
  dimension   <- dim(X)[2]
  bal = list()
  bal$bal_gri = c(1e-05,1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  
  # Set tolerances
  bal$bal_tol = 0.002
  bal$bal_std = "target"
  
  character   <- names(data.matrix)
  bal$bal_cov <- character[1:dimension]
  sbw.result <- FALSE
  
  while(sum(dim(as.matrix(sbw.result))) != 15){
    sbw.result  <- tryCatch(sbw.result <- sbw(dat = data.matrix, ind = character[1 + dimension], bal = bal, out = character[2 + dimension], par = list(par_est = "ate")), error = function(e) { skip_to_next <<- FALSE})
    bal$bal_gri <- bal$bal_gri[-1]
  }
  sbw.weight <- sbw.result$dat_weights$sbw_weights
  return(list(weight = sbw.weight))
}

#Function for Mahalanobis Balancing
MB <- function(x = x , treat = treat , group1 , group2 = NA, outcome = outcome, method = "MB", delta.space = c(100,10,1,1e-2,1e-4), iterations = 1000, dimension = dimension){
  #data 
  data.matrix         <- cbind(x , treat , outcome)
  sample.size         <- dim(x)[1]
  group1.number       <- sum(treat == group1)
  
  stopifnot(method %in% c("MB", "MB2", "kernelMB"))
  
  if (is.null(group1)){stop("group1 must be a integer")}
  
  if (method == "MB") {
    solution <- Cholesky(x = x, treat = treat, group1 = group1, group2 = NA, method = "MB", dimension = dimension)
  }
  if (method == "MB2") {
    solution <- Cholesky(x = x, treat = treat, group1 = group1, group2 = NA, method = "MB2", dimension = dimension)
  }
  if(method == "kernelMB"){
    dimension   <- dim(x)[1]
    data        <- matrix(NA,dim(x)[1],dim(x)[1])
    data.temp   <- matrix(NA,dim(x)[1],dim(x)[1])
    for(i in 1:sample.size){
      for(j in 1:sample.size){
        data.temp[i,j] <- sum((data.matrix[i,1:dim(x)[2]]-data.matrix[j,1:dim(x)[2]])^2)
      }
    }
    median.scale <- median(data.temp)
    for(i in 1:sample.size){
      for(j in 1:sample.size){
        data[i,j] <- exp( -data.temp[i,j] / median.scale )
      }
    }
    solution <- Cholesky(x = x, treat = treat, group1 = group1, group2 = NA, method = "MB2", dimension = dimension)
  }
  x <- solution
  
  #calculate the weights for group1
  w <- rep(0,group1.number)
  u <- rep(0,dimension)
  GMIM <- rep(NA,sum(delta.space > 0))
  for(i in 1:sum(delta.space > 0)){
    Opti.func <- function(u){  
      u1 <- t(x) %*% u
      sum(exp(u1-1)) + sqrt(delta.space[i]) * norm(u,type = c("2"))
    }
    solution  <- optim(par = u,fn = Opti.func, method = "BFGS", control = list(abstol=0, maxit = iterations))
    u.sol     <- t(x) %*% solution$par
    w         <- exp(u.sol-1)
    w         <- w/sum(w)
    GMIM[i]   <- t(w) %*% t(x) %*% x %*% w
  }
  
  #calculate the AT
  delta    <- delta.space[rank(GMIM) == 1]
  Opti.func <- function(u){  
    u1 <- t(x) %*% u
    sum(exp(u1-1)) + sqrt(delta) * norm(u,type = c("2"))
  }
  solution  <- optim(par = u,fn = Opti.func, method = "BFGS", control = list(abstol= 1e-8, maxit = iterations))
  u.sol     <- t(x) %*% solution$par
  w         <- exp(u.sol-1)
  w         <- w/sum(w)
  GMIM      <- t(w) %*% t(x) %*% x %*% w
  
  dimension <- dim(x)[1]
  AT        <- t(w) %*% outcome[treat == group1]
  result    <- list(AT = AT, weight = w, GMIM = GMIM, parameter = solution$par, delta = delta)
  return(result)
}

#Function for Mahalanobis Balancing
Cholesky <- function(x = x , treat = treat, group1 = 1, group2 = NA, method = "MB", dimension = dimension){
  data.matrix <- cbind(x , treat)
  group       <- data.matrix[treat == group1,1:dimension]
  group.num   <- sum(treat == group1)
  group.col   <- ncol(group)
  group.row   <- nrow(group)
  group.cov   <- matrix(0, group.col, group.col)
  sample.size <- sum(treat >= 0)
  
  if(method == "MB"){
    for(i in 1:group.col){
      meanx <- sum(x[,i]) / sample.size
      group.cov[i,i] <- sum((x[,i]-meanx)*(x[,i]-meanx)) / (sample.size-1)
    }
  }
  
  if(method == "MB2"){
    for(i in 1:group.col){
      for(j in 1:group.col){
        meanx <- sum(x[,i]) / sample.size
        meany <- sum(x[,j]) /sample.size
        group.cov[i,j] <- sum((x[,i]-meanx)*(x[,j]-meany)) / (sample.size-1)
      }
    }
  }
  
  mean.pop <- matrix(NA,1,dimension)
  mean.pop <- apply(data.matrix[,1:dimension], 2, mean)
  group.se <- data.matrix[(treat == group1),1:dimension]
  for(i in 1:sum(treat == group1)){group.se[i,] <- group.se[i,] - mean.pop}
  
  K  <- solve(group.cov)
  Q  <- chol(K)  
  x  <- Q %*% t(group.se)
  
  return(x)
}



####################



weight.esti <- function(X, Tr, delta.space = c(1,1e-2,1e-4)){
  data.matrix <- data.frame(X , factor(Tr))
  sample.size <- dim(X)[1]
  dimension   <- dim(X)[2]
  
  character   <- names(data.matrix)
  for(j in 1:(dimension+1)){character[j] <- paste(character[j])}
  myformula   <- as.formula(paste(character[1 + dimension],paste(" ~ ", paste(character[1:dimension], collapse= "+"))))
  
  ps.weight   <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "ps")$weights
  ebal.weight <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "ebal")$weights
  cbps.weight <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "cbps", over = FALSE)$weights
  energy.weight <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "energy")$weights
  MB.weight   <- rep(NA,sample.size)
  MB1.result  <- MB(x = X, treat = Tr, group1 = 1, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB", dimension = dimension)
  MB2.result  <- MB(x = X, treat = Tr, group1 = 2, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB", dimension = dimension)
  MB3.result  <- MB(x = X, treat = Tr, group1 = 3, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB", dimension = dimension)
  MB.weight[Tr == 1] <- sum(Tr == 1) * MB1.result$weight
  MB.weight[Tr == 2] <- sum(Tr == 2) * MB2.result$weight
  MB.weight[Tr == 3] <- sum(Tr == 3) * MB3.result$weight
  
  if(sum(ebal.weight[Tr == 1]) <10^(-6)){ebal.weight[Tr == 1] <- rep(1,sum(Tr == 1))}
  if(sum(ebal.weight[Tr == 2]) <10^(-6)){ebal.weight[Tr == 2] <- rep(1,sum(Tr == 2))}
  if(sum(ebal.weight[Tr == 3]) <10^(-6)){ebal.weight[Tr == 3] <- rep(1,sum(Tr == 3))}
  
  
  #kernel based balancing
  # Sobolev kernel
  Xstd <- transform.sob(X)$Xstd # standardize X to [0,1]^p
  Kern <- getGram(Xstd) # get Gram matrix using Sobolev kernel
  
  # design a grid for the tuning parameter
  nlam <- 50
  lams <- exp(seq(log(1e-8), log(1), len=nlam))
  # Tr_kernel_1: Tr==1:1;else:0
  Tr_kernel_1<-rep(0,sample.size)
  Tr_kernel_1[which(Tr==1)]<-1
  # Tr_kernel_2: Tr==2:1;else:0
  Tr_kernel_2<-rep(0,sample.size)
  Tr_kernel_2[which(Tr==2)]<-1
  # Tr_kernel_3: Tr==3:1;else:0
  Tr_kernel_3<-rep(0,sample.size)
  Tr_kernel_3[which(Tr==3)]<-1
  # compute weights for Tr=1
  fit1 <- ATE.ncb.SN(Tr_kernel_1, Kern, lam1s=lams)
  #if (sum(fit1$warns)) cat("lambda bound warning!\n")
  # compute weights for Tr=2
  fit2 <- ATE.ncb.SN(Tr_kernel_2, Kern, lam1s=lams)
  #if (sum(fit2$warns)) cat("lambda bound warning!\n")
  # compute weights for Tr=3
  fit3 <- ATE.ncb.SN(Tr_kernel_3, Kern, lam1s=lams)
  #if (sum(fit3$warns)) cat("lambda bound warning!\n")
  #kernel base balance weight
  kernel.weight<-fit1$w+fit2$w+fit3$w
  
  #MDABW method
  # Tr_MDABW_1: Tr==1:1;else:0
  Tr_MDABW_1<-rep(0,sample.size)
  Tr_MDABW_1[which(Tr==1)]<-1
  # Tr_MDABW_2: Tr==2:1;else:0
  Tr_MDABW_2<-rep(0,sample.size)
  Tr_MDABW_2[which(Tr==2)]<-1
  # Tr_MDABW_3: Tr==3:1;else:0
  Tr_MDABW_3<-rep(0,sample.size)
  Tr_MDABW_3[which(Tr==3)]<-1
  MDABW_weight1<-sbwbal(X,Tr_MDABW_1)$weight
  MDABW_weight2<-sbwbal(X,Tr_MDABW_2)$weight
  MDABW_weight3<-sbwbal(X,Tr_MDABW_3)$weight
  MDABW.weight<-rep(0,sample.size)
  MDABW.weight[Tr==1]<-MDABW_weight1[Tr==1]
  MDABW.weight[Tr==2]<-MDABW_weight2[Tr==2]
  MDABW.weight[Tr==3]<-MDABW_weight3[Tr==3]
  
  unad.imbalance   <- Imbalance(X = X, Tr = Tr, w = rep(1,sample.size))
  ps.imbalance     <- Imbalance(X = X, Tr = Tr, w = ps.weight)
  ebal.imbalance   <- Imbalance(X = X, Tr = Tr, w = ebal.weight)
  cbps.imbalance   <- Imbalance(X = X, Tr = Tr, w = cbps.weight)
  energy.imbalance <- Imbalance(X = X, Tr = Tr, w = energy.weight)
  MB.imbalance     <- Imbalance(X = X, Tr = Tr, w = MB.weight)
  
  kernel.imbalance <- Imbalance(X = X, Tr = Tr, w = kernel.weight)
  MDABW.imbalance  <- Imbalance(X = X, Tr = Tr, w = MDABW.weight)
  
  
  return(list(ps = ps.weight, ebal = ebal.weight, cbps = cbps.weight, 
              energy = energy.weight, MB = MB.weight, kernel=kernel.weight, MDABW=MDABW.weight,
              MB.parameter = list(group1 = MB1.result$parameter, 
                                  group2 = MB2.result$parameter, 
                                  group3 = MB3.result$parameter),
              Imbalance = list(unad = unad.imbalance, ps = ps.imbalance, 
                               ebal = ebal.imbalance, cbps = cbps.imbalance, 
                               energy = energy.imbalance, MB = MB.imbalance, kernel=kernel.imbalance, MDABW=MDABW.imbalance)
             
  )
  )
}


############################## weighting without kernel balancing method
###use this function when considering all variables
weight.esti.nokernel <- function(X, Tr, delta.space = c(1,1e-2,1e-4)){
  data.matrix <- data.frame(X , factor(Tr))
  sample.size <- dim(X)[1]
  dimension   <- dim(X)[2]
  
  character   <- names(data.matrix)
  for(j in 1:(dimension+1)){character[j] <- paste(character[j])}
  myformula   <- as.formula(paste(character[1 + dimension],paste(" ~ ", paste(character[1:dimension], collapse= "+"))))
  
  ps.weight   <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "ps")$weights
  ebal.weight <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "ebal")$weights
  cbps.weight <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "cbps", over = FALSE)$weights
  energy.weight <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "energy")$weights
  MB.weight   <- rep(NA,sample.size)
  MB1.result  <- MB(x = X, treat = Tr, group1 = 1, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB", dimension = dimension)
  MB2.result  <- MB(x = X, treat = Tr, group1 = 2, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB", dimension = dimension)
  MB3.result  <- MB(x = X, treat = Tr, group1 = 3, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB", dimension = dimension)
  MB.weight[Tr == 1] <- sum(Tr == 1) * MB1.result$weight
  MB.weight[Tr == 2] <- sum(Tr == 2) * MB2.result$weight
  MB.weight[Tr == 3] <- sum(Tr == 3) * MB3.result$weight
  
  if(sum(ebal.weight[Tr == 1]) <10^(-6)){ebal.weight[Tr == 1] <- rep(1,sum(Tr == 1))}
  if(sum(ebal.weight[Tr == 2]) <10^(-6)){ebal.weight[Tr == 2] <- rep(1,sum(Tr == 2))}
  if(sum(ebal.weight[Tr == 3]) <10^(-6)){ebal.weight[Tr == 3] <- rep(1,sum(Tr == 3))}
  
  
  
  #MDABW method
  # Tr_MDABW_1: Tr==1:1;else:0
  Tr_MDABW_1<-rep(0,sample.size)
  Tr_MDABW_1[which(Tr==1)]<-1
  # Tr_MDABW_2: Tr==2:1;else:0
  Tr_MDABW_2<-rep(0,sample.size)
  Tr_MDABW_2[which(Tr==2)]<-1
  # Tr_MDABW_3: Tr==3:1;else:0
  Tr_MDABW_3<-rep(0,sample.size)
  Tr_MDABW_3[which(Tr==3)]<-1
  MDABW_weight1<-sbwbal.high.dim(X,Tr_MDABW_1)$weight
  MDABW_weight2<-sbwbal.high.dim(X,Tr_MDABW_2)$weight
  MDABW_weight3<-sbwbal.high.dim(X,Tr_MDABW_3)$weight
  MDABW.weight<-rep(0,sample.size)
  MDABW.weight[Tr==1]<-MDABW_weight1[Tr==1]
  MDABW.weight[Tr==2]<-MDABW_weight2[Tr==2]
  MDABW.weight[Tr==3]<-MDABW_weight3[Tr==3]
  #gasmd: all standardized
  unad.imbalance   <- Imbalance(X = X, Tr = Tr, w = rep(1,sample.size))
  ps.imbalance     <- Imbalance(X = X, Tr = Tr, w = ps.weight)
  ebal.imbalance   <- Imbalance(X = X, Tr = Tr, w = ebal.weight)
  cbps.imbalance   <- Imbalance(X = X, Tr = Tr, w = cbps.weight)
  energy.imbalance <- Imbalance(X = X, Tr = Tr, w = energy.weight)
  MB.imbalance     <- Imbalance(X = X, Tr = Tr, w = MB.weight)
  MDABW.imbalance  <- Imbalance(X = X, Tr = Tr, w = MDABW.weight)
  
  return(list(ps = ps.weight, ebal = ebal.weight, cbps = cbps.weight, 
              energy = energy.weight, MB = MB.weight, MDABW=MDABW.weight,
              MB.parameter = list(group1 = MB1.result$parameter, 
                                  group2 = MB2.result$parameter, 
                                  group3 = MB3.result$parameter),
              Imbalance = list(unad = unad.imbalance, ps = ps.imbalance, 
                               ebal = ebal.imbalance, cbps = cbps.imbalance, 
                               energy = energy.imbalance, MB = MB.imbalance, MDABW=MDABW.imbalance)
  )
  )
}

#all variables are standardized when calculating GASMD
Imbalance <- function(X,Tr,w){
  X <- scale(X)
  abs(t(X[Tr == 1,]) %*% w[Tr == 1] / sum(w[Tr == 1])) + 
    abs(t(X[Tr == 2,]) %*% w[Tr == 2] / sum(w[Tr == 2])) + 
    abs(t(X[Tr == 3,]) %*% w[Tr == 3] / sum(w[Tr == 3]))
}

akm_rmst <- function(time, status, group, weight=NULL, tau=NULL, alpha=.05, 
                     xaxismin=0, xaxismax=max(time)){
  if(sum(time<0)>0){print("Error: times must be positive.")
  }else{
    if(sum(weight<=0)>0){print("Error: weights must be greater than 0.")
    }else{
      if(sum(status!=0 & status!=1)>0){print("Error: status must be a vector of 0s and/or 1s.")
      }else{
        
        if(is.null(weight)){weight <- rep(1, length(time))}	
        data <- data.frame(time, status, group, weight)
        data <- data[!is.na(data$group) & !is.na(data$time),]
        data <- data[order(group),] 
        
        #--- If tau not specified, use minimum tau from all groups ---
        j=length(unique(data$group))
        
        if(is.null(tau)){
          taui = rep(999, j)
          for (i in (1:j)){
            groupval <- (levels(data$group)[i])
            dat_group <- data[which(data$group==(groupval)),]
            taui[i] <- max(dat_group$time[dat_group$status==1])
          }
          tau <- min(taui)
        }
        
        #--- Calculate AKM RMST in each group ---
        rmst <- rep(999, length(1:j))
        groupval <- rep(999, length(1:j))
        rmst_var <- rep(999, length(1:j))
        rmst_se <- rep(999, length(1:j))
        plot(NULL, xlim=c(xaxismin, xaxismax), ylim=c(0,1), xlab='Time (days)',ylab='Survival Probability',cex.lab=1.3,cex.axis=0.9)
        #title(main='Adjusted Kaplan-Meier')
        
        for (i in 1:j){
          groupval[i] <- (levels(data$group)[i])
          dat_group <- data[which(data$group==(groupval[i])),]
          
          #--- AKM ---
          # Based on 'adjusted.KM' function from {IPWsurvival} package
          # Author: F. Le Borgne and Y. Foucher
          tj <- c(0,sort(unique(dat_group$time[dat_group$status==1])))
          dj <- sapply(tj, function(x){sum(dat_group$weight[dat_group$time==x & dat_group$status==1])})
          yj <- sapply(tj, function(x){sum(dat_group$weight[dat_group$time>=x])})
          st <- cumprod(1-(dj/yj))
          m <- sapply(tj, function(x){sum((dat_group$weight[dat_group$time>=x])^2)})
          mj <- ((yj^2)/m)
          #ft <- data.frame(time=tj, n_risk=yj, n_event=dj, survival=st, variable=i, m=mj)
          ft <- data.frame(tj, yj, dj, st, i, mj)
          
          #--- RMST ---
          # Based on 'rmst1 function' from {survRM2} package
          # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
          rtime <- ft$tj<=tau
          tj_r <- sort(c(ft$tj[rtime],tau))
          st_r <- ft$st[rtime]
          yj_r <- ft$yj[rtime]
          dj_r <- ft$dj[rtime]
          time_diff <- diff(c(0, tj_r))
          areas <- time_diff * c(1, st_r)
          rmst[i] <- sum(areas)
          
          #--- Variance ---
          mj_r <- ft$mj[rtime]
          var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
          #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
          var_r <- c(var_r,0)
          rmst_var[i] <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
          rmst_se[i] <- sqrt(rmst_var[i])
          
          #--- Plot AKM ---
          lines(ft$tj, ft$st,type="s", col=(i+2), lwd=2)
        }
      }
    }
  }
  
  #--- Add legend and tau to plot ---
  #abline(v=tau, col=1, lty=3, lwd=2)
  #legend( "bottomleft",paste("Treatment", groupval), lty=rep(1, j), lwd=rep(2, j), col=3:(j+2), 
  #        cex=1, bty ="n", inset = c(0, 0),seg.len = 0.2,x.intersp = 0.5)
  
  #--- Compare RMST between groups and compile output---
  results <- data.frame(groupval,rmst,rmst_var,rmst_se,tau)
  pwc <- ((j^2)-j)/2   #number of pairwise comparisons
  
  label_diff <- rep(999,(pwc))
  rmst_diff <- rep(999,(pwc))
  rmst_diff_se <- rep(999,(pwc))
  rmst_diff_low <- rep(999,(pwc))
  rmst_diff_upp <- rep(999,(pwc))
  rmst_diff_pval <- rep(999,(pwc))
  
  label_ratio <- rep(999,(pwc))
  rmst_logratio <- rep(999,(pwc))
  rmst_logratio_se <- rep(999,(pwc))
  rmst_ratio <- rep(999,(pwc))
  rmst_ratio_low <- rep(999,(pwc))
  rmst_ratio_upp <- rep(999,(pwc))
  rmst_logratio_pval <- rep(999,(pwc))
  
  output_diff <- data.frame(label_diff,rmst_diff,rmst_diff_se,rmst_diff_low,rmst_diff_upp,rmst_diff_pval)
  output_ratio <- data.frame(label_ratio,rmst_logratio,rmst_logratio_se,rmst_ratio,rmst_ratio_low,rmst_ratio_upp,rmst_logratio_pval)
  l <- 1
  
  for (i in 1:(j-1)){
    for (j in (i+1):j){
      # Based on 'rmst2 function' from {survRM2} package
      # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
      
      #--- RMST Difference ---
      output_diff[l,]$label_diff <- paste('Groups',results[j,]$groupval,'vs.',results[i,]$groupval,' ')
      output_diff[l,]$rmst_diff <- (results[j,]$rmst - results[i,]$rmst)
      output_diff[l,]$rmst_diff_se <- sqrt(results[j,]$rmst_var + results[i,]$rmst_var)
      output_diff[l,]$rmst_diff_low <- output_diff[l,]$rmst_diff - qnorm(1-alpha/2)*output_diff[l,]$rmst_diff_se
      output_diff[l,]$rmst_diff_upp <- output_diff[l,]$rmst_diff + qnorm(1-alpha/2)*output_diff[l,]$rmst_diff_se
      output_diff[l,]$rmst_diff_pval <- 2*(1-pnorm(abs(output_diff[l,]$rmst_diff)/output_diff[l,]$rmst_diff_se))
      
      #--- RMST Ratio ---
      output_ratio[l,]$label_ratio <- paste('Groups',results[j,]$groupval,'vs.',results[i,]$groupval,' ')
      output_ratio[l,]$rmst_logratio <- (log(results[j,]$rmst) - log(results[i,]$rmst))
      output_ratio[l,]$rmst_logratio_se <- sqrt(results[j,]$rmst_var/(results[j,]$rmst^2) + results[i,]$rmst_var/(results[i,]$rmst^2))
      output_ratio[l,]$rmst_ratio <- exp(output_ratio[l,]$rmst_logratio)
      output_ratio[l,]$rmst_ratio_low <- exp(output_ratio[l,]$rmst_logratio - qnorm(1-alpha/2)*output_ratio[l,]$rmst_logratio_se)
      output_ratio[l,]$rmst_ratio_upp <- exp(output_ratio[l,]$rmst_logratio + qnorm(1-alpha/2)*output_ratio[l,]$rmst_logratio_se)
      output_ratio[l,]$rmst_logratio_pval <- 2*(1-pnorm(abs(output_ratio[l,]$rmst_logratio)/output_ratio[l,]$rmst_logratio_se))
      
      l <- l+1 #move to next row
    }
  }
  
  cat("\n\n\n")
  cat(paste('RMST calculated up to tau =',round(results$tau[1],3)))
  cat("\n\n\n")
  
  cat ("Restricted Mean Survival Time (RMST) per Group \n\n")
  colnames(results) <- c("Group", "RMST", "Var", "SE", "Tau")
  rownames(results) <- c(paste("Group", results$Group,' '))
  print(round(results[c(2,4)],3))
  cat("\n\n")
  
  cat ("Restricted Mean Survival Time (RMST) Differences \n\n")
  colnames(output_diff) <- c("Groups", "Est.", "SE", "CIL", "CIU", "p")
  rownames(output_diff) <- c(output_diff$Groups)
  print(round(output_diff[c(2,3,4,5,6)],3))
  cat("\n\n")
  
  cat ("Restricted Mean Survival Time (RMST) Ratios \n\n")  
  colnames(output_ratio) <- c("Groups", "Log Est.", "SE", "Est.", "CIL", "CIU", "p")
  rownames(output_ratio) <- c(output_ratio$Groups)
  print(round(output_ratio[c(2,3,4,5,6,7)],3))
}
source('get_diff_se_ci.R')
library(simsurv)
library(survival)
library(survminer)


#data process====
library(nnet)
library(survival)
##################################load the HCC data
hcc<-read.csv('hcc.csv')
head(hcc)
group<-hcc$Group
group_cate_dum<-class.ind(hcc$Group)
group_cate<-data.frame('group1'=group_cate_dum[,1],'group2'=group_cate_dum[,2],'group3'=group_cate_dum[,3])

AFP_cate_dum<-class.ind(hcc$AFP3)
AFP_cate<-data.frame('AFP1'=AFP_cate_dum[,1],'AFP2'=AFP_cate_dum[,2],'AFP3'=AFP_cate_dum[,3]);

MELD_cate_dum<-class.ind(hcc$MELD.grading)
MELD_cate<-data.frame('MELD1'=MELD_cate_dum[,1],'MELD2'=MELD_cate_dum[,2],'MELD3'=MELD_cate_dum[,3]);

childpugh_cate_dum<-class.ind(hcc$Child.Pugh.grading)
childpugh_cate<-data.frame('childpugh1'=childpugh_cate_dum[,1],'childpugh2'=childpugh_cate_dum[,2]);

ALBI_cate_dum<-class.ind(hcc$ALBI.grading)
ALBI_cate<-data.frame('ALBI1'=ALBI_cate_dum[,1],'ALBI2'=ALBI_cate_dum[,2],'ALBI3'=ALBI_cate_dum[,3]);

APRI_cate_dum<-class.ind(hcc$APRI.grading)
APRI_cate<-data.frame('APRI1'=APRI_cate_dum[,1],'APRI2'=APRI_cate_dum[,2]);

ori_date<-hcc$OR_Date
recur_date<-hcc$Recur_Date
death_date<-hcc$Death_Date
last_date<-hcc$Last_Date

sum(!is.na(ori_date))
sum(!is.na(recur_date))
sum(!is.na(death_date))
sum(!is.na(last_date))

time_date<-rep(0,nrow(hcc))
#For individuals whose recur_date is not na, assign recur_date as its time date
time_date[!is.na(recur_date)]<-recur_date[!is.na(recur_date)]
#For individuals whose recur_date is na, assign last_date as its time date
time_date[is.na(recur_date)]<-last_date[is.na(recur_date)]
sum(time_date!=0)

time<-difftime(time_date,ori_date,units = 'days')
time

#if recur: status=1; else status=0
status<-rep(0,nrow(hcc))
status[!is.na(recur_date)]=1;
status[is.na(recur_date)]=0
sum(status)

dat<-data.frame(time=time,status=status,trt=hcc$operation,group_cate[,c(2,3)],gender=hcc$Gender,age=hcc$Age,
                AFP=hcc$AFP,MELD.SCORE=hcc$MELD.score,Childpugh.score=hcc$Child.Pugh.score,
                ALBI.score=hcc$ALBI.score,APRI.score=hcc$APRI.score,
                AFP_cate[,c(2,3)],MELD_cate[,c(2,3)],childpugh_cate[,2],
                ALBI_cate[,c(2,3)],APRI_cate[,2])

#tau=1606==================


##important variables=====
#omit (1)AFP (2)MELD.SCORE (3)Childpugh.score (4)ALBI.score (5)APRI.score
dat_discrete=dat[,-c(8:12)]
#dat_discrete$AFP=log(dat$AFP)
#names(dat_discrete)=c(names(dat_discrete)[1:7],"logAFP",names(dat_discrete)[9:16])

datcrete_con_weight<-weight.esti(dat_discrete[,-c(1:3)],dat_discrete$trt)

#save datcrete_con_weight
imp_wei_table<-cbind(ps_wei_imp=datcrete_con_weight$ps,ebal_wei_imp=datcrete_con_weight$ebal,
                     cbps_wei_imp=datcrete_con_weight$cbps,energy_wei_imp=datcrete_con_weight$energy,
                     kernel_wei_imp=datcrete_con_weight$kernel,
                     mb_wei_imp=datcrete_con_weight$MB,mdabw_wei_imp=datcrete_con_weight$MDABW
)
write.csv(file = 'imp_wei_table.csv',x=imp_wei_table[,1:7])

#warning:The optimization failed to converge. See Notes section at ?method_energy for information.
###plot KM====

png("AKM_imp.png",width = 10, height = 6, units = "in", res = 1000)

layout_matrix <- rbind(
  c(1,1,1,1),
  c(2,3,4,5),
  c(6,7,8,9)
)
layout(layout_matrix, heights = c(0.6, 3.5, 3.5))

par(mar = c(0,0,0,0))
plot.new()
legend("center",
       legend = c("LT","LR","LA"),
       col = c(3,4,5),
       lty = 1,
       lwd = 2,
       horiz = TRUE,
       cex = 1.4,
       y.intersp = 0.7,
       bty = "n")

par(mar = c(6,4,2,1))

akm_rmst(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
         group =as.factor(dat_discrete$trt))
se_ci_unad=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,tau = 1606,
                    group =as.factor(dat_discrete$trt))
title(main = "UNAD")

akm_rmst(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
         group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$ps)
se_ci_ps=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,tau = 1606,
                  group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$ps)
title(main = "IPW")


akm_rmst(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
         group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$cbps)
se_ci_cbps=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,tau = 1606,
                    group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$cbps)
title(main = "CBPS")

akm_rmst(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
         group =as.factor(dat_discrete$trt),weight = (datcrete_con_weight$ebal+min(datcrete_con_weight$ebal[which(datcrete_con_weight$ebal>10^(-7))])*(datcrete_con_weight$ebal==0)))
se_ci_ebal=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,tau = 1606,
                    group =as.factor(dat_discrete$trt),weight = (datcrete_con_weight$ebal+min(datcrete_con_weight$ebal[which(datcrete_con_weight$ebal>10^(-7))])*(datcrete_con_weight$ebal==0)))
title(main = "EBAL")

akm_rmst(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
         group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$energy)
se_ci_energy=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,tau = 1606,
                      group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$energy)
title(main = "ENERGY")

akm_rmst(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
         group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$kernel)
se_ci_kernel=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,tau = 1606,
                      group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$kernel)
title(main = "KERNEL")


akm_rmst(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
         group =as.factor(dat_discrete$trt),weight = (datcrete_con_weight$MDABW+min(datcrete_con_weight$MDABW[which(datcrete_con_weight$MDABW>10^(-7))])*(datcrete_con_weight$MDABW==0)))
se_ci_MDABW=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,tau = 1606,
                     group =as.factor(dat_discrete$trt),weight = (datcrete_con_weight$MDABW+min(datcrete_con_weight$MDABW[which(datcrete_con_weight$MDABW>10^(-7))])*(datcrete_con_weight$MDABW==0)))
title(main = "MDABW")

akm_rmst(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
         group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$MB)
se_ci_mb=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,tau = 1606,
                  group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$MB)
title(main = "MB")

dev.off()

se_ci_imp<-rbind(
  se_ci_unad[1,],
  se_ci_ps[1,],
  se_ci_ebal[1,],
  se_ci_cbps[1,],
  se_ci_energy[1,],
  se_ci_kernel[1,],
  se_ci_MDABW[1,],
  se_ci_mb[1,],
  se_ci_unad[2,],
  se_ci_ps[2,],
  se_ci_ebal[2,],
  se_ci_cbps[2,],
  se_ci_energy[2,],
  se_ci_kernel[2,],
  se_ci_MDABW[2,],
  se_ci_mb[2,],
  se_ci_unad[3,],
  se_ci_ps[3,],
  se_ci_ebal[3,],
  se_ci_cbps[3,],
  se_ci_energy[3,],
  se_ci_kernel[3,],
  se_ci_MDABW[3,],
  se_ci_mb[3,]
)
se_ci_imp[,-c(1,6)]<-round(se_ci_imp[,-c(1,6)],digits = 2)
se_ci_imp[,6]<-round(se_ci_imp[,6],digits = 3)
se_ci_imp
se_ci_imp<-cbind(weight=rep(c('Unad','IPW','EBAL','CBPS','Energy','KERNEL','mdabw','MB'),3),se_ci_imp,
                 CI=paste(
                   paste(
                     paste('(',se_ci_imp$rmst_diff_low),se_ci_imp$rmst_diff_upp,sep = ','),')'))
write.csv(file = 'se_ci_imp.1606.csv',x=se_ci_imp[,c(1,2,3,4,8,7)])

###trend with tau plot====
####get value====
ps_tau=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
cbps_tau=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
energy_tau=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
mb_tau=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
ebal_tau=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
no_tau=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
kernel_tau=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
MDABW_tau=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))

for (i in 1:1606) {
  ps_tau[i,]=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
                      group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$ps,tau=i)[,2]
  
  
  cbps_tau[i,]=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
                        group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$cbps,tau=i)[,2]
  
  
  
  energy_tau[i,]=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
                          group =as.factor(dat_discrete$trt),weight = datcrete_con_weight$energy,tau=i)[,2]
  
  
  mb_tau[i,]=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
                      group =as.factor(dat_discrete$trt),weight =datcrete_con_weight$MB,tau=i)[,2]
  
  ebal_tau[i,]=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
                        group =as.factor(dat_discrete$trt),weight =  (datcrete_con_weight$ebal+min(datcrete_con_weight$ebal[which(datcrete_con_weight$ebal>10^(-7))])*(datcrete_con_weight$ebal<10^(-7))),tau=i)[,2]
  
  no_tau[i,]=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
                      group =as.factor(dat_discrete$trt),tau=i)[,2]
  
  kernel_tau[i,]=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
                          group =as.factor(dat_discrete$trt),weight =datcrete_con_weight$kernel,tau=i)[,2]
  
  MDABW_tau[i,]=get_diff(time=as.numeric(dat_discrete[,1]),status = dat_discrete[,2],xaxismax=2000,
                         group =as.factor(dat_discrete$trt),weight =(datcrete_con_weight$MDABW+min(datcrete_con_weight$MDABW[which(datcrete_con_weight$MDABW>10^(-7))])*(datcrete_con_weight$MDABW<10^(-7))),tau=i)[,2]
  
  
}
ps_tau1=cbind(tau=1:1606,ps_tau)
cbps_tau1=cbind(tau=1:1606,cbps_tau)
energy_tau1=cbind(tau=1:1606,energy_tau)
mb_tau1=cbind(tau=1:1606,mb_tau)
ebal_tau1=cbind(tau=1:1606,ebal_tau)
no_tau1=cbind(tau=1:1606,no_tau)
kernel_tau1=cbind(tau=1:1606,kernel_tau)
MDABW_tau1=cbind(tau=1:1606,MDABW_tau)
####plot======
library(ggplot2)
library(reshape2)
library(ggpubr)
p1 <- ggplot()+
  geom_line(data = ps_tau1,aes(x = tau,y = vs12,colour = "IPW"),size=1)+
  geom_line(data = cbps_tau1,aes(x = tau,y = vs12,colour = "CBPS"),size=1)+
  geom_line(data = energy_tau1,aes(x = tau,y = vs12,colour = "ENERGY"),size=1)+
  geom_line(data = mb_tau1,aes(x = tau,y = vs12,colour = "MB"),size=1)+
  geom_line(data = ebal_tau1,aes(x = tau,y = vs12,colour = "EBAL"),size=1)+
  geom_line(data = no_tau1,aes(x = tau,y = vs12,colour = "UNAD"),size=1)+
  geom_line(data = kernel_tau1,aes(x = tau,y = vs12,colour = "KERNEL"),size=1)+
  geom_line(data = MDABW_tau1,aes(x = tau,y = vs12,colour = "MDABW"),size=1)+
  scale_y_continuous(
    limits = c(-450, 10),
    breaks = c(-450,-400,-350,-300,-250,-200,-150,-100,-50,0)) +
  scale_x_continuous(
    limits = c(0, 1650), breaks = c(0, 500, 1000, 1500))+
  ggtitle("LT vs. LR")+
  scale_color_manual(
    name='Method', 
    values=c('UNAD'='skyblue','MB'='blue','IPW'='black','EBAL'='orange',
             'CBPS'='red','ENERGY'='green','KERNEL'='#D55E00','MDABW'='#CC79A7'),
    breaks = c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')
  )+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.box.background = element_rect(color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(5,5,5,5)
  )+
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


p2 <- ggplot()+
  geom_line(data = ps_tau1,aes(x = tau,y = vs13,colour = "IPW"),size=1)+
  geom_line(data = cbps_tau1,aes(x = tau,y = vs13,colour = "CBPS"),size=1)+
  geom_line(data = energy_tau1,aes(x = tau,y = vs13,colour = "ENERGY"),size=1)+
  geom_line(data = mb_tau1,aes(x = tau,y = vs13,colour = "MB"),size=1)+
  geom_line(data = ebal_tau1,aes(x = tau,y = vs13,colour = 'EBAL'),size=1)+
  geom_line(data = no_tau1,aes(x = tau,y = vs13,colour = "UNAD"),size=1)+
  geom_line(data = kernel_tau1,aes(x = tau,y = vs13,colour = "KERNEL"),size=1)+
  geom_line(data = MDABW_tau1,aes(x = tau,y = vs13,colour = "MDABW"),size=1)+
  scale_y_continuous(
    limits = c(-750, 10),
    breaks = c(-750,-700,-650,-600,-550,-500,-450,-400,-350,-300,-250,-200,-150,-100,-50,0)) +
  scale_x_continuous(
    limits = c(0, 1650), breaks = c(0, 500, 1000, 1500))+
  ggtitle("LT vs. LR")+
  scale_color_manual(
    name='Method', 
    values=c('UNAD'='skyblue','MB'='blue','IPW'='black','EBAL'='orange',
             'CBPS'='red','ENERGY'='green','KERNEL'='#D55E00','MDABW'='#CC79A7'),
    breaks = c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')
  )+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.box.background = element_rect(color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(5,5,5,5)
  )+
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


p3 <- ggplot()+
  geom_line(data = ps_tau1,aes(x = tau,y = vs23,colour = "IPW"),size=1)+
  geom_line(data = cbps_tau1,aes(x = tau,y = vs23,colour = "CBPS"),size=1)+
  geom_line(data = energy_tau1,aes(x = tau,y = vs23,colour = "ENERGY"),size=1)+
  geom_line(data = mb_tau1,aes(x = tau,y = vs23,colour = "MB"),size=1)+
  geom_line(data = ebal_tau1,aes(x = tau,y = vs23,colour = "EBAL"),size=1)+
  geom_line(data = no_tau1,aes(x = tau,y = vs23,colour = "UNAD"),size=1)+
  geom_line(data = kernel_tau1,aes(x = tau,y = vs23,colour = "KERNEL"),size=1)+
  geom_line(data = MDABW_tau1,aes(x = tau,y = vs23,colour = "MDABW"),size=1)+
  scale_y_continuous(
    limits = c(-350, 20),
    breaks = c(-350,-300,-250,-200,-150,-100,-50,0)) +
  scale_x_continuous(
    limits = c(0, 1650), breaks = c(0, 500, 1000, 1500))+
  ggtitle("LR vs. LA")+
  scale_color_manual(
    name='Method', 
    values=c('UNAD'='skyblue','MB'='blue','IPW'='black','EBAL'='orange',
             'CBPS'='red','ENERGY'='green','KERNEL'='#D55E00','MDABW'='#CC79A7'),
    breaks = c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')
  )+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.box.background = element_rect(color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(5,5,5,5)
  )+
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


png("rmst_imp.png",width = 8, height = 5.6, units = "in", res = 1000)
ggarrange(p1,p2,p3,nrow = 1,ncol = 3, common.legend = TRUE, legend = 'top')
dev.off()





chara1<-colnames(dat_discrete)
ncolume1<-ncol(dat_discrete)
crete_formula   <- as.formula(paste(chara1[3],
                                    paste(" ~ ", paste(chara1[4:ncolume1], 
                                                       collapse= "+"))))
#names(datcrete_con_weight)<-c('IPW','EBAL','CBPS','ENERGY','MB')


library(cobalt)
plot_datcrete<-dat_discrete
plot_datcrete$trt<-factor(plot_datcrete$trt)

# obtain asmd from bal.tab function
plot_datcrete_ASMD_baltab <-  data.frame(
  Variable = rownames(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                              stats = c("mean.diffs"),
                              weights = rep(1, nrow(plot_datcrete)))$Balance.Across.Pairs),
  Unadjusted = as.numeric(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                                  stats = "mean.diffs")$Balance.Across.Pairs$Max.Diff.Un),
  PS = as.numeric(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                          stats = c("mean.diffs"), weights = datcrete_con_weight$ps)$Balance.Across.Pairs$Max.Diff.Adj),
  EBAL = as.numeric(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                            stats = c("mean.diffs"), weights = datcrete_con_weight$ebal)$Balance.Across.Pairs$Max.Diff.Adj),
  CBPS = as.numeric(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                            stats = c("mean.diffs"), weights = datcrete_con_weight$cbps)$Balance.Across.Pairs$Max.Diff.Adj),
  ENERGY = as.numeric(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                              stats = c("mean.diffs"), weights = datcrete_con_weight$energy)$Balance.Across.Pairs$Max.Diff.Adj),
  KERNEL = as.numeric(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                              stats = c("mean.diffs"), weights = datcrete_con_weight$kernel)$Balance.Across.Pairs$Max.Diff.Adj),
  MDABW = as.numeric(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                             stats = c("mean.diffs"), weights = datcrete_con_weight$MDABW)$Balance.Across.Pairs$Max.Diff.Adj),
  MB = as.numeric(bal.tab(crete_formula, data = plot_datcrete, estimand = "ATE",
                          stats = c("mean.diffs"), weights = datcrete_con_weight$MB)$Balance.Across.Pairs$Max.Diff.Adj)
)


ASMD_crete_baltab <- c(plot_datcrete_ASMD_baltab[,2],
                       plot_datcrete_ASMD_baltab[,9],
                       plot_datcrete_ASMD_baltab[,3],
                       plot_datcrete_ASMD_baltab[,4],
                       plot_datcrete_ASMD_baltab[,5],
                       plot_datcrete_ASMD_baltab[,6],
                       plot_datcrete_ASMD_baltab[,7],
                       plot_datcrete_ASMD_baltab[,8])

groupcrete_ASMD_baltab <- rep(c('UNAD','MB','IPW','EBAL','CBPS','ENERGY','KERNEL','MDABW'),
                              each = nrow(plot_datcrete_ASMD_baltab))

varnames <- plot_datcrete_ASMD_baltab$Variable

plot_datcrete_ASMD_baltab2new <- data.frame(
  id = rep(1:nrow(plot_datcrete_ASMD_baltab), 8),
  variable = rep(varnames, 8),
  ASMD_ord = as.numeric(ASMD_crete_baltab),
  asmd = pmin(1, as.numeric(ASMD_crete_baltab)),
  weight = groupcrete_ASMD_baltab,
  stringsAsFactors = FALSE
)


shape_map <- c(
  UNAD  = 15,
  MB    = 16,
  IPW   = 17,
  EBAL  = 18,
  CBPS  = 3,
  ENERGY= 4,
  KERNEL= 5,
  MDABW = 6
)

color_map <- c(
  UNAD  = 'skyblue',
  MB    = 'blue',
  IPW   = 'black',
  EBAL  = 'orange',
  CBPS  = 'red',
  ENERGY= 'green',
  KERNEL= '#D55E00',
  MDABW = '#CC79A7'
)


desired_levels <- c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')

plot_datcrete_ASMD_baltab2new$weight <- factor(plot_datcrete_ASMD_baltab2new$weight,
                                               levels = desired_levels)


order_unad <- plot_datcrete_ASMD_baltab2new %>%
  filter(weight == "UNAD") %>%
  arrange(ASMD_ord) %>%
  pull(variable)

plot_datcrete_ASMD_baltab2new$variable <- factor(plot_datcrete_ASMD_baltab2new$variable, levels = order_unad)


p_baltab <- ggplot(plot_datcrete_ASMD_baltab2new, aes(x = variable, y = asmd)) +
  geom_point(aes(color = weight, shape = weight), size = 3) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(color = 'Method', shape = 'Method',
       title = "Covariate Balance",
       subtitle = "Generalized Absolute Standardized Mean Difference",
       x = NULL, y = "Generalized Absolute Standardized Mean Difference") +
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0,0.1,0.25,0.5,0.75,1.0)) +
  coord_flip()

png("ASMD_imp_cobalt.png", width = 6.4, height = 6.4, units = "in", res = 1000)
print(p_baltab)
dev.off()


#save asmd calculated by "bal.tab" function

write.csv(file = 'asmd_cobalt_imp.csv',x=plot_datcrete_ASMD_baltab[,1:9])


plot_datcrete2 <- data.frame(datcrete_con_weight$Imbalance)


methods_order <- c('unad','mb','ipw','ebal','cbps','energy','kernel','mdabw')


if(!all(methods_order %in% tolower(colnames(plot_datcrete2)))) {
  
  imbal_crete2.1 <- c(plot_datcrete2[,1],  # UNAD
                      plot_datcrete2[,6],  # MB
                      plot_datcrete2[,2],  # IPW
                      plot_datcrete2[,3],  # EBAL
                      plot_datcrete2[,4],  # CBPS
                      plot_datcrete2[,5],  # ENERGY
                      plot_datcrete2[,7],  # KERNEL
                      plot_datcrete2[,8])  # MDABW
  
  groupcrete2.1 <- rep(c('UNAD','MB','IPW','EBAL','CBPS','ENERGY','KERNEL','MDABW'),
                       each = nrow(plot_datcrete2))
  
  plot_datcrete2new <- data.frame(
    variable = rep(rownames(plot_datcrete2), 8),
    imbal = imbal_crete2.1,
    weight = groupcrete2.1,
    stringsAsFactors = FALSE
  )
} else {
  
  plot_datcrete2_named <- plot_datcrete2 %>%
    
    setNames(tolower(names(.))) %>%
    tibble::rownames_to_column("variable")
  
  plot_datcrete2new <- plot_datcrete2_named %>%
    pivot_longer(cols = all_of(methods_order),
                 names_to = "weight", values_to = "imbal") %>%
    mutate(weight = toupper(weight)) %>%
    rename(variable = variable)
}


shape_map <- c(
  UNAD  = 15,
  MB    = 16,
  IPW   = 17,
  EBAL  = 18,
  CBPS  = 3,
  ENERGY= 4,
  KERNEL= 5,
  MDABW = 6
)

color_map <- c(
  UNAD  = 'skyblue',
  MB    = 'blue',
  IPW   = 'black',
  EBAL  = 'orange',
  CBPS  = 'red',
  ENERGY= 'green',
  KERNEL= '#D55E00',
  MDABW = '#CC79A7'
)


desired_levels <- c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')


plot_datcrete2new$weight <- factor(toupper(as.character(plot_datcrete2new$weight)),
                                   levels = desired_levels)


if("UNAD" %in% plot_datcrete2new$weight) {
  order_unad <- plot_datcrete2new %>%
    filter(weight == "UNAD") %>%
    arrange(imbal) %>%
    pull(variable) %>%
    unique()
  plot_datcrete2new$variable <- factor(plot_datcrete2new$variable, levels = order_unad)
} else {
  
  plot_datcrete2new$variable <- factor(plot_datcrete2new$variable, levels = unique(plot_datcrete2new$variable))
}


plot.GASMD_IMP <- ggplot(plot_datcrete2new, aes(x = variable, y = imbal)) + 
  geom_point(aes(color = weight, shape = weight), size = 3) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(color = 'Method', shape = 'Method',
       title = "Covariate Balance", 
       subtitle = "Generalized Absolute Standardized Mean Difference", 
       x = NULL, y = "Generalized Absolute Standardized Mean Difference") +
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  coord_flip()


png("GASMD_imp.png", width = 6.4, height = 6.4, units = "in", res = 1000)
print(plot.GASMD_IMP)
dev.off()


#save GASMD_imp 
write.csv(file = 'gasmd_imp.csv',x=plot_datcrete2[,1:8])


##all variables
###etiology-dummy variable
etiology_dum<-class.ind(hcc$Etiology)
etiology_crete<-data.frame('etiology2'=etiology_dum[,2],'etiology3'=etiology_dum[,3])
####dat is used in the important variable situation
####except BCLC
####hcc colomn=39, 35:39 denote different date (org_date,diag_date,death_date,_last_date)

#dat_all<-data.frame(dat,BMI=hcc$BMI,BMI.cat=hcc$BMI.cat,hypertension=hcc$hypertension,diabetes.type2=hcc$diabetes.type2,
#                    etiology_crete,
#                    Ascites=hcc$Ascites,Cirrhosis=hcc$Cirrhosis,BCLC=hcc$BCLC,GGT=hcc$GGT,ALB=hcc$ALB,Cr=hcc$Cr,PT=hcc$PT,PLT=hcc$PLT,
#                    portal.hypertension=hcc$portal.hypertension,AST=hcc$AST,TBIL=hcc$TBIL,ALP=hcc$ALP,Radiology.HCC.Number=hcc$Radiology.HCC.Number,
#                    Radiology.Largest.Tumor.Diameter=hcc$Radiology.Largest.Tumor.Diameter)
#use BMI.cat instead of BMI
dat_all<-data.frame(dat,BMI.cat=hcc$BMI.cat,hypertension=hcc$hypertension,diabetes.type2=hcc$diabetes.type2,
                    etiology_crete,
                    Ascites=hcc$Ascites,Cirrhosis=hcc$Cirrhosis,BCLC=hcc$BCLC,GGT=hcc$GGT,ALB=hcc$ALB,Cr=hcc$Cr,PT=hcc$PT,PLT=hcc$PLT,
                    portal.hypertension=hcc$portal.hypertension,AST=hcc$AST,TBIL=hcc$TBIL,ALP=hcc$ALP,Radiology.HCC.Number=hcc$Radiology.HCC.Number,
                    Radiology.Largest.Tumor.Diameter=hcc$Radiology.Largest.Tumor.Diameter)
dat_all_crete<-dat_all[,-(8:12)]
log_allcrete<-dat_all_crete
#log_allcrete[,c(7,16,25:29,31:33)]<-log(dat_all_crete[,c(7,16,25:29,31:33)])
log_allcrete[,c(7,24:28,30:32)]<-log(dat_all_crete[,c(7,24:28,30:32)])
dat_all_crete<-log_allcrete


###the ps score corresponding to individual 109 is 2.439892e+06
dat_all_crete<-dat_all_crete[-109,]

#names(dat_discrete)=c(names(dat_discrete)[1:7],"logAFP",names(dat_discrete)[9:16])

allcrete_con_weight<-weight.esti.nokernel(dat_all_crete[,-c(1:3)],dat_all_crete$trt)

#save allcrete_weight
all_wei_table<-cbind(ps_wei_all=allcrete_con_weight$ps,ebal_wei_all=allcrete_con_weight$ebal,
                     cbps_wei_all=allcrete_con_weight$cbps,energy_wei_all=allcrete_con_weight$energy,
                     mb_wei_all=allcrete_con_weight$MB,mdabw_wei_all=allcrete_con_weight$MDABW
)
write.csv(file = 'all_wei_table.csv',x=all_wei_table[,1:6])

png("AKM_all.png",width = 10, height = 6, units = "in", res = 1000)



layout_matrix <- rbind(
  c(1,1,1,1),
  c(2,3,4,5),
  c(6,7,8,9)
)
layout(layout_matrix, heights = c(0.6, 3.5, 3.5))

par(mar = c(0,0,0,0))
plot.new()
legend("center",
       legend = c("LT","LR","LA"),
       col = c(3,4,5),
       lty = 1,
       lwd = 2,
       horiz = TRUE,
       cex = 1.4,
       y.intersp = 0.7,
       bty = "n")

par(mar = c(6,4,2,1))

akm_rmst(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
         group =as.factor(dat_all_crete$trt))
se_ci_unad_all=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,tau = 1606,
                        group =as.factor(dat_all_crete$trt))
title(main = "UNAD")

akm_rmst(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
         group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$ps)
se_ci_ps_all=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,tau = 1606,
                      group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$ps)
title(main = "IPW")


akm_rmst(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
         group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$cbps)
se_ci_cbps_all=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,tau = 1606,
                        group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$cbps)
title(main = "CBPS")

akm_rmst(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
         group =as.factor(dat_all_crete$trt),weight = (allcrete_con_weight$ebal+min(allcrete_con_weight$ebal[which(allcrete_con_weight$ebal>10^(-7))])*(allcrete_con_weight$ebal<10^(-7))))
se_ci_ebal_all=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,tau = 1606,
                        group =as.factor(dat_all_crete$trt),weight = (allcrete_con_weight$ebal+min(allcrete_con_weight$ebal[which(allcrete_con_weight$ebal>10^(-7))])*(allcrete_con_weight$ebal<10^(-7))))
title(main = "EBAL")

akm_rmst(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
         group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$energy)
se_ci_energy_all=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,tau = 1606,
                          group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$energy)
title(main = "ENERGY")

##### we cannot obtain akm related to kernel method
plot.new()


akm_rmst(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
         group =as.factor(dat_all_crete$trt),weight = (allcrete_con_weight$MDABW+min(allcrete_con_weight$MDABW[which(allcrete_con_weight$MDABW>10^(-7))])*(allcrete_con_weight$MDABW<10^(-7))))
se_ci_MDABW_all=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,tau = 1606,
                         group =as.factor(dat_all_crete$trt),weight = (allcrete_con_weight$MDABW+min(allcrete_con_weight$MDABW[which(allcrete_con_weight$MDABW>10^(-7))])*(allcrete_con_weight$MDABW<10^(-7))))
title(main = "MDABW")

akm_rmst(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
         group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$MB)
se_ci_mb_all=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,tau = 1606,
                      group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$MB)
title(main = "MB")

dev.off()

se_ci_all<-rbind(
  se_ci_unad_all[1,],
  se_ci_ps_all[1,],
  se_ci_ebal_all[1,],
  se_ci_cbps_all[1,],
  se_ci_energy_all[1,],
  se_ci_MDABW_all[1,],
  se_ci_mb_all[1,],
  se_ci_unad_all[2,],
  se_ci_ps_all[2,],
  se_ci_ebal_all[2,],
  se_ci_cbps_all[2,],
  se_ci_energy_all[2,],
  se_ci_MDABW_all[2,],
  se_ci_mb_all[2,],
  se_ci_unad_all[3,],
  se_ci_ps_all[3,],
  se_ci_ebal_all[3,],
  se_ci_cbps_all[3,],
  se_ci_energy_all[3,],
  se_ci_MDABW_all[3,],
  se_ci_mb_all[3,]
)
se_ci_all[,-c(1,6)]<-round(se_ci_all[,-c(1,6)],digits = 2)
se_ci_all[,6]<-round(se_ci_all[,6],digits = 3)
se_ci_all
se_ci_all<-cbind(weight=rep(c('Unad','IPW','EBAL','CBPS','Energy','MDABW','MB'),3),se_ci_all,
                 CI=paste(
                   paste(
                     paste('(',se_ci_all$rmst_diff_low),se_ci_all$rmst_diff_upp,sep = ','),')'))
write.csv(file = 'se_ci_all.1606.csv',x=se_ci_all[,c(1,2,3,8,7)])

#trend with tau plot

#get value tau at each time====

ps_tau_allcrete=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
cbps_tau_allcrete=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
energy_tau_allcrete=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
mb_tau_allcrete=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
ebal_tau_allcrete=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
no_tau_allcrete=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
MDABW_tau_allcrete=data.frame(vs12=rep(0,1606),vs13=rep(0,1606),vs23=rep(0,1606))
for (i in 1:1606) {
  ps_tau_allcrete[i,]=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
                               group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$ps,tau=i)[,2]
  
  
  cbps_tau_allcrete[i,]=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
                                 group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$cbps,tau=i)[,2]
  
  
  
  energy_tau_allcrete[i,]=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
                                   group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$energy,tau=i)[,2]
  
  
  mb_tau_allcrete[i,]=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
                               group =as.factor(dat_all_crete$trt),weight = allcrete_con_weight$MB,tau=i)[,2]
  
  ebal_tau_allcrete[i,]=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
                                 group =as.factor(dat_all_crete$trt),weight = (allcrete_con_weight$ebal+min(allcrete_con_weight$ebal[which(allcrete_con_weight$ebal>10^(-7))])*(allcrete_con_weight$ebal<10^(-7))),tau=i)[,2]
  
  no_tau_allcrete[i,]=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
                               group =as.factor(dat_all_crete$trt),tau=i)[,2]
  
  MDABW_tau_allcrete[i,]=get_diff(time=as.numeric(dat_all_crete[,1]),status = dat_all_crete[,2],xaxismax=2000,
                                  group =as.factor(dat_all_crete$trt),tau=i,weight = (allcrete_con_weight$MDABW+min(allcrete_con_weight$MDABW[which(allcrete_con_weight$MDABW>10^(-7))])*(allcrete_con_weight$MDABW<10^(-7))))[,2]
}

ps_tau_allcrete1=cbind(tau=1:1606,ps_tau_allcrete)
cbps_tau_allcrete1=cbind(tau=1:1606,cbps_tau_allcrete)
energy_tau_allcrete1=cbind(tau=1:1606,energy_tau_allcrete)
mb_tau_allcrete1=cbind(tau=1:1606,mb_tau_allcrete)
ebal_tau_allcrete1=cbind(tau=1:1606,ebal_tau_allcrete)
no_tau_allcrete1=cbind(tau=1:1606,no_tau_allcrete)
MDABW_tau_allcrete1=cbind(tau=1:1606,MDABW_tau_allcrete)


p1.all <- ggplot() +
  geom_line(data = no_tau_allcrete1, aes(x = tau, y = vs12, colour = "UNAD"), size = 1) +
  geom_line(data = mb_tau_allcrete1, aes(x = tau, y = vs12, colour = "MB"), size = 1) +
  geom_line(data = ps_tau_allcrete1, aes(x = tau, y = vs12, colour = "IPW"), size = 1) +
  geom_line(data = ebal_tau_allcrete1, aes(x = tau, y = vs12, colour = "EBAL"), size = 1) +
  geom_line(data = cbps_tau_allcrete1, aes(x = tau, y = vs12, colour = "CBPS"), size = 1) +
  geom_line(data = energy_tau_allcrete1, aes(x = tau, y = vs12, colour = "ENERGY"), size = 1) +
  geom_line(data = MDABW_tau_allcrete1, aes(x = tau, y = vs12, colour = "MDABW"), size = 1) +
  scale_y_continuous(
    limits = c(-450, 60),
    breaks = seq(-450, 50, by = 50)
  ) +
  scale_x_continuous(
    limits = c(0, 1650),
    breaks = c(0, 500, 1000, 1500)
  ) +
  ggtitle("LT vs. LR") +
  scale_color_manual(
    name = "Method",
    values = c(
      UNAD = "skyblue",
      MB = "blue",
      IPW = "black",
      EBAL = "orange",
      CBPS = "red",
      ENERGY = "green",
      MDABW = "#CC79A7"
    ),
    breaks = c("UNAD", "ENERGY", "IPW", "MDABW", "CBPS", "MB", "EBAL")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


p2.all <- ggplot() +
  geom_line(data = ps_tau_allcrete1, aes(x = tau, y = vs13, colour = "IPW"), size = 1) +
  geom_line(data = cbps_tau_allcrete1, aes(x = tau, y = vs13, colour = "CBPS"), size = 1) +
  geom_line(data = energy_tau_allcrete1, aes(x = tau, y = vs13, colour = "ENERGY"), size = 1) +
  geom_line(data = mb_tau_allcrete1, aes(x = tau, y = vs13, colour = "MB"), size = 1) +
  geom_line(data = ebal_tau_allcrete1, aes(x = tau, y = vs13, colour = "EBAL"), size = 1) +
  geom_line(data = no_tau_allcrete1, aes(x = tau, y = vs13, colour = "UNAD"), size = 1) +
  geom_line(data = MDABW_tau_allcrete1, aes(x = tau, y = vs13, colour = "MDABW"), size = 1) +
  scale_y_continuous(
    limits = c(-700, 55),
    breaks = seq(-700, 50, by = 50)
  ) +
  scale_x_continuous(
    limits = c(0, 1650),
    breaks = c(0, 500, 1000, 1500)
  ) +
  ggtitle("LT vs. LA") +
  scale_color_manual(
    name = "Method",
    values = c(
      UNAD = "skyblue",
      MB = "blue",
      IPW = "black",
      EBAL = "orange",
      CBPS = "red",
      ENERGY = "green",
      MDABW = "#CC79A7"
    ),
    breaks = c("UNAD", "ENERGY", "IPW", "MDABW", "CBPS", "MB", "EBAL")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


p3.all <- ggplot() +
  geom_line(data = ps_tau_allcrete1, aes(x = tau, y = vs23, colour = "IPW"), size = 1) +
  geom_line(data = cbps_tau_allcrete1, aes(x = tau, y = vs23, colour = "CBPS"), size = 1) +
  geom_line(data = energy_tau_allcrete1, aes(x = tau, y = vs23, colour = "ENERGY"), size = 1) +
  geom_line(data = mb_tau_allcrete1, aes(x = tau, y = vs23, colour = "MB"), size = 1) +
  geom_line(data = ebal_tau_allcrete1, aes(x = tau, y = vs23, colour = "EBAL"), size = 1) +
  geom_line(data = no_tau_allcrete1, aes(x = tau, y = vs23, colour = "UNAD"), size = 1) +
  geom_line(data = MDABW_tau_allcrete1, aes(x = tau, y = vs23, colour = "MDABW"), size = 1) +
  scale_y_continuous(
    limits = c(-300, 13),
    breaks = seq(-300, 0, by = 50)
  ) +
  scale_x_continuous(
    limits = c(0, 1650),
    breaks = c(0, 500, 1000, 1500)
  ) +
  ggtitle("LR vs. LA") +
  scale_color_manual(
    name = "Method",
    values = c(
      UNAD = "skyblue",
      MB = "blue",
      IPW = "black",
      EBAL = "orange",
      CBPS = "red",
      ENERGY = "green",
      MDABW = "#CC79A7"
    ),
    breaks = c("UNAD", "ENERGY", "IPW", "MDABW", "CBPS", "MB", "EBAL")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


png("rmst_all.png", width = 8, height = 5.6, units = "in", res = 1000)
ggarrange(
  p1.all, p2.all, p3.all,
  nrow = 1, ncol = 3,
  common.legend = TRUE,
  legend = "top"
)
dev.off()



###love.plot====
chara2<-names(dat_all_crete)
ncolume2<-ncol(dat_all_crete)
allcrete_formula   <- as.formula(paste(chara2[3],
                                       paste(" ~ ", paste(chara2[4:ncolume2], 
                                                          collapse= "+"))))

# obtain asmd from bal.tab function
plot_dat2<-dat_all_crete
plot_dat2$trt<-factor(plot_dat2$trt)

plot_allcrete_ASMD_baltab <-  data.frame(
  Variable = rownames(bal.tab(allcrete_formula, data = plot_dat2, estimand = "ATE",
                              stats = c("mean.diffs"), weights = rep(1, nrow(plot_dat2)))$Balance.Across.Pairs),
  Unadjusted = as.numeric(bal.tab(allcrete_formula, data = plot_dat2, estimand = "ATE",
                                  stats = "mean.diffs")$Balance.Across.Pairs$Max.Diff.Un),
  PS = as.numeric(bal.tab(allcrete_formula, data = plot_dat2, estimand = "ATE",
                          stats = c("mean.diffs"), weights = allcrete_con_weight$ps)$Balance.Across.Pairs$Max.Diff.Adj),
  EBAL = as.numeric(bal.tab(allcrete_formula, data = plot_dat2, estimand = "ATE",
                            stats = c("mean.diffs"), weights = allcrete_con_weight$ebal)$Balance.Across.Pairs$Max.Diff.Adj),
  CBPS = as.numeric(bal.tab(allcrete_formula, data = plot_dat2, estimand = "ATE",
                            stats = c("mean.diffs"), weights = allcrete_con_weight$cbps)$Balance.Across.Pairs$Max.Diff.Adj),
  ENERGY = as.numeric(bal.tab(allcrete_formula, data = plot_dat2, estimand = "ATE",
                              stats = c("mean.diffs"), weights = allcrete_con_weight$energy)$Balance.Across.Pairs$Max.Diff.Adj),
  MDABW = as.numeric(bal.tab(allcrete_formula, data = plot_dat2, estimand = "ATE",
                             stats = c("mean.diffs"), weights = allcrete_con_weight$MDABW)$Balance.Across.Pairs$Max.Diff.Adj),
  MB = as.numeric(bal.tab(allcrete_formula, data = plot_dat2, estimand = "ATE",
                          stats = c("mean.diffs"), weights = allcrete_con_weight$MB)$Balance.Across.Pairs$Max.Diff.Adj)
)


ASMD_allcrete_baltab <- c(plot_allcrete_ASMD_baltab[,2],  # UNAD
                          plot_allcrete_ASMD_baltab[,6],  # ENERGY
                          plot_allcrete_ASMD_baltab[,3],  # IPW
                          plot_allcrete_ASMD_baltab[,7],  # MDABW
                          plot_allcrete_ASMD_baltab[,5],  # CBPS
                          plot_allcrete_ASMD_baltab[,8],  # MB
                          plot_allcrete_ASMD_baltab[,4])  # EBAL

group_allcrete_ASMD_baltab <- rep(c('UNAD','ENERGY','IPW','MDABW','CBPS','MB','EBAL'),
                                  each = nrow(plot_allcrete_ASMD_baltab))

varnames_allcrete_ASMD_baltab <- plot_allcrete_ASMD_baltab$Variable

plot_allcrete_ASMD_baltab2new <- data.frame(
  id = rep(1:nrow(plot_allcrete_ASMD_baltab), 7),
  variable = rep(varnames_allcrete_ASMD_baltab, 7),
  asmd = pmin(1, as.numeric(ASMD_allcrete_baltab)),
  ASMD_ord = as.numeric(ASMD_allcrete_baltab),
  weight = group_allcrete_ASMD_baltab,
  stringsAsFactors = FALSE
)


desired_levels <- c('UNAD','ENERGY','IPW','MDABW','CBPS','MB','EBAL')

color_map <- c(
  UNAD  = 'skyblue',
  ENERGY= 'green',
  IPW   = 'black',
  MDABW = '#CC79A7',
  CBPS  = 'red',
  MB    = 'blue',
  EBAL  = 'orange'
)

shape_map <- c(
  UNAD  = 15,
  ENERGY= 4,
  IPW   = 17,
  MDABW = 6,
  CBPS  = 3,
  MB    = 16,
  EBAL  = 18
)

plot_allcrete_ASMD_baltab2new$weight <- factor(plot_allcrete_ASMD_baltab2new$weight,
                                               levels = desired_levels)


order_unad <- plot_allcrete_ASMD_baltab2new %>%
  filter(weight == "UNAD") %>%
  arrange(ASMD_ord) %>%
  pull(variable)

plot_allcrete_ASMD_baltab2new$variable <- factor(plot_allcrete_ASMD_baltab2new$variable,
                                                 levels = order_unad)


p_baltab_all <- ggplot(plot_allcrete_ASMD_baltab2new, aes(x = variable, y = asmd)) +
  geom_point(aes(color = weight, shape = weight), size = 3) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(color = 'Method', shape = 'Method',
       title = "Covariate Balance",
       subtitle = "Generalized Absolute Standardized Mean Difference",
       x = NULL, y = "Generalized Absolute Standardized Mean Difference") +
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1.0)) +
  coord_flip()


png("ASMD_all_cobalt.png", width = 6.4, height = 6.4, units = "in", res = 1000)
print(p_baltab_all)
dev.off()


#save asmd calculated by "bal.tab" function

write.csv(file = 'asmd_cobalt_all.csv',x=plot_allcrete_ASMD_baltab[,1:8])


plot_dat2.1 <- data.frame(allcrete_con_weight$Imbalance)

imbal_2.1 <- c(plot_dat2.1[,1],
               plot_dat2.1[,6],
               plot_dat2.1[,2],
               plot_dat2.1[,3],
               plot_dat2.1[,4],
               plot_dat2.1[,5],
               plot_dat2.1[,7])

group2.1 <- rep(c("UNAD", "MB", "IPW", "EBAL", "CBPS", "ENERGY","MDABW"),
                each = nrow(plot_dat2.1))

plot_dat2.1new <- data.frame(
  variable = rep(rownames(plot_dat2.1), 7),
  id = rep(1:nrow(plot_dat2.1), 7),
  imbal = imbal_2.1,
  weight = group2.1,
  stringsAsFactors = FALSE
)


desired_levels <- c('UNAD','ENERGY','IPW','MDABW','CBPS','MB','EBAL')

color_map <- c(
  UNAD  = 'skyblue',
  ENERGY= 'green',
  IPW   = 'black',
  MDABW = '#CC79A7',
  CBPS  = 'red',
  MB    = 'blue',
  EBAL  = 'orange'
)

shape_map <- c(
  UNAD  = 15,
  ENERGY= 4,
  IPW   = 17,
  MDABW = 6,
  CBPS  = 3,
  MB    = 16,
  EBAL  = 18
)


plot_dat2.1new$weight <- factor(plot_dat2.1new$weight, levels = desired_levels)


order_unad <- plot_dat2.1 %>%
  arrange(unad) %>%
  rownames()

plot_dat2.1new$variable <- factor(plot_dat2.1new$variable, levels = order_unad)


plot.GASMD_all <- ggplot(plot_dat2.1new, aes(x = variable, y = imbal)) + 
  geom_point(aes(color = weight, shape = weight), size = 3) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(color = 'Method', shape = 'Method',
       title = "Covariate Balance", 
       subtitle = "Generalized Absolute Standardized Mean Difference", 
       x = NULL, y = "Generalized Absolute Standardized Mean Difference") +
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  coord_flip()


png("GASMD_all.png", width = 6.4, height = 6.4, units = "in", res = 1000)
print(plot.GASMD_all)
dev.off()


# save gasmd_all
write.csv(file = 'gasmd_all.csv',x=plot_dat2.1[,1:7])



