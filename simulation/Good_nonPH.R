rm(list = ls())
n_total=500000
#n_total is the whole simulation num
n=100000
#n is the sample size used to estimate true rmst
#best include n_total
h=500
#h is a simulation contain h samples
s=n_total/h
#s is the times of simulation
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
################MDABW method
sbwbal <- function(X, Tr){
  data.matrix <- data.frame(X , factor(Tr))
  dimension   <- dim(X)[2]
  bal = list()
  bal$bal_gri = c(1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  character   <- names(data.matrix)
  bal$bal_cov <- character[1:dimension]
  sbw.result <- FALSE
  ######sbw的output有15项，若未输出正确结果则将bal_grid的搜索范围减少1项
  while(sum(dim(as.matrix(sbw.result))) != 15){
    sbw.result  <- tryCatch(sbw.result <- sbw(dat = data.matrix, ind = character[1 + dimension], bal = bal, out = character[2 + dimension], par = list(par_est = "ate")), error = function(e) { skip_to_next <<- FALSE})
    bal$bal_gri <- bal$bal_gri[-1]
  }
  sbw.weight <- sbw.result$dat_weights$sbw_weights
  return(list(weight = sbw.weight))
}
#two functions for Mahalanobis balancing
#1 Function for Mahalanobis Balancing
library(cobalt)
library(ggplot2)


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
  MB1.result  <- MB(x = X, treat = Tr, group1 = 1, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB")
  MB2.result  <- MB(x = X, treat = Tr, group1 = 2, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB")
  MB3.result  <- MB(x = X, treat = Tr, group1 = 3, outcome = rep(0,sample.size), delta.space = delta.space, method = "MB")
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
  if (sum(fit1$warns)) cat("lambda bound warning!\n")
  # compute weights for Tr=2
  fit2 <- ATE.ncb.SN(Tr_kernel_2, Kern, lam1s=lams)
  if (sum(fit2$warns)) cat("lambda bound warning!\n")
  # compute weights for Tr=3
  fit3 <- ATE.ncb.SN(Tr_kernel_3, Kern, lam1s=lams)
  if (sum(fit3$warns)) cat("lambda bound warning!\n")
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
        plot(NULL, xlim=c(xaxismin, xaxismax), ylim=c(0,1), xlab='Time',ylab='Adjusted Survival Probability')
        title(main='Adjusted Kaplan-Meier')
        
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
  abline(v=tau, col=1, lty=3, lwd=2)
  legend('bottomleft', paste("Group", groupval), lty=rep(1, j), lwd=rep(2, j), col=3:(j+2), 
         cex=.75, bty ="n", inset = c(0, 0))
  
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

get_akm_rmst<-function(time, status, group, weight=NULL, tau=NULL, alpha=.05, 
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
  return(data.frame('1_2'=output_diff$rmst_diff[1],'1_3'=output_diff$rmst_diff[2],'2_3'=output_diff$rmst_diff[3]))
}


library(simsurv)
library(survival)
library(survminer)
si.data <- function(sample.size = 500, dimension = 10){
  set.seed(12345)
  X   <- matrix(rnorm(sample.size * dimension), nrow = sample.size, ncol = dimension)
  ps1 <- rep(1,sample.size)
  ps2 <- exp(1 * apply(X[,1:5],1,sum))
  ps3 <- exp(1 * apply(X[,6:10],1,sum)) 
  ps  <- cbind(ps1 / (ps1 + ps2 + ps3),ps2 / (ps1 + ps2 + ps3),ps3 / (ps1 + ps2 + ps3))
  uniform <- runif(sample.size)
  index   <- cbind(ps[,1] < uniform, ps[,1] + ps[,2] < uniform)
  Tr      <- apply(index,1,sum) + 1
  return(list(X = X, Tr = Tr))
}


#non ph with weibull=======
library(mlogit)

data   <- si.data(sample.size = n_total)


#get true rmst of simulation=====

set.seed(2022)
true_sample=sample(1:n_total,n,replace = FALSE)

dat_true<-data$X[true_sample,];trt<-data$Tr[true_sample]


ncov1<-data.frame(id=1:n,trt=rep((1-1),n),X=dat_true[,c(1,2,3,4,5)])
ncov2<-data.frame(id=1:n,trt=rep((2-1),n),X=dat_true[,c(1,2,3,4,5)])
ncov3<-data.frame(id=1:n,trt=rep((3-1),n),X=dat_true[,c(1,2,3,4,5)])
set.seed(2022);
nsdat1<-simsurv(dist = 'weibull',lambdas=0.03,gammas = 3,
                betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=ncov1,
                maxt=10)

set.seed(2022);
nsdat2<-simsurv(dist = 'weibull',lambdas=0.2,gammas = 1.5,
                betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=ncov2,
                maxt=10)

set.seed(2022);
nsdat3<-simsurv(dist = 'weibull',lambdas=0.1,gammas = 0.9,
                betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=ncov3,
                maxt=10)

pro_cen=(sum(nsdat1$status==0)+sum(nsdat2$status==0)+sum(nsdat3$status==0))/(3*n)
#percent of censored
akm_rmst(time = c(nsdat1$eventtime,nsdat2$eventtime,nsdat3$eventtime)
         ,c(nsdat1$status,nsdat2$status,nsdat3$status),tau=10,
         group = as.factor(c(rep(0,n),rep(1,n),rep(2,n))))

theta_true_non=get_akm_rmst(time = c(nsdat1$eventtime,nsdat2$eventtime,nsdat3$eventtime)
                            ,c(nsdat1$status,nsdat2$status,nsdat3$status),tau=10,
                            group = as.factor(c(rep(0,n),rep(1,n),rep(2,n))))

#plot non ph weighted======
# 
# X=dat_true
# Tr=trt
# cov<-data.frame(id=1:nrow(X),trt=(Tr-1),X=X[,c(1,2,3,4,5)])
# #Change the order of X, tr, cov, and survival time according to tr
# cov1<-cov[cov$trt==0,]
# sdat1<-simsurv(dist = 'weibull',lambdas=0.03,gammas = 3,
#                betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=cov1,
#                maxt=10)
# 
# cov2<-cov[cov$trt==1,]
# sdat2<-simsurv(dist = 'weibull',lambdas=0.2,gammas = 1.5,
#                betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=cov2,
#                maxt=10)
# 
# cov3<-cov[cov$trt==2,]
# sdat3<-simsurv(dist = 'weibull',lambdas=0.1,gammas = 0.9,
#                betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=cov3,
#                maxt=10)
# 
# n1=sum(cov$trt==0);n2=sum(cov$trt==1);n3=sum(cov$trt==2)
# strt<-c(rep(0,n1),rep(1,n2),rep(2,n3))
# sdat<-rbind(sdat1,sdat2,sdat3)
# time<-data.frame(sdat,trt=strt)
# X1=X[Tr==1,];X2=X[Tr==2,];X3=X[Tr==3,]
# 
# X_new=rbind(X1,X2,X3)
# result <- weight.esti(X_new,strt+1)
# 
# par(mfrow=c(2,3))
# akm_rmst(time=time$eventtime, status=time$status, 
#                               group=as.factor(time$trt), weight=result$ps, tau=10)
# title(sub = "weighted:ps")
# 
# akm_rmst(time=time$eventtime, status=time$status, 
#                               group=as.factor(time$trt), weight=result$cbps, tau=10)
# title(sub = "weighted:cbps")
# 
# akm_rmst(time=time$eventtime, status=time$status, 
#                               group=as.factor(time$trt), weight=result$energy, tau=10)
# title(sub = "weighted:energy")
# 
# akm_rmst(time=time$eventtime, status=time$status, 
#                                 group=as.factor(time$trt), weight=result$MB, tau=10)
# title(sub = "weighted:MB")
# 
# akm_rmst(time=time$eventtime, status=time$status, 
#                                 group=as.factor(time$trt),weight=(result$ebal+(result$ebal==0)*1),  tau=10)
# title(sub = "weighted:ebal")
# 
# akm_rmst(time=time$eventtime, status=time$status, 
#                                 group=as.factor(time$trt),  tau=10)
# title(sub = "weighted:no")
# 
# #all ebal weighted are NA

#weight estimate========

theta_non=data.frame( ps_theta_1_2=rep(0,s),ps_theta_1_3=rep(0,s),ps_theta_2_3=rep(0,s),
                      cbps_theta_1_2=rep(0,s),cbps_theta_1_3=rep(0,s),cbps_theta_2_3=rep(0,s),
                      energy_theta_1_2=rep(0,s),energy_theta_1_3=rep(0,s),energy_theta_2_3=rep(0,s),
                      mb_theta_1_2=rep(0,s),mb_theta_1_3=rep(0,s),mb_theta_2_3=rep(0,s),
                      ebal_theta_1_2=rep(0,s),ebal_theta_1_3=rep(0,s),ebal_theta_2_3=rep(0,s),
                      no_theta_1_2=rep(0,s),no_theta_1_3=rep(0,s),no_theta_2_3=rep(0,s),
                      kernel_theta_1_2=rep(0,s),kernel_theta_1_3=rep(0,s),kernel_theta_2_3=rep(0,s),
                      MDABW_theta_1_2=rep(0,s),MDABW_theta_1_3=rep(0,s),MDABW_theta_2_3=rep(0,s)) 

Im_ps_non_1_2=matrix(rep(0,10*s),ncol = 10);Im_ps_non_1_3=matrix(rep(0,10*s),ncol = 10);Im_ps_non_2_3=matrix(rep(0,10*s),ncol = 10)
Im_cbps_non_1_2=matrix(rep(0,10*s),ncol = 10);Im_cbps_non_1_3=matrix(rep(0,10*s),ncol = 10);Im_cbps_non_2_3=matrix(rep(0,10*s),ncol = 10)
Im_energy_non_1_2=matrix(rep(0,10*s),ncol = 10);Im_energy_non_1_3=matrix(rep(0,10*s),ncol = 10);Im_energy_non_2_3=matrix(rep(0,10*s),ncol = 10)
Im_MB_non_1_2=matrix(rep(0,10*s),ncol = 10);Im_MB_non_1_3=matrix(rep(0,10*s),ncol = 10);Im_MB_non_2_3=matrix(rep(0,10*s),ncol = 10)
Im_ebal_non_1_2=matrix(rep(0,10*s),ncol = 10);Im_ebal_non_1_3=matrix(rep(0,10*s),ncol = 10);Im_ebal_non_2_3=matrix(rep(0,10*s),ncol = 10)
Im_no_non_1_2=matrix(rep(0,10*s),ncol = 10);Im_no_non_1_3=matrix(rep(0,10*s),ncol = 10);Im_no_non_2_3=matrix(rep(0,10*s),ncol = 10)

Im_kernel_non_1_2=matrix(rep(0,10*s),ncol = 10);Im_kernel_non_1_3=matrix(rep(0,10*s),ncol = 10);Im_kernel_non_2_3=matrix(rep(0,10*s),ncol = 10)
Im_MDABW_non_1_2=matrix(rep(0,10*s),ncol = 10);Im_MDABW_non_1_3=matrix(rep(0,10*s),ncol = 10);Im_MDABW_non_2_3=matrix(rep(0,10*s),ncol = 10)

Im_GASMD <- array(NA,c(8,s,10))

k=matrix(1:n_total,nrow = h)

for (i in 1:s) {
  X=data$X[k[,i],]
  Tr=data$Tr[k[,i]]
  cov<-data.frame(id=1:nrow(X),trt=(Tr-1),X=X[,c(1,2,3,4,5)])
  #Change the order of X, tr, cov, and survival time according to tr
  set.seed(i+2025)   
  cov1<-cov[cov$trt==0,]
  sdat1<-simsurv(dist = 'weibull',lambdas=0.03,gammas = 3,
                 betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=cov1,
                 maxt=10)
  
  cov2<-cov[cov$trt==1,]
  sdat2<-simsurv(dist = 'weibull',lambdas=0.2,gammas = 1.5,
                 betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=cov2,
                 maxt=10)
  
  cov3<-cov[cov$trt==2,]
  sdat3<-simsurv(dist = 'weibull',lambdas=0.1,gammas = 0.9,
                 betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),x=cov3,
                 maxt=10)
  
  n1=sum(cov$trt==0);n2=sum(cov$trt==1);n3=sum(cov$trt==2)
  strt<-c(rep(0,n1),rep(1,n2),rep(2,n3))
  sdat<-rbind(sdat1,sdat2,sdat3)
  time<-data.frame(sdat,trt=strt)
  X1=X[Tr==1,];X2=X[Tr==2,];X3=X[Tr==3,]
  
  X_new=rbind(X1,X2,X3)
  result <- weight.esti(X_new,strt+1)
  
  theta_non[i,1:3]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                group=as.factor(time$trt), weight=result$ps, tau=10)
  theta_non[i,4:6]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                group=as.factor(time$trt), weight=result$cbps, tau=10)
  theta_non[i,7:9]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                group=as.factor(time$trt), weight=result$energy, tau=10)
  theta_non[i,10:12]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                  group=as.factor(time$trt), weight=result$MB, tau=10)
  theta_non[i,13:15]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                  group=as.factor(time$trt),weight=(result$ebal+min(result$ebal[which(result$ebal>10^(-7))])*(result$ebal==0)),  tau=10)
  theta_non[i,16:18]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                  group=as.factor(time$trt),  tau=10)
  
  theta_non[i,19:21]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                  group=as.factor(time$trt), weight=result$kernel,tau=10)
  theta_non[i,22:24]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                  group=as.factor(time$trt), weight=(result$MDABW+min(result$MDABW[which(result$MDABW>10^(-7))])*(result$MDABW==0)),tau=10)
  
  #calculate imbalance
  sample.size <- dim(X_new)[1]
  dimension   <- dim(X_new)[2]
  data.matrix <- data.frame(X_new , factor(strt+1))
  character   <- names(data.matrix)
  for(j in 1:(dimension+1)){character[j] <- paste(character[j])}
  myformula   <- as.formula(paste(character[1 + dimension],paste(" ~ ", paste(character[1:dimension], collapse= "+"))))
  
  #ASMD
  ps.bal     <- bal.tab(X_new, treat = character[1 + dimension], data = data.matrix, weights = result$ps)
  ebal.bal   <- bal.tab(X_new, treat = character[1 + dimension], data = data.matrix, weights = result$ebal+min(result$ebal[which(result$ebal>10^(-7))])*(result$ebal==0)) 
  cbps.bal   <- bal.tab(X_new, treat = character[1 + dimension], data = data.matrix, weights = result$cbps)
  energy.bal <- bal.tab(X_new, treat = character[1 + dimension], data = data.matrix, weights = result$energy)
  MB.bal     <- bal.tab(X_new, treat = character[1 + dimension], data = data.matrix, weights = result$MB)
  no.bal     <- bal.tab(X_new, treat = character[1 + dimension], data = data.matrix, weights = rep(1,h))
  
  kernel.bal <- bal.tab(X_new, treat = character[1 + dimension], data = data.matrix, weights = result$kernel)
  MDABW.bal <- bal.tab(X_new, treat = character[1 + dimension], data = data.matrix, weights = (result$MDABW+min(result$MDABW[which(result$MDABW>10^(-7))])*(result$MDABW==0)))
  
  Im_ps_non_1_2[i,]=ps.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_ps_non_1_3[i,]=ps.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_ps_non_2_3[i,]=ps.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_ebal_non_1_2[i,]=ebal.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_ebal_non_1_3[i,]=ebal.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_ebal_non_2_3[i,]=ebal.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_cbps_non_1_2[i,]=cbps.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_cbps_non_1_3[i,]=cbps.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_cbps_non_2_3[i,]=cbps.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_energy_non_1_2[i,]=energy.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_energy_non_1_3[i,]=energy.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_energy_non_2_3[i,]=energy.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_MB_non_1_2[i,]=MB.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_MB_non_1_3[i,]=MB.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_MB_non_2_3[i,]=MB.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_no_non_1_2[i,]=no.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_no_non_1_3[i,]=no.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_no_non_2_3[i,]=no.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  
  Im_kernel_non_1_2[i,]=kernel.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_kernel_non_1_3[i,]=kernel.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_kernel_non_2_3[i,]=kernel.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_MDABW_non_1_2[i,]=MDABW.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_MDABW_non_1_3[i,]=MDABW.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_MDABW_non_2_3[i,]=MDABW.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_GASMD[1,i,] <- result$Imbalance$unad
  Im_GASMD[2,i,] <- result$Imbalance$ps
  Im_GASMD[3,i,] <- result$Imbalance$ebal
  Im_GASMD[4,i,] <- result$Imbalance$cbps
  Im_GASMD[5,i,] <- result$Imbalance$energy
  Im_GASMD[6,i,] <- result$Imbalance$MB
  
  Im_GASMD[7,i,] <- result$Imbalance$kernel
  Im_GASMD[8,i,] <- result$Imbalance$MDABW
}

##GASMD==========
Im_GASMD_mean = list(unad = apply(Im_GASMD[1,,],2,mean),
                     ps = apply(Im_GASMD[2,,],2,mean),
                     ebal = apply(Im_GASMD[3,,],2,mean),
                     cbps = apply(Im_GASMD[4,,],2,mean),
                     energy = apply(Im_GASMD[5,,],2,mean),
                     MB = apply(Im_GASMD[6,,],2,mean),
                     kernel=apply(Im_GASMD[7,,],2,mean),
                     MDABW=apply(Im_GASMD[8,,],2,mean)
)
max(Im_GASMD_mean$unad);max(Im_GASMD_mean$ps)
max(Im_GASMD_mean$ebal);max(Im_GASMD_mean$cbps)
max(Im_GASMD_mean$energy);max(Im_GASMD_mean$MB)

max(Im_GASMD_mean$kernel)
max(Im_GASMD_mean$MDABW)
#ASMD==================
#mean
Im_ps_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ps_non_1_2), 2,mean),apply(abs(Im_ps_non_1_3), 2,mean),apply(abs(Im_ps_non_2_3), 2,mean)))
Im_cbps_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_cbps_non_1_2), 2,mean),apply(abs(Im_cbps_non_1_3), 2,mean),apply(abs(Im_cbps_non_2_3), 2,mean)))
Im_energy_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_energy_non_1_2), 2,mean),apply(abs(Im_energy_non_1_3), 2,mean),apply(abs(Im_energy_non_2_3), 2,mean)))
Im_MB_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MB_non_1_2), 2,mean),apply(abs(Im_MB_non_1_3), 2,mean),apply(abs(Im_MB_non_2_3), 2,mean)))
Im_ebal_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ebal_non_1_2), 2,mean),apply(abs(Im_ebal_non_1_3), 2,mean),apply(abs(Im_ebal_non_2_3), 2,mean)))
Im_no_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_no_non_1_2), 2,mean),apply(abs(Im_no_non_1_3), 2,mean),apply(abs(Im_no_non_2_3), 2,mean)))


Im_kernel_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_kernel_non_1_2), 2,mean),apply(abs(Im_kernel_non_1_3), 2,mean),apply(abs(Im_kernel_non_2_3), 2,mean)))
Im_MDABW_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MDABW_non_1_2), 2,mean),apply(abs(Im_MDABW_non_1_3), 2,mean),apply(abs(Im_MDABW_non_2_3), 2,mean)))

Im_non_mean=t(data.frame(ps=t(Im_ps_mean),
                         cbps=t(Im_cbps_mean),
                         energy=t(Im_energy_mean),
                         MB=t(Im_MB_mean),
                         ebal=t(Im_ebal_mean),
                         no=t(Im_no_mean),
                         kernel=t(Im_kernel_mean),
                         MDABW=t(Im_MDABW_mean))) 

Im_non_mean_max=data.frame(
  no=apply(Im_non_mean[16:18,], 2, max),
  ps=apply(Im_non_mean[1:3,], 2, max),
  ebal=apply(Im_non_mean[13:15,], 2, max),
  cbps=apply(Im_non_mean[4:6,], 2, max),
  energy=apply(Im_non_mean[7:9,], 2, max),
  MB=apply(Im_non_mean[10:12,], 2, max),
  kernel=apply(Im_non_mean[19:21,], 2, max),
  MDABW=apply(Im_non_mean[22:24,], 2, max)
)


#sum
Im_ps_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ps_non_1_2), 2,sum),apply(abs(Im_ps_non_1_3), 2,sum),apply(abs(Im_ps_non_2_3), 2,sum)))
Im_cbps_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_cbps_non_1_2), 2,sum),apply(abs(Im_cbps_non_1_3), 2,sum),apply(abs(Im_cbps_non_2_3), 2,sum)))
Im_energy_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_energy_non_1_2), 2,sum),apply(abs(Im_energy_non_1_3), 2,sum),apply(abs(Im_energy_non_2_3), 2,sum)))
Im_MB_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MB_non_1_2), 2,sum),apply(abs(Im_MB_non_1_3), 2,sum),apply(abs(Im_MB_non_2_3), 2,sum)))
Im_ebal_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ebal_non_1_2), 2,sum),apply(abs(Im_ebal_non_1_3), 2,sum),apply(abs(Im_ebal_non_2_3), 2,sum)))
Im_no_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_no_non_1_2), 2,sum),apply(abs(Im_no_non_1_3), 2,sum),apply(abs(Im_no_non_2_3), 2,sum)))


Im_kernel_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_kernel_non_1_2), 2,sum),apply(abs(Im_kernel_non_1_3), 2,sum),apply(abs(Im_kernel_non_2_3), 2,sum)))
Im_MDABW_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MDABW_non_1_2), 2,sum),apply(abs(Im_MDABW_non_1_3), 2,sum),apply(abs(Im_MDABW_non_2_3), 2,sum)))

Im_non_sum=t(data.frame(ps=t(Im_ps_sum),
                        cbps=t(Im_cbps_sum),
                        energy=t(Im_energy_sum),
                        MB=t(Im_MB_sum),
                        ebal=t(Im_ebal_sum),
                        no=t(Im_no_sum),
                        kernel=t(Im_kernel_sum),
                        MDABW=t(Im_MDABW_sum))) 

Im_non_all_sum=data.frame(ps=apply(Im_non_sum[1:3,], 2, sum),
                          cbps=apply(Im_non_sum[4:6,], 2, sum),
                          energy=apply(Im_non_sum[7:9,], 2, sum),
                          MB=apply(Im_non_sum[10:12,], 2, sum),
                          ebal=apply(Im_non_sum[13:15,], 2, sum),
                          no=apply(Im_non_sum[16:18,], 2, sum),
                          kernel=apply(Im_non_sum[19:21,], 2, sum),
                          MDABW=apply(Im_non_sum[22:24,], 2, sum)) 
#max
Im_ps_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ps_non_1_2), 2,max),apply(abs(Im_ps_non_1_3), 2,max),apply(abs(Im_ps_non_2_3), 2,max)))
Im_cbps_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_cbps_non_1_2), 2,max),apply(abs(Im_cbps_non_1_3), 2,max),apply(abs(Im_cbps_non_2_3), 2,max)))
Im_energy_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_energy_non_1_2), 2,max),apply(abs(Im_energy_non_1_3), 2,max),apply(abs(Im_energy_non_2_3), 2,max)))
Im_MB_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MB_non_1_2), 2,max),apply(abs(Im_MB_non_1_3), 2,max),apply(abs(Im_MB_non_2_3), 2,max)))
Im_ebal_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ebal_non_1_2), 2,max),apply(abs(Im_ebal_non_1_3), 2,max),apply(abs(Im_ebal_non_2_3), 2,max)))
Im_no_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_no_non_1_2), 2,max),apply(abs(Im_no_non_1_3), 2,max),apply(abs(Im_no_non_2_3), 2,max)))


Im_kernel_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_kernel_non_1_2), 2,max),apply(abs(Im_kernel_non_1_3), 2,max),apply(abs(Im_kernel_non_2_3), 2,max)))
Im_MDABW_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MDABW_non_1_2), 2,max),apply(abs(Im_MDABW_non_1_3), 2,max),apply(abs(Im_MDABW_non_2_3), 2,max)))

Im_non_max=t(data.frame(ps=t(Im_ps_max),
                        cbps=t(Im_cbps_max),
                        energy=t(Im_energy_max),
                        MB=t(Im_MB_max),
                        ebal=t(Im_ebal_max),
                        no=t(Im_no_max),
                        kernel=t(Im_kernel_max),
                        MDABW=t(Im_MDABW_max)))
#mean============
ps_theta_non_1_2_hat=mean(as.numeric(theta_non$ps_theta_1_2))
cbps_theta_non_1_2_hat=mean(as.numeric(theta_non$cbps_theta_1_2))
energy_theta_non_1_2_hat=mean(as.numeric(theta_non$energy_theta_1_2))
mb_theta_non_1_2_hat=mean(as.numeric(theta_non$mb_theta_1_2))
ebal_theta_non_1_2_hat=mean(as.numeric(theta_non$ebal_theta_1_2))
no_theta_non_1_2_hat=mean(as.numeric(theta_non$no_theta_1_2))

kernel_theta_non_1_2_hat=mean(as.numeric(theta_non$kernel_theta_1_2))
MDABW_theta_non_1_2_hat=mean(as.numeric(theta_non$MDABW_theta_1_2))

ps_theta_non_1_3_hat=mean(as.numeric(theta_non$ps_theta_1_3))
cbps_theta_non_1_3_hat=mean(as.numeric(theta_non$cbps_theta_1_3))
energy_theta_non_1_3_hat=mean(as.numeric(theta_non$energy_theta_1_3))
mb_theta_non_1_3_hat=mean(as.numeric(theta_non$mb_theta_1_3))
ebal_theta_non_1_3_hat=mean(as.numeric(theta_non$ebal_theta_1_3))
no_theta_non_1_3_hat=mean(as.numeric(theta_non$no_theta_1_3))

kernel_theta_non_1_3_hat=mean(as.numeric(theta_non$kernel_theta_1_3))
MDABW_theta_non_1_3_hat=mean(as.numeric(theta_non$MDABW_theta_1_3))

ps_theta_non_2_3_hat=mean(as.numeric(theta_non$ps_theta_2_3))
cbps_theta_non_2_3_hat=mean(as.numeric(theta_non$cbps_theta_2_3))
energy_theta_non_2_3_hat=mean(as.numeric(theta_non$energy_theta_2_3))
mb_theta_non_2_3_hat=mean(as.numeric(theta_non$mb_theta_2_3))
ebal_theta_non_2_3_hat=mean(as.numeric(theta_non$ebal_theta_2_3))
no_theta_non_2_3_hat=mean(as.numeric(theta_non$no_theta_2_3))

kernel_theta_non_2_3_hat=mean(as.numeric(theta_non$kernel_theta_2_3))
MDABW_theta_non_2_3_hat=mean(as.numeric(theta_non$MDABW_theta_2_3))

non_ph_estimate=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"), 
                           true=c(theta_true_non[1,1],theta_true_non[1,2],theta_true_non[1,3]),
                           ps=c(ps_theta_non_1_2_hat,ps_theta_non_1_3_hat,ps_theta_non_2_3_hat),
                           cbps=c(cbps_theta_non_1_2_hat,cbps_theta_non_1_3_hat,cbps_theta_non_2_3_hat),
                           energy=c(energy_theta_non_1_2_hat,energy_theta_non_1_3_hat,energy_theta_non_2_3_hat),
                           mb=c(mb_theta_non_1_2_hat,mb_theta_non_1_3_hat,mb_theta_non_2_3_hat),
                           ebal=c(ebal_theta_non_1_2_hat,ebal_theta_non_1_3_hat,ebal_theta_non_2_3_hat),
                           no=c(no_theta_non_1_2_hat,no_theta_non_1_3_hat,no_theta_non_2_3_hat),
                           kernel=c(kernel_theta_non_1_2_hat,kernel_theta_non_1_3_hat,kernel_theta_non_2_3_hat),
                           MDABW=c(MDABW_theta_non_1_2_hat,MDABW_theta_non_1_3_hat,MDABW_theta_non_2_3_hat))

imbalance_of_non_ph=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"))



#bias==================
n_bias=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),
                  ps=c(((non_ph_estimate[1,2]-non_ph_estimate[1,1])/non_ph_estimate[1,1]),((non_ph_estimate[2,2]-non_ph_estimate[2,1])/non_ph_estimate[2,1]),((non_ph_estimate[3,2]-non_ph_estimate[3,1])/non_ph_estimate[3,1])),
                  cbps=c(((non_ph_estimate[1,3]-non_ph_estimate[1,1])/non_ph_estimate[1,1]),((non_ph_estimate[2,3]-non_ph_estimate[2,1])/non_ph_estimate[2,1]),((non_ph_estimate[3,3]-non_ph_estimate[3,1])/non_ph_estimate[3,1])),
                  energy=c(((non_ph_estimate[1,4]-non_ph_estimate[1,1])/non_ph_estimate[1,1]),((non_ph_estimate[2,4]-non_ph_estimate[2,1])/non_ph_estimate[2,1]),((non_ph_estimate[3,4]-non_ph_estimate[3,1])/non_ph_estimate[3,1])),
                  mb=c(((non_ph_estimate[1,5]-non_ph_estimate[1,1])/non_ph_estimate[1,1]),((non_ph_estimate[2,5]-non_ph_estimate[2,1])/non_ph_estimate[2,1]),((non_ph_estimate[3,5]-non_ph_estimate[3,1])/non_ph_estimate[3,1])),
                  ebal=c(((non_ph_estimate[1,6]-non_ph_estimate[1,1])/non_ph_estimate[1,1]),((non_ph_estimate[2,6]-non_ph_estimate[2,1])/non_ph_estimate[2,1]),((non_ph_estimate[3,6]-non_ph_estimate[3,1])/non_ph_estimate[3,1])),
                  no=c(((non_ph_estimate[1,7]-non_ph_estimate[1,1])/non_ph_estimate[1,1]),((non_ph_estimate[2,7]-non_ph_estimate[2,1])/non_ph_estimate[2,1]),((non_ph_estimate[3,7]-non_ph_estimate[3,1])/non_ph_estimate[3,1])),
                  kernel=c(((non_ph_estimate[1,8]-non_ph_estimate[1,1])/non_ph_estimate[1,1]),((non_ph_estimate[2,8]-non_ph_estimate[2,1])/non_ph_estimate[2,1]),((non_ph_estimate[3,8]-non_ph_estimate[3,1])/non_ph_estimate[3,1])),
                  MDABW=c(((non_ph_estimate[1,9]-non_ph_estimate[1,1])/non_ph_estimate[1,1]),((non_ph_estimate[2,9]-non_ph_estimate[2,1])/non_ph_estimate[2,1]),((non_ph_estimate[3,9]-non_ph_estimate[3,1])/non_ph_estimate[3,1])))

n_bias_mean=apply(abs(n_bias), 2, mean)
#mse======
n_mse_ps_1_2<-sum((as.numeric(theta_non$ps_theta_1_2)-as.numeric(theta_true_non[1]))^2)/s
n_mse_cbps_1_2<-sum((as.numeric(theta_non$cbps_theta_1_2)-as.numeric(theta_true_non[1]))^2)/s
n_mse_energy_1_2<-sum((as.numeric(theta_non$energy_theta_1_2)-as.numeric(theta_true_non[1]))^2)/s
n_mse_mb_1_2<-sum((as.numeric(theta_non$mb_theta_1_2)-as.numeric(theta_true_non[1]))^2)/s
n_mse_ebal_1_2<-sum((as.numeric(theta_non$ebal_theta_1_2)-as.numeric(theta_true_non[1]))^2)/s
n_mse_no_1_2<-sum((as.numeric(theta_non$no_theta_1_2)-as.numeric(theta_true_non[1]))^2)/s

n_mse_kernel_1_2<-sum((as.numeric(theta_non$kernel_theta_1_2)-as.numeric(theta_true_non[1]))^2)/s
n_mse_MDABW_1_2<-sum((as.numeric(theta_non$MDABW_theta_1_2)-as.numeric(theta_true_non[1]))^2)/s

n_mse_ps_1_3<-sum((as.numeric(theta_non$ps_theta_1_3)-as.numeric(theta_true_non[2]))^2)/s
n_mse_cbps_1_3<-sum((as.numeric(theta_non$cbps_theta_1_3)-as.numeric(theta_true_non[2]))^2)/s
n_mse_energy_1_3<-sum((as.numeric(theta_non$energy_theta_1_3)-as.numeric(theta_true_non[2]))^2)/s
n_mse_mb_1_3<-sum((as.numeric(theta_non$mb_theta_1_3)-as.numeric(theta_true_non[2]))^2)/s
n_mse_ebal_1_3<-sum((as.numeric(theta_non$ebal_theta_1_3)-as.numeric(theta_true_non[2]))^2)/s
n_mse_no_1_3<-sum((as.numeric(theta_non$no_theta_1_3)-as.numeric(theta_true_non[2]))^2)/s

n_mse_kernel_1_3<-sum((as.numeric(theta_non$kernel_theta_1_3)-as.numeric(theta_true_non[2]))^2)/s
n_mse_MDABW_1_3<-sum((as.numeric(theta_non$MDABW_theta_1_3)-as.numeric(theta_true_non[2]))^2)/s

n_mse_ps_2_3<-sum((as.numeric(theta_non$ps_theta_2_3)-as.numeric(theta_true_non[3]))^2)/s
n_mse_cbps_2_3<-sum((as.numeric(theta_non$cbps_theta_2_3)-as.numeric(theta_true_non[3]))^2)/s
n_mse_energy_2_3<-sum((as.numeric(theta_non$energy_theta_2_3)-as.numeric(theta_true_non[3]))^2)/s
n_mse_mb_2_3<-sum((as.numeric(theta_non$mb_theta_2_3)-as.numeric(theta_true_non[3]))^2)/s
n_mse_ebal_2_3<-sum((as.numeric(theta_non$ebal_theta_2_3)-as.numeric(theta_true_non[3]))^2)/s
n_mse_no_2_3<-sum((as.numeric(theta_non$no_theta_2_3)-as.numeric(theta_true_non[3]))^2)/s

n_mse_kernel_2_3<-sum((as.numeric(theta_non$kernel_theta_2_3)-as.numeric(theta_true_non[3]))^2)/s
n_mse_MDABW_2_3<-sum((as.numeric(theta_non$MDABW_theta_2_3)-as.numeric(theta_true_non[3]))^2)/s

n_mse=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),
                 ps=c(n_mse_ps_1_2,n_mse_ps_1_3,n_mse_ps_2_3),
                 cbps=c(n_mse_cbps_1_2,n_mse_cbps_1_3,n_mse_cbps_2_3),
                 energy=c(n_mse_energy_1_2,n_mse_energy_1_3,n_mse_energy_2_3),
                 mb=c(n_mse_mb_1_2,n_mse_mb_1_3,n_mse_mb_2_3),
                 ebal=c(n_mse_ebal_1_2,n_mse_ebal_1_3,n_mse_ebal_2_3),
                 no=c(n_mse_no_1_2,n_mse_no_1_3,n_mse_no_2_3),
                 kernel=c(n_mse_kernel_1_2,n_mse_kernel_1_3,n_mse_kernel_2_3),
                 MDABW=c(n_mse_MDABW_1_2,n_mse_MDABW_1_3,n_mse_MDABW_2_3))

n_mse_mean=apply(n_mse, 2, mean)


n_sam_ps_1_2<-sum((as.numeric(theta_non$ps_theta_1_2)-mean(as.numeric(theta_non$ps_theta_1_2)))^2)/s
n_sam_cbps_1_2<-sum((as.numeric(theta_non$cbps_theta_1_2)-mean(as.numeric(theta_non$cbps_theta_1_2)))^2)/s
n_sam_energy_1_2<-sum((as.numeric(theta_non$energy_theta_1_2)-mean(as.numeric(theta_non$energy_theta_1_2)))^2)/s
n_sam_mb_1_2<-sum((as.numeric(theta_non$mb_theta_1_2)-mean(as.numeric(theta_non$mb_theta_1_2)))^2)/s
n_sam_ebal_1_2<-sum((as.numeric(theta_non$ebal_theta_1_2)-mean(as.numeric(theta_non$ebal_theta_1_2)))^2)/s
n_sam_no_1_2<-sum((as.numeric(theta_non$no_theta_1_2)-mean(as.numeric(theta_non$no_theta_1_2)))^2)/s

n_sam_kernel_1_2<-sum((as.numeric(theta_non$kernel_theta_1_2)-mean(as.numeric(theta_non$kernel_theta_1_2)))^2)/s
n_sam_MDABW_1_2<-sum((as.numeric(theta_non$MDABW_theta_1_2)-mean(as.numeric(theta_non$MDABW_theta_1_2)))^2)/s

n_sam_ps_1_3<-sum((as.numeric(theta_non$ps_theta_1_3)-mean(as.numeric(theta_non$ps_theta_1_3)))^2)/s
n_sam_cbps_1_3<-sum((as.numeric(theta_non$cbps_theta_1_3)-mean(as.numeric(theta_non$cbps_theta_1_3)))^2)/s
n_sam_energy_1_3<-sum((as.numeric(theta_non$energy_theta_1_3)-mean(as.numeric(theta_non$energy_theta_1_3)))^2)/s
n_sam_mb_1_3<-sum((as.numeric(theta_non$mb_theta_1_3)-mean(as.numeric(theta_non$mb_theta_1_3)))^2)/s
n_sam_ebal_1_3<-sum((as.numeric(theta_non$ebal_theta_1_3)-mean(as.numeric(theta_non$ebal_theta_1_3)))^2)/s
n_sam_no_1_3<-sum((as.numeric(theta_non$no_theta_1_3)-mean(as.numeric(theta_non$no_theta_1_3)))^2)/s
 
n_sam_kernel_1_3<-sum((as.numeric(theta_non$kernel_theta_1_3)-mean(as.numeric(theta_non$kernel_theta_1_3)))^2)/s
n_sam_MDABW_1_3<-sum((as.numeric(theta_non$MDABW_theta_1_3)-mean(as.numeric(theta_non$MDABW_theta_1_3)))^2)/s

n_sam_ps_2_3<-sum((as.numeric(theta_non$ps_theta_2_3)-mean(as.numeric(theta_non$ps_theta_2_3)))^2)/s
n_sam_cbps_2_3<-sum((as.numeric(theta_non$cbps_theta_2_3)-mean(as.numeric(theta_non$cbps_theta_2_3)))^2)/s
n_sam_energy_2_3<-sum((as.numeric(theta_non$energy_theta_2_3)-mean(as.numeric(theta_non$energy_theta_2_3)))^2)/s
n_sam_mb_2_3<-sum((as.numeric(theta_non$mb_theta_2_3)-mean(as.numeric(theta_non$mb_theta_2_3)))^2)/s
n_sam_ebal_2_3<-sum((as.numeric(theta_non$ebal_theta_2_3)-mean(as.numeric(theta_non$ebal_theta_2_3)))^2)/s
n_sam_no_2_3<-sum((as.numeric(theta_non$no_theta_2_3)-mean(as.numeric(theta_non$no_theta_2_3)))^2)/s

n_sam_kernel_2_3<-sum((as.numeric(theta_non$kernel_theta_2_3)-mean(as.numeric(theta_non$kernel_theta_2_3)))^2)/s
n_sam_MDABW_2_3<-sum((as.numeric(theta_non$MDABW_theta_2_3)-mean(as.numeric(theta_non$MDABW_theta_2_3)))^2)/s

n_sam=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),
                 ps=c(n_sam_ps_1_2,n_sam_ps_1_3,n_sam_ps_2_3),
                 cbps=c(n_sam_cbps_1_2,n_sam_cbps_1_3,n_sam_cbps_2_3),
                 energy=c(n_sam_energy_1_2,n_sam_energy_1_3,n_sam_energy_2_3),
                 mb=c(n_sam_mb_1_2,n_sam_mb_1_3,n_sam_mb_2_3),
                 ebal=c(n_sam_ebal_1_2,n_sam_ebal_1_3,n_sam_ebal_2_3),
                 no=c(n_sam_no_1_2,n_sam_no_1_3,n_sam_no_2_3),
                 kernel=c(n_sam_kernel_1_2,n_sam_kernel_1_3,n_sam_kernel_2_3),
                 MDABW=c(n_sam_MDABW_1_2,n_sam_MDABW_1_3,n_sam_MDABW_2_3))

n_sam[1,]+(n_bias[1,]*as.numeric(theta_true_non[1]))^2
n_mse[1,]

#write table======
###### under a specific method, for a specific covariate, the maximum mean(average of 1000 simulations) asmd between different pair(1vs2,1vs3,2vs3)
wASMD=data.frame(
  Unad=Im_non_mean_max$no,IPW=Im_non_mean_max$ps,
  EBAL=Im_non_mean_max$ebal,CBPS=Im_non_mean_max$cbps,
  Energy=Im_non_mean_max$energy,MB=Im_non_mean_max$MB,kernel=Im_non_mean_max$kernel,MDABW=Im_non_mean_max$MDABW)

rownames(wASMD)=c(paste("Max-ASMD",c(1:10)))

wSD=data.frame(Unad=mean(sqrt(n_sam[,6])),IPW=mean(sqrt(n_sam[,1])),
               EBAL=mean(sqrt(n_sam[,5])),CBPS=mean(sqrt(n_sam[,2])),
               Energy=mean(sqrt(n_sam[,3])),MB=mean(sqrt(n_sam[,4])),
               kernel=mean(sqrt(n_sam[,7])),MDABW=mean(sqrt(n_sam[,8])))
rownames(wSD)=c("SD-RMST")

wRMSE=data.frame(Unad=sqrt(n_mse_mean[6]),IPW=sqrt(n_mse_mean[1]),
                 EBAL=sqrt(n_mse_mean[5]),CBPS=sqrt(n_mse_mean[2]),
                 Energy=sqrt(n_mse_mean[3]),MB=sqrt(n_mse_mean[4]),
                 kernel=sqrt(n_mse_mean[7]),MDABW=sqrt(n_mse_mean[8]))
rownames(wRMSE)=c("RMSE-RMST")

wBias=data.frame(Unad=n_bias_mean[6],IPW=n_bias_mean[1],
                 EBAL=n_bias_mean[5],CBPS=n_bias_mean[2],
                 Energy=n_bias_mean[3],MB=n_bias_mean[4],kernel= n_bias_mean[7], MDABW= n_bias_mean[8])
rownames(wBias)=c("Relative Bias-RMST")

###### sumgasmd independent of treatment groups, reflecting overall covariate balance of a specific covariate
###### under a specific method, maximum of mean(over 1000 simulations) of sumgasmd of different covariates

wGASMD=data.frame(Unad=max(Im_GASMD_mean$unad),IPW=max(Im_GASMD_mean$ps),
                  EBAL=max(Im_GASMD_mean$ebal),CBPS=max(Im_GASMD_mean$cbps),
                  Energy=max(Im_GASMD_mean$energy),MB=max(Im_GASMD_mean$MB),kernel=max(Im_GASMD_mean$kernel),MDABW=max(Im_GASMD_mean$MDABW))
rownames(wGASMD)=c("Max-GASMD")

all_table=rbind(wBias,wRMSE,wSD,wASMD,wGASMD)
all_table=round(all_table,digits = 3)


write.csv(file = 'Good_Non_result.csv',x=all_table)

#appendix
#save bias\sd\rmse for all pairwise comparison rather than mean of bias\sd\rmse for all pairwise comparison
#save RMST estimate
write.csv(file = 'Good_NonPH_rmst_est.csv',x=non_ph_estimate)
#save bias
write.csv(file = 'Good_NonPH_bias.csv',x=n_bias)
#save sd
write.csv(file = 'Good_NonPH_sd.csv',x=n_sam)
#save mse
write.csv(file = 'Good_NonPH_mse.csv',x=n_mse)
#save GASMD for all covariates
Im_GASMD_mean_df<-data.frame(Im_GASMD_mean)
write.csv(file = 'Good_NonPH_gasmd.csv',x=Im_GASMD_mean_df)
