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
library(dplyr)

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
  
  while(sum(dim(as.matrix(sbw.result))) != 15){
    sbw.result  <- tryCatch(sbw.result <- sbw(dat = data.matrix, ind = character[1 + dimension], bal = bal, out = character[2 + dimension], par = list(par_est = "ate")), error = function(e) { skip_to_next <<- FALSE})
    bal$bal_gri <- bal$bal_gri[-1]
  }
  sbw.weight <- sbw.result$dat_weights$sbw_weights
  return(list(weight = sbw.weight))
}




####################



weight.esti <- function(X, Tr, delta.space = c(1,1e-2,1e-4)){
  options(error = function() traceback(2))
  data.matrix <- data.frame(X , factor(Tr))
  sample.size <- dim(X)[1]
  dimension   <- dim(X)[2]
  
  character   <- names(data.matrix)
  for(j in 1:(dimension+1)){character[j] <- paste(character[j])}
  myformula   <- as.formula(paste(character[1 + dimension],paste(" ~ ", paste(character[1:dimension], collapse= "+"))))
  
  ps.weight   <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "ps")$weights
  ebal.weight <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "ebal",method.options = list(maxit = 100000000))$weights
  cbps.weight <- weightit(myformula , data = data.matrix, estimand = "ATE", method = "cbps", over = FALSE,method.options = list(maxit = 100000000))$weights
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



#generate data.
#X denotes observed covariate.
#Tr denotes treatment indicator, where 1,2,3 represent group1,group2,group3.
#Propensity score model:
#Pr(Tr = 1 \mid X) = 1 / (1 + \exp(\beta^{\top} X[1:5]) + \exp(\beta^{\top} X[6:10])).
#Pr(Tr = 2 \mid X) = \exp(\beta^{\top} X[1:5]) / (1 + \exp(\beta^{\top} X[1:5]) + \exp(\beta^{\top} X[6:10])).
#Pr(Tr = 3 \mid X) = \exp(\beta^{\top} X[6:10]) / (1 + \exp(\beta^{\top} X[1:5]) + \exp(\beta^{\top} X[6:10])).
#\beta = (1,1,1,1,1)&{\top}


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
          #lines(ft$tj, ft$st,type="s", col=(i+2), lwd=2)
          #--- Plot AKM (extend to xaxismax) ---
          if (max(ft$tj) < xaxismax) {
            ft <- rbind(
              ft,
              data.frame(
                tj = xaxismax,
                yj = NA,
                dj = NA,
                st = ft$st[length(ft$st)],
                i = i,
                mj = NA
              )
            )
          }
          
          lines(ft$tj, ft$st, type = "s", col = (i+2), lwd = 2)
          
        }
      }
    }
  }
  
  #--- Add legend and tau to plot ---
  #abline(v=tau, col=1, lty=3, lwd=2)
  #legend( "bottomleft",paste("Treatment", groupval), lty=rep(1, j), lwd=rep(2, j), col=3:(j+2), 
  #       cex=1, bty ="n", inset = c(0, 0),seg.len = 0.2,x.intersp = 0.5)
  
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

#source('get_diff_se_ci.R')
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

library(simsurv)
library(survival)
library(survminer)
library(ggplot2)
library(reshape2)
library(cobalt)



#data process====
library(nnet)
library(mice)
View(colon)
#md.pattern(colon)
notna.colon<-na.omit(colon)

notna.colon<-notna.colon[which(notna.colon$etype==1),]

notna.colon$age<-log(notna.colon$age)

trt<-notna.colon$rx#?õ?treatment
extent_cate_dum<-class.ind(notna.colon$extent)
extent_cate<-data.frame('extent2'=extent_cate_dum[,2],'extent3'=extent_cate_dum[,3],'extent4'=extent_cate_dum[,4])

differ_cate_dum<-class.ind(notna.colon$differ)
differ_cate<-data.frame('differ2'=differ_cate_dum[,2],'differ3'=differ_cate_dum[,3]);

notna.colon$nodes<-log(notna.colon$nodes+1)

dat<-data.frame(time=notna.colon$time,status=notna.colon$status,rx=notna.colon$rx,notna.colon[,c(4:9,13)],
                extent_cate,differ_cate,age_square=notna.colon$age^2,nodes_square=notna.colon$nodes^2)
dat <- dat %>% mutate(
  # nodes_cross
  nodes_age = nodes * age,
  nodes_sex = nodes * sex,
  nodes_obstruct = nodes * obstruct,
  nodes_perfor = nodes * perfor,
  nodes_adhere = nodes * adhere,
  nodes_surg = nodes * surg,
  nodes_extent2 = nodes * extent2,
  nodes_extent3 = nodes * extent3,
  nodes_extent4 = nodes * extent4,
  nodes_differ2 = nodes * differ2,
  nodes_differ3 = nodes * differ3,
  
  
  
  extent4_sex       = extent4 * sex,
  extent4_age       = extent4 * age,
  extent4_obstruct  = extent4 * obstruct,
  
  extent4_adhere    = extent4 * adhere,
  extent4_nodes     = extent4 * nodes,
  extent4_surg      = extent4 * surg,
  #extent4_differ2 = extent4 * differ2,
  #extent4_differ3 = extent4 * differ3,
  
  surg_sex        = surg * sex,
  surg_age        = surg * age,
  surg_obstruct   = surg * obstruct,
  surg_perfor     = surg * perfor,
  surg_adhere     = surg * adhere,
  surg_extent2    = surg * extent2,
  surg_extent3    = surg * extent3,
  
  #surg_differ2    = surg * differ2,
  #surg_differ3    = surg * differ3

  #sex_nodes_extent4 = sex * nodes * extent4,
  #sex_nodes_surg    = sex * nodes * surg,
  
  
)


#obs-1,lev-2,lev+5fu-3
levels(dat$rx)<-c(1,2,3)
dat$rx<-as.numeric(dat$rx)
#dat$rx[dat$rx=='Lev+5FU']<-1
#$rx[dat$rx=='Obs']<-2
#dat$rx[dat$rx=='Lev']<-3
###????log?????ᱨ??
#dat$nodes<-log(dat$nodes)###log????
#weighted=====
plot_dat<-dat
plot_dat$rx<-factor(plot_dat$rx)
chara1<-colnames(dat)
ncolume1<-ncol(dat)
notna_formula   <- as.formula(paste(chara1[3],
                                    paste(" ~ ", paste(chara1[4:ncolume1], 
                                                       collapse= "+"))))
bal.tab(notna_formula, data = plot_dat, estimand = "ATE",
        stats = c("mean.diffs"))
wei_notna<-weight.esti(dat[,-c(1:3)],dat$rx)

#warning:The optimization failed to converge. See Notes section at ?method_energy for information.
png("D:/covariates balancing/MB_new_latex/colon/crossmore/AKM_colon_cross_etype1.png",width = 10, height = 6, units = "in", res = 1000)

layout_matrix <- rbind(
  c(1,1,1,1),
  c(2,3,4,5),
  c(6,7,8,9)
)
layout(layout_matrix, heights = c(0.6, 3.5, 3.5))

par(mar = c(0,0,0,0))
plot.new()
legend("center",
       legend = c("Obs","Lev","Lev+5Fu"),
       col = c(3,4,5),
       lty = 1,
       lwd = 2,
       horiz = TRUE,
       cex = 1.4,
       y.intersp = 0.7,
       bty = "n")

par(mar = c(6,4,2,1))

akm_rmst(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
         group =as.factor(dat$rx))

title(main = "UNAD")

akm_rmst(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
         group =as.factor(dat$rx),weight = wei_notna$ps)

title(main = "IPW")


akm_rmst(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
         group =as.factor(dat$rx),weight = wei_notna$cbps)

title(main = "CBPS")

akm_rmst(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
         group =as.factor(dat$rx),weight = (wei_notna$ebal+min(wei_notna$ebal[which(wei_notna$ebal>10^(-7))])*(wei_notna$ebal<10^(-7))))

title(main = "EBAL")



akm_rmst(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
         group =as.factor(dat$rx),weight = wei_notna$energy)

title(main = "ENERGY")


akm_rmst(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
         group =as.factor(dat$rx),weight = wei_notna$kernel)

title(main = "KERNEL")

akm_rmst(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
         group =as.factor(dat$rx),weight = (wei_notna$MDABW+min(wei_notna$MDABW[which(wei_notna$MDABW>10^(-7))])*(wei_notna$MDABW<10^(-7))))

title(main = "MDABW")

akm_rmst(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
         group =as.factor(dat$rx),weight = (wei_notna$MB+min(wei_notna$MB[which(wei_notna$MB>10^(-7))])*(wei_notna$MB<10^(-7))))

title(main = "MB")

dev.off()


se_ci_unad=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,tau=2725,
                    group =as.factor(dat$rx))
se_ci_ebal=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,tau=2725,
                    group =as.factor(dat$rx),weight = (wei_notna$ebal+min(wei_notna$ebal[which(wei_notna$ebal>10^(-7))])*(wei_notna$ebal<10^(-7))))
se_ci_mb=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,tau=2725,
                  group =as.factor(dat$rx),weight = (wei_notna$MB+min(wei_notna$MB[which(wei_notna$MB>10^(-7))])*(wei_notna$MB<10^(-7))))
se_ci_energy=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,tau=2725,
                      group =as.factor(dat$rx),weight = wei_notna$energy)
se_ci_cbps=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,tau=2725,
                    group =as.factor(dat$rx),weight = wei_notna$cbps)
se_ci_ps=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,tau=2725,
                  group =as.factor(dat$rx),weight = wei_notna$ps)
se_ci_kernel=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,tau=2725,
                      group =as.factor(dat$rx),weight = wei_notna$kernel)
se_ci_MDABW=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,tau=2725,
                     group =as.factor(dat$rx),weight = (wei_notna$MDABW+min(wei_notna$MDABW[which(wei_notna$MDABW>10^(-7))])*(wei_notna$MDABW<10^(-7))))

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
se_ci_imp<-cbind(weight=rep(c('Unad','IPW','EBAL','CBPS','Energy','Kernel','MDABW','MB'),3),se_ci_imp,
                 CI=paste(
                   paste(
                     paste('(',se_ci_imp$rmst_diff_low),se_ci_imp$rmst_diff_upp,sep = ','),')'))
write.csv(file = 'se_ci_colon_cross_etype1.csv',x=se_ci_imp)

#trend with tau plot
#get value====
ps_tau=data.frame(vs12=rep(0,2725),vs13=rep(0,2725),vs23=rep(0,2725))
cbps_tau=data.frame(vs12=rep(0,2725),vs13=rep(0,2725),vs23=rep(0,2725))
energy_tau=data.frame(vs12=rep(0,2725),vs13=rep(0,2725),vs23=rep(0,2725))
mb_tau=data.frame(vs12=rep(0,2725),vs13=rep(0,2725),vs23=rep(0,2725))
ebal_tau=data.frame(vs12=rep(0,2725),vs13=rep(0,2725),vs23=rep(0,2725))
no_tau=data.frame(vs12=rep(0,2725),vs13=rep(0,2725),vs23=rep(0,2725))
kernel_tau=data.frame(vs12=rep(0,2725),vs13=rep(0,2725),vs23=rep(0,2725))
MDABW_tau=data.frame(vs12=rep(0,2725),vs13=rep(0,2725),vs23=rep(0,2725))
#plot x:tau; y: RMST difference
for (i in 1:2725) {
  ps_tau[i,]=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
                      group =as.factor(dat$rx),weight = wei_notna$ps,tau=i)[,2]
  
  
  cbps_tau[i,]=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
                        group =as.factor(dat$rx),weight = wei_notna$cbps,tau=i)[,2]
  
  
  
  energy_tau[i,]=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
                          group =as.factor(dat$rx),weight = wei_notna$energy,tau=i)[,2]
  
  
  mb_tau[i,]=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
                      group =as.factor(dat$rx),weight =(wei_notna$MB+min(wei_notna$MB[which(wei_notna$MB>10^(-7))])*(wei_notna$MB<10^(-7))),tau=i)[,2]
  
  ebal_tau[i,]=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
                        group =as.factor(dat$rx),weight = (wei_notna$ebal+min(wei_notna$ebal[which(wei_notna$ebal>10^(-7))])*(wei_notna$ebal<10^(-7))),tau=i)[,2]
  
  no_tau[i,]=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
                      group =as.factor(dat$rx),tau=i)[,2]
  
  kernel_tau[i,]=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
                          group =as.factor(dat$rx),weight = wei_notna$kernel,tau=i)[,2]
  
  MDABW_tau[i,]=get_diff(time=as.numeric(dat[,1]),status = dat[,2],xaxismax=2725,
                         group =as.factor(dat$rx),weight = (wei_notna$MDABW+min(wei_notna$MDABW[which(wei_notna$MDABW>10^(-7))])*(wei_notna$MDABW<10^(-7))),tau=i)[,2]
  
}
ps_tau1=cbind(tau=1:2725,ps_tau)
cbps_tau1=cbind(tau=1:2725,cbps_tau)
energy_tau1=cbind(tau=1:2725,energy_tau)
mb_tau1=cbind(tau=1:2725,mb_tau)
ebal_tau1=cbind(tau=1:2725,ebal_tau)
no_tau1=cbind(tau=1:2725,no_tau)
kernel_tau1=cbind(tau=1:2725,kernel_tau)
MDABW_tau1=cbind(tau=1:2725,MDABW_tau)

#plot======
library(ggplot2)
library(reshape2)
library(ggpubr)
p1 <- ggplot() +
  geom_line(data = ps_tau1, aes(x = tau, y = vs12, colour = "IPW"), size = 1) +
  geom_line(data = cbps_tau1, aes(x = tau, y = vs12, colour = "CBPS"), size = 1) +
  geom_line(data = energy_tau1, aes(x = tau, y = vs12, colour = "ENERGY"), size = 1) +
  geom_line(data = mb_tau1, aes(x = tau, y = vs12, colour = "MB"), size = 1) +
  geom_line(data = ebal_tau1, aes(x = tau, y = vs12, colour = "EBAL"), size = 1) +
  geom_line(data = no_tau1, aes(x = tau, y = vs12, colour = "UNAD"), size = 1) +
  geom_line(data = kernel_tau1, aes(x = tau, y = vs12, colour = "KERNEL"), size = 1) +
  geom_line(data = MDABW_tau1, aes(x = tau, y = vs12, colour = "MDABW"), size = 1) +
  scale_y_continuous(limits = c(-20, 100), breaks = seq(-20, 100, 20)) +
  scale_x_continuous(limits = c(0, 2725), breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  ggtitle("Lev vs. Obs") +
  scale_color_manual(
    name = "Method",
    values = c(
      'UNAD'='skyblue','MB'='blue','IPW'='black','EBAL'='orange',
      'CBPS'='red','ENERGY'='green','KERNEL'='#D55E00','MDABW'='#CC79A7'
    ),
    breaks = c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


p2 <- ggplot() +
  geom_line(data = ps_tau1, aes(x = tau, y = vs13, colour = "IPW"), size = 1) +
  geom_line(data = cbps_tau1, aes(x = tau, y = vs13, colour = "CBPS"), size = 1) +
  geom_line(data = energy_tau1, aes(x = tau, y = vs13, colour = "ENERGY"), size = 1) +
  geom_line(data = mb_tau1, aes(x = tau, y = vs13, colour = "MB"), size = 1) +
  geom_line(data = ebal_tau1, aes(x = tau, y = vs13, colour = "EBAL"), size = 1) +
  geom_line(data = no_tau1, aes(x = tau, y = vs13, colour = "UNAD"), size = 1) +
  geom_line(data = kernel_tau1, aes(x = tau, y = vs13, colour = "KERNEL"), size = 1) +
  geom_line(data = MDABW_tau1, aes(x = tau, y = vs13, colour = "MDABW"), size = 1) +
  scale_y_continuous(limits = c(-5, 420), breaks = seq(0, 420, 20)) +
  scale_x_continuous(limits = c(0, 2725), breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  ggtitle("Lev+5Fu vs. Obs") +
  scale_color_manual(
    name = "Method",
    values = c(
      'UNAD'='skyblue','MB'='blue','IPW'='black','EBAL'='orange',
      'CBPS'='red','ENERGY'='green','KERNEL'='#D55E00','MDABW'='#CC79A7'
    ),
    breaks = c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


p3 <- ggplot() +
  geom_line(data = ps_tau1, aes(x = tau, y = vs23, colour = "IPW"), size = 1) +
  geom_line(data = cbps_tau1, aes(x = tau, y = vs23, colour = "CBPS"), size = 1) +
  geom_line(data = energy_tau1, aes(x = tau, y = vs23, colour = "ENERGY"), size = 1) +
  geom_line(data = mb_tau1, aes(x = tau, y = vs23, colour = "MB"), size = 1) +
  geom_line(data = ebal_tau1, aes(x = tau, y = vs23, colour = "EBAL"), size = 1) +
  geom_line(data = no_tau1, aes(x = tau, y = vs23, colour = "UNAD"), size = 1) +
  geom_line(data = kernel_tau1, aes(x = tau, y = vs23, colour = "KERNEL"), size = 1) +
  geom_line(data = MDABW_tau1, aes(x = tau, y = vs23, colour = "MDABW"), size = 1) +
  scale_y_continuous(limits = c(-5, 420), breaks = seq(0, 420, 20)) +
  scale_x_continuous(limits = c(0, 2725), breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  ggtitle("Lev+5Fu vs. Lev") +
  scale_color_manual(
    name = "Method",
    values = c(
      'UNAD'='skyblue','MB'='blue','IPW'='black','EBAL'='orange',
      'CBPS'='red','ENERGY'='green','KERNEL'='#D55E00','MDABW'='#CC79A7'
    ),
    breaks = c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Time (days)", y = "Restricted Mean Survival Time Difference")


png("rmst_colon_cross_etype1.png",
    width = 8, height = 5.6, units = "in", res = 1000)

ggarrange(
  p1, p2, p3,
  nrow = 1, ncol = 3,
  common.legend = TRUE,
  legend = "top"
)

dev.off()


library(cobalt)
library(dplyr)
library(ggplot2)

chara1 <- colnames(dat)
ncolume1 <- ncol(dat)

notna_formula <- as.formula(
  paste(chara1[3],
        paste(" ~ ", paste(chara1[4:ncolume1], collapse = "+"))
  )
)

plot_dat <- dat
plot_dat$rx <- factor(plot_dat$rx)

plot_dat_ASMD <- data.frame(
  Variable = rownames(
    bal.tab(notna_formula, data = plot_dat, estimand = "ATE",
            stats = "mean.diffs",
            weights = rep(1, nrow(plot_dat)),
            binary = "std")$Balance.Across.Pairs
  ),
  UNAD  = as.numeric(bal.tab(notna_formula, data = plot_dat,
                             estimand = "ATE", stats = "mean.diffs",
                             binary = "std")$Balance.Across.Pairs$Max.Diff.Un),
  MB    = as.numeric(bal.tab(notna_formula, data = plot_dat,
                             estimand = "ATE", stats = "mean.diffs",
                             weights = wei_notna$MB,
                             binary = "std")$Balance.Across.Pairs$Max.Diff.Adj),
  IPW   = as.numeric(bal.tab(notna_formula, data = plot_dat,
                             estimand = "ATE", stats = "mean.diffs",
                             weights = wei_notna$ps,
                             binary = "std")$Balance.Across.Pairs$Max.Diff.Adj),
  EBAL  = as.numeric(bal.tab(notna_formula, data = plot_dat,
                             estimand = "ATE", stats = "mean.diffs",
                             weights = wei_notna$ebal,
                             binary = "std")$Balance.Across.Pairs$Max.Diff.Adj),
  CBPS  = as.numeric(bal.tab(notna_formula, data = plot_dat,
                             estimand = "ATE", stats = "mean.diffs",
                             weights = wei_notna$cbps,
                             binary = "std")$Balance.Across.Pairs$Max.Diff.Adj),
  ENERGY= as.numeric(bal.tab(notna_formula, data = plot_dat,
                             estimand = "ATE", stats = "mean.diffs",
                             weights = wei_notna$energy,
                             binary = "std")$Balance.Across.Pairs$Max.Diff.Adj),
  KERNEL= as.numeric(bal.tab(notna_formula, data = plot_dat,
                             estimand = "ATE", stats = "mean.diffs",
                             weights = wei_notna$kernel,
                             binary = "std")$Balance.Across.Pairs$Max.Diff.Adj),
  MDABW = as.numeric(bal.tab(notna_formula, data = plot_dat,
                             estimand = "ATE", stats = "mean.diffs",
                             weights = wei_notna$MDABW,
                             binary = "std")$Balance.Across.Pairs$Max.Diff.Adj)
)

ASMD_all <- c(plot_dat_ASMD$UNAD,
              plot_dat_ASMD$MB,
              plot_dat_ASMD$IPW,
              plot_dat_ASMD$EBAL,
              plot_dat_ASMD$CBPS,
              plot_dat_ASMD$ENERGY,
              plot_dat_ASMD$KERNEL,
              plot_dat_ASMD$MDABW)

group_all <- rep(c("UNAD","MB","IPW","EBAL","CBPS","ENERGY","KERNEL","MDABW"),
                 each = nrow(plot_dat_ASMD))

varnames <- plot_dat_ASMD$Variable

plot_dat_ASMD_long <- data.frame(
  id = rep(1:nrow(plot_dat_ASMD), 8),
  variable = rep(varnames, 8),
  ASMD_ord = ASMD_all,
  asmd = pmin(0.3, ASMD_all),
  weight = group_all,
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

desired_levels <- c("UNAD","MB","IPW","EBAL","CBPS","ENERGY","KERNEL","MDABW")

plot_dat_ASMD_long$weight <- factor(plot_dat_ASMD_long$weight,
                                    levels = desired_levels)

order_unad <- plot_dat_ASMD_long %>%
  filter(weight == "UNAD") %>%
  arrange(ASMD_ord) %>%
  pull(variable)

plot_dat_ASMD_long$variable <- factor(plot_dat_ASMD_long$variable,
                                      levels = order_unad)

p_ASMD <- ggplot(plot_dat_ASMD_long, aes(x = variable, y = asmd)) +
  geom_point(aes(color = weight, shape = weight), size = 3) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(color = 'Method', shape = 'Method',
       #title = "Covariate Balance",
       #subtitle = "Absolute Standardized Mean Difference",
       x = NULL, y = "ASMD") +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  ) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "solid") +
  scale_y_continuous(limits = c(0,0.3),
                     breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3)) +
  coord_flip()

png("ASMD_colon_cross_etype1.png", width = 6.4, height = 6.4, units = "in", res = 1000)
print(p_ASMD)
dev.off()


####====
plot_dat2 <- data.frame(wei_notna$Imbalance)
names(plot_dat2) <- tolower(names(plot_dat2))
plot_dat2$variable <- rownames(plot_dat2)

methods_order <- c("unad","mb","ps","ebal","cbps","energy","kernel","mdabw")
display_names <- c("UNAD","MB","IPW","EBAL","CBPS","ENERGY","KERNEL","MDABW")

imbal_list <- lapply(seq_along(methods_order), function(i){
  m <- methods_order[i]
  data.frame(
    variable = plot_dat2$variable,
    weight = display_names[i],
    imbal = plot_dat2[[m]],
    stringsAsFactors = FALSE
  )
})

plot_dat2_long <- do.call(rbind, imbal_list)

order_unad <- subset(plot_dat2_long, weight == "UNAD")
order_unad <- order_unad[order(order_unad$imbal), "variable"]
order_unad <- unique(order_unad)

plot_dat2_long$variable <- factor(plot_dat2_long$variable, levels = order_unad)

desired_levels <- c('UNAD','ENERGY','IPW','KERNEL','CBPS','MDABW','EBAL','MB')
plot_dat2_long$weight <- factor(plot_dat2_long$weight, levels = desired_levels)

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

library(ggplot2)

plot.GASMD_IMP <- ggplot(plot_dat2_long, aes(x = variable, y = imbal)) +
  geom_point(aes(color = weight, shape = weight), size = 3) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map, breaks = desired_levels) +
  scale_y_continuous(
    limits = c(0, 0.32),
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
  ) +
  labs(color = 'Method', shape = 'Method',
       #title = "Covariate Balance",
       #subtitle = "Generalized Absolute Standardized Mean Difference",
       x = NULL, y = "GASMD") +
  coord_flip() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black")

png("GASMD_notna_cross_etype1.png",
    width = 6.4, height = 6.4, units = "in", res = 1000)
print(plot.GASMD_IMP)
dev.off()

