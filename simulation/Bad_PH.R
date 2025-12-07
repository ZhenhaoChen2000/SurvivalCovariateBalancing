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
  ############修改了
  if(sum(ebal.weight[Tr == 1]) < 10^(-6)){ebal.weight[Tr == 1] <- rep(1,sum(Tr == 1))}
  if(sum(ebal.weight[Tr == 2]) < 10^(-6)){ebal.weight[Tr == 2] <- rep(1,sum(Tr == 2))}
  if(sum(ebal.weight[Tr == 3]) < 10^(-6)){ebal.weight[Tr == 3] <- rep(1,sum(Tr == 3))}
  
  
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


#ph data generation======
si.data <- function(sample.size = 500, dimension = 10){
  set.seed(12345)
  X   <- matrix(rnorm(sample.size * dimension), nrow = sample.size, ncol = dimension)
  ps1 <- rep(1,sample.size)
  ps2 <- exp(4 * apply(X[,1:5],1,sum))
  ps3 <- exp(2 * apply(X[,6:10],1,sum)) 
  ps  <- cbind(ps1 / (ps1 + ps2 + ps3),ps2 / (ps1 + ps2 + ps3),ps3 / (ps1 + ps2 + ps3))
  uniform <- runif(sample.size)
  index   <- cbind(ps[,1] < uniform, ps[,1] + ps[,2] < uniform)
  Tr      <- apply(index,1,sum) + 1
  return(list(X = X, Tr = Tr))
}

data3   <- si.data(sample.size = n_total)


#estimate true theta of ph=====
n_true<-n
set.seed(2022)
sam_ph_true=sample(1:n_total,n)


dat3<-data3$X[sam_ph_true,];trt3<-data3$Tr[sam_ph_true]
cov3<-data.frame(id=1:n,trt=(trt3-1),X=dat3[,c(1,2,3,4,5)])

sdat.event3<-simsurv(dist = 'weibull',lambdas=0.3,gammas = 1.5,
                     betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),
                     x=cov3,maxt=10)
pcov1<-data.frame(id=1:n,trt3=rep((1-1),n),X=dat3[,c(1,2,3,4,5)])
pcov2<-data.frame(id=1:n,trt3=rep((2-1),n),X=dat3[,c(1,2,3,4,5)])
pcov3<-data.frame(id=1:n,trt3=rep((3-1),n),X=dat3[,c(1,2,3,4,5)])
set.seed(2022);
cpsdat1<-simsurv(dist = 'weibull',lambdas=0.3,gammas = 1.5,
                 betas = c(trt3=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),
                 x=pcov1,maxt=10)
set.seed(2022);
cpsdat2<-simsurv(dist = 'weibull',lambdas=0.3,gammas = 1.5,
                 betas = c(trt3=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),
                 x=pcov2,maxt=10)
set.seed(2022);
cpsdat3<-simsurv(dist = 'weibull',lambdas=0.3,gammas = 1.5,
                 betas = c(trt3=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),
                 x=pcov3,maxt=10)

pro_cen=(sum(cpsdat1$status==0)+sum(cpsdat2$status==0)+sum(cpsdat3$status==0))/(3*n)
#percent of censored

akm_rmst(time = c(cpsdat1$eventtime,cpsdat2$eventtime,cpsdat3$eventtime)
         ,c(cpsdat1$status,cpsdat2$status,cpsdat3$status),tau=10,
         group = as.factor(c(rep(0,n),rep(1,n),rep(2,n))))


theta_ph_true=get_akm_rmst(time = c(cpsdat1$eventtime,cpsdat2$eventtime,cpsdat3$eventtime)
                         ,c(cpsdat1$status,cpsdat2$status,cpsdat3$status),tau=10,
                         group = as.factor(c(rep(0,n),rep(1,n),rep(2,n))))





#weighted ph======

theta_ph=data.frame( ps_theta_1_2=rep(0,s),ps_theta_1_3=rep(0,s),ps_theta_2_3=rep(0,s),
                   cbps_theta_1_2=rep(0,s),cbps_theta_1_3=rep(0,s),cbps_theta_2_3=rep(0,s),
                   energy_theta_1_2=rep(0,s),energy_theta_1_3=rep(0,s),energy_theta_2_3=rep(0,s),
                   mb_theta_1_2=rep(0,s),mb_theta_1_3=rep(0,s),mb_theta_2_3=rep(0,s),
                   ebal_theta_1_2=rep(0,s),ebal_theta_1_3=rep(0,s),ebal_theta_2_3=rep(0,s),
                   no_theta_1_2=rep(0,s),no_theta_1_3=rep(0,s),no_theta_2_3=rep(0,s),
                   kernel_theta_1_2=rep(0,s),kernel_theta_1_3=rep(0,s),kernel_theta_2_3=rep(0,s),
                   MDABW_theta_1_2=rep(0,s),MDABW_theta_1_3=rep(0,s),MDABW_theta_2_3=rep(0,s))

k=matrix(1:n_total,nrow = h)

Im_ps_ph_1_2=matrix(rep(0,10*s),ncol = 10);Im_ps_ph_1_3=matrix(rep(0,10*s),ncol = 10);Im_ps_ph_2_3=matrix(rep(0,10*s),ncol = 10)
Im_cbps_ph_1_2=matrix(rep(0,10*s),ncol = 10);Im_cbps_ph_1_3=matrix(rep(0,10*s),ncol = 10);Im_cbps_ph_2_3=matrix(rep(0,10*s),ncol = 10)
Im_energy_ph_1_2=matrix(rep(0,10*s),ncol = 10);Im_energy_ph_1_3=matrix(rep(0,10*s),ncol = 10);Im_energy_ph_2_3=matrix(rep(0,10*s),ncol = 10)
Im_MB_ph_1_2=matrix(rep(0,10*s),ncol = 10);Im_MB_ph_1_3=matrix(rep(0,10*s),ncol = 10);Im_MB_ph_2_3=matrix(rep(0,10*s),ncol = 10)
Im_ebal_ph_1_2=matrix(rep(0,10*s),ncol = 10);Im_ebal_ph_1_3=matrix(rep(0,10*s),ncol = 10);Im_ebal_ph_2_3=matrix(rep(0,10*s),ncol = 10)
Im_no_ph_1_2=matrix(rep(0,10*s),ncol = 10);Im_no_ph_1_3=matrix(rep(0,10*s),ncol = 10);Im_no_ph_2_3=matrix(rep(0,10*s),ncol = 10)

Im_kernel_ph_1_2=matrix(rep(0,10*s),ncol = 10);Im_kernel_ph_1_3=matrix(rep(0,10*s),ncol = 10);Im_kernel_ph_2_3=matrix(rep(0,10*s),ncol = 10)
Im_MDABW_ph_1_2=matrix(rep(0,10*s),ncol = 10);Im_MDABW_ph_1_3=matrix(rep(0,10*s),ncol = 10);Im_MDABW_ph_2_3=matrix(rep(0,10*s),ncol = 10)

Im_GASMD <- array(NA,c(8,s,10))

for (i in 1:s) {
  X=data3$X[k[,i],]
  Tr=data3$Tr[k[,i]]
  result <- weight.esti(X , Tr)
  cov<-data.frame(id=1:nrow(X),trt=(Tr-1),X=X[,c(1,2,3,4,5)])
  set.seed(i+2025)  
  sdat<-simsurv(dist = 'weibull',lambdas=0.3,gammas = 1.5,
                betas = c(trt=-0.8,X.1=-1,X.2=1,X.3=0.5,X.4=-0.5,X.5=log(1.5)),
                x=cov,maxt=10)
  
  time<-cbind(sdat,Tr)
  theta_ph[i,1:3]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                      group=as.factor(time$Tr), weight=result$ps, tau=10)
  theta_ph[i,4:6]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                      group=as.factor(time$Tr), weight=result$cbps, tau=10)
  theta_ph[i,7:9]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                        group=as.factor(time$Tr), weight=result$energy, tau=10)
  theta_ph[i,10:12]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                        group=as.factor(time$Tr), weight=result$MB, tau=10)
  theta_ph[i,13:15]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                          group=as.factor(time$Tr), weight=(result$ebal+min(result$ebal[which(result$ebal>10^(-7))])*(result$ebal==0)), tau=10)
  theta_ph[i,16:18]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                          group=as.factor(time$Tr), tau=10)
  
  theta_ph[i,19:21]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                 group=as.factor(time$Tr), weight=result$kernel,tau=10)
  theta_ph[i,22:24]=get_akm_rmst(time=time$eventtime, status=time$status, 
                                 group=as.factor(time$Tr), weight=(result$MDABW+min(result$MDABW[which(result$MDABW>10^(-7))])*(result$MDABW==0)),tau=10)
  #calculate imbalance
  sample.size <- dim(X)[1]
  dimension   <- dim(X)[2]
  data.matrix <- data.frame(X , factor(Tr))
  character   <- names(data.matrix)
  for(j in 1:(dimension+1)){character[j] <- paste(character[j])}
  myformula   <- as.formula(paste(character[1 + dimension],paste(" ~ ", paste(character[1:dimension], collapse= "+"))))
  
  #ASMD
  ps.bal     <- bal.tab(X, treat = character[1 + dimension], data = data.matrix, weights = result$ps)
  ebal.bal   <- bal.tab(X, treat = character[1 + dimension], data = data.matrix, weights = (result$ebal+min(result$ebal[which(result$ebal>10^(-7))])*(result$ebal==0)))
  cbps.bal   <- bal.tab(X, treat = character[1 + dimension], data = data.matrix, weights = result$cbps)
  energy.bal <- bal.tab(X, treat = character[1 + dimension], data = data.matrix, weights = result$energy)
  MB.bal     <- bal.tab(X, treat = character[1 + dimension], data = data.matrix, weights = result$MB)
  no.bal     <- bal.tab(X, treat = character[1 + dimension], data = data.matrix, weights = rep(1,h))
  
  kernel.bal <- bal.tab(X, treat = character[1 + dimension], data = data.matrix, weights = result$kernel)
  MDABW.bal <- bal.tab(X, treat = character[1 + dimension], data = data.matrix, weights = (result$MDABW+min(result$MDABW[which(result$MDABW>10^(-7))])*(result$MDABW==0)))
  
  Im_ps_ph_1_2[i,]=ps.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_ps_ph_1_3[i,]=ps.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_ps_ph_2_3[i,]=ps.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_ebal_ph_1_2[i,]=ebal.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_ebal_ph_1_3[i,]=ebal.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_ebal_ph_2_3[i,]=ebal.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_cbps_ph_1_2[i,]=cbps.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_cbps_ph_1_3[i,]=cbps.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_cbps_ph_2_3[i,]=cbps.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_energy_ph_1_2[i,]=energy.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_energy_ph_1_3[i,]=energy.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_energy_ph_2_3[i,]=energy.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_MB_ph_1_2[i,]=MB.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_MB_ph_1_3[i,]=MB.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_MB_ph_2_3[i,]=MB.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_no_ph_1_2[i,]=no.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_no_ph_1_3[i,]=no.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_no_ph_2_3[i,]=no.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  
  Im_kernel_ph_1_2[i,]=kernel.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_kernel_ph_1_3[i,]=kernel.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_kernel_ph_2_3[i,]=kernel.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_MDABW_ph_1_2[i,]=MDABW.bal$Pair.Balance$`2 vs. 1`$Balance$Diff.Adj
  Im_MDABW_ph_1_3[i,]=MDABW.bal$Pair.Balance$`3 vs. 1`$Balance$Diff.Adj
  Im_MDABW_ph_2_3[i,]=MDABW.bal$Pair.Balance$`3 vs. 2`$Balance$Diff.Adj
  
  Im_GASMD[1,i,] <- result$Imbalance$unad
  Im_GASMD[2,i,] <- result$Imbalance$ps
  Im_GASMD[3,i,] <- result$Imbalance$ebal
  Im_GASMD[4,i,] <- result$Imbalance$cbps
  Im_GASMD[5,i,] <- result$Imbalance$energy
  Im_GASMD[6,i,] <- result$Imbalance$MB
  
  Im_GASMD[7,i,] <- result$Imbalance$kernel
  Im_GASMD[8,i,] <- result$Imbalance$MDABW
}
#ASMD==================

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
#mean
Im_ps_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ps_ph_1_2), 2,mean),apply(abs(Im_ps_ph_1_3), 2,mean),apply(abs(Im_ps_ph_2_3), 2,mean)))
Im_cbps_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_cbps_ph_1_2), 2,mean),apply(abs(Im_cbps_ph_1_3), 2,mean),apply(abs(Im_cbps_ph_2_3), 2,mean)))
Im_energy_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_energy_ph_1_2), 2,mean),apply(abs(Im_energy_ph_1_3), 2,mean),apply(abs(Im_energy_ph_2_3), 2,mean)))
Im_MB_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MB_ph_1_2), 2,mean),apply(abs(Im_MB_ph_1_3), 2,mean),apply(abs(Im_MB_ph_2_3), 2,mean)))
Im_ebal_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ebal_ph_1_2), 2,mean),apply(abs(Im_ebal_ph_1_3), 2,mean),apply(abs(Im_ebal_ph_2_3), 2,mean)))
Im_no_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_no_ph_1_2), 2,mean),apply(abs(Im_no_ph_1_3), 2,mean),apply(abs(Im_no_ph_2_3), 2,mean)))

Im_kernel_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_kernel_ph_1_2), 2,mean),apply(abs(Im_kernel_ph_1_3), 2,mean),apply(abs(Im_kernel_ph_2_3), 2,mean)))
Im_MDABW_mean=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MDABW_ph_1_2), 2,mean),apply(abs(Im_MDABW_ph_1_3), 2,mean),apply(abs(Im_MDABW_ph_2_3), 2,mean)))

Im_ph_mean=t(data.frame(ps=t(Im_ps_mean),
                      cbps=t(Im_cbps_mean),
                      energy=t(Im_energy_mean),
                      MB=t(Im_MB_mean),
                      ebal=t(Im_ebal_mean),
                      no=t(Im_no_mean),
                      kernel=t(Im_kernel_mean),
                      MDABW=t(Im_MDABW_mean)))

Im_ph_mean_max=data.frame(ps=apply(Im_ph_mean[1:3,], 2, max),
                         cbps=apply(Im_ph_mean[4:6,], 2, max),
                         energy=apply(Im_ph_mean[7:9,], 2, max),
                         MB=apply(Im_ph_mean[10:12,], 2, max),
                         ebal=apply(Im_ph_mean[13:15,], 2, max),
                         no=apply(Im_ph_mean[16:18,], 2, max),
                         kernel=apply(Im_ph_mean[19:21,], 2, max),
                         MDABW=apply(Im_ph_mean[22:24,], 2, max))



#sum
Im_ps_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ps_ph_1_2), 2,sum),apply(abs(Im_ps_ph_1_3), 2,sum),apply(abs(Im_ps_ph_2_3), 2,sum)))
Im_cbps_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_cbps_ph_1_2), 2,sum),apply(abs(Im_cbps_ph_1_3), 2,sum),apply(abs(Im_cbps_ph_2_3), 2,sum)))
Im_energy_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_energy_ph_1_2), 2,sum),apply(abs(Im_energy_ph_1_3), 2,sum),apply(abs(Im_energy_ph_2_3), 2,sum)))
Im_MB_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MB_ph_1_2), 2,sum),apply(abs(Im_MB_ph_1_3), 2,sum),apply(abs(Im_MB_ph_2_3), 2,sum)))
Im_ebal_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ebal_ph_1_2), 2,sum),apply(abs(Im_ebal_ph_1_3), 2,sum),apply(abs(Im_ebal_ph_2_3), 2,sum)))
Im_no_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_no_ph_1_2), 2,sum),apply(abs(Im_no_ph_1_3), 2,sum),apply(abs(Im_no_ph_2_3), 2,sum)))

Im_kernel_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_kernel_ph_1_2), 2,sum),apply(abs(Im_kernel_ph_1_3), 2,sum),apply(abs(Im_kernel_ph_2_3), 2,sum)))
Im_MDABW_sum=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MDABW_ph_1_2), 2,sum),apply(abs(Im_MDABW_ph_1_3), 2,sum),apply(abs(Im_MDABW_ph_2_3), 2,sum)))

Im_ph_sum=t(data.frame(ps=t(Im_ps_sum),
                        cbps=t(Im_cbps_sum),
                        energy=t(Im_energy_sum),
                        MB=t(Im_MB_sum),
                        ebal=t(Im_ebal_sum),
                        no=t(Im_no_sum),
                       kernel=t(Im_kernel_sum),
                       MDABW=t(Im_MDABW_sum)))

Im_ph_all_sum=data.frame(ps=apply(Im_ph_sum[1:3,], 2, sum),
                          cbps=apply(Im_ph_sum[4:6,], 2, sum),
                          energy=apply(Im_ph_sum[7:9,], 2, sum),
                          MB=apply(Im_ph_sum[10:12,], 2, sum),
                          ebal=apply(Im_ph_sum[13:15,], 2, sum),
                          no=apply(Im_ph_sum[16:18,], 2, sum),
                         kernel=apply(Im_ph_sum[19:21,], 2, sum),
                         MDABW=apply(Im_ph_sum[22:24,], 2, sum))


#max
Im_ps_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ps_ph_1_2), 2,max),apply(abs(Im_ps_ph_1_3), 2,max),apply(abs(Im_ps_ph_2_3), 2,max)))
Im_cbps_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_cbps_ph_1_2), 2,max),apply(abs(Im_cbps_ph_1_3), 2,max),apply(abs(Im_cbps_ph_2_3), 2,max)))
Im_energy_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_energy_ph_1_2), 2,max),apply(abs(Im_energy_ph_1_3), 2,max),apply(abs(Im_energy_ph_2_3), 2,max)))
Im_MB_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MB_ph_1_2), 2,max),apply(abs(Im_MB_ph_1_3), 2,max),apply(abs(Im_MB_ph_2_3), 2,max)))
Im_ebal_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_ebal_ph_1_2), 2,max),apply(abs(Im_ebal_ph_1_3), 2,max),apply(abs(Im_ebal_ph_2_3), 2,max)))
Im_no_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_no_ph_1_2), 2,max),apply(abs(Im_no_ph_1_3), 2,max),apply(abs(Im_no_ph_2_3), 2,max)))

Im_kernel_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_kernel_ph_1_2), 2,max),apply(abs(Im_kernel_ph_1_3), 2,max),apply(abs(Im_kernel_ph_2_3), 2,max)))
Im_MDABW_max=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),rbind(apply(abs(Im_MDABW_ph_1_2), 2,max),apply(abs(Im_MDABW_ph_1_3), 2,max),apply(abs(Im_MDABW_ph_2_3), 2,max)))
Im_ph_max=t(data.frame(ps=t(Im_ps_max),
                        cbps=t(Im_cbps_max),
                        energy=t(Im_energy_max),
                        MB=t(Im_MB_max),
                        ebal=t(Im_ebal_max),
                        no=t(Im_no_max),
                       kernel=t(Im_kernel_max),
                       MDABW=t(Im_MDABW_max)))

Im_ph_all_max=data.frame(ps=apply(Im_ph_max[1:3,], 2, max),
                         cbps=apply(Im_ph_max[4:6,], 2, max),
                         energy=apply(Im_ph_max[7:9,], 2, max),
                         MB=apply(Im_ph_max[10:12,], 2, max),
                         ebal=apply(Im_ph_max[13:15,], 2, max),
                         no=apply(Im_ph_max[16:18,], 2, max),
                         kernel=apply(Im_ph_max[19:21,], 2, max),
                         MDABW=apply(Im_ph_max[22:24,], 2, max))



#calculate the mean,bias and mse===========
###mean
ps_theta_ph_1_2_hat=mean(as.numeric(theta_ph$ps_theta_1_2))
cbps_theta_ph_1_2_hat=mean(as.numeric(theta_ph$cbps_theta_1_2))
energy_theta_ph_1_2_hat=mean(as.numeric(theta_ph$energy_theta_1_2))
mb_theta_ph_1_2_hat=mean(as.numeric(theta_ph$mb_theta_1_2))
ebal_theta_ph_1_2_hat=mean(as.numeric(theta_ph$ebal_theta_1_2))
no_theta_ph_1_2_hat=mean(as.numeric(theta_ph$no_theta_1_2))

kernel_theta_ph_1_2_hat=mean(as.numeric(theta_ph$kernel_theta_1_2))
MDABW_theta_ph_1_2_hat=mean(as.numeric(theta_ph$MDABW_theta_1_2))

ps_theta_ph_1_3_hat=mean(as.numeric(theta_ph$ps_theta_1_3))
cbps_theta_ph_1_3_hat=mean(as.numeric(theta_ph$cbps_theta_1_3))
energy_theta_ph_1_3_hat=mean(as.numeric(theta_ph$energy_theta_1_3))
mb_theta_ph_1_3_hat=mean(as.numeric(theta_ph$mb_theta_1_3))
ebal_theta_ph_1_3_hat=mean(as.numeric(theta_ph$ebal_theta_1_3))
no_theta_ph_1_3_hat=mean(as.numeric(theta_ph$no_theta_1_3))

kernel_theta_ph_1_3_hat=mean(as.numeric(theta_ph$kernel_theta_1_3))
MDABW_theta_ph_1_3_hat=mean(as.numeric(theta_ph$MDABW_theta_1_3))

ps_theta_ph_2_3_hat=mean(as.numeric(theta_ph$ps_theta_2_3))
cbps_theta_ph_2_3_hat=mean(as.numeric(theta_ph$cbps_theta_2_3))
energy_theta_ph_2_3_hat=mean(as.numeric(theta_ph$energy_theta_2_3))
mb_theta_ph_2_3_hat=mean(as.numeric(theta_ph$mb_theta_2_3))
ebal_theta_ph_2_3_hat=mean(as.numeric(theta_ph$ebal_theta_2_3))
no_theta_ph_2_3_hat=mean(as.numeric(theta_ph$no_theta_2_3))

kernel_theta_ph_2_3_hat=mean(as.numeric(theta_ph$kernel_theta_2_3))
MDABW_theta_ph_2_3_hat=mean(as.numeric(theta_ph$MDABW_theta_2_3))

ph_estimate=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"), 
                           true=c(theta_ph_true[1,1],theta_ph_true[1,2],theta_ph_true[1,3]),
                           ps=c(ps_theta_ph_1_2_hat,ps_theta_ph_1_3_hat,ps_theta_ph_2_3_hat),
                           cbps=c(cbps_theta_ph_1_2_hat,cbps_theta_ph_1_3_hat,cbps_theta_ph_2_3_hat),
                           energy=c(energy_theta_ph_1_2_hat,energy_theta_ph_1_3_hat,energy_theta_ph_2_3_hat),
                           mb=c(mb_theta_ph_1_2_hat,mb_theta_ph_1_3_hat,mb_theta_ph_2_3_hat),
                           ebal=c(ebal_theta_ph_1_2_hat,ebal_theta_ph_1_3_hat,ebal_theta_ph_2_3_hat),
                           no=c(no_theta_ph_1_2_hat,no_theta_ph_1_3_hat,no_theta_ph_2_3_hat),
                       kernel=c(kernel_theta_ph_1_2_hat,kernel_theta_ph_1_3_hat,kernel_theta_ph_2_3_hat),
                       MDABW=c(MDABW_theta_ph_1_2_hat,MDABW_theta_ph_1_3_hat,MDABW_theta_ph_2_3_hat))

#bias==========
p_bias=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),
                  ps=c(((ph_estimate[1,2]-ph_estimate[1,1])/ph_estimate[1,1]),((ph_estimate[2,2]-ph_estimate[2,1])/ph_estimate[2,1]),((ph_estimate[3,2]-ph_estimate[3,1])/ph_estimate[3,1])),
                  cbps=c(((ph_estimate[1,3]-ph_estimate[1,1])/ph_estimate[1,1]),((ph_estimate[2,3]-ph_estimate[2,1])/ph_estimate[2,1]),((ph_estimate[3,3]-ph_estimate[3,1])/ph_estimate[3,1])),
                  energy=c(((ph_estimate[1,4]-ph_estimate[1,1])/ph_estimate[1,1]),((ph_estimate[2,4]-ph_estimate[2,1])/ph_estimate[2,1]),((ph_estimate[3,4]-ph_estimate[3,1])/ph_estimate[3,1])),
                  mb=c(((ph_estimate[1,5]-ph_estimate[1,1])/ph_estimate[1,1]),((ph_estimate[2,5]-ph_estimate[2,1])/ph_estimate[2,1]),((ph_estimate[3,5]-ph_estimate[3,1])/ph_estimate[3,1])),
                  ebal=c(((ph_estimate[1,6]-ph_estimate[1,1])/ph_estimate[1,1]),((ph_estimate[2,6]-ph_estimate[2,1])/ph_estimate[2,1]),((ph_estimate[3,6]-ph_estimate[3,1])/ph_estimate[3,1])),
                  no=c(((ph_estimate[1,7]-ph_estimate[1,1])/ph_estimate[1,1]),((ph_estimate[2,7]-ph_estimate[2,1])/ph_estimate[2,1]),((ph_estimate[3,7]-ph_estimate[3,1])/ph_estimate[3,1])),
                  kernel=c(((ph_estimate[1,8]-ph_estimate[1,1])/ph_estimate[1,1]),((ph_estimate[2,8]-ph_estimate[2,1])/ph_estimate[2,1]),((ph_estimate[3,8]-ph_estimate[3,1])/ph_estimate[3,1])),
                  MDABW=c(((ph_estimate[1,9]-ph_estimate[1,1])/ph_estimate[1,1]),((ph_estimate[2,9]-ph_estimate[2,1])/ph_estimate[2,1]),((ph_estimate[3,9]-ph_estimate[3,1])/ph_estimate[3,1])))
p_bias_mean=apply(abs(p_bias), 2, mean)

#p_mse=============
p_mse_ps_1_2<-sum((as.numeric(theta_ph$ps_theta_1_2)-as.numeric(theta_ph_true[1]))^2)/s
p_mse_cbps_1_2<-sum((as.numeric(theta_ph$cbps_theta_1_2)-as.numeric(theta_ph_true[1]))^2)/s
p_mse_energy_1_2<-sum((as.numeric(theta_ph$energy_theta_1_2)-as.numeric(theta_ph_true[1]))^2)/s
p_mse_mb_1_2<-sum((as.numeric(theta_ph$mb_theta_1_2)-as.numeric(theta_ph_true[1]))^2)/s
p_mse_ebal_1_2<-sum((as.numeric(theta_ph$ebal_theta_1_2)-as.numeric(theta_ph_true[1]))^2)/s
p_mse_no_1_2<-sum((as.numeric(theta_ph$no_theta_1_2)-as.numeric(theta_ph_true[1]))^2)/s

p_mse_kernel_1_2<-sum((as.numeric(theta_ph$kernel_theta_1_2)-as.numeric(theta_ph_true[1]))^2)/s
p_mse_MDABW_1_2<-sum((as.numeric(theta_ph$MDABW_theta_1_2)-as.numeric(theta_ph_true[1]))^2)/s

p_mse_ps_1_3<-sum((as.numeric(theta_ph$ps_theta_1_3)-as.numeric(theta_ph_true[2]))^2)/s
p_mse_cbps_1_3<-sum((as.numeric(theta_ph$cbps_theta_1_3)-as.numeric(theta_ph_true[2]))^2)/s
p_mse_energy_1_3<-sum((as.numeric(theta_ph$energy_theta_1_3)-as.numeric(theta_ph_true[2]))^2)/s
p_mse_mb_1_3<-sum((as.numeric(theta_ph$mb_theta_1_3)-as.numeric(theta_ph_true[2]))^2)/s
p_mse_ebal_1_3<-sum((as.numeric(theta_ph$ebal_theta_1_3)-as.numeric(theta_ph_true[2]))^2)/s
p_mse_no_1_3<-sum((as.numeric(theta_ph$no_theta_1_3)-as.numeric(theta_ph_true[2]))^2)/s

p_mse_kernel_1_3<-sum((as.numeric(theta_ph$kernel_theta_1_3)-as.numeric(theta_ph_true[2]))^2)/s
p_mse_MDABW_1_3<-sum((as.numeric(theta_ph$MDABW_theta_1_3)-as.numeric(theta_ph_true[2]))^2)/s

p_mse_ps_2_3<-sum((as.numeric(theta_ph$ps_theta_2_3)-as.numeric(theta_ph_true[3]))^2)/s
p_mse_cbps_2_3<-sum((as.numeric(theta_ph$cbps_theta_2_3)-as.numeric(theta_ph_true[3]))^2)/s
p_mse_energy_2_3<-sum((as.numeric(theta_ph$energy_theta_2_3)-as.numeric(theta_ph_true[3]))^2)/s
p_mse_mb_2_3<-sum((as.numeric(theta_ph$mb_theta_2_3)-as.numeric(theta_ph_true[3]))^2)/s
p_mse_ebal_2_3<-sum((as.numeric(theta_ph$ebal_theta_2_3)-as.numeric(theta_ph_true[3]))^2)/s
p_mse_no_2_3<-sum((as.numeric(theta_ph$no_theta_2_3)-as.numeric(theta_ph_true[3]))^2)/s

p_mse_kernel_2_3<-sum((as.numeric(theta_ph$kernel_theta_2_3)-as.numeric(theta_ph_true[3]))^2)/s
p_mse_MDABW_2_3<-sum((as.numeric(theta_ph$MDABW_theta_2_3)-as.numeric(theta_ph_true[3]))^2)/s

p_mse=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),
                 ps=c(p_mse_ps_1_2,p_mse_ps_1_3,p_mse_ps_2_3),
                 cbps=c(p_mse_cbps_1_2,p_mse_cbps_1_3,p_mse_cbps_2_3),
                 energy=c(p_mse_energy_1_2,p_mse_energy_1_3,p_mse_energy_2_3),
                 mb=c(p_mse_mb_1_2,p_mse_mb_1_3,p_mse_mb_2_3),
                 ebal=c(p_mse_ebal_1_2,p_mse_ebal_1_3,p_mse_ebal_2_3),
                 no=c(p_mse_no_1_2,p_mse_no_1_3,p_mse_no_2_3),
                 kernel=c(p_mse_kernel_1_2,p_mse_kernel_1_3,p_mse_kernel_2_3),
                 MDABW=c(p_mse_MDABW_1_2,p_mse_MDABW_1_3,p_mse_MDABW_2_3))

p_mse_mean=apply(p_mse, 2, mean)

############
p_sam_ps_1_2<-sum((as.numeric(theta_ph$ps_theta_1_2)-mean(as.numeric(theta_ph$ps_theta_1_2)))^2)/s
p_sam_cbps_1_2<-sum((as.numeric(theta_ph$cbps_theta_1_2)-mean(as.numeric(theta_ph$cbps_theta_1_2)))^2)/s
p_sam_energy_1_2<-sum((as.numeric(theta_ph$energy_theta_1_2)-mean(as.numeric(theta_ph$energy_theta_1_2)))^2)/s
p_sam_mb_1_2<-sum((as.numeric(theta_ph$mb_theta_1_2)-mean(as.numeric(theta_ph$mb_theta_1_2)))^2)/s
p_sam_ebal_1_2<-sum((as.numeric(theta_ph$ebal_theta_1_2)-mean(as.numeric(theta_ph$ebal_theta_1_2)))^2)/s
p_sam_no_1_2<-sum((as.numeric(theta_ph$no_theta_1_2)-mean(as.numeric(theta_ph$no_theta_1_2)))^2)/s

p_sam_kernel_1_2<-sum((as.numeric(theta_ph$kernel_theta_1_2)-mean(as.numeric(theta_ph$kernel_theta_1_2)))^2)/s
p_sam_MDABW_1_2<-sum((as.numeric(theta_ph$MDABW_theta_1_2)-mean(as.numeric(theta_ph$MDABW_theta_1_2)))^2)/s

p_sam_ps_1_3<-sum((as.numeric(theta_ph$ps_theta_1_3)-mean(as.numeric(theta_ph$ps_theta_1_3)))^2)/s
p_sam_cbps_1_3<-sum((as.numeric(theta_ph$cbps_theta_1_3)-mean(as.numeric(theta_ph$cbps_theta_1_3)))^2)/s
p_sam_energy_1_3<-sum((as.numeric(theta_ph$energy_theta_1_3)-mean(as.numeric(theta_ph$energy_theta_1_3)))^2)/s
p_sam_mb_1_3<-sum((as.numeric(theta_ph$mb_theta_1_3)-mean(as.numeric(theta_ph$mb_theta_1_3)))^2)/s
p_sam_ebal_1_3<-sum((as.numeric(theta_ph$ebal_theta_1_3)-mean(as.numeric(theta_ph$ebal_theta_1_3)))^2)/s
p_sam_no_1_3<-sum((as.numeric(theta_ph$no_theta_1_3)-mean(as.numeric(theta_ph$no_theta_1_3)))^2)/s

p_sam_kernel_1_3<-sum((as.numeric(theta_ph$kernel_theta_1_3)-mean(as.numeric(theta_ph$kernel_theta_1_3)))^2)/s
p_sam_MDABW_1_3<-sum((as.numeric(theta_ph$MDABW_theta_1_3)-mean(as.numeric(theta_ph$MDABW_theta_1_3)))^2)/s

p_sam_ps_2_3<-sum((as.numeric(theta_ph$ps_theta_2_3)-mean(as.numeric(theta_ph$ps_theta_2_3)))^2)/s
p_sam_cbps_2_3<-sum((as.numeric(theta_ph$cbps_theta_2_3)-mean(as.numeric(theta_ph$cbps_theta_2_3)))^2)/s
p_sam_energy_2_3<-sum((as.numeric(theta_ph$energy_theta_2_3)-mean(as.numeric(theta_ph$energy_theta_2_3)))^2)/s
p_sam_mb_2_3<-sum((as.numeric(theta_ph$mb_theta_2_3)-mean(as.numeric(theta_ph$mb_theta_2_3)))^2)/s
p_sam_ebal_2_3<-sum((as.numeric(theta_ph$ebal_theta_2_3)-mean(as.numeric(theta_ph$ebal_theta_2_3)))^2)/s
p_sam_no_2_3<-sum((as.numeric(theta_ph$no_theta_2_3)-mean(as.numeric(theta_ph$no_theta_2_3)))^2)/s
 
p_sam_kernel_2_3<-sum((as.numeric(theta_ph$kernel_theta_2_3)-mean(as.numeric(theta_ph$kernel_theta_2_3)))^2)/s
p_sam_MDABW_2_3<-sum((as.numeric(theta_ph$MDABW_theta_2_3)-mean(as.numeric(theta_ph$MDABW_theta_2_3)))^2)/s

p_sam=data.frame(row.names = c("1 vs 2","1 vs 3","2 vs 3"),
                 ps=c(p_sam_ps_1_2,p_sam_ps_1_3,p_sam_ps_2_3),
                 cbps=c(p_sam_cbps_1_2,p_sam_cbps_1_3,p_sam_cbps_2_3),
                 energy=c(p_sam_energy_1_2,p_sam_energy_1_3,p_sam_energy_2_3),
                 mb=c(p_sam_mb_1_2,p_sam_mb_1_3,p_sam_mb_2_3),
                 ebal=c(p_sam_ebal_1_2,p_sam_ebal_1_3,p_sam_ebal_2_3),
                 no=c(p_sam_no_1_2,p_sam_no_1_3,p_sam_no_2_3),
                 kernel=c(p_sam_kernel_1_2,p_sam_kernel_1_3,p_sam_kernel_2_3),
                 MDABW=c(p_sam_MDABW_1_2,p_sam_MDABW_1_3,p_sam_MDABW_2_3))


#plot PH==========

#result3 <- weight.esti(data3$X[sam_ph_true,] , data3$Tr[sam_ph_true])

#time3<-cbind(sdat.event3,trt3)

#par(mfrow=c(2,3))
#akm_rmst(time=time3$eventtime, status=time3$status, 
#         group=as.factor(time3$trt3), weight=result3$ps, tau=10)
#title(sub = "weighted:ps")


#akm_rmst(time=time3$eventtime, status=time3$status, 
#         group=as.factor(time3$trt3), weight=result3$cbps, tau=10)
#title(sub = "weighted:CBPS")


#akm_rmst(time=time3$eventtime, status=time3$status, 
#         group=as.factor(time3$trt3), weight=result3$energy, tau=10)
#title(sub = "weighted:energy")


#akm_rmst(time=time3$eventtime, status=time3$status, 
#         group=as.factor(time3$trt3), weight=result3$MB, tau=10)
#title(sub = "weighted:MB")

#akm_rmst(time=time3$eventtime, status=time3$status, 
#         group=as.factor(time3$trt3), weight=(result3$ebal+1*(result3$ebal==0)), tau=10)
#title(sub = "weighted:ebal")

#akm_rmst(time=time3$eventtime, status=time3$status, 
#         group=as.factor(time3$trt3), tau=10)
#title(sub = "No other weight")

#write table======
wASMD=data.frame(
  Unad=Im_ph_mean_max$no,IPW=Im_ph_mean_max$ps,
  EBAL=Im_ph_mean_max$ebal,CBPS=Im_ph_mean_max$cbps,
  Energy=Im_ph_mean_max$energy,MB=Im_ph_mean_max$MB, kernel= Im_ph_mean_max$kernel, MDABW= Im_ph_mean_max$MDABW)

rownames(wASMD)=c(paste("Max-ASMD",c(1:10)))

wRMSE=data.frame(Unad=sqrt(p_mse_mean[6]),IPW=sqrt(p_mse_mean[1]),
                 EBAL=sqrt(p_mse_mean[5]),CBPS=sqrt(p_mse_mean[2]),
                 Energy=sqrt(p_mse_mean[3]),MB=sqrt(p_mse_mean[4]),
                 kernel=sqrt(p_mse_mean[7]),MDABW=sqrt(p_mse_mean[8]))
rownames(wRMSE)=c("RMSE-RMST")

wBias=data.frame(Unad=p_bias_mean[6],IPW=p_bias_mean[1],
                 EBAL=p_bias_mean[5],CBPS=p_bias_mean[2],
                 Energy=p_bias_mean[3],MB=p_bias_mean[4], kernel=p_bias_mean[7],MDABW=p_bias_mean[8])
rownames(wBias)=c("Relative Bias-RMST")

wSD=data.frame(Unad=mean(sqrt(p_sam[,6])),IPW=mean(sqrt(p_sam[,1])),
               EBAL=mean(sqrt(p_sam[,5])),CBPS=mean(sqrt(p_sam[,2])),
               Energy=mean(sqrt(p_sam[,3])),MB=mean(sqrt(p_sam[,4])), kernel=mean(sqrt(p_sam[,7])),MDABW=mean(sqrt(p_sam[,8])))
rownames(wSD)=c("SD-RMST")

wGASMD=data.frame(Unad=max(Im_GASMD_mean$unad),IPW=max(Im_GASMD_mean$ps),
                  EBAL=max(Im_GASMD_mean$ebal),CBPS=max(Im_GASMD_mean$cbps),
                  Energy=max(Im_GASMD_mean$energy),MB=max(Im_GASMD_mean$MB), kernel=max(Im_GASMD_mean$kernel),MDABW=max(Im_GASMD_mean$MDABW))
rownames(wGASMD)=c("Max-GASMD")

all_table=rbind(wBias,wRMSE,wSD,wASMD,wGASMD)
all_table=round(all_table,digits = 3)

write.csv(file = 'Bad_PH_result.csv',x=all_table)

#appendix
#save bias\sd\rmse for all pairwise comparison rather than mean of bias\sd\rmse for all pairwise comparison
#save RMST estimate
write.csv(file = 'Bad_PH_rmst_est.csv',x=ph_estimate)
#save bias
write.csv(file = 'Bad_PH_bias.csv',x=p_bias)
#save sd
write.csv(file = 'Bad_PH_sd.csv',x=p_sam)
#save mse
write.csv(file = 'Bad_PH_mse.csv',x=p_mse)
#save GASMD for all covariates
Im_GASMD_mean_df<-data.frame(Im_GASMD_mean)
write.csv(file = 'Bad_PH_gasmd.csv',x=Im_GASMD_mean_df)
