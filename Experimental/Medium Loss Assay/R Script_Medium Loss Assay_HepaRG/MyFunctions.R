TransData <- function(x) {
  
  x_hat <- log(x+0.00000001)
  
  return(x_hat)
  
}

DeTransData <- function(x_hat) {
  
  x <- exp(x_hat)-0.00000001
  
  return(x)
  
}



GetPred3 <- function(est,dataset,t) {
  
  est["sd_w"]     <- est[1]
  est["k_0.1_w"]  <- est[2]
  est["sd_c"]     <- est[3]
  est["k_0.1_c"]  <- est[4]
  
  
  # Create event table
  ev <- eventTable(amount.units ="percent", time.units="hours")
  
  # Add time span to event table
  ev$add.sampling(t)
  
  # Define initial conditions
  inits <- c(S_0.1_w=1,P_0.1_w=0,
             S_0.1_c=1,P_0.1_c=0)
  
  #solve ODE system
  out <- solve(mMonoexponential,est,ev,inits)
  
  #transform predicitons
  # create new rows with predicitons
  
  
  # 0.1 warm
  new_rows_0.1_w <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_0.1_w)              <- colnames(dataset)
  new_rows_0.1_w$Species                <- "warm"
  new_rows_0.1_w$Time                   <- t
  new_rows_0.1_w$Concentration          <- "0.1"
  new_rows_0.1_w$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="warm"&dataset$Concentration=="0.1")])
  new_rows_0.1_w$DepletionRate1         <- est["k_0.1_w"]
  new_rows_0.1_w$ModelSD                <- est["sd_w"]
  new_rows_0.1_w$Fit                    <- TransData(out$S_0.1_w) + rnorm(n=length(t),mean=0,sd=est["sd_w"])
  
  # 0.1 cold
  new_rows_0.1_c <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_0.1_c)              <- colnames(dataset)
  new_rows_0.1_c$Species                <- "cold"
  new_rows_0.1_c$Time                   <- t
  new_rows_0.1_c$Concentration          <- "0.1"
  new_rows_0.1_c$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="cold"&dataset$Concentration=="0.1")])
  new_rows_0.1_c$DepletionRate1         <- est["k_0.1_c"]
  new_rows_0.1_c$ModelSD                <- est["sd_c"]
  new_rows_0.1_c$Fit                    <- TransData(out$S_0.1_c) + rnorm(n=length(t),mean=0,sd=est["sd_c"])
  
  
  new_rows <- rbind(new_rows_0.1_w,new_rows_0.1_c)
  
  return(new_rows)
}

GetPred2 <- function(est,dataset,t) {
  
  est["sd_w"]     <- est[1]
  est["k_0.1_w"]  <- est[2]
  est["sd_c"]     <- est[3]
  est["k_0.1_c"]  <- est[4]
  
  
  # Create event table
  ev <- eventTable(amount.units ="percent", time.units="hours")
  
  # Add time span to event table
  ev$add.sampling(t)
  
  # Define initial conditions
  inits <- c(S_0.1_w=1,P_0.1_w=0,
             S_0.1_c=1,P_0.1_c=0)
  
  #solve ODE system
  out <- solve(mMonoexponential,est,ev,inits)
  
  #transform predicitons
  # create new rows with predicitons
  
   # 0.1 warm
  new_rows_0.1_w <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_0.1_w)              <- colnames(dataset)
  new_rows_0.1_w$Species                <- "warm"
  new_rows_0.1_w$Time                   <- t
  new_rows_0.1_w$Concentration          <- "0.1"
  new_rows_0.1_w$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="warm"&dataset$Concentration=="0.1")])
  new_rows_0.1_w$DepletionRate1         <- est["k_0.1_w"]
  new_rows_0.1_w$ModelSD                <- est["sd_w"]
  new_rows_0.1_w$Fit                    <- TransData(out$S_0.1_w)
  
  # 0.1 cold
  new_rows_0.1_c <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_0.1_c)              <- colnames(dataset)
  new_rows_0.1_c$Species                <- "cold"
  new_rows_0.1_c$Time                   <- t
  new_rows_0.1_c$Concentration          <- "0.1"
  new_rows_0.1_c$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="cold"&dataset$Concentration=="0.1")])
  new_rows_0.1_c$DepletionRate1         <- est["k_0.1_c"]
  new_rows_0.1_c$ModelSD                <- est["sd_c"]
  new_rows_0.1_c$Fit                    <- TransData(out$S_0.1_c)
  
 
  new_rows <- rbind(new_rows_0.1_w,new_rows_0.1_c)
  
  return(new_rows)
}


GetPred <- function(est,dataset,t) {
  
  est["sd_w"]     <- est[1]
  est["k_0.1_w"]  <- est[2]
  est["sd_c"]     <- est[3]
  est["k_0.1_c"]  <- est[4]
  
 
  # Create event table
  ev <- eventTable(amount.units ="percent", time.units="hours")
  
  # Add time span to event table
  ev$add.sampling(t)
  
  
  # Define initial conditions
  inits <- c(S_0.1_w=1,P_0.1_w=0,
             S_0.1_c=1,P_0.1_c=0)
  
  #solve ODE system
  out <- solve(mMonoexponential,est,ev,inits)
  
  #transform predicitons
  #transfer predictions for substrate into dataset
  dataset$Fit[which(dataset$Concentration==0.1 & dataset$Species=="warm")]  <- rep(TransData(out$S_0.1_w),each=2)
  dataset$Fit[which(dataset$Concentration==0.1 & dataset$Species=="cold")]  <- rep(TransData(out$S_0.1_c),each=2)
  
  return(dataset)
}





GetMinus2LL <- function(est,dataset,t)  {
  
  est["sd_w"]     <- est[1]
  est["k_0.1_w"]  <- est[2]
  est["sd_c"]     <- est[3]
  est["k_0.1_c"]  <- est[4]
  
  #--- SOLVE ODE SYSTEM ---------------------------------------------------------
  dataset_pred <- GetPred(est,dataset,t)
  
  #--- CALCULATE SSR ------------------------------------------------------------
 
   # Remove columns with NormalizedConcentration is NA
  dataset_pred_narm <- dataset_pred[complete.cases(dataset_pred[ ,"NormalizedConcentration"]),]
  
  SSR_w <-  sum((dataset_pred_narm$NormalizedConcentration[which(dataset_pred_narm$Species=="warm")] - dataset_pred_narm$Fit[which(dataset_pred_narm$Species == "warm")])^2)
  SSR_c <-  sum((dataset_pred_narm$NormalizedConcentration[which(dataset_pred_narm$Species=="cold")] - dataset_pred_narm$Fit[which(dataset_pred_narm$Species == "cold")])^2)
  
  SSR <- SSR_w+SSR_c
  
  #--- DERIVE -2 log(L) WITH a^2 AS VARIANCE FOR IID (INDEPENDENT IDENTICALLY DISTRUBUTED) NORMAL DISTRUBUTED ERROR
  
 
  # constant
  const_w <- length(dataset_pred$Time[which(dataset_pred$Species=="warm")]) * log(2*pi*est["sd_w"]^2)
  const_c <- length(dataset_pred$Time[which(dataset_pred$Species=="cold")]) * log(2*pi*est["sd_c"]^2)

  
  # minus 2 log likelihood
  minus2LL_w <- 1/est["sd_w"]^2 * SSR_w + const_w 
  minus2LL_c <- 1/est["sd_c"]^2 * SSR_c + const_c
  
  minus2LL <- minus2LL_w + minus2LL_c  
  
  
  # rv <- sample(1:100,1,replace=TRUE)
  # if(rv>95){
  #  print(minus2LL)
  #  print(est)
  # }
  
  if (is.na(minus2LL)) {
    minus2LL <- Inf
    
    print("Minus2LL returned NaN for est:")
    print(est)
    print(SSR)
    print(dataset_pred)
    
    break
  }
  
  return(minus2LL)
}







