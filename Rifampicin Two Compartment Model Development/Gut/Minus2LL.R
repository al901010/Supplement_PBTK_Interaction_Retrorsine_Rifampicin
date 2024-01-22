
GetMinus2LL <- function(est,
                        t,
                        dataset, 
                        runtable) {
  
  est["a_RIF"]  <- est[1]
  est["ka"]     <- est[2]
  est["k12"]    <- est[3]
  est["k21"]    <- est[4]
  est["ke"]     <- est[5]
  
  
  #--- SOLVE ODE SYSTEM ---------------------------------------------------------
  
  dataset_pred <- PredRxODE(est,
                            t,
                            dataset=dataset_rifampicin, 
                            runtable=runtable_rifampicin, 
                            add_timepoints = FALSE)
  
  #--- DERIVE -2 log(L) WITH a^2 AS VARIANCE FOR IID (INDEPENDENT IDENTICALLY DISTRUBUTED) NORMAL DISTRUBUTED ERROR
  
  dataset_pred$SQUAREDRESIDUALS <- (dataset_pred$Y - dataset_pred$YPRED)^2
  
  # SSR <- sum(dataset_pred$SQUAREDRESIDUALS[which(
  #                    (!is.na(dataset_pred$SQUAREDRESIDUALS)) & 
  #                       ((dataset_pred$ID=="Rifampicin_human_gut_600mg_po_SD")|
  #                        (dataset_pred$ID=="Rifampicin_human_gut_600mg_po_MD"))   
  #                    )])
  # 
  # minus2LL <- 1/est["a_RIF"]^2 * SSR + length(dataset_pred$XORIGINAL[which((!(is.na(dataset_pred$Y)))&
  #                                            ((dataset_pred$ID=="Rifampicin_human_gut_600mg_po_SD")|
  #                                             (dataset_pred$ID=="Rifampicin_human_gut_600mg_po_MD")))]) *
  #                                          log(2*pi*est["a_RIF"]^2)
  
  SSR <- sum(dataset_pred$SQUAREDRESIDUALS[which(
    (!is.na(dataset_pred$SQUAREDRESIDUALS)) & 
      ((dataset_pred$ID=="Rifampicin_human_gut_600mg_po_MD"))   
  )])
  
  minus2LL <- 1/est["a_RIF"]^2 * SSR + length(dataset_pred$XORIGINAL[which((!(is.na(dataset_pred$Y)))&
                                                                             ((dataset_pred$ID=="Rifampicin_human_gut_600mg_po_MD")))]) *
    log(2*pi*est["a_RIF"]^2)
  
  # rv <- sample(1:100,1,replace=TRUE)
  # if(rv>90){
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


