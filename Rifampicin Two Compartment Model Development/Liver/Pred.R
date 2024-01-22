
###############################################################################
# PREDICTION FOR ESTIMATION INCLUDES ONLY OBSERVED TIMEPOINTS
###############################################################################

PredRxODE <- function(est,
                      t,
                      dataset, 
                      runtable, 
                      add_timepoints) {
  
  est["a_RIF"]  <- est[1]
  est["ka"]     <- est[2]
  est["k12"]    <- est[3]
  est["k21"]    <- est[4]
  est["ke"]     <- est[5]

#--- SOLVE ODE SYSTEM FOR RUNTABLE AND PROCESS SOLUTION -----------------------
  
  # Initialise empty lists of length = nrow(runtable)
  results_ev   <- vector("list",nrow(runtable))
  results_wide <- vector("list",nrow(runtable))
  results_long <- vector("list",nrow(runtable))
  results_sum  <- vector("list",nrow(runtable))

  #--- LOOP THROUGH RUNTABLE
  
  for(row in 1:nrow(runtable_rifampicin)) {
    
    # event table
    ev <- eventTable(amount.units = "mol", time.units="h")
    
    # Add time span to event table
    ev$add.sampling(t)
    
    # Dose [mg] is converted to [mol] 
    dosing_value <- runtable_rifampicin[row,"DOSEORIGINAL"]/1000/822.9 #molar mass rifampicin (Pubchem)
   
     ev$add.dosing(
        dose=dosing_value,
        nbr.doses=runtable[row,]$NODOSES,
        dosing.interval=runtable[row,]$DOSINGINTERVAL,
        start.time=24*4,
        dosing.to=1)
      
  
    # Define initial conditions
    inits <- c(RIFdep=0,RIFcen=0,RIFper=0,RIFbag=0)
  
    # Solve ODE system 
    wide <- mRifampicin %>% rxSolve(est,ev,inits,method="lsoda",atol=1e-12,rtol=1e-10) #,maxsteps=500000,atol=1e-12,rtol=1e-10) 
  
    # Check balance
    sum <- rowSums(wide) - parse_number(as.character(wide$time)) 
  
    # Convert units [mol] in [umol]
    wide[,2: ncol(wide)] <-  wide[,2: ncol(wide)]*10^6
  
  
    #--- TRANSFORM PREDICTIONS ------------------------------------------------
    
      wide$RIFdep   <- TransData(wide$RIFdep)
      wide$RIFcen  <- TransData(wide$RIFcen)
      wide$RIFper  <- TransData(wide$RIFper)
      wide$RIFbag  <- TransData(wide$RIFbag)
      
    # Convert wide to long with tidyr package with new key column "compartment" and new value column "value"
    long <- gather(wide, compartment, value, 2:ncol(wide), factor_key=TRUE)
    
    # Save results in list
    results_ev[[row]]   <- ev
    results_sum[[row]]  <- sum
    results_wide[[row]] <- wide
    results_long[[row]] <- long
  }
  
  
  # # View dosing in event tables
  # print("Dosing")
  # print(results_ev[[1]]$get.dosing())
  # print(results_ev[[2]]$get.dosing())

  # # View head of balance check and results in wide format
  # print("Balance check")
  # print(head(results_sum))
  # print("Predictions wide-format")
  # print(head(results_wide))
  
  # Assign variable in the global environment as it is needed for plotting
  results_long <<- results_long 
  
  
  #--- INCLUDE PREDICTIONS IN DATASET -------------------------------------------
  
  # Initialize new columns 
  dataset$YPRED <- NA
  dataset$YPREDUNIT <- NA

  # Extract ids
  ids <- unique(dataset$ID)
  
  # Loop through ids
  for (row in ids) {
    
    # Identify runid for id
    runid <-  unique(dataset$RUNID[which(dataset$ID==row)])
    
    # Identify rownumbers of timepoints from results_wide[[runid]]$time that match timepoints of observed data
    is_timepoints <- match(dataset$XORIGINAL[which(dataset$ID == row & 
                                                     is.na(dataset$DOSEORIGINAL))],
                           results_wide[[runid]]$time)
    
    # Transfer predictions for specific compartment at these timepoints to dataset column YPRED
    YPRED <- rep(NA, length(is_timepoints))
    
    #if(row=="Rifampicin_human_liver_600mg_po_SD")      {YPRED         <- results_wide[[runid]]$RIFcen[is_timepoints]    }
    if(row=="Rifampicin_human_liver_600mg_po_MD")      {YPRED         <- results_wide[[runid]]$RIFcen[is_timepoints]    }

    dataset$YPRED[which(dataset$ID==row & !(is.na(dataset$Y)))] <- YPRED
  }
  
  #add units
  dataset$YPREDUNIT[which(is.na(dataset$DOSEORIGINAL))] <- "umol"
  
  
  #--- INCLUDE MORE TIMEPOINTS FOR PLOTTING OF YPRED in dataset
  
  if(add_timepoints==TRUE){
    
    #loop over ids
    for (row in ids) {
      # New rows to be inserted in dataframe
      new_rows             <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
      colnames(new_rows)   <- colnames(dataset)
      new_rows$ID          <- row
      new_rows$XORIGINAL   <- t
      new_rows$XORIGINALUNIT  <- unique(dataset$XORIGINALUNIT[which((dataset$ID==row)&(is.na(dataset[["DOSEORIGINAL"]])))])
      new_rows$YPREDUNIT   <- unique(dataset$YPREDUNIT[which((dataset$ID==row)&(is.na(dataset[["DOSEORIGINAL"]])))])
      
      # Identify runid for id
      runid <-  unique(dataset$RUNID[which(dataset$ID==row)])
      
      # Identify rownumbers of timepoints from results_wide[[runid]]$time that match timepoints of t
      is_timepoints <- match(t,results_wide[[runid]]$time)
      
      # Transfer predictions for specific compartment at these timepoints to new_rows columns YPRED
      if(row=="Rifampicin_human_liver_600mg_po_MD")       {new_rows$YPRED         <- results_wide[[runid]]$RIFcen[is_timepoints]      }
      
      #add new rows to dataset
      dataset <- rbind(dataset, new_rows)
    }}
  
    
    
  
  return(dataset)
  
}



