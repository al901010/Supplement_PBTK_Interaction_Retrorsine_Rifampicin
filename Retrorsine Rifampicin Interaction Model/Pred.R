

###############################################################################
# PREDICTION FOR ESTIMATION INCLUDES ONLY OBSERVED TIMEPOINTS
###############################################################################

PredRxODE <- function(t,
                      runtable) {
  
  
  
  #--- SOLVE ODE SYSTEM FOR RUNTABLE AND PROCESS SOLUTION -----------------------
  
  
  # Initialise empty lists of length = nrow(runtable)
  results_ev   <- vector("list",nrow(runtable))
  results_wide <- vector("list",nrow(runtable))
  results_sum  <- vector("list",nrow(runtable))
  
  
  # FOR LOOP THROUGH ROWS OF RUNTABLE
  for(row in 1:nrow(runtable)) {
    
    # Speciestype related data
    speciestype <- runtable[row,"SPECIESTYPE"]
    theta       <- GetParameterVector(drug="Retrorsine",speciestype,method_partcoeff="RodgersAndRowland")
    
    # Event table
    RETdose   <- runtable[row,"RETDOSE"]/10^6/351.4*get(paste("theta",speciestype,sep="_"))["bw"] #ugperkgbw to mol
    RETdoseno <- runtable[row,"RETDOSENO"]
    RETstart  <- runtable[row,"RETSTART"]
    
    RIFdose   <- runtable[row,"RIFDOSE"]/10^3/822.9 #mg to mol
    RIFdoseno <- runtable[row,"RIFDOSENO"]
    RIFstart  <- runtable[row,"RIFSTART"]
    
    
    ev <- eventTable(amount.units="mol", time.units="hours") %>%
      add.dosing(dose=RETdose, nbr.doses=RETdoseno, start.time=RETstart,dosing.interval=24,dosing.to="RETlum") %>%
      add.dosing(dose=RIFdose, nbr.doses=RIFdoseno, start.time=RIFstart,dosing.interval=24,dosing.to = "RIFdeplivc") %>%
      add.dosing(dose=RIFdose, nbr.doses=RIFdoseno, start.time=RIFstart,dosing.interval=24,dosing.to = "RIFdepgut") %>%
      add.sampling(t);
    
    CYP3A4baseline <- as.numeric(get(paste("theta",speciestype,sep="_"))["cCYP3A40"]*get(paste("theta",speciestype,sep="_"))["Vlivc"])
    CYP3A4gutbaseline <- as.numeric(get(paste("theta",speciestype,sep="_"))["cCYP3A40gut"]*get(paste("theta",speciestype,sep="_"))["Vgut"])
    
    inits <- c(RETlum=0,
               RETper=0,
               RETven=0,
               RETart=0,
               RETadi=0,
               RETbon=0,
               RETbra=0,
               RETgut=0,
               REThea=0,
               RETkid=0,
               RETlivvi=0,
               RETlivc=0,
               RETlun=0,
               RETmus=0,
               RETski=0,
               RETspl=0,
               RETMetlivnonCYP3A4=0,
               RETMetlivCYP3A4=0,
               RETMetgutnonCYP3A4=0,
               RETMetgutCYP3A4=0,
               RETuri=0,
               RETbil=0,
               CYP3A4livc=CYP3A4baseline,
               RIFdeplivc=0,
               RIFcenlivc=0,
               RIFperlivc=0,
               RIFbag=0,
               CYP3A4gut=CYP3A4gutbaseline,
               RIFdepgut=0,
               RIFcengut=0,
               RIFpergut=0,
               RIFbaggut=0
               )
    
    # Solve ODE system 
    wide <- mInteraction %>% rxSolve(theta,ev,inits,method="lsoda",atol=1e-12,rtol=1e-10) #,maxsteps=500000,atol=1e-12,rtol=1e-10) 
    
    # Check balance
    sum <- rowSums(wide)-parse_number(as.character(wide$time))
    #- parse_number(as.character(wide$CYP3A4livc)) - parse_number(as.character(wide$CYP3A4gut)) 
    

    
    # Determine RETliv, RETpla, RETMetliv
    wide$RETliv <- wide$RETlivc + wide$RETlivvi
    wide$RETpla <- (wide$RETven + wide$RETart)*(1-get(paste("theta",speciestype,sep="_"))["hct"])/get(paste("theta",speciestype,sep="_"))["BP"]
    wide$RETMetliv <- wide$RETMetlivnonCYP3A4 + wide$RETMetlivCYP3A4
    wide$RETMetgut <- wide$RETMetgutnonCYP3A4 + wide$RETMetgutCYP3A4
    
    
    #optional FOD: Convert units [mol] in [fraction of dose]
    #wide[,2: ncol(wide)] <-  wide[,2: ncol(wide)]*100/RETdose
    
    # Convert units [mol] in [?mol]
    wide[,2: ncol(wide)] <-  wide[,2: ncol(wide)]*10^6


    #Convert units [?mol] to [?mol/L]
    wide$RETlum   #no conversion
    wide$RETper   #no conversion
    wide$RETven   <- wide$RETven/get(paste("theta",speciestype,sep="_"))["Vven"]
    wide$RETart   <- wide$RETart/get(paste("theta",speciestype,sep="_"))["Vart"]
    wide$RETadi   <- wide$RETadi/get(paste("theta",speciestype,sep="_"))["Vadi"]
    wide$RETbon   <- wide$RETbon/get(paste("theta",speciestype,sep="_"))["Vbon"]
    wide$RETbra   <- wide$RETbra/get(paste("theta",speciestype,sep="_"))["Vbra"]
    wide$RETgut   <- wide$RETgut/get(paste("theta",speciestype,sep="_"))["Vgut"]
    wide$REThea   <- wide$REThea/get(paste("theta",speciestype,sep="_"))["Vhea"]
    wide$RETkid   <- wide$RETkid/get(paste("theta",speciestype,sep="_"))["Vkid"]
    wide$RETlivvi <- wide$RETlivvi/get(paste("theta",speciestype,sep="_"))["Vlivvi"]
    wide$RETlivc  <- wide$RETlivc/get(paste("theta",speciestype,sep="_"))["Vlivc"]
    wide$RETlun   <- wide$RETlun/get(paste("theta",speciestype,sep="_"))["Vlun"]
    wide$RETmus   <- wide$RETmus/get(paste("theta",speciestype,sep="_"))["Vmus"]
    wide$RETski   <- wide$RETski/get(paste("theta",speciestype,sep="_"))["Vski"]
    wide$RETspl   <- wide$RETspl/get(paste("theta",speciestype,sep="_"))["Vspl"]

    wide$RETMetlivnonCYP3A4 <- wide$RETMetlivnonCYP3A4 #/get(paste("theta",speciestype,sep="_"))["Vlivc"]
    wide$RETMetlivCYP3A4    <- wide$RETMetlivCYP3A4    #/get(paste("theta",speciestype,sep="_"))["Vlivc"]
    wide$RETMetliv          <- wide$RETMetliv          #/get(paste("theta",speciestype,sep="_"))["Vlivc"]
    wide$RETMetgutnonCYP3A4 <- wide$RETMetgutnonCYP3A4 #/get(paste("theta",speciestype,sep="_"))["Vgut"]
    wide$RETMetgutCYP3A4    <- wide$RETMetgutCYP3A4    #/get(paste("theta",speciestype,sep="_"))["Vgut"]
    wide$RETMetgut          <- wide$RETMetgut          #/get(paste("theta",speciestype,sep="_"))["Vgut"]

    wide$RETuri      #no conversion
    wide$RETbil      #no conversion
    wide$RETliv    <- wide$RETliv/get(paste("theta",speciestype,sep="_"))["Vliv"]
    wide$RETpla    <- wide$RETpla/get(paste("theta",speciestype,sep="_"))["Vpla"]

    wide$CYP3A4livc   <- wide$CYP3A4livc/get(paste("theta",speciestype,sep="_"))["Vlivc"]
    wide$RIFdeplivc  #no conversion
    wide$RIFcenlivc  <- wide$RIFcenlivc/get(paste("theta",speciestype,sep="_"))["Vlivc"]
    wide$RIFperlivc  #no conversion
    wide$RIFbag      #no conversion

    wide$CYP3A4gut   <- wide$CYP3A4gut/get(paste("theta",speciestype,sep="_"))["Vgut"]
    wide$RIFdepgut   #no conversion
    wide$RIFcengut   <- wide$RIFcengut/get(paste("theta",speciestype,sep="_"))["Vgut"]
    wide$RIFpergut    #no conversion
    wide$RIFbaggut      #no conversion

    # Save results in list
    results_ev[[row]]   <- ev
    results_sum[[row]]  <- sum
    results_wide[[row]] <- wide
  }
  
  
  # # View dosing in event tables
  # print("Dosing")
  # print(results_ev[[1]]$get.dosing())
  # print(results_ev[[2]]$get.dosing())
  # print(results_ev[[3]]$get.dosing())
  # print(results_ev[[4]]$get.dosing())
  # # print(results_ev[[5]]$get.dosing())
  # 
  # # View head of balance check and results in wide format
  # print("Balance check")
  # print(head(results_sum))
  # print("Predictions wide-format")
  # print(head(results_wide))
  
  return(results_wide)
  
}



