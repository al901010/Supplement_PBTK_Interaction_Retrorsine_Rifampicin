
GetParameterVector <- function (drug, 
                                speciestype, 
                                method_partcoeff) {
  

  #############################################################################
  #    SPECIES-SPECIFIC DATA                                                    
  #############################################################################
  
  SpeciesData <- GetSpeciesData(speciestype)
  
  # BODYWEIGHT ################################################################
  bw      <- SpeciesData$BW$bodyweight # body weight [kg]
  
  # Hematorcrit ###############################################################
  hct     <- SpeciesData$hct
  
  # VOLUMES ###################################################################
  Vvas    <- SpeciesData$fV$fV_vas_tot*SpeciesData$OW_and_V$V                   # Volume of tissues that is vascular space [L]
  SumVvas <- sum(Vvas["adi"],Vvas["bon"],Vvas["bra"],Vvas["gut"],Vvas["hea"], 
                 Vvas["kid"],Vvas["liv"],Vvas["lun"],Vvas["mus"],Vvas["ski"],
                 Vvas["spl"])                                                   # Sum of vascular volume across all tissues [L]
  Vblo    <- as.vector(SpeciesData$OW_and_V$V["blo"]-SumVvas)                   # Blood volume reduced by vascular volume of all tissues [L]
  
  Vven    <- as.vector(SpeciesData$fVvenart["ven"]*Vblo)   # Volume of venous blood [L]
  Vart    <- as.vector(SpeciesData$fVvenart["art"]*Vblo)   # Volume of arterial blood [L]
  
  Vpla    <- Vblo*(1-SpeciesData$hct)
  Very    <- Vblo - Vpla
  
  Vadi    <- as.vector(SpeciesData$OW_and_V$V["adi"]) # Tissue volumes V [L]
  Vbon    <- as.vector(SpeciesData$OW_and_V$V["bon"])
  Vbra    <- as.vector(SpeciesData$OW_and_V$V["bra"])
  Vgut    <- as.vector(SpeciesData$OW_and_V$V["gut"])
  Vhea    <- as.vector(SpeciesData$OW_and_V$V["hea"])
  Vkid    <- as.vector(SpeciesData$OW_and_V$V["kid"])
  Vliv    <- as.vector(SpeciesData$OW_and_V$V["liv"])
  Vlun    <- as.vector(SpeciesData$OW_and_V$V["lun"])
  Vmus    <- as.vector(SpeciesData$OW_and_V$V["mus"])
  Vski    <- as.vector(SpeciesData$OW_and_V$V["ski"])
  Vspl    <- as.vector(SpeciesData$OW_and_V$V["spl"]) 
  Vlivv   <- as.vector(SpeciesData$fV$fV_vas_tot["liv"]*SpeciesData$OW_and_V$V["liv"]) # Liver vascular volume [L]
  Vlivi   <- as.vector(SpeciesData$fV$fV_int_tot["liv"]*SpeciesData$OW_and_V$V["liv"]) # Liver interstitital volume  [L]
  Vlivvi  <- Vlivv + Vlivi                                                             # Liver volume of lumped vascular and interstitial space [L]
  Vlivc   <- as.vector(SpeciesData$fV$fV_cel_tot["liv"]*SpeciesData$OW_and_V$V["liv"]) # Liver intracellular volume [L]
  V       <- c(Vblo=Vblo,Vven=Vven, Vart=Vart, Vpla=Vpla,Very=Very,Vadi=Vadi,
               Vbon=Vbon,Vbra=Vbra,Vgut=Vgut,Vhea=Vhea,Vkid=Vkid,Vliv=Vliv,
               Vlun=Vlun,Vmus=Vmus,Vski=Vski,Vspl=Vspl,Vlivv=Vlivv,Vlivi=Vlivi,
               Vlivvi=Vlivvi,Vlivc=Vlivc)
  
  # BLOOD FLOWS ###############################################################
  Qadi    <- as.vector(SpeciesData$co_and_Q$Q["adi"]) # Tisse blood flows Q [L/min] and cardiac output co [L/min]
  Qbon    <- as.vector(SpeciesData$co_and_Q$Q["bon"])
  Qbra    <- as.vector(SpeciesData$co_and_Q$Q["bra"])
  Qgut    <- as.vector(SpeciesData$co_and_Q$Q["gut"])
  Qhea    <- as.vector(SpeciesData$co_and_Q$Q["hea"])
  Qkid    <- as.vector(SpeciesData$co_and_Q$Q["kid"])
  Qliv    <- as.vector(SpeciesData$co_and_Q$Q["liv"])
  Qmus    <- as.vector(SpeciesData$co_and_Q$Q["mus"])
  Qski    <- as.vector(SpeciesData$co_and_Q$Q["ski"])
  Qspl    <- as.vector(SpeciesData$co_and_Q$Q["spl"])
  Qc      <- sum(Qadi,Qbon,Qbra,Qhea,Qkid,Qliv,Qmus,Qski)     #update cardiac output: sum of tissue blood blows that go into venous blood
  Q       <- c(Qadi=Qadi, Qbon=Qbon, Qbra=Qbra, Qgut=Qgut, 
               Qhea=Qhea, Qkid=Qkid, Qliv=Qliv, Qmus=Qmus, 
               Qski=Qski, Qspl=Qspl, Qc=Qc)
  
  # GLOMERULAR FILTRATION RATE ################################################
  
  GFR     <- as.vector(SpeciesData$GFR) # Glomerular filtration rate [L/h]
  
  # GLUTATHIONE ###############################################################
  GSH0liv   <- as.vector(SpeciesData$GSH0liv)   # GSH amount in liver cellular compartment at steady state [mol]
  kdegGSH   <- as.vector(SpeciesData$kdegGSH)   # 1st-order degradation rate constant for glutathione [1/h]
  ksynGSH   <- as.vector(GSH0liv/Vlivc*kdegGSH)
  # ksynGSH <- 0
  # kdegGSH   <- as.vector(ksynGSH*Vlivc/GSH0liv)
  # kdegGSH <-0
  GSH       <- c(GSH0liv=GSH0liv,kdegGSH=kdegGSH, ksynGSH=ksynGSH)
  
  # CYTOCHROMEP450 (CYP) ######################################################
  # human: CYP3A4
  CYPss   <- NA                                     # CYP concentration in liver cellular compartment at steady-state [mol/L]
  ksynCYP <- NA                                     # 0th-order synthesis rate constant for CYP3A4 [mol/L/h]
  kdegCYP <- NA                                     # 1st-order degradation rate constant for CYP3A4 [1/h]
  CYP     <- c(CYPss=CYPss, ksynCYP=ksynCYP, kdegCYP=kdegCYP)

  
  #############################################################################
  #    DRUG-SPECIFIC PARAMETERS                                                 
  #############################################################################
  DrugData <- GetDrugData(drug,speciestype)
  

  # CLEARANCE RATES ###########################################################
  CLuri        <- DrugData$CLuri                                 # Renal clearance [L/h]
  CLcanef      <- DrugData$CLcanef                               # Canicular efflux clearance [L/h]
  CLactin      <- DrugData$CLactin_perkgliver*Vliv               # Active transport into liver cell [L/h]
  CLactef      <- DrugData$CLactef                               # Active transport out of liver cell [L/h]
  PSdiff       <- DrugData$PSdiff_perkgliver*Vliv                # Passive diffusion into liver cell [L/h]
  Km           <- DrugData$Km                                    # Michaelis Menten constant [M]
  Vmaxliv      <- DrugData$Vmax_perkgliver*Vliv                  # Maximum reaction velocity liver [M*L/h]
  Fsi          <- 0.882*(0.72/(0.72+0.41))                       # Volume fraction that is intracellular of
                                                                 # the small intestine to confine retrorsine 
                                                                 # metabolism in the gut tissue
  fmcyp3a4     <- 0.527                                          # fraction of retrorsine that is metabolized
                                                                 # by CYP3A4

  CL           <- c(CLuri=CLuri, CLcanef=CLcanef,CLactin=CLactin,CLactef=CLactef,
                    PSdiff=PSdiff,Vmaxliv=Vmaxliv, Km=Km, Fsi=Fsi, fmcyp3a4=fmcyp3a4)
  

  # PERITONEAL ABSORPTION RATE CONSTANT [1/h] 
  kip <- DrugData$kip
  
  #############################################################################
  #    PARTITION COEFFICIENTS AND RELATED PARAMETERS                                                   
  #############################################################################
  
  RodgersAndRowland <- UseMethodRodgersAndRowland(DrugData,SpeciesData) 
  Schmitt <- UseMethodSchmitt(DrugData,SpeciesData)

  # Check in dataframe ########################################################
  KCheck         <- data.frame(RodgersAndRowland$Kpu, Schmitt$Kpu)
  KCheck$RoverS  <- KCheck$RodgersAndRowland.Kpu/KCheck$Schmitt.Kpu
  
  # SWITCH METHOD #############################################################
  switch(method_partcoeff,
         "RodgersAndRowland"={
           Kblo <- as.vector(RodgersAndRowland$Kpu["blo"])
           Kpla <- as.vector(RodgersAndRowland$Kpu["pla"])
           Kery <- as.vector(RodgersAndRowland$Kpu["ery"])
           Kadi <- as.vector(RodgersAndRowland$Kpu["adi"])
           Kbon <- as.vector(RodgersAndRowland$Kpu["bon"])
           Kbra <- as.vector(RodgersAndRowland$Kpu["bra"])
           Kgut <- as.vector(RodgersAndRowland$Kpu["gut"])
           Khea <- as.vector(RodgersAndRowland$Kpu["hea"])
           Kkid <- as.vector(RodgersAndRowland$Kpu["kid"])
           Kliv <- as.vector(RodgersAndRowland$Kpu["liv"])
           Klun <- as.vector(RodgersAndRowland$Kpu["lun"])
           Kmus <- as.vector(RodgersAndRowland$Kpu["mus"])
           Kski <- as.vector(RodgersAndRowland$Kpu["ski"])
           Kspl <- as.vector(RodgersAndRowland$Kpu["spl"])
           fnblo <- as.vector(RodgersAndRowland$fn["blo"])
           fnpla <- as.vector(RodgersAndRowland$fn["pla"])
           fnery <- as.vector(RodgersAndRowland$fn["ery"])
           fnadi <- as.vector(RodgersAndRowland$fn["adi"])
           fnbon <- as.vector(RodgersAndRowland$fn["bon"])
           fnbra <- as.vector(RodgersAndRowland$fn["bra"])
           fngut <- as.vector(RodgersAndRowland$fn["gut"])
           fnhea <- as.vector(RodgersAndRowland$fn["hea"])
           fnkid <- as.vector(RodgersAndRowland$fn["kid"])
           fnliv <- as.vector(RodgersAndRowland$fn["liv"])
           fnlun <- as.vector(RodgersAndRowland$fn["lun"])
           fnmus <- as.vector(RodgersAndRowland$fn["mus"])
           fnski <- as.vector(RodgersAndRowland$fn["ski"])
           fnspl <- as.vector(RodgersAndRowland$fn["spl"])
           },
         "Schmitt"          ={
           Kblo <- as.vector(Schmitt$Kpu["blo"])
           Kpla <- as.vector(Schmitt$Kpu["pla"])
           Kery <- as.vector(Schmitt$Kpu["ery"])
           Kadi <- as.vector(Schmitt$Kpu["adi"])
           Kbon <- as.vector(Schmitt$Kpu["bon"])
           Kbra <- as.vector(Schmitt$Kpu["bra"])
           Kgut <- as.vector(Schmitt$Kpu["gut"])
           Khea <- as.vector(Schmitt$Kpu["hea"])
           Kkid <- as.vector(Schmitt$Kpu["kid"])
           Kliv <- as.vector(Schmitt$Kpu["liv"])
           Klun <- as.vector(Schmitt$Kpu["lun"])
           Kmus <- as.vector(Schmitt$Kpu["mus"])
           Kski <- as.vector(Schmitt$Kpu["ski"])
           Kspl <- as.vector(Schmitt$Kpu["spl"])
           fnblo <- as.vector(Schmitt$fn["blo"])
           fnpla <- as.vector(Schmitt$fn["pla"])
           fnery <- as.vector(Schmitt$fn["ery"])
           fnadi <- as.vector(Schmitt$fn["adi"])
           fnbon <- as.vector(Schmitt$fn["bon"])
           fnbra <- as.vector(Schmitt$fn["bra"])
           fngut <- as.vector(Schmitt$fn["gut"])
           fnhea <- as.vector(Schmitt$fn["hea"])
           fnkid <- as.vector(Schmitt$fn["kid"])
           fnliv <- as.vector(Schmitt$fn["liv"])
           fnlun <- as.vector(Schmitt$fn["lun"])
           fnmus <- as.vector(Schmitt$fn["mus"])
           fnski <- as.vector(Schmitt$fn["ski"])
           fnspl <- as.vector(Schmitt$fn["spl"])
         }
  )
  
  
  # BP ratio [-] 
  BP           <- DrugData$BP  
  # Fraction unbound in plasma 
  fuP          <- DrugData$fuP 
  # Fraction unbound in interstitial space (assumed to be identical to fuP)  
  fuInt        <- fuP                    
  
  # Fraction unbound in liver cellular space (assumed to be identical to fuP) 
  fuliv <- fuP
  
  #############################################################################
  #DATA:   Kvasvi (liver vascular-to-lumped compartment partition coefficient)                    
  #UNIT:   []                                                                 
  #SOURCE: Schweinoch 2014                                                       
  #NOTE:   vi: lumped compartments vascular and interstitial space       
  #############################################################################  
  
  Kvasvi <- BP/fuP/((Vlivv/Vlivvi*BP/fuP)+(Vlivi/Vlivvi*1/fuInt))
  
  #############################################################################
  #DATA:   Kintuvi (liver unbound interstital-to-lumped compartment partition 
  #                 coefficient)                    
  #UNIT:   []                                                                 
  #SOURCE: Schweinoch 2014                                                       
  #NOTE:   vi: lumped compartments vascular and interstitial space       
  #############################################################################  
  
  Kintuvi <- fuInt/((Vlivv/Vlivvi*BP/fuP)+(Vlivi/Vlivvi*1/fuInt))
  
  
  
  partitioning  <- c(Kblo=Kblo,Kpla=Kpla,Kery=Kery,Kadi=Kadi,Kbon=Kbon,Kbra=Kbra,
                     Kgut=Kgut,Khea=Khea,Kkid=Kkid,Kliv=Kliv,Klun=Klun,Kmus=Kmus,
                     Kski=Kski,Kspl=Kspl,
                     Kvasvi=Kvasvi,Kintuvi=Kintuvi,
                     fnblo=fnblo,fnpla=fnpla,fnery=fnery,fnadi=fnadi,fnbon=fnbon,
                     fnbra=fnbra,fngut=fngut,fnhea=fnhea,fnkid=fnkid,fnliv=fnliv,
                     fnlun=fnlun,fnmus=fnmus,fnski=fnski,fnspl=fnspl,
                     BP=BP,fuP=fuP,fuInt=fuInt,fuliv=fuliv)
  

    
  #############################################################################
  # ORAL ABSOPRITON RELATED PARAMETERS
  #############################################################################
  
  #Intestinal absorption rate constant [1/h]
  ka <- DrugData$ka
  
  #Fraction absorbed from gut tissue [-], according to Skolnik et al. 2010
  Fa <- (0.01+(1-0.01))/(1+(exp((-5.74-log10(DrugData$Papp))/0.39)))
  
  #############################################################################
  # RIFAMPICIN MODEL PARAMETERS                                                                       
  #############################################################################

  kaRIFlivc  <- 0.3565737  
  k12RIFlivc <- 0.2785765  
  k21RIFlivc <- 0.1550562 
  keRIFlivc  <- 0.5483147 
  kaRIFgut   <- 0.1870029  
  k12RIFgut  <- 0.7159209  
  k21RIFgut  <- 0.02581702 
  keRIFgut   <- 5.765528 
  EC50       <- 0.34/10^6     #µmol/L to mol/L
  Emax       <- 9
  kdeg       <- log(2)/36     #1/h
  kdeggut    <- log(2)/23     #1/h
  cCYP3A40   <- 6.48/10^6     #µmol/L to mol/L
  cCYP3A40gut<- 0.356/10^6    #µmol/L to mol/L
  kcatCYP3A4 <- 5.59*60       #1/min to 1/h
  kcatCYP3A4gut <- 13.0*60    #1/min to 1/h
  VmaxlivnonCYP3A4 <- Vmaxliv #mol/h
  Ki         <- 18.5/10^6     #µmol/L to mol/L
  
  RIF <- c(kaRIFlivc=kaRIFlivc,
           k12RIFlivc=k12RIFlivc,
           k21RIFlivc=k21RIFlivc,
           keRIFlivc=keRIFlivc,
           kaRIFgut=kaRIFgut,
           k12RIFgut=k12RIFgut,
           k21RIFgut=k21RIFgut,
           keRIFgut=keRIFgut,
           EC50=EC50,
           Emax=Emax,
           kdeg=kdeg,
           kdeggut=kdeggut,
           cCYP3A40=cCYP3A40,
           cCYP3A40gut=cCYP3A40gut,
           kcatCYP3A4=kcatCYP3A4,
           kcatCYP3A4gut=kcatCYP3A4gut,
           VmaxlivnonCYP3A4=VmaxlivnonCYP3A4,
           Ki=Ki)
  
  #############################################################################
  # THETA                                                                       
  #############################################################################
  
  theta <- c(bw=bw,hct=hct,V,Q,GFR=GFR,GSH,CYP,CL,kip=kip,partitioning,
             Fa=Fa,ka=ka,RIF)
  
  
  return(theta)
}

