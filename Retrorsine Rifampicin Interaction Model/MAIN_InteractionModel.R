###############################################################################
# AUTHOR:  Anja Lehmann                                                       
# VERSION: 22 January 2024  
# CITATION: The functions GetSpeciesData() and UseMethodRodgersAndRowland() 
# are based on a MATBLAB script originally written by                             
# W. Huisinga, A. Solms, L. Fronton, S. Pilari, Modeling Interindividual       
# Variability in Physiologically Based Pharmacokinetics and Its Link to       
# Mechanistic Covariate Modeling, CPT: Pharmacometrics & Systems Pharmacology 
# (2012) 1, e5; doi:10.1038/psp.2012.3                                        
###############################################################################

#--- INITIALIZATION 

# Remove all objects from current workspace
rm(list=ls())

# Reset graphics
graphics.off()
cat("\014")

# Set working directory to current path
PATH <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(PATH)

# Load dependencies
source("DrugDatabase.R")
source("SpeciesDatabase.R")
source("ParameterVector.R")
source("PartitionCoefficients.R")
source("Pred.R")


# Load packages
library(readr)
library(crayon)
library(rxode2)
library(deSolve)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(dfoptim)
library(cowplot)
library(profvis)
library(FME)
library(nloptr)
library(matlib)
library(bayestestR)
library(NonCompart)


# Define color-blind friendly palette (Okabe and Ito 2002)
palette_OkabeIto <- c("#E69F00", "#000000", "#009E73", "#56B4E9", "#D55E00", 
                      "#CC79A7", "#999999", "#F0E442") 


###############################################################################
# PART: INITIALIZATION OF PREDICTION                                              
###############################################################################

#--- DEFINE SPECIESTYPES AND DOSING REGIMEN

# speciestype is a string defining species, sex (only human 15 and 35 years), and age, e.g. "rat10weeks"
# Speciestypes, for which data is available in SpeciesData.R: 
# "rat10weeks",   "rat70weeks", 
# "mouse9weeks",  "mouse70weeks",
# "humannewborn", "human1year", "human5years",  "human10years", 
# "humanm15years","humanf15years","humanm35years","humanf35years"


#--- CREATE RUNTABLE AND RUNID 

runtable <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("RETDOSE", 
                                                                "RETDOSEUNIT", 
                                                                "RETDOSENO",
                                                                "RETSTART",
                                                                "RIFDOSE",
                                                                "RIFDOSEUNIT",
                                                                "RIFDOSENO",
                                                                "RIFSTART",
                                                                "SPECIESTYPE",
                                                                "ADM",
                                                                "RUNID"))

runtable[1,] <- c(0.019,"ugperkgbw",35,24*0,    0,"mg",14, 24*4,"humanm35years","po",1)
runtable[2,] <- c(0.019,"ugperkgbw",35,24*0,  600,"mg",14, 24*4,"humanm35years","po",2) 


runtable$RETDOSE   <- as.numeric(runtable$RETDOSE)
runtable$RETDOSENO <- as.numeric(runtable$RETDOSENO)
runtable$RETSTART  <- as.numeric(runtable$RETSTART)
runtable$RIFDOSE   <- as.numeric(runtable$RIFDOSE)
runtable$RIFDOSENO <- as.numeric(runtable$RIFDOSENO)
runtable$RIFSTART  <- as.numeric(runtable$RIFSTART)
runtable$RUNID     <- as.numeric(runtable$RUNID)


#--- LOAD THETA FOR EACH SPECIESTYPE and PRINT PLAUSIBILITY CHECKS 

speciestypes <- unique(runtable$SPECIESTYPE)

for(i in seq_along(speciestypes)) {
  speciestype       <- speciestypes[i]
  drug <- "Retrorsine"
  
  theta             <- GetParameterVector(drug,speciestype,method_partcoeff="RodgersAndRowland")
  # "Retrorsine" has to be included in DrugDatabase.R
  # method_partcoeff choose between "RodgersAndRowland" / "Schmitt", NOTE: also change in Pred.R
  
  SpeciesData       <- GetSpeciesData(speciestype)
  DrugData          <- GetDrugData(drug="Retrorsine",speciestype)  
  RodgersAndRowland <- UseMethodRodgersAndRowland(DrugData,SpeciesData)
  Schmitt           <- UseMethodSchmitt(DrugData,SpeciesData)
  
  #Assign variable names xxx_speciestype
  assign(paste("theta",speciestype,sep="_"), theta)
 
  #Print PlausibilityCheck
  print("Plausibilitycheck")
  cat(cyan(speciestype),
      "\n RodgersAndRowland$PlausibilityCheckKA_PR:",  
      yellow(RodgersAndRowland$PlausibilityCheckKA_PR),
      "\n RodgersAndRowland$PlausibilityCheckKA_AP:",  
      yellow(RodgersAndRowland$PlausibilityCheckKA_AP),
      "\n Schmitt$PlausibilityCheckfuC:",            
      yellow(Schmitt$PlausibilityCheckfuC),"\n")
  

}

#  if (BP<(1-hct)) {cat(red("Blood-to-plasma ratio smaller than (1-hematocrit)!"))}

#--- DEFINE ODE SYSTEM FOR RxODE PACKAGE 

mInteraction <-RxODE({
  d/dt (RETlum)            = -ka*Fa*RETlum;
  d/dt (RETper)            = -kip*RETper;
  d/dt (RETven)            =  Qadi*RETadi/Vadi/Kadi
  +Qbon*RETbon/Vbon/Kbon
  +Qbra*RETbra/Vbra/Kbra
  +Qhea*REThea/Vhea/Khea
  +Qkid*RETkid/Vkid/Kkid
  +Qliv*Kvasvi*RETlivvi/Vlivvi
  +Qmus*RETmus/Vmus/Kmus
  +Qski*RETski/Vski/Kski
  -Qc*RETven/Vven;
  d/dt (RETart)            =  Qc*RETlun/Vlun/Klun
  -Qc*RETart/Vart;
  d/dt (RETadi)            =  Qadi*RETart/Vart
  -Qadi*RETadi/Vadi/Kadi;
  d/dt (RETbon)            =  Qbon*RETart/Vart
  -Qbon*RETbon/Vbon/Kbon;
  d/dt (RETbra)            =  Qbra*RETart/Vart
  -Qbra*RETbra/Vbra/Kbra;
  d/dt (RETgut)            =  Qgut*RETart/Vart
  -Qgut*RETgut/Vgut/Kgut
  -Fsi*(1-fmcyp3a4)*((VmaxlivnonCYP3A4/10)/(Km+fuP*(RETgut/Vgut)))*fuP*RETgut/Vgut
  -Fsi*fmcyp3a4*((kcatCYP3A4gut*CYP3A4gut)/(Km*(1+(RIFcengut/Vgut/Ki))+fuP*(RETgut/Vgut)))*fuP*RETgut/Vgut
  +ka*Fa*RETlum;
  d/dt (REThea)            =  Qhea*RETart/Vart
  -Qhea*REThea/Vhea/Khea;
  d/dt (RETkid)            =  Qkid*RETart/Vart
  -Qkid*RETkid/Vkid/Kkid
  -fuP*GFR*RETkid/Vkid;
  d/dt (RETlivvi)          =  kip*RETper
  +(Qliv-Qspl-Qgut)*RETart/Vart
  +Qspl*RETspl/Vspl/Kspl
  +Qgut*RETgut/Vgut/Kgut
  -Qliv*Kvasvi*RETlivvi/Vlivvi
  +CLactef*fuliv*RETlivc/Vlivc
  +PSdiff*fnliv/fnpla*fuliv*RETlivc/Vlivc
  -CLactin*Kintuvi*RETlivvi/Vlivvi
  -PSdiff*Kintuvi*RETlivvi/Vlivvi;
  d/dt (RETlivc)           =  CLactin*Kintuvi*RETlivvi/Vlivvi
  +PSdiff*Kintuvi*RETlivvi/Vlivvi
  -CLactef*fuliv*RETlivc/Vlivc
  -PSdiff*fnliv/fnpla*fuliv*RETlivc/Vlivc
  -CLcanef*fuliv*RETlivc/Vlivc
  -(1-fmcyp3a4)*(VmaxlivnonCYP3A4/(Km+fuliv*(RETlivc/Vlivc)))*fuliv*RETlivc/Vlivc
  -fmcyp3a4*(kcatCYP3A4*CYP3A4livc/(Km*(1+(RIFcenlivc/Vlivc/Ki))+fuliv*(RETlivc/Vlivc)))*fuliv*RETlivc/Vlivc;
  d/dt (RETlun)            =  Qc*RETven/Vven
  -Qc*RETlun/Vlun/Klun;
  d/dt (RETmus)            =  Qmus*RETart/Vart
  -Qmus*RETmus/Vmus/Kmus;
  d/dt (RETski)            =  Qski*RETart/Vart
  -Qski*RETski/Vski/Kski;
  d/dt (RETspl)            =  Qspl*RETart/Vart
  -Qspl*RETspl/Vspl/Kspl;
  d/dt (RETMetlivnonCYP3A4)= (1-fmcyp3a4)*(VmaxlivnonCYP3A4/(Km+fuliv*(RETlivc/Vlivc)))*fuliv*RETlivc/Vlivc;
  d/dt (RETMetlivCYP3A4)   = fmcyp3a4*(kcatCYP3A4*CYP3A4livc/(Km*(1+(RIFcenlivc/Vlivc/Ki))+fuliv*(RETlivc/Vlivc)))*fuliv*RETlivc/Vlivc;
  d/dt (RETMetgutnonCYP3A4)= Fsi*(1-fmcyp3a4)*((VmaxlivnonCYP3A4/10)/(Km+fuP*(RETgut/Vgut)))*fuP*RETgut/Vgut;
  d/dt (RETMetgutCYP3A4)   = Fsi*fmcyp3a4*(kcatCYP3A4gut*CYP3A4gut/(Km*(1+(RIFcengut/Vgut/Ki))+fuP*(RETgut/Vgut)))*fuP*RETgut/Vgut;
  d/dt (RETuri)            =  fuP*GFR*RETkid/Vkid;
  d/dt (RETbil)            =  CLcanef*fuliv*RETlivc/Vlivc;
  
  
  
  
  d/dt (CYP3A4livc)        =  kdeg*cCYP3A40*Vlivc*(1+((Emax*RIFcenlivc/Vlivc)/(EC50+RIFcenlivc/Vlivc))) 
  -kdeg*CYP3A4livc;
  d/dt (CYP3A4gut)         =  kdeggut*cCYP3A40gut*Vgut*(1+((Emax*RIFcengut/Vgut)/(EC50+RIFcengut/Vgut)))
  -kdeggut*CYP3A4gut;
  d/dt (RIFdeplivc)        = -kaRIFlivc*RIFdeplivc;
  d/dt (RIFcenlivc)        =  kaRIFlivc*RIFdeplivc 
  +k21RIFlivc*RIFperlivc 
  -k12RIFlivc*RIFcenlivc 
  -keRIFlivc*RIFcenlivc;
  d/dt (RIFperlivc)        =  k12RIFlivc*RIFcenlivc 
  -k21RIFlivc*RIFperlivc;
  d/dt (RIFbag)            =  keRIFlivc*RIFcenlivc;
  d/dt (RIFdepgut)         = -kaRIFgut*RIFdepgut;
  d/dt (RIFcengut)         =  kaRIFgut*RIFdepgut 
  +k21RIFgut*RIFpergut 
  -k12RIFgut*RIFcengut 
  -keRIFgut*RIFcengut;
  d/dt (RIFpergut)         =  k12RIFgut*RIFcengut 
  -k21RIFgut*RIFpergut;
  d/dt (RIFbaggut)         =  keRIFgut*RIFcengut;
  
})




################################################################################
# PREDICTION 
################################################################################

pred_wide <- PredRxODE(t=seq(0,24*35,by=0.01), 
                        runtable=runtable)



###############################################################################
# Calculate CL and F
###############################################################################

#liver weight [g]
Vliv <- theta_humanm35years["Vliv"]*1000

#liver clearance [L/h], extended clearance model of the liver
CLliv <- theta_humanm35years["Qliv"]*theta_humanm35years["fuP"]*theta_humanm35years["CLactin"]*((theta_humanm35years["Vmaxliv"]/theta_humanm35years["Km"])+0)/
  (theta_humanm35years["Qliv"]*((theta_humanm35years["Vmaxliv"]/theta_humanm35years["Km"])+0) + theta_humanm35years["fuP"]*theta_humanm35years["CLactin"]*((theta_humanm35years["Vmaxliv"]/theta_humanm35years["Km"])+0 ))

#liver clearance [mL/min/kg bodyweight]
CLliv_2 <- CLliv*1000/60/theta_humanm35years["bw"]

#liver clearance [L/h], well-stirred model of the liver
CLliv_ws <- (theta_humanm35years["Qliv"]*theta_humanm35years["fuP"]*theta_humanm35years["Vmaxliv"]/theta_humanm35years["Km"])/
  (theta_humanm35years["Qliv"]+theta_humanm35years["fuP"]*theta_humanm35years["Vmaxliv"]/theta_humanm35years["Km"])

#fraction ecaping liver first pass metabolism
Fh <- 1-(CLliv/theta_humanm35years["Qliv"])

#gut clearance [L/h], well-stirred model of the gut
CLgut <- (theta_humanm35years["Qgut"]*theta_humanm35years["fuP"]*0.1*theta_humanm35years["Vmaxliv"]/theta_humanm35years["Km"])/
  (theta_humanm35years["Qgut"]+theta_humanm35years["fuP"]*0.1*theta_humanm35years["Vmaxliv"]/theta_humanm35years["Km"])

#gut clearance [mL/min/kg bodyweight]
CLgut_2 <- CLgut*1000/60/theta_humanm35years["bw"]

#fraction escaping gut first pass metabolism
Fg <- 1-(theta_humanm35years["Fsi"]*CLgut/theta_humanm35years["Qgut"])

#bioavailability  
F <- theta_humanm35years["Fa"]*Fh*Fg

#renal clearance [L/h]
CLuri <- theta_humanm35years["fuP"]*theta_humanm35years["GFR"]

#renal clearance [mL/min/kg bodyweight]
CLuri_2 <- CLuri*1000/60/theta_humanm35years["bw"]

#total clearance [L/h]
CLtot <- CLliv + CLgut + theta_humanm35years["fuP"]*theta_humanm35years["GFR"]

#Fraction of CLuri that is total retrorsine CL
fCLuri <- 100*CLuri/CLtot
#Fraction of CLgut that is total retrorsine CL
fCLgut <- 100*CLgut/CLtot
#Fraction of CLliv that is total retrorsine CL
fCLliv <- 100*CLliv/CLtot


###############################################################################
# PART: Plots                                                            
###############################################################################

# Set working directory to subfolder
setwd(paste0(PATH,"/Figures"))

# Convert wide to long with tidyr package with new key column "compartment" and 
# new value column "value"
pred_long   <- vector("list",nrow(runtable))

for(row in 1:nrow(runtable)) {
pred_long[[row]] <- gather(pred_wide[[row]], 
                           compartment, 
                           value, 
                           2:ncol(pred_wide[[row]]),                           
                           factor_key=TRUE)
}

###############################################################################
# optional: fraction of dose
###############################################################################
# 
# #-- COLOR
# c25 <- c(
#   "dodgerblue2", #"#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", "gold1",
#   "skyblue2",
#   #"#FB9A99", # lt pink
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", #"khaki2",
#   "maroon", "orchid1",
#   #"deeppink1",
#   "blue1", "steelblue4",
#   "darkturquoise", #"darkorange4",
#   "yellow4", "yellow3",
#   "brown", "green1"
# )
# 
# 
# 
# 
#   pred_long[[1]] <- subset(pred_long[[1]],
#                                compartment!="RETbil"
#                               & compartment!="RETven"
#                               & compartment!="RETart"
#                               & compartment!="RETper"
#                               & compartment!="RETliv"
#                               & compartment!="RETMetlivnonCYP3A4"
#                               & compartment!="RETMetlivCYP3A4"
#                               & compartment!="CYP3A4livc"
#                               & compartment!="RIFdeplivc"
#                               & compartment!="RIFcenlivc"
#                               & compartment!="RIFperlivc"
#                               & compartment!="RIFbag"
#                               & compartment!="CYP3A4gut"
#                               & compartment!="RIFdepgut"
#                               & compartment!="RIFcengut"
#                               & compartment!="RIFpergut"
#                               & compartment!="RIFbaggut"
# 
# 
#                            )
# 
# 
#   pred_long[[1]]$compartment <- factor(pred_long[[1]]$compartment ,
#                                           levels = c("RETadi",#
#                                                      "RETbil",
#                                                      "RETbon",
#                                                      "RETbra",
#                                                      "RETgut",
#                                                      "REThea",
#                                                      "RETkid",
#                                                      "RETlivc",
#                                                      "RETlivvi",
#                                                      "RETlum",#
#                                                      "RETlun",
#                                                      "RETmus",
#                                                      "RETpla",
#                                                      "RETski",
#                                                      "RETspl",
#                                                      "RETuri",
#                                                      "RETMetgut",
#                                                      "RETMetliv"))
# 
#   dataset          <- pred_long[[1]]
#   xlimits          <- c(0,24)
#   xbreaks          <- seq(0,24,by = 4)
#   ylimits_log      <- c(1e-2,100)
#   pdfname_log      <- paste("RUNID",i,".pdf", sep = "")
#   pngname_log      <- paste("RUNID",i,".png", sep = "")
# 
# 
#   fod_plot <-
#     ggplot(data=dataset,mapping=aes(color=factor(compartment)))+
#     geom_line(aes(x=time, y=value), size=1.5)+
#     labs(x="Time (h)",y="Fraction of dose (%)")+
#     #ggtitle(title)+
#     theme(
#       plot.background   = element_rect(fill = NA),
#       panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
#       panel.background  = element_rect(fill=NA),
#       plot.title        = element_text(size=18, hjust = 0.5),
#       axis.line         = element_line(size = 0.5),
#       axis.text         = element_text(colour = "black",size = 26),
#       axis.title        = element_text(size = 28),
#       # axis.title.y      = element_text(margin=margin(r=15)),
#       # axis.title.x      = element_text(margin=margin(t=12)),
#       axis.ticks        = element_line(colour = "black"),
#       axis.ticks.length = unit(2.5,"mm"),
#       legend.text       = element_text(size=18),
#       legend.title      = element_text(size =18),
#       legend.position   = "none",
#       legend.background = element_rect(fill = NA),
#       aspect.ratio=0.66,
#       legend.key.width = unit(2, "cm"),
#       legend.key=element_blank(),
#       legend.direction = "vertical"
# 
#     )+
#     guides(col = guide_legend(nrow = 6,title.position = "top")) +
#     scale_y_log10(
#       # breaks = scales::trans_breaks("log10", function(x) 10^x),
#       # labels = scales::trans_format("log10", scales::math_format(10^.x)),
#       limits = ylimits_log,
#       labels = function(x) sprintf("%g", x)
#     )+
#     scale_x_continuous(
#       limits = xlimits,
#       breaks = xbreaks)+
#     annotation_logticks(sides="l")+
#     scale_color_manual(values = c25)+
#     guides(color=guide_legend(ncol=3,override.aes = list(size = 3))) #,byrow=TRUE, direction="vertical"
# 
# 
# fod_plot
# ggsave("fractionofdose.png", bg="white", width = 8, height =8)
# 


###############################################################################
# RET                                                      
###############################################################################

#Without RIF
dataset <- subset(pred_long[[1]],compartment=="RETlivc")
#with RIF
dataset2 <- subset(pred_long[[2]],compartment=="RETlivc")


xlimits          <- c(0,35)
xbreaks          <- seq(0,35,by = 7)
ylimits_log      <- c(0.01,0.1)


retrorsine <-
  ggplot(data=dataset)+
  #geom_line(aes(x=time/24, y=value*10^3, color="#D55E00"), size=1.5)+
  geom_line(data=dataset2,aes(x=time/24, y=value*10^3, color="#56B4E9"), 
            linewidth=1.5)+
  
  labs(x="Time (days)",y="Concentration (nmol/L)")+
  #ggtitle(title)+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     linewidth = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(linewidth = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
   scale_y_continuous(
    limits = c(-0.0015,0.025),
    breaks = seq(0,0.025,by=0.005))+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  scale_color_manual(values = c("#85c4c9", "#000000") )+
  annotate("text", x = Inf, y = Inf, label = "Retrorsine liver ",size=11, hjust=1.05, vjust=1.5)+
  geom_segment(aes(x = 4, y = -0.0015, xend = 18, yend = -0.0015),size=8,color="#85c4c9")

retrorsine
ggsave("retrorsine.png", bg="white", width = 9, height =9)



# AUC and Cmax
#w/o RIF
dataset <- subset(pred_long[[2]],compartment=="RETlivc" & time>=24*0 & time<=24*1)

#1 day after 1st dose of RIF
dataset2 <- subset(pred_long[[2]],compartment=="RETlivc" & time>=24*4 & time<=24*5)

#3 days after 1st dose of RIF
dataset3 <- subset(pred_long[[2]],compartment=="RETlivc" & time>=24*6 & time<=24*7)

#14 days after 1st dose of RIF
dataset4 <- subset(pred_long[[2]],compartment=="RETlivc" & time>=24*17 & time<=24*18)

#2 days after last dose of RIF
dataset5 <- subset(pred_long[[2]],compartment=="RETlivc" & time>=24*19 & time<=24*20)

#6 days after last dose of RIF
dataset6 <- subset(pred_long[[2]],compartment=="RETlivc" & time>=24*23 & time<=24*24)

#14 days after last dose of RIF
dataset7 <- subset(pred_long[[2]],compartment=="RETlivc" & time>=24*31 & time<=24*32)


#?mol/L*h
AUC <- AUC(dataset$time, dataset$value, down = "Linear")
AUC2 <- AUC(dataset2$time, dataset2$value, down = "Linear")
AUC3 <- AUC(dataset3$time, dataset3$value, down = "Linear")
AUC4 <- AUC(dataset4$time, dataset4$value, down = "Linear")
AUC5 <- AUC(dataset5$time, dataset5$value, down = "Linear")
AUC6 <- AUC(dataset6$time, dataset6$value, down = "Linear")
AUC7 <- AUC(dataset7$time, dataset7$value, down = "Linear")

# nmol/L*h
max(AUC[,"AUC"])*10^3
max(AUC2[,"AUC"])*10^3
max(AUC3[,"AUC"])*10^3
max(AUC4[,"AUC"])*10^3
max(AUC5[,"AUC"])*10^3
max(AUC6[,"AUC"])*10^3
max(AUC7[,"AUC"])*10^3

# ?mol/L
Cmax <- max(dataset$value)
Cmax2 <- max(dataset2$value)
Cmax3 <- max(dataset3$value)
Cmax4 <- max(dataset4$value)
Cmax5 <- max(dataset5$value)
Cmax6 <- max(dataset6$value)
Cmax7 <- max(dataset7$value)

# nmol/L
Cmax*10^3
Cmax2*10^3
Cmax3*10^3
Cmax4*10^3
Cmax5*10^3
Cmax6*10^3
Cmax7*10^3



###############################################################################
# RET gut                                                       
###############################################################################

#Without RIF
dataset <- subset(pred_long[[1]],compartment=="RETgut")
#with RIF
dataset2 <- subset(pred_long[[2]],compartment=="RETgut")


xlimits          <- c(0,35)
xbreaks          <- seq(0,35,by = 7)
ylimits_log      <- c(0.1,0.3)


retrorsinegut <-
  ggplot(data=dataset)+
  geom_line(data=dataset2,aes(x=time/24, y=value*10^3, color="#56B4E9"), 
            size=1.5)+
   labs(x="Time (days)",y="Concentration (nmol/L)")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_continuous(
    limits = c(-0.003,0.07),
    breaks = seq(0,0.07,by=0.01))+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  scale_color_manual(values = c("#d39c83", "#000000") )+
  annotate("text", x = Inf, y = Inf, label = "Retrorsine gut ",size=11, 
           hjust=1.05, vjust=1.5)+
  geom_segment(aes(x = 4, y = -0.003, xend = 18, yend = -0.003),size=8,
               color="#d39c83")

retrorsinegut
ggsave("retrorsinegut.png", bg="white", width = 9, height =9)




# AUC and Cmax
#w/o RIF
dataset <- subset(pred_long[[2]],compartment=="RETgut" & time>=24*0 & time<=24*1)

#0 days after 1st dose of RIF
dataset2 <- subset(pred_long[[2]],compartment=="RETgut" & time>=24*4 & time<=24*5)

#2 days after 1st dose of RIF
dataset3 <- subset(pred_long[[2]],compartment=="RETgut" & time>=24*6 & time<=24*7)

#12 days after 1st dose of RIF
dataset4 <- subset(pred_long[[2]],compartment=="RETgut" & time>=24*17 & time<=24*18)

#2 days after last dose of RIF
dataset5 <- subset(pred_long[[2]],compartment=="RETgut" & time>=24*19 & time<=24*20)

#6 days after last dose of RIF
dataset6 <- subset(pred_long[[2]],compartment=="RETgut" & time>=24*23 & time<=24*24)

#12 days after last dose of RIF
dataset7 <- subset(pred_long[[2]],compartment=="RETgut" & time>=24*31 & time<=24*32)


#?mol/L*h
AUC <- AUC(dataset$time, dataset$value, down = "Linear")
AUC2 <- AUC(dataset2$time, dataset2$value, down = "Linear")
AUC3 <- AUC(dataset3$time, dataset3$value, down = "Linear")
AUC4 <- AUC(dataset4$time, dataset4$value, down = "Linear")
AUC5 <- AUC(dataset5$time, dataset5$value, down = "Linear")
AUC6 <- AUC(dataset6$time, dataset6$value, down = "Linear")
AUC7 <- AUC(dataset7$time, dataset7$value, down = "Linear")

# nmol/L*h
max(AUC[,"AUC"])*10^3
max(AUC2[,"AUC"])*10^3
max(AUC3[,"AUC"])*10^3
max(AUC4[,"AUC"])*10^3
max(AUC5[,"AUC"])*10^3
max(AUC6[,"AUC"])*10^3
max(AUC7[,"AUC"])*10^3

# ?mol/L
Cmax <- max(dataset$value)
Cmax2 <- max(dataset2$value)
Cmax3 <- max(dataset3$value)
Cmax4 <- max(dataset4$value)
Cmax5 <- max(dataset5$value)
Cmax6 <- max(dataset6$value)
Cmax7 <- max(dataset7$value)

# nmol/L
Cmax*10^3
Cmax2*10^3
Cmax3*10^3
Cmax4*10^3
Cmax5*10^3
Cmax6*10^3
Cmax7*10^3




###############################################################################
# RET pla AUC and Cmax                                                    
###############################################################################

# AUC and Cmax
#w/o RIF
dataset <- subset(pred_long[[2]],compartment=="RETpla" & time>=24*0 & time<=24*1)

#0 days after 1st dose of RIF
dataset2 <- subset(pred_long[[2]],compartment=="RETpla" & time>=24*4 & time<=24*5)

#2 days after 1st dose of RIF
dataset3 <- subset(pred_long[[2]],compartment=="RETpla" & time>=24*6 & time<=24*7)

#12 days after 1st dose of RIF
dataset4 <- subset(pred_long[[2]],compartment=="RETpla" & time>=24*17 & time<=24*18)

#2 days after last dose of RIF
dataset5 <- subset(pred_long[[2]],compartment=="RETpla" & time>=24*19 & time<=24*20)

#6 days after last dose of RIF
dataset6 <- subset(pred_long[[2]],compartment=="RETpla" & time>=24*23 & time<=24*24)

#12 days after last dose of RIF
dataset7 <- subset(pred_long[[2]],compartment=="RETpla" & time>=24*31 & time<=24*32)


#µmol/L*h
AUC <- AUC(dataset$time, dataset$value, down = "Linear")
AUC2 <- AUC(dataset2$time, dataset2$value, down = "Linear")
AUC3 <- AUC(dataset3$time, dataset3$value, down = "Linear")
AUC4 <- AUC(dataset4$time, dataset4$value, down = "Linear")
AUC5 <- AUC(dataset5$time, dataset5$value, down = "Linear")
AUC6 <- AUC(dataset6$time, dataset6$value, down = "Linear")
AUC7 <- AUC(dataset7$time, dataset7$value, down = "Linear")

# nmol/L*h
max(AUC[,"AUC"])*10^3
max(AUC2[,"AUC"])*10^3
max(AUC3[,"AUC"])*10^3
max(AUC4[,"AUC"])*10^3
max(AUC5[,"AUC"])*10^3
max(AUC6[,"AUC"])*10^3
max(AUC7[,"AUC"])*10^3

# µmol/L
Cmax <- max(dataset$value)
Cmax2 <- max(dataset2$value)
Cmax3 <- max(dataset3$value)
Cmax4 <- max(dataset4$value)
Cmax5 <- max(dataset5$value)
Cmax6 <- max(dataset6$value)
Cmax7 <- max(dataset7$value)

# nmol/L
Cmax*10^3
Cmax2*10^3
Cmax3*10^3
Cmax4*10^3
Cmax5*10^3
Cmax6*10^3
Cmax7*10^3



###############################################################################
# RIF                                                    
###############################################################################

dataset <- subset(pred_long[[2]],compartment=="RIFcenlivc"|compartment=="RIFcengut")
xlimits     <- c(0,35)
xbreaks     <- seq(0,35,by = 7)
ylimits_log <- c(0.1,400)
ylimits          <- c(0,300)
ybreaks          <- seq(0,300,by = 100)

rifampicin <-
  ggplot(data=dataset,mapping=aes(color=factor(compartment)))+
  geom_line(aes(x=time/24, y=value), size=1.5)+
  labs(x="Time (days)",y="Concentration (µmol/L)")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    limits = ylimits_log,
    labels = function(x) sprintf("%g", x)
  )+
  annotation_logticks(sides="l")+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  scale_color_manual(values = c( "#85c4c9","#d39c83") )+
annotate("text", x = Inf, y = Inf, label = "Rifampicin ",size=11, hjust=1.05, 
         vjust=1.5)
  

rifampicin
ggsave("rifampicin.png", bg="white", width = 9, height =9)



###############################################################################
# CYP3A4 conc.                                              
###############################################################################

dataset <- subset(pred_long[[2]],compartment=="CYP3A4livc" | compartment == "CYP3A4gut")
xlimits          <- c(0,35)
xbreaks          <- seq(0,35,by = 7)
ylimits_log      <- c(5,100)
ylimits <- c(0,12)
ybreaks <- seq(0,12,by=2)
#fold
baseline <- dataset$value[which(dataset$time==0)]

CYP3A4 <-
  ggplot(data=dataset,mapping=aes(color=factor(compartment)))+
  geom_line(aes(x=time/24, y=value), size=1.5)+
  labs(x="Time (days)",y="Concentration (µmol/L)")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",
                                     linetype="solid",size = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.title.y.right = element_text( angle = 90),
     axis.title.y      = element_text(margin=margin(r=15)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +

  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  scale_color_manual(values = c( "#85c4c9","#d39c83") )+
  annotate("text", x = Inf, y = Inf, label = "CYP3A4 ",size=11, hjust=1.05, vjust=1.5)

CYP3A4
ggsave("CYP3A4.png", bg="white", width = 9, height =9)




###############################################################################
# CYP3A4 fold change                                            
###############################################################################

dataset <- subset(pred_long[[2]],compartment=="CYP3A4livc" | compartment == "CYP3A4gut")
xlimits          <- c(0,35)
xbreaks          <- seq(0,35,by = 7)
ylimits <- c(-1,12)
ybreaks <- seq(0,12,by=2)

#fold
baseline_livc <- dataset$value[which(dataset$time==0 & dataset$compartment == "CYP3A4livc")]
baseline_gut  <- dataset$value[which(dataset$time==0 & dataset$compartment == "CYP3A4gut")]
dataset$foldchange <- NA
dataset$foldchange[which(dataset$compartment == "CYP3A4livc")] <- dataset$value[which(dataset$compartment == "CYP3A4livc")]/baseline_livc
dataset$foldchange[which(dataset$compartment == "CYP3A4gut")] <- dataset$value[which(dataset$compartment == "CYP3A4gut")]/baseline_gut


CYP3A42 <-
  ggplot(data=dataset,mapping=aes(color=factor(compartment)))+
  geom_line(aes(x=time/24, y=foldchange), size=1.5)+
  labs(x="Time (days)",y="Relative conc. (fold)")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.title.y.right = element_text( angle = 90),
    axis.title.y      = element_text(margin=margin(r=15)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    limits = c(0.7,10),
    labels = function(x) sprintf("%g", x),
  )+
   annotation_logticks(sides="l")+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  scale_color_manual(values = c( "#85c4c9","#d39c83") )+
  annotate("text", x = Inf, y = Inf, label = "CYP3A4 ",size=11, hjust=1.05, 
           vjust=1.5)+
  geom_segment(aes(x = 4, y = 0.8, xend = 18, yend = 0.8),size=8,color="grey")

CYP3A42
ggsave("CYP3A42.png", bg="white", width = 10, height =10)




###############################################################################
# Metliv Cumulative amount
###############################################################################

#nonCYP3A4
#Without RIF
dataset <- subset(pred_long[[1]],compartment=="RETMetlivnonCYP3A4")
dataset$compartment <- "nonCYP3A4woRIF"
#with RIF
dataset2 <- subset(pred_long[[2]],compartment=="RETMetlivnonCYP3A4")
dataset2$compartment <- "nonCYP3A4wRIF"

#CYP3A4
#Without RIF
dataset3 <- subset(pred_long[[1]],compartment=="RETMetlivCYP3A4")
dataset3$compartment <- "CYP3A4woRIF"

#with RIF
dataset4 <- subset(pred_long[[2]],compartment=="RETMetlivCYP3A4")
dataset4$compartment <- "CYP3A4wRIF"


dataset_plot <- rbind(dataset,dataset2,dataset3,dataset4)


xlimits          <- c(0,35)
xbreaks          <- seq(0,35,by = 7)
ylimits_log      <- c(1,1000)


Metliv <-
  ggplot(data=dataset_plot, aes(color=factor(compartment)))+
  geom_line(aes(x=time/24, y=value*10^3),size=1.5)+
  labs(x="Time (days)",y="Cumulative amount (nmol)")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"

  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_continuous(
    limits = c(-5,70),
    breaks = seq(0,70,by=10))+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  scale_color_manual(values = c("80000000", "#cee7e9","#000000", "#85c4c9") )+
  annotate("text", x = Inf, y = Inf, label = "Metabolites liver ",size=11, 
           hjust=1.05, vjust=1.5)+
  geom_segment(aes(x = 4, y = -5, xend = 18, yend = -5),size=8,color="#85c4c9")


Metliv
ggsave("Metliv2.png", bg="white", width = 9, height =9)


#nonCYP3A4

#day 5
n_diff_1wo_5 <- (dataset$value[which(dataset$time==5*24)]-dataset$value[which(dataset$time==4*24)])*10^3
n_diff_2w_5 <- (dataset2$value[which(dataset2$time==5*24)]-dataset2$value[which(dataset2$time==4*24)])*10^3
n_diff_1wo_5
n_diff_2w_5

#day 7
n_diff_1wo_7 <- (dataset$value[which(dataset$time==7*24)]-dataset$value[which(dataset$time==6*24)])*10^3
n_diff_2w_7 <- (dataset2$value[which(dataset2$time==7*24)]-dataset2$value[which(dataset2$time==6*24)])*10^3
n_diff_1wo_7
n_diff_2w_7

#day 18
n_diff_1wo_18 <- (dataset$value[which(dataset$time==18*24)]-dataset$value[which(dataset$time==17*24)])*10^3
n_diff_2w_18 <- (dataset2$value[which(dataset2$time==18*24)]-dataset2$value[which(dataset2$time==17*24)])*10^3
n_diff_1wo_18
n_diff_2w_18

#day 20
n_diff_1wo_20 <- (dataset$value[which(dataset$time==20*24)]-dataset$value[which(dataset$time==19*24)])*10^3
n_diff_2w_20 <- (dataset2$value[which(dataset2$time==20*24)]-dataset2$value[which(dataset2$time==19*24)])*10^3
n_diff_1wo_20
n_diff_2w_20

#day 24
n_diff_1wo_24 <- (dataset$value[which(dataset$time==24*24)]-dataset$value[which(dataset$time==23*24)])*10^3
n_diff_2w_24 <- (dataset2$value[which(dataset2$time==24*24)]-dataset2$value[which(dataset2$time==23*24)])*10^3
n_diff_1wo_24
n_diff_2w_24

#day 32
n_diff_1wo_32 <- (dataset$value[which(dataset$time==32*24)]-dataset$value[which(dataset$time==31*24)])*10^3
n_diff_2w_32 <- (dataset2$value[which(dataset2$time==32*24)]-dataset2$value[which(dataset2$time==31*24)])*10^3
n_diff_1wo_32
n_diff_2w_32



#CYP3A4

#day 5
n_diff_3wo_5 <- (dataset3$value[which(dataset3$time==5*24)]-dataset3$value[which(dataset3$time==4*24)])*10^3
n_diff_4w_5 <- (dataset4$value[which(dataset4$time==5*24)]-dataset4$value[which(dataset4$time==4*24)])*10^3
n_diff_3wo_5
n_diff_4w_5

#day 7
n_diff_3wo_7 <- (dataset3$value[which(dataset3$time==7*24)]-dataset3$value[which(dataset3$time==6*24)])*10^3
n_diff_4w_7 <- (dataset4$value[which(dataset4$time==7*24)]-dataset4$value[which(dataset4$time==6*24)])*10^3
n_diff_3wo_7
n_diff_4w_7

#day 18
n_diff_3wo_18 <- (dataset3$value[which(dataset3$time==18*24)]-dataset3$value[which(dataset3$time==17*24)])*10^3
n_diff_4w_18 <- (dataset4$value[which(dataset4$time==18*24)]-dataset4$value[which(dataset4$time==17*24)])*10^3
n_diff_3wo_18
n_diff_4w_18

#day 20
n_diff_3wo_20 <- (dataset3$value[which(dataset3$time==20*24)]-dataset3$value[which(dataset3$time==19*24)])*10^3
n_diff_4w_20 <- (dataset4$value[which(dataset4$time==20*24)]-dataset4$value[which(dataset4$time==19*24)])*10^3
n_diff_3wo_20
n_diff_4w_20

#day 24
n_diff_3wo_24 <- (dataset3$value[which(dataset3$time==24*24)]-dataset3$value[which(dataset3$time==23*24)])*10^3
n_diff_4w_24 <- (dataset4$value[which(dataset4$time==24*24)]-dataset4$value[which(dataset4$time==23*24)])*10^3
n_diff_3wo_24
n_diff_4w_24

#day 32
n_diff_3wo_32 <- (dataset3$value[which(dataset3$time==32*24)]-dataset3$value[which(dataset3$time==31*24)])*10^3
n_diff_4w_32 <- (dataset4$value[which(dataset4$time==32*24)]-dataset4$value[which(dataset4$time==31*24)])*10^3
n_diff_3wo_32
n_diff_4w_32




###############################################################################
# Barplot liver metabolites - Daily formation of metabolite (nmol/day)
###############################################################################

daily_CYP3A4_liv <- data.frame(Compartment = c("RETMetlivCYP3A4", "RETMetlivCYP3A4", "RETMetlivCYP3A4",  "RETMetlivCYP3A4", "RETMetlivCYP3A4", "RETMetlivCYP3A4","RETMetlivCYP3A4"),
                            Value =       c(n_diff_3wo_5,      n_diff_4w_5,        n_diff_4w_7,       n_diff_4w_18,      n_diff_4w_20,      n_diff_4w_24,      n_diff_4w_32),
                            Unit =        c("nmol/day",        "nmol/day",         "nmol/day",        "nmol/day",        "nmol/day",        "nmol/day",        "nmol/day"),
                            RIFcondition= c("reference",         "*1",               "*3",              "*14",             "#2",              "#6",              "#14"),
                            CYPcondition= c("CYP3A4",          "CYP3A4",           "CYP3A4",          "CYP3A4",          "CYP3A4",          "CYP3A4",          "CYP3A4")
)

daily_CYP3A4_liv


daily_nonCYP3A4_liv <- data.frame(Compartment = c("RETMetlivnonCYP3A4", "RETMetlivnonCYP3A4", "RETMetlivnonCYP3A4",  "RETMetlivnonCYP3A4", "RETMetlivnonCYP3A4", "RETMetlivnonCYP3A4","RETMetlivnonCYP3A4"),
                                  Value =       c(n_diff_1wo_5,          n_diff_2w_5,         n_diff_2w_7,            n_diff_2w_18,         n_diff_2w_20,        n_diff_2w_24,        n_diff_2w_32),
                                  Unit =        c("nmol/day",           "nmol/day",           "nmol/day",             "nmol/day",           "nmol/day",          "nmol/day",          "nmol/day"),
                                  RIFcondition= c("reference",            "*1",                 "*3",                   "*14",                "#2",                "#6",                "#14"),
                                  CYPcondition= c("nonCYP3A4",          "nonCYP3A4",          "nonCYP3A4",            "nonCYP3A4",          "nonCYP3A4",         "nonCYP3A4",         "nonCYP3A4")
)

daily_nonCYP3A4_liv


dataset_barplot <- rbind(daily_CYP3A4_liv,daily_nonCYP3A4_liv)


dataset_barplot$Value <- signif(dataset_barplot$Value, digits = 3)
dataset_barplot$lab_ypos <- NA
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="reference")] <-  0.205/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="reference")]    <-  1.470+1.640/2 
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="*1")]  <-  0.205/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="*1")]     <-  2.450+0.652/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="*3")]  <-  0.205/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="*3")]     <-  1.210+1.380/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="*14")] <-  0.205/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="*14")]    <-  0.908+1.550/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="#2")]  <-  0.158/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="#2")]     <-  0.205+2.030/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="#6")]  <-  0.205/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="#6")]     <-  0.716+2.100/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="#14")] <-  0.205/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="#14")]    <-  1.440+2.100/2

dataset_barplot$RIFcondition <- factor(dataset_barplot$RIFcondition ,levels = c("reference","*1","*3","*14","#2", "#6", "#14"))

#stack, fill, dodge
barplot_liv <-
  ggplot(dataset_barplot, aes(fill=RIFcondition,x = RIFcondition, y = Value))+
  geom_bar(position="stack", stat="identity", width = 0.7,aes(alpha=CYPcondition))+
  geom_text(aes(y = lab_ypos, label = Value, group =CYPcondition), color = "white",size=7)+
  labs(x=" ",y="Daily formation (nmol/day)")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"
  )+
  coord_cartesian(ylim=c(0,3.5))+
  scale_fill_manual(values = c("#000000","#85c4c9", "#85c4c9", "#85c4c9", "#85c4c9",
                               "#85c4c9", "#85c4c9", "#85c4c9") )+
  scale_alpha_discrete(range=c(0.5,1))+
  annotate("text", x = Inf, y = Inf, label = "Metabolites liver ",size=11, 
           hjust=1.05, vjust=1.5)

barplot_liv
ggsave("barplot_daily_liv.png", bg="white", width = 9, height =9)



###############################################################################
# Metgut Cumulative amount
###############################################################################

#nonCYP3A4
#Without RIF
dataset <- subset(pred_long[[1]],compartment=="RETMetgutnonCYP3A4")
dataset$compartment <- "nonCYP3A4woRIF"
#with RIF
dataset2 <- subset(pred_long[[2]],compartment=="RETMetgutnonCYP3A4")
dataset2$compartment <- "nonCYP3A4wRIF"

#CYP3A4
#Without RIF
dataset3 <- subset(pred_long[[1]],compartment=="RETMetgutCYP3A4")
dataset3$compartment <- "CYP3A4woRIF"

#with RIF
dataset4 <- subset(pred_long[[2]],compartment=="RETMetgutCYP3A4")
dataset4$compartment <- "CYP3A4wRIF"


dataset_plot <- rbind(dataset,dataset2,dataset3,dataset4)

xlimits          <- c(0,35)
xbreaks          <- seq(0,35,by = 7)
ylimits_log      <- c(1,1000)
mapping=aes(color=ID)

Metgut <-
  ggplot(data=dataset_plot, aes(color=factor(compartment)))+
  geom_line(aes(x=time/24, y=value*10^3),size=1.5)+
  labs(x="Time (days)",y="Cumulative amount (nmol)")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"

  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_continuous(
    limits = c(-2,30),
    breaks = seq(0,30,by=5))+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  scale_color_manual(values = c("80000000", "#edd7cd","#000000", "#d39c83") )+
  annotate("text", x = Inf, y = Inf, label = "Metabolites gut ",size=11, 
           hjust=1.05, vjust=1.5)+
  geom_segment(aes(x = 4, y = -2, xend = 18, yend = -2),size=8,color="#d39c83")

Metgut
ggsave("Metgut2.png", bg="white", width = 9, height =9)


#nonCYP3A4

#day 5
n_diff_1wo_5 <- (dataset$value[which(dataset$time==5*24)]-dataset$value[which(dataset$time==4*24)])*10^3
n_diff_2w_5 <- (dataset2$value[which(dataset2$time==5*24)]-dataset2$value[which(dataset2$time==4*24)])*10^3
n_diff_1wo_5
n_diff_2w_5

#day 7
n_diff_1wo_7 <- (dataset$value[which(dataset$time==7*24)]-dataset$value[which(dataset$time==6*24)])*10^3
n_diff_2w_7 <- (dataset2$value[which(dataset2$time==7*24)]-dataset2$value[which(dataset2$time==6*24)])*10^3
n_diff_1wo_7
n_diff_2w_7

#day 18
n_diff_1wo_18 <- (dataset$value[which(dataset$time==18*24)]-dataset$value[which(dataset$time==17*24)])*10^3
n_diff_2w_18 <- (dataset2$value[which(dataset2$time==18*24)]-dataset2$value[which(dataset2$time==17*24)])*10^3
n_diff_1wo_18
n_diff_2w_18

#day 20
n_diff_1wo_20 <- (dataset$value[which(dataset$time==20*24)]-dataset$value[which(dataset$time==19*24)])*10^3
n_diff_2w_20 <- (dataset2$value[which(dataset2$time==20*24)]-dataset2$value[which(dataset2$time==19*24)])*10^3
n_diff_1wo_20
n_diff_2w_20

#day 24
n_diff_1wo_24 <- (dataset$value[which(dataset$time==24*24)]-dataset$value[which(dataset$time==23*24)])*10^3
n_diff_2w_24 <- (dataset2$value[which(dataset2$time==24*24)]-dataset2$value[which(dataset2$time==23*24)])*10^3
n_diff_1wo_24
n_diff_2w_24

#day 32
n_diff_1wo_32 <- (dataset$value[which(dataset$time==32*24)]-dataset$value[which(dataset$time==31*24)])*10^3
n_diff_2w_32 <- (dataset2$value[which(dataset2$time==32*24)]-dataset2$value[which(dataset2$time==31*24)])*10^3
n_diff_1wo_32
n_diff_2w_32



#CYP3A4

#day 5
n_diff_3wo_5 <- (dataset3$value[which(dataset3$time==5*24)]-dataset3$value[which(dataset3$time==4*24)])*10^3
n_diff_4w_5 <- (dataset4$value[which(dataset4$time==5*24)]-dataset4$value[which(dataset4$time==4*24)])*10^3
n_diff_3wo_5
n_diff_4w_5

#day 7
n_diff_3wo_7 <- (dataset3$value[which(dataset3$time==7*24)]-dataset3$value[which(dataset3$time==6*24)])*10^3
n_diff_4w_7 <- (dataset4$value[which(dataset4$time==7*24)]-dataset4$value[which(dataset4$time==6*24)])*10^3
n_diff_3wo_7
n_diff_4w_7

#day 18
n_diff_3wo_18 <- (dataset3$value[which(dataset3$time==18*24)]-dataset3$value[which(dataset3$time==17*24)])*10^3
n_diff_4w_18 <- (dataset4$value[which(dataset4$time==18*24)]-dataset4$value[which(dataset4$time==17*24)])*10^3
n_diff_3wo_18
n_diff_4w_18

#day 20
n_diff_3wo_20 <- (dataset3$value[which(dataset3$time==20*24)]-dataset3$value[which(dataset3$time==19*24)])*10^3
n_diff_4w_20 <- (dataset4$value[which(dataset4$time==20*24)]-dataset4$value[which(dataset4$time==19*24)])*10^3
n_diff_3wo_20
n_diff_4w_20

#day 24
n_diff_3wo_24 <- (dataset3$value[which(dataset3$time==24*24)]-dataset3$value[which(dataset3$time==23*24)])*10^3
n_diff_4w_24 <- (dataset4$value[which(dataset4$time==24*24)]-dataset4$value[which(dataset4$time==23*24)])*10^3
n_diff_3wo_24
n_diff_4w_24

#day 32
n_diff_3wo_32 <- (dataset3$value[which(dataset3$time==32*24)]-dataset3$value[which(dataset3$time==31*24)])*10^3
n_diff_4w_32 <- (dataset4$value[which(dataset4$time==32*24)]-dataset4$value[which(dataset4$time==31*24)])*10^3
n_diff_3wo_32
n_diff_4w_32




###############################################################################
# Barplot gut metabolites - Daily formation of metabolite (nmol/day)
###############################################################################

daily_CYP3A4_gut <- data.frame(Compartment = c("RETMetgutCYP3A4", "RETMetgutCYP3A4", "RETMetgutCYP3A4",  "RETMetgutCYP3A4", "RETMetgutCYP3A4", "RETMetgutCYP3A4","RETMetgutCYP3A4"),
                               Value =       c(n_diff_3wo_5,      n_diff_4w_5,        n_diff_4w_7,       n_diff_4w_18,      n_diff_4w_20,      n_diff_4w_24,      n_diff_4w_32),
                               Unit =        c("nmol/day",        "nmol/day",         "nmol/day",        "nmol/day",        "nmol/day",        "nmol/day",        "nmol/day"),
                               RIFcondition= c("reference",         "*1",               "*3",              "*14",             "#2",              "#6",              "#14"),
                               CYPcondition= c("CYP3A4",          "CYP3A4",           "CYP3A4",          "CYP3A4",          "CYP3A4",          "CYP3A4",          "CYP3A4")
)

daily_CYP3A4_gut


daily_nonCYP3A4_gut <- data.frame(Compartment = c("RETMetgutnonCYP3A4", "RETMetgutnonCYP3A4", "RETMetgutnonCYP3A4",  "RETMetgutnonCYP3A4", "RETMetgutnonCYP3A4", "RETMetgutnonCYP3A4","RETMetgutnonCYP3A4"),
                                  Value =       c(n_diff_1wo_5,          n_diff_2w_5,         n_diff_2w_7,            n_diff_2w_18,         n_diff_2w_20,        n_diff_2w_24,        n_diff_2w_32),
                                  Unit =        c("nmol/day",           "nmol/day",           "nmol/day",             "nmol/day",           "nmol/day",          "nmol/day",          "nmol/day"),
                                  RIFcondition= c("reference",            "*1",                 "*3",                   "*14",                "#2",                "#6",                "#14"),
                                  CYPcondition= c("nonCYP3A4",          "nonCYP3A4",          "nonCYP3A4",            "nonCYP3A4",          "nonCYP3A4",         "nonCYP3A4",         "nonCYP3A4")
)

daily_nonCYP3A4_gut


dataset_barplot <- rbind(daily_CYP3A4_gut,daily_nonCYP3A4_gut)


dataset_barplot$Value <- signif(dataset_barplot$Value, digits = 3)
dataset_barplot$lab_ypos <- NA
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="reference")] <-  0.171/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="reference")]    <-  0.247+0.276/2 
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="*1")]  <-  0.171/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="*1")]     <-  0.255+0.261/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="*3")]  <-  0.171/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="*3")]     <-  0.206+0.893/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="*14")] <-  0.171/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="*14")]    <-  0.193+1.060/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="#2")]  <-  0.171/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="#2")]     <-  0.171+1.830/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="#6")]  <-  0.171/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="#6")]     <-  0.218+0.644/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="nonCYP3A4" & dataset_barplot$RIFcondition=="#14")] <-  0.171/2
dataset_barplot$lab_ypos[which(dataset_barplot$CYPcondition=="CYP3A4" & dataset_barplot$RIFcondition=="#14")]    <-  0.246+0.285/2

dataset_barplot$RIFcondition <- factor(dataset_barplot$RIFcondition ,levels = c("reference","*1","*3","*14","#2", "#6", "#14"))

#stack, fill, dodge
barplot_gut <-
  ggplot(dataset_barplot, aes(fill=RIFcondition,x = RIFcondition, y = Value))+
  geom_bar(position="stack", stat="identity", width = 0.7,aes(alpha=CYPcondition))+
  geom_text(aes(y = lab_ypos, label = Value, group =CYPcondition), color = "white",size=7)+
  labs(x=" ",y="Daily formation (nmol/day)")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 2),
    panel.background  = element_rect(fill=NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=18),
    legend.title      = element_text(size =18),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66,
    legend.key.width = unit(2, "cm"),
    legend.key=element_blank(),
    legend.direction = "vertical"
  )+
  coord_cartesian(ylim=c(0,2.0))+
  scale_fill_manual(values = c("#000000","#d39c83", "#d39c83", "#d39c83", "#d39c83",
                               "#d39c83", "#d39c83", "#d39c83") )+
  scale_alpha_discrete(range=c(0.5,1))+
  annotate("text", x = Inf, y = Inf, label = "Metabolites gut ",size=11, hjust=1.05, vjust=1.5)

barplot_gut
ggsave("barplot_daily_gut.png", bg="white", width = 9, height =9)








