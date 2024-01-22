
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
source("MyFunctions.R")
source("Pred.R")
source("Minus2LL.R")


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


# Define color-blind friendly palette (Okate and Ito 2002)
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", 
                      "#CC79A7", "#999999", "#F0E442") 


###############################################################################
# PART: OBSERVED DATA                                                       
###############################################################################

#--- READ DATASET 

library(readxl)
dataset_rifampicin <- read_excel("ObservedDataRifampicin.xlsx", 
                                     col_types = c("text", "numeric", "text", 
                                                   "numeric", "text", "text", 
                                                   "numeric", "numeric", 
                                                   "numeric", "text", "numeric", 
                                                   "numeric", "text", "text",
                                                   "text", "text", "text", 
                                                   "text", "text", "text", 
                                                   "text", "text"))
View(dataset_rifampicin)

#---Subset
dataset_rifampicin <- subset(dataset_rifampicin,
                             ID=="Rifampicin_human_liver_600mg_po_MD")


#--- TRANSFORM UNITS 

dataset_rifampicin$Y <- NA
dataset_rifampicin$Y <- dataset_rifampicin$YORIGINAL*1.2996  # fraction of liver intracellular space * Vliv [L]
dataset_rifampicin$YUNIT <- "umol"

# Initialize 
dataset_rifampicin$SQUAREDRESIDUALS <- NA
dataset_rifampicin$RESIDUALS        <- NA

#--- TRANSFORM OBSERVED DATA 
dataset_rifampicin$Y<- TransData(dataset_rifampicin$Y)


###############################################################################
# PART: INITIALIZATION OF PREDICTION                                              
###############################################################################

#--- CREATE RUNTABLE AND RUNID 

# Create runtable and give runid for each run
runtable_rifampicin <- data.frame(unique(subset(dataset_rifampicin,!(is.na(dataset_rifampicin$DOSEORIGINAL)),
                                     select = c("DOSEORIGINAL","DOSEORIGINALUNIT","NODOSES","DOSINGINTERVAL","DOSINGINTERVALUNIT","ADM"))))
runtable_rifampicin$RUNID <- seq.int(nrow(runtable_rifampicin))
View(runtable_rifampicin)
  
#Transfer RUNID into dataset
dataset_rifampicin$RUNID <- NA

for(row in 1:nrow(dataset_rifampicin)){
  
  if(dataset_rifampicin$ID[row]=="Rifampicin_human_liver_600mg_po_MD"){
    dataset_rifampicin$RUNID[row] <- 1
  }
}


#--- DEFINE ODE SYSTEM FOR RxODE PACKAGE 

mRifampicin <-RxODE({
  d/dt (RIFdep)  = -ka*RIFdep;
  d/dt (RIFcen) = +ka*RIFdep + k21*RIFper - k12*RIFcen - ke*RIFcen;
  d/dt (RIFper) = +k12*RIFcen - k21*RIFper;
  d/dt (RIFbag) = +ke*RIFcen;
  
})


#################################################################################
# PART: MAXIMUM LIKELIHOOD ESTIMATION                                       
#################################################################################

################################################################################# 
# CAUTION: There is a bug in function nmkb of package dfoptim. To solve the bug: 
# the position of the elements in vector 'est' have to be assigned to 
# the names of the elements of 'est' in functions GetMinus2LL and GetPred 
#################################################################################

# timepoints
timepoints <- unique(subset(subset(dataset_rifampicin,is.na(dataset_rifampicin[["DOSEORIGINAL"]])), select = XORIGINAL))
timepoints <- sort(as.vector(timepoints$XORIGINAL))

#--- DEFINE START VALUES FOR PARAMETERS TO BE ESTIMATED 


est0 <- c(a_RIF=0.12, 
          ka= 0.4, #0.21, 
          k12=0.2, #0.12,
          k21=0.2, #0.03,  
          ke=0.4) #0.57 )



dataset_rifampicin_pred <- PredRxODE(est=est0,
                                     t=timepoints,
                                     dataset=dataset_rifampicin, 
                                     runtable=runtable_rifampicin, 
                                     add_timepoints = FALSE)

View(dataset_rifampicin_pred)

#--- DEFINE BOUNDARIES FOR PARAMETERS TO BE ESTIMATED
lboundary <- c(1e-10,
               0,
               0,
               0,
               0)
uboundary <- c(Inf,
               Inf,
               Inf,
               Inf,
               Inf)

#--- Hooke Jeeves direct search from Package dfoptim

set.seed(10)

# out_hjkb <- hjkb(par=est0,
#                  t=timepoints,
#                  fn=GetMinus2LL,
#                  lower = lboundary,
#                  upper = uboundary,
#                  control = list(tol=1e-10), # default 1e-6
#                  dataset=dataset_rifampicin,
#                  runtable=runtable_rifampicin)
# out_hjkb
# 
# save(out_hjkb, file="out_hjkb.rda")

load(file = "out_hjkb.rda")
est_hat_hjkb <- out_hjkb$par

#--- AKAIKE INFORMATION CRITERION

minus2LL <- out_hjkb$value
minus2LL
AIC <- 2*length(est0) + minus2LL
AIC


#################################################################################
# PART: MCMC
#################################################################################

#--- PERFORM MCMC (PACKAGE FME) TO ESTIMATE PARAMETER DISTRIBUTION

numberiter <- 10000 
burnin     <- 0.1*numberiter
update     <- 0.1*numberiter
dr         <- 3
dispersion <- sample(seq(0.7,1.3,by=0.1),length(est_hat_hjkb),replace=TRUE)
dispersion2 <- sample(seq(0.7,1.3,by=0.1),length(est_hat_hjkb),replace=TRUE)

# # Generate first chain
# out_MCMC1 <- modMCMC(f=GetMinus2LL,
#                      p=est_hat_hjkb,   # initial values for parameters
#                      jump=NULL,          # SD of proposal normal distribution; if NULL 10% of p
#                      # it can be efficient to use covar from model fit
#                      lower=lboundary,
#                      upper=uboundary,
#                      prior=NULL,         # -2log(parameter prior probability); NULL: non-informative prior, all parameters are equally likely
#                      var0=NULL,          # NULL: it is assumned that model variance is 1 and the return element from f is -2logL
#                      wvar0=NULL,         # "weight" for initial model variance, NULL: error variances are assumed to be fixed
#                      n0=NULL,            # parameter used for weihghing initial model variance, NULL: n0=wvar0*n
#                      niter=numberiter,   # no of iterations for the MCMC
#                      updatecov=update,   # setting updatecov smaller than niter will trigger adaptive MH,
#                      # proposal distribution is only updated during burnin when burninlenth is positive
#                      burninlength=burnin,# no of initial iterations to be removed, about 10% of niter
#                      ntrydr=dr,          # max no of tries for delayed rejection procedure
#                      dataset=dataset_retrorsine,
#                      runtable=runtable_retrorsine,
#                      t=timepoints)
# 
# save(out_MCMC1, file="out_MCMC1.rda")
# 
# 
# 
# # Generate second chain
# out_MCMC2 <- modMCMC(f=GetMinus2LL,
#                      p=dispersion*est_hat_hjkb,
#                      jump=NULL,
#                      lower=lboundary,
#                      upper=uboundary,
#                      prior=NULL,
#                      var0=NULL,
#                      wvar0=NULL,
#                      n0=NULL,
#                      niter=numberiter,
#                      updatecov=update,
#                      burninlength=burnin,
#                      ntrydr=dr,
#                      dataset=dataset_retrorsine,
#                      runtable=runtable_retrorsine,
#                      t=timepoints)
# 
# save(out_MCMC2, file="out_MCMC2.rda")
# 
# 
# # Generate third chain
# out_MCMC3 <- modMCMC(f=GetMinus2LL,
#                      p=dispersion2*est_hat_hjkb,
#                      jump=NULL,
#                      lower=lboundary,
#                      upper=uboundary,
#                      prior=NULL,
#                      var0=NULL,
#                      wvar0=NULL,
#                      n0=NULL,
#                      niter=numberiter,
#                      updatecov=update,
#                      burninlength=burnin,
#                      ntrydr=dr,
#                      dataset=dataset_retrorsine,
#                      runtable=runtable_retrorsine,
#                      t=timepoints)
# 
# save(out_MCMC3, file="out_MCMC3.rda")
# 

# load files
load(file = "out_MCMC1.rda")
load(file = "out_MCMC2.rda")
load(file = "out_MCMC3.rda")

# Join  chains
out_MCMC <- rbind(out_MCMC1$pars,
                  out_MCMC2$pars,
                  out_MCMC3$pars)

# MCMC output can be used as functions from the coda package
MC1 <- as.mcmc(out_MCMC1$pars)
MC2 <- as.mcmc(out_MCMC2$pars)
MC3 <- as.mcmc(out_MCMC3$pars)


# Gelman-Rubin convergence dignostc
combinedchains <- mcmc.list(MC1,MC2,MC3)
gelman.plot(combinedchains)

gelmandiag     <- gelman.diag(combinedchains, confidence = 0.95,autoburnin = TRUE, multivariate = TRUE)
gelmandiag

#trace plot
plot(combinedchains)

# pairs plot
sample_params <- out_MCMC[sample(nrow(out_MCMC),size=1000,replace = TRUE),]

sample_params_plot <- subset(sample_params, select = c(ka,k12,k21,ke))

sample_params_plot   <- sample_params_plot[seq(1, NROW(sample_params_plot), 1), ]
pairs(sample_params_plot, diag.panel = panel.hist,upper.panel = panel.cor,panel=panel.smooth,
      cex.labels = 1.5, font.labels = 2, cex = 1.5,pch = 21)


# summary of posterior distribution

# Initialize
summary<- data.frame(matrix(NA,nrow = 1, ncol = 1))

# loop through posterior 
for(i in colnames(out_MCMC)) {
  print(i)
  
  #summary <- summary(out_MCMC[,i])
  mode_i    <- map_estimate(out_MCMC[,i])
  # mode_i    <- estimate_mode(out_MCMC[,i])
  mean_i    <- mean(out_MCMC[,i])
  median_i  <- median(out_MCMC[,i])
  ci_i      <- ci(out_MCMC[,i], ci=0.95,method = "HDI")
  summary_i <- rbind(median=median_i,mean=mean_i,mode=mode_i,ci=ci_i$CI,
                     ci_low=ci_i$CI_low,ci_high=ci_i$CI_high)
  colnames(summary_i) <- i
  
  print(summary_i)
  summary <- cbind(summary, summary_i)
}

# remove initialisation column
summary <- summary[colSums(!is.na(summary)) > 0]

View(summary)

#mode
est_mode <- unlist(summary[3,])

# minus 2LL out_MCMC1
minus2LL_out_MCMC1 <- out_MCMC1$bestfunp
minus2LL_out_MCMC1
AIC <- 2*length(est0) + minus2LL_out_MCMC1
AIC



#################################################################################
# PART: PREDICTION
#################################################################################

#RxODE

#--- SOLVE ODE SYSTEM USING OPTIMIZED PARAMETERS

dataset_rifampicin_pred <- PredRxODE(est=est_mode,
                                     t=c(timepoints,seq(0,24*42,by=0.1)),
                                     dataset=dataset_rifampicin,
                                     runtable=runtable_rifampicin,
                                     add_timepoints = TRUE)

View(dataset_rifampicin_pred)



#################################################################################
# PART: SSR                                   
#################################################################################
  
dataset_rifampicin_pred$RESIDUALS        <- (dataset_rifampicin_pred$Y - dataset_rifampicin_pred$YPRED)
dataset_rifampicin_pred$SQUAREDRESIDUALS <- (dataset_rifampicin_pred$Y - dataset_rifampicin_pred$YPRED)^2


SSR_RIF <- sum(dataset_rifampicin_pred$SQUAREDRESIDUALS[which(
  (!is.na(dataset_rifampicin_pred$SQUAREDRESIDUALS)) & 
    ( (dataset_rifampicin_pred$ID=="Rifampicin_human_liver_600mg_po_MD") )
)]
)
  

################################################################################  
# PART: POSTPROCESSING
################################################################################

#--- DETRANSFORM PREDICTIONS 

for(row in 1:nrow(runtable_rifampicin)) {
  results_long[[row]]$value <- DeTransData(results_long[[row]]$value)
}  

dataset_rifampicin_pred$Y     <- DeTransData(dataset_rifampicin_pred$Y)
dataset_rifampicin_pred$YPRED <- DeTransData(dataset_rifampicin_pred$YPRED)

# dataset_ppd$Y     <- DeTransData(dataset_ppd$Y)
# dataset_ppd$YPRED <- DeTransData(dataset_ppd$YPRED)


#--- CONVERT Y AND YPRED INTO Y_PERCENT AND YPRED_PERCENT IN DATASET 

# Initialize new columns
dataset_rifampicin_pred$Y_PERCENT <- NA
dataset_rifampicin_pred$Y_PERCENTUNIT <- "%"
dataset_rifampicin_pred$YPRED_PERCENT <- NA
dataset_rifampicin_pred$YPRED_PERCENTUNIT <- "%"

###############################################################################
# PART: Plots                                                            
###############################################################################

# Subset of datasets
dataset_rifampicin_plot     <- subset(dataset_rifampicin_pred,
                                      is.na(dataset_rifampicin_pred$DOSEORIGINAL))

# Dataset for plotting in long format (TYPE YPRED or Y)
dataset_rifampicin_long_to_be   <- subset(dataset_rifampicin_plot, 
                                          select = c(ID,COMPARTMENT,SPECIES,XORIGINAL,Y,YPRED))
dataset_rifampicin_long_to_be_2 <- gather(dataset_rifampicin_long_to_be, TYPE, Y, 5:6, factor_key=TRUE)                                  
dataset_rifampicin_long         <- subset(dataset_rifampicin_long_to_be_2, (!(is.na(Y))))



###############################################################################
# OBSERVED AND PREDICTED log scale
###############################################################################

dataset_rifampicin_long_RIF <-  subset(dataset_rifampicin_long, 
                                       ID=="Rifampicin_human_liver_600mg_po_MD")


dataset1    <- subset(dataset_rifampicin_long_RIF, dataset_rifampicin_long_RIF$TYPE=="Y")
dataset2    <- subset(dataset_rifampicin_long_RIF, dataset_rifampicin_long_RIF$TYPE=="YPRED")
aest        <- aes(x=X, y=Y)
xlabel      <- "Time (days)"
ylabel      <- "Concentration (µmol/L)"
title       <- ""
ylimits_log <- c(0.1,1000)
ylimits     <- c(0,100)
ybreaks     <- seq(0,100,by = 10)
xlimits     <- c(0,21)
xbreaks     <- seq(0,21,by = 7)

rifampicin <- 
  ggplot(data=dataset1,mapping=aes(color=TYPE,linetype=TYPE))+
  geom_line(aes(x=XORIGINAL/24, y=Y/1.2996),size=1.5)+ #divided by liver volume
  geom_line(data=dataset2,aes(x=XORIGINAL/24, y=Y/1.2996),size=1.5)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=22),
    legend.title      = element_text(size =22),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    #breaks = scales::trans_breaks("log10", function(x) 10^x),
    #labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log,
    labels = function(x) sprintf("%g", x)
  )+
  annotation_logticks(sides="l")+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  scale_color_manual(values=c("grey","black"))+
  annotate("text", x = Inf, y = Inf, label = "Rifampicin (liver intracellular) ",size=11, hjust=1.05, vjust=1.5)


rifampicin
ggsave("rifampicin.png", bg="white", width = 9, height =9)


