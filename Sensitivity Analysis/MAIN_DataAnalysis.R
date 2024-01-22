
#--- INITIALIZE ---------------------------------------------------------------

# Remove all objects from current workspace
rm(list=ls())

# Reset graphics
graphics.off()
cat("\014")

# Set working directory to current path
PATH <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(PATH)

# install packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load("rstudioapi","ggplot2","readxl","readr")

# Load packages
library(ggplot2)
library(readxl)
library(readr)


###############################################################################
# LOAD AND PROCESS DATA 
###############################################################################

# read file
dataset <- read_excel("SensitivityAnalysis.xlsx", 
                      sheet = "plot",
                      col_names = c("ModelParameter",
                                    "SensitivityCoefficient",
                                    "Tissue"),skip=1)
View(dataset)


###############################################################################
# Plot
###############################################################################

dataset$ModelParameter <- factor(dataset$ModelParameter,
                                 levels = c("EC50·1.1","EC50·0.9",
                                            "Ki·1.1","Ki·0.9",
                                            "Emax·1.1","Emax·0.9",
                                            "fm,CYP3A4·1.1","fm,CYP3A4·0.9" ))


plot <-
  ggplot(dataset, aes(x = ModelParameter, y = SensitivityCoefficient)) +
  geom_col( aes(fill = Tissue),position = "dodge", width=0.8) + 
  coord_flip() +
  labs(x="Model parameter variation",y="Sensitivity coefficient")+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 2),
    panel.background  = element_rect(fill=NA),
    panel.grid.major.x = element_line(color = "gray",
                                      size = 0.5,
                                      linetype = 1),
    panel.grid.minor.x = element_line(color = "gray",
                                      size = 0.5,
                                      linetype = 1),
    panel.grid.major.y = element_line(color = NA),
    panel.grid.minor.y = element_line(color = NA),
    plot.title        = element_text(size=18, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 22),
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
  scale_y_continuous(
    limits = c(-1,1),
    breaks = seq(-1,1,by=0.5))+
  geom_hline(yintercept = 0,size=1)+
  scale_fill_manual(values = c("#d39c83","#85c4c9"))

plot

ggsave("sensitivity.png", bg="white", width = 9, height =9)

