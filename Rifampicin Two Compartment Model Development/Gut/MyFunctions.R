###############################################################################
# Function to log transform data
###############################################################################
TransData <- function(x) {
  
  x_hat <- log(x+ 1)
  
  return(x_hat)
  
}

###############################################################################
# Function to reverse log transformation of data
###############################################################################

DeTransData <- function(x_hat) {
  
  x <- exp(x_hat) - 1
  
  return(x)
  
}


###############################################################################
# Functions to create histogram and pairs plot for MCMC results
###############################################################################

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}

###############################################################################
# Function to calculate mode 
###############################################################################

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}


###############################################################################
# Function to measure runtime
###############################################################################

tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}

###############################################################################
# function returns head of a list 
###############################################################################

head.list <- function(obj, n = 6L, ...)
{
  stopifnot(length(n) == 1L)
  origN <- n
  n <- if (n < 0L)
    max(length(obj) + n, 0L)
  else min(n, length(obj))
  lapply(obj[seq_len(n)], function(x)
  {
    tryCatch({
      head(x, origN, ...)
    }, error = function(e) {
      x
    })
  })
}


###############################################################################
# HISTOGRAM
###############################################################################
GetPlotHistogram <- function(dataset,aest,xlabel,ylabel,title,bw) {

  
  
plot <- 
  
  ggplot(dataset,aest)+
  geom_histogram(binwidth = bw,colour="black", fill="white",aes(y=..density..))+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+
  theme(
    plot.background   = element_rect(fill = NA),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
    panel.background  = element_rect(fill=NA),
    plot.title       = element_text(size=10, hjust = 0.5),
    axis.line        = element_line(size = 0.5),
    axis.text        = element_text(colour = "black",size = 10),
    axis.title       = element_text(size = 10),
    axis.title.y     = element_text(margin=margin(r=15)),
    axis.title.x     = element_text(margin=margin(t=15)),
    axis.ticks       = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=10),
    legend.title      = element_text(size =10),
    legend.position   = "bottom",
    legend.background = element_rect(fill = NA))+
  stat_function(fun = dnorm, color="blue",args = list(mean = mean(dataset$RESIDUALS), sd = sd(dataset$RESIDUALS)))


return(plot)
}



###############################################################################
# GEOM-POINT (OBSERVED VS. PREDICTED)                     
# LINE AT INTERCEPT 0
###############################################################################

#plot data normal scale 
GetPlotResiduals2 <- function(dataset,aest,xlabel,ylabel,title,xlimits,ylimits) {
  library(ggplot2)
  plot <- 
    ggplot(dataset,aest)+
    geom_point(shape=16,size=2)+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title       = element_text(size=10, hjust = 0.5),
      axis.line        = element_line(size = 0.5),
      axis.text        = element_text(colour = "black",size = 10),
      axis.title       = element_text(size = 10),
      axis.title.y     = element_text(margin=margin(r=15)),
      axis.title.x     = element_text(margin=margin(t=15)),
      axis.ticks       = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=10),
      legend.title      = element_text(size =10),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA)
    )+
    geom_abline(slope=1, intercept=0, linetype="dashed", color = "black",show.legend = TRUE)+
    scale_x_continuous(limits = xlimits)+
    scale_y_continuous(limits = ylimits)+
    stat_smooth(formula = y ~ x, method = 'lm',se = FALSE,size = 0.25,fullrange=TRUE,show.legend =TRUE)

  
  return(plot)
}




###############################################################################
# GEOM-POINT (RESIDUALS VS. PREDICTED)                     
# LINE AT INTERCEPT 0
###############################################################################

#plot data normal scale 
GetPlotResiduals <- function(dataset,aest,xlabel,ylabel,title) {
  library(ggplot2)
  plot <- 
    ggplot(dataset,aest)+
    geom_point(shape=16,size=2)+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=10, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 10),
      axis.title        = element_text(size = 10),
      axis.title.y      = element_text(margin=margin(r=15)),
      axis.title.x      = element_text(margin=margin(t=15)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=10),
      legend.title      = element_text(size =10),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA),
      aspect.ratio = 0.5
    )+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    scale_y_continuous(
           limits = c(-3,3),
           breaks = seq(-3,3,by = 1))

  return(plot)
}


###############################################################################
# GEOM_LINE (Prediction) and GEOM-POINT (Observed data)                       
# LOG SCALE Y AXIS  
# MAP BOTH BY COLOUR
# INSERT HORIZONTAL LINES
###############################################################################

GetPlotPointLineMappingLogScaleY2 <- function(dataset1,dataset2,aest,xlabel,ylabel,title,pdfname_log,pngname_log, ylimits_log,xlimits,xbreaks) {
  library(ggplot2)
  plot <- 
    ggplot(data=dataset1,mapping=aes(color=ID))+
    geom_point(aest,shape=16,size=4)+
    geom_line(data=dataset2,aest, size=1)+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=18, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 18),
      axis.title        = element_text(size = 18),
      # axis.title.y      = element_text(margin=margin(r=15)),
      # axis.title.x      = element_text(margin=margin(t=12)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA),
      aspect.ratio=1
      
    )+
    guides(col = guide_legend(nrow = 6,title.position = "top")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = ylimits_log
    )+
    scale_x_continuous(
      limits = xlimits,
      breaks = xbreaks)+
    annotation_logticks(sides="l")+
    scale_color_manual(values=c("#CEAB07","#899DA4","#D8A499","#35274A","#0B775E","#DD8D29","#972D15","#D8B70A")) +
    geom_hline(yintercept=0, linetype="dashed", color = "black")
  
  
  ggsave(pdfname_log, width = 7.8, height = 6)
  ggsave(pngname_log, bg="transparent", width = 7.8, height = 6)
  
  return(plot)
}



###############################################################################
# GEOM_LINE (Prediction) and GEOM-POINT (Observed data)                       
# LOG SCALE Y AXIS  
#AND MAP BOTH BY COLOUR
###############################################################################

GetPlotPointLineMappingLogScaleY <- function(dataset1,dataset2,aest,xlabel,ylabel,title,pdfname_log,pngname_log, ylimits_log,xlimits,xbreaks) {
  library(ggplot2)
  plot <- 
    ggplot(data=dataset1,mapping=aes(color=ID))+
    geom_point(aest,shape=16,size=4)+
    geom_line(data=dataset2,aest, size=1)+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=18, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 18),
      axis.title        = element_text(size = 18),
      # axis.title.y      = element_text(margin=margin(r=15)),
      # axis.title.x      = element_text(margin=margin(t=12)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA),
      aspect.ratio=1
      
    )+
    guides(col = guide_legend(nrow = 6,title.position = "top")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = ylimits_log
    )+
    scale_x_continuous(
      limits = xlimits,
      breaks = xbreaks)+
    annotation_logticks(sides="l")+
    scale_color_manual(values=c("#0B775E","#DD8D29","#972D15","#D8B70A")) # "#35274A","#CEAB07","#899DA4","#D8A499",
  
  ggsave(pdfname_log, width = 7.8, height = 6)
  ggsave(pngname_log, bg="transparent", width = 7.8, height = 6)
  
  return(plot)
}

  
###############################################################################
# GEOM_LINE (Prediction) and GEOM-POINT (Observed data)                       
# LOG SCALE Y AXIS                                                            
###############################################################################

GetPlotPointLineLogScaleY <- function(dataset,aest,aest2,xlabel,ylabel,title,pdfname_log,pngname_log, ylimits_log) {
  library(ggplot2)
  plot <- 
    ggplot(data=dataset,aest)+
    geom_point(shape=16,size=2)+
    geom_line(data=dataset,aest2)+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=11, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 11),
      axis.title        = element_text(size = 11),
      axis.title.y      = element_text(margin=margin(r=15)),
      axis.title.x      = element_text(margin=margin(t=15)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=11),
      legend.title      = element_text(size =11),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA)
    )+
    guides(col = guide_legend(nrow = 6,title.position = "top")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = ylimits_log
      )+
    annotation_logticks(sides="l")
  
  # ggsave(pdfname_log)
  # ggsave(pngname_log, bg="transparent")
  # 
  return(plot)
}


###############################################################################
# GEOM_LINE                                                                   
# MAPPING SPECIFIED IN MAIN (FACET WRAP)                                      
# LOG SCALE Y AXIS                                                            
###############################################################################

GetPlotLineMappingLogScaleY <- function(dataset,aest,xlabel,ylabel,title,pdfname_log,pngname_log,mapping,xlimits,xbreaks,ylimits_log) {
  library(ggplot2)
  plot <- 
    ggplot(dataset,aest)+
    geom_line()+
    mapping+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=18, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 10),
      axis.title        = element_text(size = 18),
      axis.title.y      = element_text(margin=margin(r=15)),
      axis.title.x      = element_text(margin=margin(t=15)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA)
    )+
    guides(col = guide_legend(nrow = 6,title.position = "top")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = ylimits_log
      )+
    annotation_logticks(sides="l")+
    scale_x_continuous(
     limits = xlimits,
     breaks = xbreaks)

  # ggsave(pdfname_log)
  # ggsave(pngname_log, bg="transparent")
  
  return(plot)
}

###############################################################################
# GEOM_POINT                                                                  
# MAPPING BY COLOR (ID)                                                       
# LOG SCALE Y AXIS                                                            
###############################################################################

# plot geom_point by color mapping log scale y axis
GetPlotMappingLogScaleY <- function(dataset,aest,xlabel,ylabel,title,xlimits,xbreaks,ylimits,ybreaks,pdfname,pngname) {
  library(ggplot2)
  plot <- 
    ggplot(dataset,aest)+
    geom_point(shape=16,size=5,mapping=aes(color=ID))+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=23, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 23),
      axis.title        = element_text(size = 23),
      axis.title.y      = element_text(margin=margin(r=15)),
      axis.title.x      = element_text(margin=margin(t=15)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA)
      
      # plot.background   = element_rect(fill="transparent",colour=NA), #transparent background
      )+
    guides(col = guide_legend(nrow = 6,title.position = "top")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = ylimits_log)+
    scale_x_continuous(
      limits = xlimits,
      breaks = xbreaks)+
    annotation_logticks(sides="l")+
    scale_color_manual(values=c("#D8A499","#0B775E","#DD8D29","#972D15","#D8B70A")) # "#899DA4","#CEAB07","#35274A",
  
  ggsave(pdfname_log)
  ggsave(pngname_log, bg="transparent")
  
  return(plot)
}


###############################################################################
# plot by shape mapping log scale y axis
###############################################################################

GetPlotShapeMappingLogScaleY <- function(dataset,aest,xlabel,ylabel,title,xlimits,xbreaks,ylimits,ybreaks,pdfname,pngname) {
  library(ggplot2)
  plot <- 
    ggplot(dataset,aest)+
    geom_point(size=4,mapping=aes(shape=compartment))+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=18, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 18),
      axis.title        = element_text(size = 18),
      axis.title.y      = element_text(margin=margin(r=15)),
      axis.title.x      = element_text(margin=margin(t=15)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA)
      )+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = ylimits_log)+
    scale_x_continuous(
      limits = xlimits,
      breaks = xbreaks)+
    annotation_logticks(sides="l")
    
  
  ggsave(pdfname_log)
  ggsave(pngname_log, bg="transparent")
  
  return(plot)
}

###############################################################################
# Plot data normal scale
###############################################################################

GetPlotNormalScale <- function(dataset,aest,xlabel,ylabel,title,xlimits,xbreaks,ylimits,ybreaks,pdfname,pngname) {
 library(ggplot2)
  plot <- 
    ggplot(dataset,aest)+
    geom_point(shape=16,size=4)+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title       = element_text(size=18, hjust = 0.5),
      axis.line        = element_line(size = 0.5),
      axis.text        = element_text(colour = "black",size = 18),
      axis.title       = element_text(size = 18),
      axis.title.y     = element_text(margin=margin(r=15)),
      axis.title.x     = element_text(margin=margin(t=15)),
      axis.ticks       = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA)
      )+
    scale_y_continuous(
      limits = ylimits,
      breaks = ybreaks)+
    scale_x_continuous(
      limits = xlimits,
      breaks = xbreaks)
  
  ggsave(pdfname)
  ggsave(pngname, bg="transparent")
  
  return(plot)
}



###############################################################################
# Plot data log scale Y axis
###############################################################################

GetPlotLogScaleY <- function(dataset,aest,xlabel,ylabel,title,xlimits,xbreaks,ylimits_log,pdfname_log,pngname_log) {
  library(scales, ggplot2)
  
  plot <- 
    ggplot(dataset,aest)+
    geom_point(shape=16,size=3.5,color="black")+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=18, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 18),
      axis.title        = element_text(size = 18),
      axis.title.y      = element_text(margin=margin(r=15)),
      axis.title.x      = element_text(margin=margin(t=15)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA)
     )+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = ylimits_log)+
    scale_x_continuous(
      limits = xlimits,
      breaks = xbreaks)+
     annotation_logticks(sides="l")
  
  ggsave(pdfname_log)
  ggsave(pngname_log, bg="transparent")
  
  return(plot)
}


###############################################################################
# Plot data log scale X axis
###############################################################################

GetPlotLogScaleX <- function(dataset,aest,xlabel,ylabel,title,xlimits_log,ylimits,ybreaks,pdfname_log,pngname_log) {
  library(scales, ggplot2)
  
  plot <- 
    ggplot(dataset,aest)+
    geom_point(shape=1,size=4)+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=18, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 18),
      axis.title        = element_text(size = 18),
      axis.title.y      = element_text(margin=margin(r=15)),
      axis.title.x      = element_text(margin=margin(t=15)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA)
      )+
    scale_y_continuous(
      limits = ylimits,
      breaks = ybreaks)+
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = xlimits_log)+
    annotation_logticks(sides="b")
      
  
  ggsave(pdfname_log)
  ggsave(pngname_log, bg="transparent")
  
  return(plot )
}



