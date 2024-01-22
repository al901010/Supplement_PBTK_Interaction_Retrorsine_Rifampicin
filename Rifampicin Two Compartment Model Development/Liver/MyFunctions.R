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

