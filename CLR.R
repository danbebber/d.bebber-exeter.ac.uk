#~~~~~~~~~~~~~~~~~~~~~~~~#
# COFFEE LEAF RUST MODEL #
#~~~~~~~~~~~~~~~~~~~~~~~~#

# LIBRARIES
library(raster)
library(fields)
require(wq)
require(accelerometry)
require(ggplot2)
require(RColorBrewer)

# DIRECTORY
setwd("~/Dropbox/Coffee/CLR/xFarm")

# GRAPHICS DEFAULTS
theme_set(theme_bw(base_size = 10))

# FUNCTIONS
# Temperature-dependent acceleration factor
rT <- function(Temp, card){ # Temp = vector of temperatures, card = cardinal temperatures
  tmin <- card[1]; topt <- card[2]; tmax <- card[3]
  r <- ((tmax-Temp)/(tmax-topt))*((Temp-tmin)/(topt-tmin))^
    ((topt-tmin)/(tmax-topt))
  r[Temp <= tmin | Temp >= tmax] <- 0 #set values outside temperature range to zero
  r}
# Cumulative distribution function for accelerated failure time model
Ft <- function(u, pars){ # u = effective elapsed time; pars = c(alpha, beta) of Weibull
  1 - exp(-(u*pars[1])^pars[2])
}

# Find wet periods in timeseries
library(accelerometry)
wetper <- function(wet, min = 6){ # wet = binary vector, # limit = min LWD
  ind <- rle2(wet, indices = TRUE) # run lengths
  ind <- ind[ind[,1] == 1,] # wet periods
  ind <- ind[ind[,4] >= min, 2:3]
  if(nrow(ind)==0) warning(paste("No wet periods with length >= ", min, sep = "")) else ind
  }

# Calculate infections from a timeseries
infect <- function(Temp, par.G, par.A, int = 1){ 
  # Temp = vector (T(0),...,T(W))
  # par.G = c(tmin, topt, tmax, alpha, beta)
  # int = time interval of timeseries
  W <- length(Temp) - 1
  # Germination
  rG <- rT(Temp, par.G[1:3]) # rate at each hour
  mG <- rG[1:W] + diff(rG)/2 # mean over hour
  mGm <- matrix(mG, nr = W, nc = W) # Repeat gi for each cohort
  mGm[upper.tri(mGm)] <- 0;  # Set to zero for non-valid periods
  eGm <- int*(apply(mGm, 2, cumsum)) # Effective elapsed time
  FGm <- apply(eGm, 2, Ft, pars = par.G[4:5]) # F(t) per cohort
  FG <- rowSums(FGm) # F(t) sum over all cohorts
  #Appressorium formation
  rA <- rT(Temp, par.A[1:3]) # rate at each hour
  mA <- rA[1:W] + diff(rA)/2 # mean over hour
  mAm <- matrix(mA, nr = W, nc = W) # Repeat gi for each cohort
  mAm[upper.tri(mAm)] <- 0;  # Set to zero for non-valid periods
  eAm <- int*(apply(mAm, 2, cumsum)) # Effective elapsed time
  FAm <- apply(eAm, 2, Ft, pars = par.G[4:5]) # F(t) per cohort
  FGc <- matrix(c(0, FG[-W]), nr = W, nc = W, byrow = TRUE) # Scaling for germinated cohorts
  FA <- rowSums(FAm*FGc) # F(t) sum over all cohorts
  c(0,FA) # Leading zero for 0th hour measurement
}

# DATA
# CLR parameters (Tmin, Topt, Tmax, alpha, gamma)
par.G <- c(12.909, 21.410, 30.943, 1/13.363,  1.287) # Germination parameters
par.A <- c(11.000, 11.500, 32.109, 1/19.121,  2.141) # Infection parameters
names(par.G) <- names(par.A) <- c("tmin", "topt", "tmax", "alpha", "beta")
# Meterological timeseries
dat <- read.csv("site1.csv", stringsAsFactors=FALSE)

# DATA PROCESSING
# Convert ISO 8601 date and time to POSIXlt format
dat$time <- gsub("T", " ", dat$time) # remove T
dat$time <- as.POSIXct(strptime(dat$time, format = "%F %T", tz = "UTC"))

# Determine wet periods
# RH = 100 (high humidity alone is insufficient, according to Rayner (1961))
# Dewpoint (assume wet if air temp < dew point)
# Rain > 0 (assume rain wets leaves)
dat$Wet <- as.numeric(dat$AH >= 100 | dat$AT <= dat$ATD | dat$RM > 0)

# Extract start and end of wet periods
wet <- wetper(dat$Wet, min = 6)

# Run infection model over those periods
results <- list()
for(i in 1:nrow(wet)){
  Temp <- dat$AT[wet[i,1]:wet[i,2]] # extract temperatures
  res <- infect(Temp, par.G, par.A, int = 0.5)
  names(res) <- wet[i,1]:wet[i,2] # name with indices
  results[[i]] <- res
}



