	#' Aerodynamic resistance terms
	#'
	#' @param u2 Wind speed at 2m height [m/s]
	#' @param lai Leaf Area Index [m2/m2]
	#' @param fvc Fractional vegetation (canopy) cover
	#' @param h Canopy height [m]
	#' @param lw Typical leaf width [m]
	#' @return Returns a vector of 3 aerodynamic resistance terms: above-canopy (raa), within-canopy (rca) and sub-canopy (rsa)
	#' @export
eraeroCan <- function(u2,lai,fvc,h,lw){


#################################
# Aerodynamic resistance for canopy. Simplified after Shuttleworth-Wallace (1985) and Guan and Wilson (2009)
#
# References:
#
# Shuttleworth, W.J. and Wallace, J.S. (1985): Evaporation from sparse crops-an energy combination theory.
# Q. J. Roy. Meteor. Soc. 111(469): 839-855
# 10.1002/qj.49711146910
#
# Guan, H. and Wilson, J.L. (2009): A hybrid dual-source model for potential evaporation and transpiration partitioning.
# J. Hydrol. 377(3-4): 405-416
# 10.1016/j.jhydrol.2009.08.037
#
# Iritz, Z. et al. (1999): Test of a modified Shuttleworth–Wallace estimate of boreal forest evaporation.
# Agric. For. Meteorol. 98-99: 605-619
# 10.1016/S0168-1923(99)00127-6
#
# Zhou et al. (2006): Estimating potential evapotranspiration using Shuttleworth–Wallace model
# and NOAA-AVHRR NDVI data to feed a distributed hydrological model over the Mekong River basin.
# J. Hydrol. 317(1-2): 151-173
# 10.1016/j.jhydrol.2005.11.013

# Inputs:
# - u2: Wind speed at 2 m above ground [m s-1]
# - lai: Leaf Area Index [m2 m-2]
# - fvc: Fractional tree canopy cover [m2 m-2]
# - h: Canopy height [m]
# - lw: Characteristic leaf width [m]

# Outputs:
# - raa: above-canopy aerodynamic resistance [s m-1]
# - rca: within-canopy aerodynamic resistance [s m-1]
# (below-canopy aerodynamic resistance is not calculated, as it is not needed for modeling soil evaporation in this implementation)

# Set constants
	k <- 0.4	# von Kármán's constant
	critLai <- 4	# LAI at which canopy is considered fully closed (Shuttleworth-Wallace)
	zmS <- 0.01		# Soil roughness length for momentum transfer [m]
	zhS <- 0.01		# Soil roughness length for heat/water vapor transfer [m]
	zb <- 296.97	# Height of internal boundary layer (for wind speed conversion) [m]
	zow <- 0.005	# Assumed roughness length at weather stations [m]
	b <- 0.01		# Proportionality constant for rca
	
  # Set aerodynamic parameters as a function of canopy height (all in m)
  	z <- h+2.0			# Reference height above canopy
 	zm <- 0.1*h			# Roughness length for momentum transfer
  	zh <- 0.1*zm		# Roughness length for heat/water vapor transfer
  	d <- 0.6*h			# Zero plane displacement height
  	
 # Get wind at reference level (uz) and canopy level (uh)
 # See Zhou et al. 2006, eqs. 27 and 28
  	uz <- u2*(log(zb/zow)/log(zb/zm))*(log((z-d)/zm)/log(2/zow))
  	uhA <- u2*(log(zb/zow)/log(zb/zm))*(log((h-d)/zm)/log(2/zow))
  	uh <- fvc*uhA + (1-fvc)*uz
  	
# Above-canopy aerodynamic resistance
  	raaF <- (log((z-d)/zm)*log((z-d)/zh))/((k^2)*uz)
  	raaZ <- (log(z/zmS)*log(z/zhS))/((k^2)*uz) - ((log(h/zmS)*log(h/zhS))/((k^2)*uz))
  	
# Sub-canopy aerodynamic resistance
	rsaF <- (log(h/zmS)*log(h/zhS))/((k^2)*uhA)
	rsaZ <- ((log(h/zmS)*log(h/zhS))/((k^2)*uz))
	
  	if(lai >= 4){  # Interpolation following Shuttleworth&Wallace
    	raa = raaF
    	rsa <- rsaF
  	}else if (lai == 0){
    	raa = raaZ
    	rsa <- rsaZ
  	}else{
    	raa = 0.25*lai*raaF + 0.25*(4-lai)*raaZ
    	rsa = 0.25*lai*rsaF + 0.25*(4-lai)*rsaZ
    }
    
# Intra-canopy aerodynamic resistance
# See Iritz et al. 1999
  	rca = 1/(b*sqrt(uh/lw))
  	
  	return(c(raa,rca,rsa))
  	

  	
 }
 
#==================================================================================================================