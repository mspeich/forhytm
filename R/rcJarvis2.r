#' Calculation of canopy resistance
#'
#' @param laic Clumped LAI (LAI times FVC)
#' @param rsmin Minimum stomatal resistance [s/m]
#' @param rsmax Maximum stomatal resistance [s/m]
#' @param rg Global radiation [W/m2]
#' @param sm Current plant available soil moisture [mm]
#' @param sfc Size of plant available soil moisture reservoir [mm]
#' @param vp Atmospheric water vapor pressure [hPa]
#' @param vpsat Saturated vapor pressure [hPa]
#' @param tAir Air temperature [°C]
#' @param jR Parameter indicating stomatal sensitivity to radiation [W/m2]
#' @param jvpd Parameter indicating stomatal sensitivity to VPD [hPa-1]
#' @param jsm Parameter indicating stomatal sensitivity to relative extractable water [-]
#' @param tl Minimum temperature for photosynthesis [°C]
#' @param th Maximum temperature for photosynthesis [°C]
#' @param topt Optimal temperature for photosynthesis [°C]
#' @param pa Control parameter: Calculate potential transpiration if pa==1, and actual transpiration of pa==2
#' @export

rcJarvis2 <- function(laic,rsmin,rsmax,rg,sm,sfc,vp,vpsat,tAir,jR,jvpd,jsm,tl,th,topt,pa){
#################################
# Canopy surface resistance following Jarvis-Stewart.
#
# References:
#
# Stewart, J.B. (1988): Modelling surface conductance of pine forest.
# Agric. For. Meteorol. 43:19-35
# 10.1016/0168-1923(88)90003-2
#
# Task Committee on Hydrology Handbook of Management Group D of ASCE (1996): Hydrology Handbook. Second edition.
# Chapter 4: Evaporation and Transpiration. American Society of Civil Engineers, New York.


# Inputs:
# - laic: clumped LAI
# - rsmin: minimum surface resistance
# - rg: global radiation [W m-2]
# - sm: plant available soil moisture [mm]
# - sfc: size of plant available soil moisture reservoir [mm]
# - vp: Vapor pressure [hPa]
# - vpsat: Saturated vapor pressure [hPa]
# - tAir: Air temperature [°C]
# - pa: Control parameter indicating whether rc for potential (pa==1) or actual transpiration (pa==2) should be calculated
# - land: land cover type (3: coniferous forest; 4: broadleaved forest; 5: mixed forest)

# Output:
# - rs: Canopy resistance [s m-1]

# Set constants (Jarvis-Stewart parameters)
	# jR <- 100
	# rsmax <- 5000
	# jv1 <- 0.05
	# jv2 <- 15
	# jsm <- 6.7
	# tl <- 0
	# topt <- 18
	# th <- 40
	
# 1st reduction function: global radiation
  	f1 <- (1000/(1000+jR)) * ((rg+jR)/rg)
  
# 2nd reduction function: Air temperature
  	if(tAir <= tl){
    	f2 <- 20  # ~ the value it would take with tAir==tl in Stewart parameterization (tl = 0)
  	}else if(tAir >= th){
  		f2 <- 20
  	}else{
    	a <- (th-topt)/(topt-tl)
    	f2 <- 1/(((tAir-tl)*((th-tAir)^a))/((topt-tl)*((th-topt)^a)))
  	}

# 3rd reduction function: vapor pressure deficit
  	if (pa==1){
    	f3 <- 1
  	}else{
    	vpd <- vpsat-vp
      	gvpd <- max(0.001,(1-(jvpd*vpd)))
      	f3 <- 1/gvpd
	}
	
	# 4th reduction function: soil moisture
  	if(pa==1){
    	f4 <- 1
  	}else{
    	#f4 <- 1/(1-exp(-jsm*(sm/sfc)))
    	if(sm/sfc >= jsm){
    		f4 <- 1	# No stress above REWc
    	} else {
    		f4 <- 1/(sm/(jsm*sfc))
    	}
    }
    

  	
   	rs <- (rsmin/laic) * (f1*f2*f3*f4)
  	rs <- max(rs,rsmin/laic)
  	rs <- min(rs,rsmax)
  	
  	return(rs)
}
 #==================================================================================================================