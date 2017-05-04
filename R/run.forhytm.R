#' Model core of FORHYTM
#'
#' @param t Vector containing hourly values of air temperature [°C]
#' @param p Vector containing hourly values of precipitation [mm]
#' @param u Vector containing hourly values of wind speed [m/s]
#' @param rh Vector containing hourly values of relative humidity (0...1)
#' @param rg Vector containing hourly values of global radiation [W/m2]
#' @param rssd Vector containing daily values of relative sunshine duration (0...1)
#' @param elv Elevation [m asl]
#' @param sfc Size of plant available soil moisture reservoir [mm]
#' @param kappa e-folding of soil evaporation reduction function [days]
#' @param laimax Leaf Area Index at full foliage cover [m2/m2]
#' @param gsLen Length of leafed season, from bud burst to leaf fall [days]
#' @param cBeta Shape parameter of soil moisture recharge [-]
#' @param rsmin Minimum stomatal resistance [s/m]
#' @param jvpd Parameter indicating stomatal sensitivity to VPD [hPa-1]
#' @param topt Optimal temperature for photosynthesis [°C]
#' @param ff Parameter linking fractional canopy cover to LAI
#' @param z Canopy height [m]
#' @param ksimax Parameter linking canopy interception capacity to LAI
#' @param jsm Parameter indicating stomatal sensitivity to relative extractable water [-]
#' @param tl Minimum temperature for photosynthesis [°C]
#' @param th Maximum temperature for photosynthesis [°C]
#' @param extCoeff Canopy light extinction coefficient
#' @param Parameter indicating stomatal sensitivity to radiation [W/m2]
#' @param rsmax Maximum stomatal resistance [s/m]
#' @param alb albedo (0...1)
#' @param lw Typical leaf width [m]
#' @return Returns a data.frame containing daily values of water fluxes (potential transpiration, actual transpiration, soil evaporation, interception evaporation, water yield) and storage (canopy storage, soil moisture, snowpack).
#' @export

run.forhytm <- function(t,p,u,rh,rg,rssd,elv=1190,sfc=150,kappa=15,laimax=4,gsLen=180,cBeta=2,rsmin=180,jvpd=0.025,topt=25,ff=0.386,z=25,ksimax=2,jsm=0.4,tl=0,th=40,extCoeff=0.5,jR=100,rsmax=5000,alb=0.1,lw=0.01){
	
	### Set constants 
	
	# Minimum LAI and fractional cover
	laimin <- 0.2
	vfmin <- 0.2
	
	# Auxiliary variables for the calculation of vapor pressure and transpiration/evaporation
	cp <- 1004
	at <- 0.95
	sigma <- 5.6703E-8
	a <- 3600000
	zeroK <- 273.15
	hundredK <- 373.15
	pZero <- 1013.25
	gammaT <- 0.005
	g <- 9.81
	rAir <- 287.04
	rWat <- 461.5
	rhoWat <- 999.941
	H <- c(13.3185,-1.976,-0.6445,-0.1299)
	pp <- c(-137.0, -75.0,  30.0, 167.0, 236.0, 252.0, 213.0,  69.0, -85.0,-206.0,-256.0,-206.0)

	# Parameters for HBV-light snow model
	tt <- 0
	cfmax <- 0.125
	scf <- 0.9
	cwh <- 0.05
	cfr <- 0.05
	
	
	### Prepare vectors for output
	ndays <- length(rssd$year)	# Number of days, taken from daily input file (RSSD)
	pt <- vector("numeric",ndays)
	at <- vector("numeric",ndays)
	pes <- vector("numeric",ndays)
	aes <- vector("numeric",ndays)
	ei <- vector("numeric",ndays)
	smV <- vector("numeric",ndays)
	siV <- vector("numeric",ndays)
	qV <- vector("numeric",ndays)
	rgV <- vector("numeric",ndays)
	rnV <- vector("numeric",ndays)
	snowV <- vector("numeric",ndays)
	
	### Initialization
	sm <- sfc
	si <- 0
	sso <- 0
	sliq <- 0
	ndrydays <- 0
	hourcount <- 1
	p16 <- rep(1,16)
	e16 <- rep(1,16)
	flp <- 1
	
	# Calculation of phenology dates
	center <- 218
	doyBB <- round(center-(gsLen/2))
	doySen <- round(center+(gsLen/2))
	doyMax <- doyBB+20
	doyFall <- doySen + 20	
	
	### Daily loop
	for(iday in 1:ndays){

		# Calculate day of year
		doy <- as.numeric((strftime(as.Date(paste(padstr(rssd$mon[iday]),padstr(rssd$day[iday]),rssd$year[iday],sep="/"),format="%m/%d/%Y"),format="%j"))) 
	
	## Daily leaf phenology(simplified)
	  	if(doy >= doyMax & doy <= doySen){
	  		lai <- laimax
	  	} else if(doy <= doyBB | doy >= doyFall){
	  		lai <- laimin
	  	} else if(doy > doyBB & doy < doyMax){
	  		lai <- laimin + ((laimax-laimin)/(doyMax-doyBB))*(doy-doyBB)
	  	} else if(doy > doySen & doy < doyFall){
	  		lai <- laimax + ((laimin-laimax)/(doyFall-doySen))*(doy-doySen)
	  	}
	  	vfrac <- min(1,(ff+(lai*0.054)))
	  	vfrac <- max(vfmin,vfrac)
	  	
		#Size of interception storage
		simax <- (ksimax*log10(1+lai))#*vfrac
		
		ssd <- rssd$rssd[iday]
		
		# Temporary vectors for hourly values
		siTemp <- vector("numeric",24)
		smTemp <- vector("numeric",24)
		
		# Initialize daily sums
		cumulrain <- 0
		ptDay <- 0
		atDay <- 0
		pesDay <- 0
		aesDay <- 0
		eiDay <- 0
		qDay <- 0
			  	
	  	for(hour in 1:24){
	### Hourly loop
		# Get meteorological variables for current hour
			ta <- t$t[hourcount]
			prc <- p$p[hourcount]
			u2 <- max(0.01,u$u[hourcount])
			h <- rh$h[hourcount]
			rglob <- rg$r[hourcount]
			
			# Current state of storage elements (interception and soil moisture)
			if(hour==1){
				siH <- siV[iday-1]
				smH <- smV[iday-1]
			} else{
				siH <- siTemp[hour-1]
				smH <- smTemp[hour-1]
			}
			
			if(iday==1 & hour==1){
				siH <- 0
				smH <- sfc
			}
				
	## Basic snow routine (HBV light, Seibert and Vis (2012)) - snow accumulation
	if(ta < tt){
		sso <- sso + (prc*scf)
		prl <- 0
	}
	
	if((sso >= 5) | (ta<tt)) {
		# Snowpack is present - no other source of evaporation/transpiration
		peTrans <- 0
		aTrans <- 0
		peInt <- 0
		aeSoil <- 0
		peSoil <- 0
		
		if(ta >= tt){
			# Temperature above snow threshold - melt
			m <- min(sso,(cfmax*(ta-tt)))
			sso <- sso-m
			sliq <- sliq + m + prc
			sliq.out <- max(0,(sliq-(cwh*sso)))
			sliq <- sliq-sliq.out
			prl <- sliq.out
		} else {
			# Temperature below snow threshold - refreezing
			refr <- min(sliq,(cfr*cfmax*(tt-ta)))
			sliq <- sliq-refr
			sso <- sso+refr
		}
	} else {
		
		# No snow - actual model core
		sso <- 0
		sliq <- 0
		
	## Auxiliary calculations for evaporation - get current and saturated water vapor pressure
    # as well as delta and gamma parameters for the Penman-Monteith equation
      		tk <- ta + zeroK
      		tr <- 1-(hundredK/tk)
      		tp <- tk+gammaT*elv/2
      		PT <- pZero * exp(-g*elv/rAir/tp)
      		l <- 2500800*(zeroK/tk)
      		gamma <- PT*cp/l/(rAir/rWat)
      		rhoAir <- 100*PT/rAir/tk
      		x<-0
      		dx<-0
      		for (iter in 1:4){
        		x <- x + H[iter] * tr^iter
        		dx <- dx + iter * H[iter] * tr^(iter-1)
      		}
      		elx <- pZero * exp(x)
      		del <- hundredK*elx/tk/tk*dx
      		el <- elx*h
      		      		
	## Parameterization of SW/LW radiation
			rs <- rglob*(1-alb)
			rl <- -sigma*((ta+zeroK)^4.)*(0.52-(0.065*sqrt(el)))* (0.23+(0.77*ssd))  # Schulla 1997
			Rn <- rs+rl		   # Net radiation [W/m2]
			ANet <- max(Rn,0)  # Available energy for evaporation and transpiration - split between canopy and soil:
			ACan <- ANet*(1.0-exp(-extCoeff*lai))
			ASoil <- ANet*exp(-extCoeff*lai)
			
      			
    ## Aerodynamic resistance
    		raero <- eraeroCan(u2,lai,vfrac,z,lw)
    		raa <- raero[1]
    		rca <- raero[2]
    		rsa <- raero[3]
    		ra <- raa+rca
	
	## Surface resistance
			laic <- lai/vfrac
			rcPot <- rcJarvis2(laic=laic,rsmin=rsmin,rsmax=rsmax,rg=rglob,sm=smH,sfc=sfc,vp=el,vpsat=elx,tAir=ta,jR=jR,jvpd=jvpd,jsm=jsm,tl=tl,th=th,topt=topt,pa=1)
			rcAct <- rcJarvis2(laic=laic,rsmin=rsmin,rsmax=rsmax,rg=rglob,sm=smH,sfc=sfc,vp=el,vpsat=elx,tAir=ta,jR=jR,jvpd=jvpd,jsm=jsm,tl=tl,th=th,topt=topt,pa=2)

	## Interception filling
	# For very small storages, use a linear storage, otherwise simulate filling following Menzel (1996)
      		SIold <- siH
      		if(simax < 0.1){
        		siH <- SIold + prc*0.8
        		if( siH >= simax){
          			siH <- simax
          			prl <- prc-(simax-SIold)
        		}else{
          			prl <- prc-prc*0.8
          		}
      		}else{
        		siH <- SIold + (simax-SIold)*(1-exp(-1/simax*prc))
        		if( siH >= simax){
          			siH <- simax
          			prl <- prc-(simax-SIold)
        		}else{
          			prl <- prc - (simax-SIold)*(1-exp(-1/simax*prc))
        		}
        	}
        	
# Wet fraction of foliage (Deardorff 1978)
      		if(siH == 0){
        		wetFrac = 0
      		}else{
        		wetFrac = min(1.0,(siH/simax)^0.6666)
      		}
      		
      		# Cumulative effective rainfall over the whole day.
      		# Used to decide whether the current day counts as a rainy day or not
      		cumulrain <- cumulrain + prl
      		if(cumulrain > 0.5){ndrydays <- 0}
	
	## Interception evaporation
			num <- del*ACan + vfrac*rhoAir*cp*(elx-el) / ra
			denom <- del+gamma  # No surface resistance here
			energy <- num/denom								# Expressed in terms of energy
        	peiZero <- (energy * a/rhoWat/l) * wetFrac		# Expressed as a water flux
        	peInt <- max(0.0,min(siH,peiZero))				# Minimum of stored water and potential EInt
        	siH <- siH-peInt								# Update canopy storage
        # Calculate energy used for interception evaporation
        	AInt <- peInt * (rhoWat*l)/a					# Express EInt in terms of energy
        	ACan <- max(0.0, ACan - AInt)					# Energy used for interception evaporation is subtracted from energy available for transpiration
	
	## Potential and actual transpiration
        	num <- del*(max(0.0,ACan)) + vfrac*rhoAir*cp*(elx-el) / ra
        	denom <- del + gamma * (1.0 + rcPot/ra)
        	energy <- num/denom
        	peTrans <- energy * a/rhoWat/l * (1.0-wetFrac)
        	peTrans <- max(0.0,peTrans)	

        	num <- del*(max(0.0,ACan)) + vfrac*rhoAir*cp*(elx-el) / ra
        	denom <- del + gamma * (1.0 + rcAct/ra)
        	energy <- num/denom
        	aTrans <- energy * a/rhoWat/l * (1.0-wetFrac)
        	aTrans <- max(0.0,aTrans)	
        	        	
	## Soil evap
			
			# Update state variables
			if(hour==1){
				p16[2:16] <- p16[1:15]
				e16[2:16] <- e16[1:15]
				p16[1] <- 0
				e16[1] <- 0
			}
			ra <- raa + rsa
			ASoil <- 0.7*ASoil	# 30% of ASoil go to ground heat flux
			num <- del*max(0.0,ASoil) + (1-vfrac)*rhoAir*cp*(elx-el) / ra
			denom <- del + gamma
			energy <- num/denom
			peSoil <- energy * a/rhoWat/l	# Potential soil evaporation (during rainy days)
			# Reduction of soil evaporation as a function of antecedent precipitation and soil evaporation, following Morillas et al. (2013)
			if(cumulrain <= 0.5){
				aeSoil <- peSoil *flp* exp(-ndrydays/kappa)
			} else {
				fZhang <- min(1,(sum(p16)/sum(e16)))
				aeSoil <- peSoil * fZhang
			}
			p16[1] <- p16[1] + prl
			e16[1] <- e16[1] + peSoil
			
			}  # Snow / no snow
			
			
	## Soil moisture refilling/depletion
	# Explicit Euler iteration
      		IM <- max(3,((prl*10)/30)+1.5)	# Number of iterations as a function of precipitation (more iterations needed if precipitation is high)
      		ts <- 1/IM
      		sme <- smH
      		eSoil <- 0
      		eTrans <- 0
      		q <- 0  # Runoff
      		prl <- prl*ts
      		for (tstep in 1:IM){
        		esIter <- aeSoil*ts
        		etIter <- aTrans*ts
        		dsuzIter <- prl*((sme/sfc)^cBeta)
        		sme <- sme - etIter - esIter + (prl-dsuzIter)
        		if( sme > sfc){
        			dsuzIter = dsuzIter + (sme-sfc)
          			sme <- sfc
          		}
        		if( sme < 0){
          			etIter <- etIter + esIter + sme
          			esIter <- 0.0
          			sme <- 0.0
          		}
        		eSoil = eSoil + esIter
        		eTrans = eTrans + etIter
        		q = q + dsuzIter
      		}  # End of sub-hourly iteration 
      		smTemp[hour] <- sme
      		siTemp[hour] <- siH     	
      		# Update daily totals
      		ptDay <- ptDay + peTrans
      		atDay <- atDay + eTrans
      		pesDay <- pesDay + peSoil
      		aesDay <- aesDay + aeSoil
      		eiDay <- eiDay + peInt
      		qDay <- qDay + q	
      		
			hourcount <- hourcount+1
		} # End of hourly loop
		
		# For soil evaporation drying function
		if(cumulrain <= 0.5){
			ndrydays <- ndrydays + 1
		} else{
			ndrydays <- 0
			flp <- fZhang
		}
		
		# Vectors of daily values (will be written as output)
		pt[iday] <- ptDay
		at[iday] <- atDay
		pes[iday] <- pesDay
		aes[iday] <- aesDay
		ei[iday] <- eiDay
		qV[iday] <- qDay
		smV[iday] <- sme
		siV[iday] <- siH
		snowV[iday] <- sso

		
	} # iday (Daily loop)
	
	out <- data.frame(rssd$year,rssd$mon,rssd$day,pt,at,aes,ei,qV,smV,siV,snowV)
	names(out) <- c("year","mon","day","pt","at","esoil","eint","wy","ssm","si","snow")
	return(out)
}
