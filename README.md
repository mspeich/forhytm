# forhytm
A dynamic water balance model for numerical experiments

This repository contains the source code of the toy model FORHYTM, as described in a paper in preparation by Speich et al. In a first step, the aim is to make the code available to the reviewers of this paper. Later, this repository shall provide the source code to anyone interested. Along with the source code, a sample dataset to run the model with is provided.

# Package installation

The easiest way to install this package is a direct installation from Github, which can be done by executing the following code in R:

```R
install.packages("devtools")
library(devtools)
install_github("mspeich/forhytm")
library(forhytm)
```
This way, all functions and the sample data are loaded directly into the R workspace.

Alternatively, the source code files may be downloaded from the folder ```R```. The sample data may be downloaded from the folder ```data``` and read into an R session using the function ```load```.

# Sample data

The sample data contains five years (2002-2007) of meteorological data measured at Erlenbach, Switzerland. The aim is to provide some data with which the model can be run and tested. This dataset contains the following variables:

Hourly data:
- Precipitation [mm]
- Air temperature [°C]
- Relative humidity (0 ... 1)
- Global radiation [W/m2]
- Wind speed [m/s]

Daily data:
- Relative sunshine duration (0...1) (obtained from actual SSD, expressed as a fraction of day length as calculated with the R package sirad)

If you wish to use the Erlenbach data for your own research, please see: http://www.wsl.ch/fe/gebirgshydrologie/testgebiet_alptal/data/index_EN

This is not the same dataset that was used in the paper, as this data may not be shared publicly. However, it is available from the Swiss Federal Office of Meteorology and Climatology MeteoSwiss: http://www.meteoswiss.admin.ch/home/services-and-publications/advice-and-service/datenportal-fuer-lehre-und-forschung.html

For a quick-and-dirty variation of climatic conditions, the meteorological variables may be modified. For example, temperature may be increased by a given offset ```deltaT``` by entering:

```R
erl.temp$t <- erl.temp$t+deltaT
```

Precipitation may be scaled by a given factor ```fPrec``` by typing:
```R
erl.prec$p <- erl.prec$p+fPrec
```

# Running the model
Once the package is installed and loaded, the model may be run with all parameters set to a standard value with the following line:
```R
sim <- run.forhytm(t=erl.temp,p=erl.prec,rh=erl.rh,rg=erl.rg,u=erl.wind,rssd=erl.rssd)
```
This returns an object ```sim``` ```(data.frame)``` containing the daily simulation results. To vary parameter values, specify them as arguments in the function call, e.g.:
```R
sim <- run.forhytm(t=erl.temp,p=erl.prec,rh=erl.rh,rg=erl.rg,u=erl.wind,rssd=erl.rssd,sfc=100,laimax=5)
```
A list of all parameters can be found on the function description page, accessible by typing ```?run.forhytm``` in the R console.
