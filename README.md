buoyancyModifiedTurbulenceModels
================================

## General

The buoyancy-modified turbulence models are developed to simulate offshore 
and coastal engineering processes. The buoyancy-modified turbulence 
models not only result in a stable wave propagation model without wave 
damping but they also predict the turbulence level inside the flow field 
more accurately in the surf zone.

## References and citing

The buoyancy-modified turbulence models have been developed within the 
PhD thesis of Brecht DEVOLDER at the Department of Civil Engineering at 
Ghent University and KU Leuven, funded by the Research Foundation – 
Flanders (FWO), Belgium (Ph.D. fellowship 1133817N).

If you want to reference the model in your publications, you can use the 
following references in which the implementation and validation details 
are published:  
- Devolder, B., Rauwoens, P., & Troch, P. (2017). Application of a buoyancy-modified *k-ω SST* turbulence model to simulate wave run-up around a monopile subjected to regular waves using OpenFOAM<sup>®</sup>. Coastal Engineering, 125, 81–94. [doi:10.1016/j.coastaleng.2017.04.004](https://doi.org/10.1016/j.coastaleng.2017.04.004).
- Devolder, B., Troch, P., & Rauwoens, P. (2018). Performance of a buoyancy-modified *k-ω* and *k-ω SST* turbulence model for simulating wave breaking under regular waves using OpenFOAM<sup>®</sup>. Coastal Engineering, 138, 49–65. [doi:10.1016/j.coastaleng.2018.04.011](https://doi.org/10.1016/j.coastaleng.2018.04.011).

## Installation

- In a linux terminal, download the package using git:

      git clone https://github.com/BrechtDevolder-UGent-KULeuven/buoyancyModifiedTurbulenceModels.git
      cd buoyancyModifiedTurbulenceModels

- Source your OpenFOAM environment, e.g. for OpenFOAM-3.0.1: 

      source $HOME/OpenFOAM/OpenFOAM-3.0.1/etc/bashrc        

- Go into the folder for the specific OpenFOAM installation (e.g. OF301 for OpenFOAM-3.0.1):

      cd OF301
        
- Compile the source code to build a shared library:

      wmake libso

## How to use

- Include the buoyancyModifiedTurbulenceModels library in system/controlDict:

      libs
      (
          "libbuoyancyModifiedTurbulenceModels.so"
      );
		
- Add the correct turbulence model in constant/turbulenceProperties:

      simulationType  RAS;

      RAS
      {
          RASModel        kOmegaSSTBuoyancy; //kOmegaBuoyancy;
          turbulence      on;
          printCoeffs     on;
      }

- Modify system/fvSchemes:

      ddt(k)
      ddt(omega)
      ...
      div(phi,k)
      div(phi,omega)
		
    to

      ddt(rho,k)
      ddt(rho,omega)
      ...
      div((interpolate(rho)*phi),k)
      div((interpolate(rho)*phi),omega)

## Tutorials

*coming soon...*


## Contributors

- Brecht DEVOLDER  

	> Department of Civil Engineering, Ghent University  
	> Technologiepark 904, B-9052 Zwijnaarde (GENT), BELGIUM  
	> <brecht.devolder@ugent.be>  
	
	> Construction Technology Cluster, Campus Bruges, Department of Civil Engineering, KU Leuven  
	> Spoorwegstraat 12, B-8200 Bruges, BELGIUM  
	> <brecht.devolder@kuleuven.be>  

	
- prof. Peter TROCH  

	> Department of Civil Engineering, Ghent University  
	> Technologiepark 904, B-9052 Zwijnaarde (GENT), BELGIUM  
	> <peter.troch@ugent.be>  

	
- prof. Pieter RAUWOENS  

	> Construction Technology Cluster, Campus Bruges, Department of Civil Engineering, KU Leuven  
	> Spoorwegstraat 12, B-8200 Bruges, BELGIUM  
	> <pieter.rauwoens@kuleuven.be>  
