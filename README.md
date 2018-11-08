# README
Branch for cleaned model codes for pyAPES development

## Migration to Python3
https://docs.python.org/3.0/whatsnew/3.0.html
Code is mainly written so that it works in both python 2.7 and 3.6.6
Code should run python version >3.5
* main things:
  - in python3, print is a function; brackets are needed
  - iterating a list: if a index is needed use enumarate() then list(range(len(foo))) is not needed. e.g. 'for index, item in foo:' instead of 'for index in list(range(len(foo))):'
  - iterating a dict: if you are iterating through dict.items(), dict.keys(), and dict.values() AND adding/deleting item from dict wrap it inside a list(). e.g. 'for key in list(dict.keys()): del dict[key]'.   
  - absolute import vs. relative import: from .foo import spam if you want to import from same package.

### SoilProfile
Water flow:
* Equilibrium within vertical column during each timestep 
* OR Richards 1D equation 
! problem with infiltration (pF curves and unsaturated hydraulic conductivity)
! description of preferential flow?

Heat flow:
* solves heat conduction
-> soil temperature and ice content 
! neglects heat convection (heat transfer with water flow)

### Canopy
Multilayer canopy description
* Radiation model: canopy SW (PAR&NIR) and LW including multiple scattering in horizontally homogenous porous media Zhao & Qualls (2005, 2006), sunlit/shade leaves
! range of zenith angle?
* Interception model: interception of rainfall and snow following approach by Tanaka (2002, Ecol. Mod.) 
-> solves wet leaf temperature
! restrict to snow surface level? sublimation of snow from canopy should be checked
* Leaf gas-exchange: Photosynthesis calculated based on biochemical model of Farquhar et al. (1980) coupled with various stomatal control schemes (Medlyn, Ball-Woodrow-Berry, Hari, Katul-Vico et al.)
-> solves dry leaf temperature
* Momentum, H2O, CO2 and T within canopy: 1st-order closure model (sources/sinks: evaporation, transpiration, photosynthesis, respiration, sensible heat etc.)

Forest floor and snowpack:
* Moss cover (present during snow free periods): Interceps rainfall and evaporates interception storage, CO2 exchange (respiration and photo?)
! capillary flux not included
* Bare soil surface energy balance
-> solves soil surface temperature 
* Soil respiration
! simplified?
* Snow model: Temperature-based snow accumulation and melt
In future two layer energy balance snow scheme + soil freezing thawing (FEMMA)?
		
### Forcing
Lettosuo (2010-2018)
Hyytiälä (1997-2016)
see tools/dataprocessing_scripts

### Todos..

* Documentation!
* Updating bryophyte model (energy and water)
* Description of soil respiration? branches etc..
* Marklund biomass functions ok for drained peatlands?
* Feedbacks from soil to canopy
* Parallelize model run and result writing (netCDF4) as in Climoss
* Running sensitivity analysis easily?
