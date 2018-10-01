# CCFPeat project's README
kokeilu

### Soil_model
Water flow:
* Equilibrium within vertical column during each timestep 
* OR Richards 1D equation 
! problem with infiltration (pF curves and unsaturated hydraulic conductivity)
! description of preferential flow?

Heat flow:
* solves heat conduction
-> soil temperature and ice content 
! neglects heat convection (heat transfer with water flow)

### Canopy_model
Simple canopy description (based on SpaFHy)
* Big leaf interception 
! Interception by shrubs separately?
* Interception evaporation and tranpiration besad on Penman-Montieth
* Wind speed amd net radiation at ground and canopy + boundary layer conductances

OR Multilayer canopy description
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

* Updating bryophyte model (energy and water)
* Description of soil respiration? branches etc..
* Marklund biomass functions ok for drained peatlands?
* Feedbacks from soil to canopy
* Parallelize model run and result writing (netCDF4) as in Climoss
* Running sensitivity analysis easily?
