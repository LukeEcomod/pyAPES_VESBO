# CCFPeat project's README

Soil_model
	Water flow 
		Equilibrium within vertical column during each timestep
			OR
		Richards 1D equation 
			! problem with infiltration (pF curves and unsaturated hydraulic conductivity)
			! description of preferential flow?
	Heat flow
		As in Climoss, not tested yet
Canopy_model
	SpatHy canopy component 
		Simple schemes for computing water flows and storages within vegetation canopy and snowpack
		! Feedback from soil to canopy not yet included (moisture restriction)
		
Forcing included for Lettosuo (from FMI Jokioinen) and Hyytiälä + evaluation data 

! Writing results to file (netCDF?)
! Checking units in "documentation"