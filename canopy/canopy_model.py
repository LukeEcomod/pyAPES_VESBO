# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:01:50 2017

@author: slauniai

******************************************************************************
CanopyModel:

Gridded canopy and snow hydrology model for SpatHy -integration
Uses simple schemes for computing water flows and storages within vegetation
canopy and snowpack at daily or sub-daily timesteps.

(C) Samuli Launiainen, 2016-2017
last edit: 1.11.2017: Added CO2-response to dry_canopy_et
******************************************************************************

"""
import numpy as np
import canopy_utils as cu
eps = np.finfo(float).eps


class CanopyModel():
    def __init__(self, cpara, lai_conif=None, lai_decid=None, cf=None, hc=None,  cmask=np.ones(1), outputs=False):
        """
        initializes CanopyModel -object

        Args:
            cpara - parameter dict:
                {'spatial': no, 'wmaxsnow': 1.0, 'lon': 24.21, 'cf': 0.7,
                'wmax': 0.5, 'albedo': 0.1, 'swe': 40.0,pmax': 20.0,
                'q50': 30.0, 'lai': 4.0, 'kfreeze': 5.79e-06, 'elev': 100.0,
                'lat': 60.38, 'hc': 22.0, 'dt': 86400.0, 'kmelt': 2.8934e-05,
                'ka': 0.6, 'f': 0.4, 'emi': 0.98, 'kp': 0.6, 'gsref': 0.0018,
                'w': 0.0, 'clump': 0.7}
            lai - (conifer) leaf area index grid; if spatial: yes
            lai_decid - (deciduous) leaf area index grid
            cf - canopy closure grid
            hc - canopy height grid
            cmask - catchment mask (when lai=None)
            outputs - True saves output grids to list at each timestep

        Returns:
            self - object
        self.
            LAI - leaf area index [m2 m-2]
            hc - canopy heigth [m]
            cf - closure [-]
            physpara - physiologic parameters dict
            phenopara - phenology param. dict
            cpara - copy of inputs, for testing snow model
            Wmax - maximum water storage capacity [m]
            WmaxSnow -  maximum snow storage capacity [m]
            Kmelt - degree day snowmelt factor [m degC-1 s-1]
            Kfreeze - degree day freezing factor [m degC-1 s-1]
            R - snow water retention coeff [-]

            State variables:
            X - 'phenological state', i.e. delayed temperature [degC]
            W - canopy water or snow storage [m]
            SWE - snow water equivalent [m]
            SWEi - snow water as ice [m]
            SWEl - snow water as liquid [m]
            ddsum - degree-day sum [degC]
        """
        epsi = eps  #0.01  --- Miksi olit näin suuri?

        self.Lat = cpara['lat']
        self.Lon = cpara['lon']

        # physiology: transpi + floor evap
        self.physpara = {'q50': cpara['q50'],
                         'gsref_conif': cpara['gsref_conif'],
                         'gsref_decid': cpara['gsref_decid'],
                         'kp': cpara['kp'],
                         'f': cpara['f'], 'rw': cpara['rw'],
                         'rwmin': cpara['rwmin'], 'ga': cpara['ga']}

        # phenology
        self.phenopara = {'Smax': cpara['smax'], 'Xo': cpara['xo'],
                          'tau': cpara['tau'], 'fmin': cpara['fmin'],
                          'ddo': cpara['ddo'], 'ddur': cpara['ddur'],
                          'sdl': cpara['sdl'], 'sdur': cpara['sdur'],
                          'lai_decid_min': cpara['lai_min']}

        if lai_conif is not None:
            # print '**** CanopyModel - stands from lai data *******'
            self._LAIconif = lai_conif + epsi  # m2m-2

            if lai_decid is not None:
                self._LAIdecid_max = lai_decid + epsi
            else:
                self._LAIdecid_max = np.zeros(np.shape(lai_conif)) + epsi

            self.hc = hc + epsi
            self.cf = cf + epsi
            # self.cf = 0.1939 * ba / (0.1939 * ba + 1.69) + epsi
            # canopy closure [-] as function of basal area ba m2ha-1;
            # fitted to Korhonen et al. 2007 Silva Fennica Fig.2

        else:  # spatially constant properties are used, given in cpara
            # print '**** CanopyModel - stands with constant values *******'
            self._LAIconif = cpara['lai_conif'] * cmask
            self._LAIdecid_max = cpara['lai_decid'] * cmask
            self.hc = cpara['hc'] * cmask
            self.cf = cpara['cf'] * cmask

        # current deciduous and total LAI
        self._LAIdecid = self.phenopara['lai_decid_min'] * self._LAIdecid_max
        self.LAI = self._LAIconif + self._LAIdecid

        gridshape = np.shape(self.LAI)
        # self.cf = np.ones(gridshape)*0.9
        # deciduous leaf growth stage
        self._growth_stage = 0.0
        self._senesc_stage = 0.0

        # senescence starts at first doy when daylength < self.phenopara['sdl']
        doy = np.arange(1, 366)
        dl = cu.daylength(self.Lat, self.Lon, doy)

        ix = np.max(np.where(dl > self.phenopara['sdl']))
        self.phenopara['sso'] = doy[ix]  # this is onset date for senescence
        del ix

        # self.cpara = cpara  # added new parameters self.cpara['kmt'],
        # self.cpara['kmr'] here for testing radiation-based snow melt model
        # snow
        self.snowpara = {'wmax': cpara['wmax'], 'wmaxsnow': cpara['wmaxsnow'],
                         'kmelt': cpara['kmelt'], 'kfreeze': cpara['kfreeze'],
                         'retention': cpara['retention'], 'Tmin': cpara['Tmin'],
                         'Tmax': cpara['Tmax'], 'Tmelt': cpara['Tmelt'],
                         'kp': cpara['kp']}

        # --- for computing aerodynamic resistances
        self.zmeas = 2.0
        self.zground = 0.5
        self.zo_ground = 0.01
        self.soilrp = 300.0

        # --- state variables
        self.X = 0.0
        self.W = np.minimum(cpara['w'], cpara['wmax'] * self.LAI)
        self.SWE = {'SWE': cpara['swe'] * np.zeros(gridshape),
                    'SWEi': cpara['swe'] * np.zeros(gridshape),
                    'SWEl': np.zeros(gridshape)}
        self.DDsum = 0.0

        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'PotInf': [], 'Trfall': [], 'Interc': [], 'Evap': [], 'ET': [],
                            'Transpi': [], 'Efloor': [], 'SWE': [], 'LAI': [],
                            'LAIfract': [], 'Mbe': [], 'LAIdecid': [], 'erate': [], 'Unload': [],
                            'fact': []}

    def _run(self, doy, dt, Ta, Prec, Rg, Par, VPD, U=2.0, CO2=380.0, Rew=1.0, beta=1.0, P=101300.0):  # Tänne forcing niinkuin bryotypessä
        """
        Runs CanopyModel instance for one timestep
        IN:
            doy - day of year
            dt - timestep [s]
            Ta - air temperature  [degC], scalar or (n x m) -matrix
            prec - precipitatation rate [m s-1]
            Rg - global radiation [W m-2], scalar or matrix
            Par - photos. act. radiation [W m-2], scalar or matrix
            VPD - vapor pressure deficit [kPa], scalar or matrix
            U - mean wind speed at ref. height above canopy top [m s-1], scalar or matrix
            CO2 - atm. CO2 mixing ratio [ppm]
            Rew - relative extractable water [-], scalar or matrix
            beta - term for soil evaporation resistance (Wliq/FC) [-]
            P - pressure [Pa], scalar or matrix
        OUT:
            updated CanopyModel instance state variables
            flux grids PotInf, Trfall, Interc, Evap, ET, MBE [m]
        """

        # Rn = 0.7 * Rg #net radiation
        Rn = np.maximum(2.57 * self.LAI / (2.57 * self.LAI + 0.57) - 0.2,
                        0.55) * Rg  # Launiainen et al. 2016 GCB, fit to Fig 2a


        f = self.physpara['f'] * np.ones(np.shape(self.SWE['SWE']))
        f[self.SWE['SWE'] > 0] = eps  # in case of snow, neglect forest floor evap

        """ --- update phenology: self.ddsum & self.X ---"""
        self._degreeDays(Ta, doy)
        fPheno = self._photoacclim(Ta)

        """ --- update deciduous leaf area index --- """
        laifract = self._lai_dynamics(doy)

        """ --- aerodynamic conductances --- """
        Ra, Rb, Ras, ustar, Uh, Ug = cu.aerodynamics(self.LAI,
                                                     self.hc,
                                                     U,
                                                     w=0.01,
                                                     zm=self.zmeas,
                                                     zg=self.zground,
                                                     zos=self.zo_ground)

        """ --- interception, evaporation and snowpack --- """
        self.W, self.SWE, PotInf, Trfall, Evap, Interc, MBE, erate, unload, fact = \
            cu.canopy_water_snow(W=self.W,
                                 SWE=self.SWE,
                                 LAI=self.LAI,
                                 cf=self.cf,
                                 snowpara=self.snowpara,
                                 dt=dt,
                                 T=Ta,
                                 Prec=Prec,
                                 AE=Rn,
                                 VPD=VPD,
                                 Ra=Ra,
                                 U=2.0)

        """--- dry-canopy evapotranspiration [mm s-1] --- """
        Transpi, Efloor, Gc = cu.dry_canopy_et(LAI=self.LAI,
                                               LAIconif=self._LAIconif,
                                               LAIdecid=self._LAIdecid,
                                               physpara=self.physpara,
                                               soilrp=self.soilrp,
                                               D=VPD,
                                               Qp=Par,
                                               AE=Rn,
                                               Ta=Ta,
                                               Ra=Ra,
                                               Ras=Ras,
                                               CO2=CO2,
                                               f=f,
                                               Rew=Rew,
                                               beta=beta,
                                               fPheno=fPheno)

        Transpi = Transpi * dt
        Efloor = Efloor * dt
        ET = Transpi + Efloor

        # append results to lists; use only for testing small grids!
        if hasattr(self, 'results'):
            self.results['PotInf'].append(PotInf)
            self.results['Trfall'].append(Trfall)
            self.results['Interc'].append(Interc)
            self.results['Evap'].append(Evap)
            self.results['ET'].append(ET)
            self.results['Transpi'].append(Transpi)
            self.results['Efloor'].append(Efloor)
            self.results['SWE'].append(self.SWE['SWE'])
            self.results['LAI'].append(self.LAI)
            self.results['LAIfract'].append(laifract)
            self.results['Mbe'].append(np.nanmax(MBE))
            self.results['LAIdecid'].append(self._LAIdecid)
            self.results['erate'].append(erate)
            self.results['Unload'].append(unload)
            self.results['fact'].append(fact)

        # return state and fluxes in dictionary
        state = {"SWE": self.SWE['SWE'],
                 "LAI": self.LAI,
                 "LAIdecid": self._LAIdecid,
                 }
        fluxes = {'PotInf': PotInf,
                  'Trfall': Trfall,
                  'Interc': Interc,
                  'CanEvap': Evap,
                  'Transp': Transpi, 
                  'Efloor': Efloor,
                  'MBE': np.nanmax(MBE)
                  }

        return fluxes, state

    def _degreeDays(self, T, doy):
        """
        Calculates and updates degree-day sum from the current mean Tair.
        INPUT:
            T - daily mean temperature (degC)
            doy - day of year 1...366 (integer)
        """
        To = 5.0  # threshold temperature
        if doy == 1:  # reset in the beginning of the year
            self.DDsum = 0.
        else:
            self.DDsum += np.maximum(0.0, T - To)

    def _photoacclim(self, T):
        """
        computes new stage of temperature acclimation and phenology modifier.
        Peltoniemi et al. 2015 Bor.Env.Res.
        IN: object, T = daily mean air temperature
        OUT: fPheno - phenology modifier [0...1], updates object state
        """

        self.X = self.X + 1.0 / self.phenopara['tau'] * (T - self.X)  # degC
        S = np.maximum(self.X - self.phenopara['Xo'], 0.0)
        fPheno = np.maximum(self.phenopara['fmin'],
                            np.minimum(S / self.phenopara['Smax'], 1.0))
        return fPheno

    def _lai_dynamics(self, doy):
        """
        Seasonal cycle of deciduous leaf area

        Args:
            self - object
            doy - day of year

        Returns:
            none, updates state variables self.LAIdecid, self._growth_stage,
            self._senec_stage
        """
        lai_min = self.phenopara['lai_decid_min']
        ddo = self.phenopara['ddo']
        ddur = self.phenopara['ddur']
        sso = self.phenopara['sso']
        sdur = self.phenopara['sdur']

        # growth phase
        if self.DDsum <= ddo:
            f = lai_min
            self._growth_stage = 0.
            self._senesc_stage = 0.
        elif self.DDsum > ddo:
            self._growth_stage += 1.0 / ddur
            f = np. minimum(1.0, lai_min + (1.0 - lai_min) * self._growth_stage)

        # senescence phase
        if doy > sso:
            self._growth_stage = 0.
            self._senesc_stage += 1.0 / sdur
            f = 1.0 - (1.0 - lai_min) * np.minimum(1.0, self._senesc_stage)

        # update self.LAIdecid and total LAI
        self._LAIdecid = self._LAIdecid_max * f
        self.LAI = self._LAIconif + self._LAIdecid
        return f





