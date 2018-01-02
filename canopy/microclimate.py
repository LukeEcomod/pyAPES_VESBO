# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:24:41 2017

@author: slauniai
"""
import os
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import matplotlib.pyplot as plt

from radiation import solar_dec_azimuth, SIGMA, NT  # , RAD_TO_DEG, DEG_TO_RAD
from radiation import canopy_sw_ZhaoQualls as sw_radi
from radiation import canopy_lw_ZhaoQualls as lw_radi

from micromet import wind_profile_canopy as wind_profile
from micromet import bulk_aerodynamic_conductance, penman_monteith, vpd_from_rh

eps = np.finfo(float).eps  # machine epsilon

class Model():
    def __init__(self, properties, dt, canopy_grid):
    # def __init__(self, dt, loc, LAI, grid, radi_para, phys_para, beta, wmax):
        """
        initialize model for simple estimates of microclimate within canopy
        
        Args:
            dt - timestep (s)
            loc [dict] - site location
                'Lat' (deg)
                'Lon' (deg)
                'Elev' (m) - not used now
            LAI [float] - leaf-area index, 1-sided (m2m-2)
            grid [dict] - computation grid
                'Nlayers'
                'zmax' - height of upper boundary (m)
                'hc' - canopy height (m)
                
            radi_para - radiation parameters
            'Clump' - clumping index
            'leaf_angle' - leaf angle distribution parameter x
            'Par_albedo' - effective leaf/shoot Par albedo
            'Ni_albedo' - effective leaf/shoot Nir albedo

            phys_para - physiologic parameters
                'gsref' - max. leaf-scale stomatal conductance of CO2 (ms-1)
                'q50' - half-saturation rate of leaf light response (Wm-2 of PAR)
                'kp' - attenuation coefficient for par within canopy (-)
                'rw', 'rwmin' - drought-response parameters
            
            beta - parameter of exponential wind profile
            wmax - maximum canopy water storage per unit LAI, mm
        """
        self.properties = properties
        self.dt = dt

        loc = properties['loc']
        LAI = properties['LAI']
        grid = properties['grid']
        radi_para = properties['radiation']

        physpara = properties['ecophysiology']

        # --- site location
        self.Lat = loc['Lat']
        self.Lon = loc['Lon']
        self.Elev = loc['Elev']

        # --- flow field parameter
        self.beta = properties['beta']
        # --- canopy water storage
        self.W = 0.0  # water storage of canopy (mm)
        self.wmax = properties['wmax']  # mm / unit of LAI

        # --- above-ground computation grid
        # self.Nlayers = grid['Nlayers']
        self.z = np.linspace(canopy_grid[0], canopy_grid[1], canopy_grid[2])  # grid, m
        self.Nlayers = len(self.z)
        self.dz = self.z[1] - self.z[0]  # gridsize, m
        self.hc = grid['hc']  # canopy height, m
        self.ones = np.ones(self.Nlayers)

        self.LAI = LAI
        self.LAIz = np.zeros(self.Nlayers)
        self.LAIz[0:-1] = LAI / (self.Nlayers - 1)  # leaf area index of each layer

        # --- radiation parameters
        self.Clump = radi_para['Clump']
        self.leaf_angle = radi_para['leaf_angle']
        self.Par_albedo = radi_para['Par_alb']
        self.Nir_albedo = radi_para['Nir_alb']
        # these are not coming from parameters, but assigned in an event loop
        # self.soil_Par_albedo = radi_para['soil_Par_alb']
        # self.soil_Nir_albedo = radi_para['soil_Nir_alb']

        # --- for  computing aerodynamic resistance (from Leuning et al 2008 WRR)
        self.zm = self.hc + 8.0  # wind speed measurement height at hyytiälä ~23m
        self.d = 2. / 3. * self.hc  # displacement height (m)
        self.zom = 0.123 * self.hc  # roughness lengths for momentum and H2O
        self.zov = 0.1 * self.zom

        # physiologic parameters (transpiration rate), uses big-leaf approximation
        self.physpara = physpara


    def transpiration_rate(self, Qp, AE, D, T, Ga=0.04, f=0.0, Rew=1.0, fPheno=1.0):
        """
        Computes transpiration rate using Penman-Monteith equation and simple
        formulation for canopy conductance.
        IN:
           self - object
           Qp - PAR in Wm-2
           AE - available energy in Wm-2
           D - vpd in kPa 
           T - air temperature (degC)
           Ga - aerodynamic conductance (ms-1)
           f - forest floor evaporation ratio [-] (f=1.0 for unlimited rate, 
                                                   0 sets forest floor component to zero and 
                                                   computes only vegetation transpiration)
           Rew - relative extractable water [-]
           fPheno - phenology modifier [-]
        OUT:
           Transpi - transpiration rate (mm/s = kg m-2 s-1)
           Gs - surface conductance (m/s)
           Gc - canopy conductance (m/s)
        NOTE:  Gc is canopy condutance (integrated stomatal conductance) and
            Gs bulk ecosystem conductance (canopy + soil/forest floor)
        
        SOURCES:
        Launiainen et al. (2016). Do the energy fluxes and surface conductance
        of boreal coniferous forests in Europe scale with leaf area?
        Global Change Biol.
        Modified from: Leuning et al. 2008. A Simple surface conductance model
        to estimate regional evaporation using MODIS leaf area index and the
        Penman-Montheith equation. Water. Resources. Res., 44, W10419
        Original idea Kelliher et al. (1995). Maximum conductances for
        evaporation from global vegetation types. Agric. For. Met 85, 135-147

        Samuli Launiainen, Luke 6-11/2015. Converted to Python 21.9.2016
        """

        rho = 1.25  # kgm-3; air density
        cp = 1000.4  # Jkg-1K-1
        gamma = 66.0  # Pa K^-1, psycrometric const. at 20degC
        delta = 145.0  # Pa K-1, slope of sat. vapor pressure curve at 20degC
        epsi = delta / gamma

        # ---- model parameters:
        gsref = self.physpara['gsref']  # maximum (leaf-scale) gs (m/s)
        kp = self.physpara['kp']  # (-) attenuation coefficient for PAR
        q50 = self.physpara['q50']  # Wm-2, half-sat. of leaf light response
        rw = self.physpara['rw']  # rew parameter
        rwmin = self.physpara['rwmin']  # rew parameter

        tau = np.exp(-kp * self.LAI)

        # soil moisture modifier [0..1] """
        # ... use one based on Hyytiälä pine EC-data, 2006 & 2009 drought
        fRew = np.minimum(1.0, np.maximum(Rew / rw, rwmin))

        """--- compute conductances ----- """
        # isothermal conductance
        Gi = AE / (rho * cp / gamma * 1e3 * D)  # m/s
        # canopy conductance (integrated stomatal)
        fQ = np.log((Qp + q50) / (Qp * np.exp(-kp * self.LAI) + q50) + eps)
        fD = 1.0 / np.sqrt(D)
        # fD=1.0

        Gc = gsref / kp * fQ * fD * fRew * fPheno

        NN = 1.0 + (tau * Ga / ((epsi + 1.0) * Gc)) \
            * (f - (epsi + 1.0) * (1.0 - f) * Gc / Ga) + Ga / (epsi * Gi)

        DN = 1.0 - tau * (f - (epsi + 1.0) * (1.0 - f) * Gc / Ga) \
            + Ga / (epsi * Gi)

        Gs = Gc * NN / DN

        if not np.isfinite(Gs):
            Gs = 1e-6

        # transpiration rate
        Transpi = penman_monteith(AE, 1e3*D, T, Gs, Ga, units='mm')

        return Transpi, Gs, Gc

    def interception(self, dt, T, Prec, erate):
        """
        Calculates canopy rainfall water interception and canopy water storage changes
        Args: 
            self - object
            dt - timestep [s]
            T - air temperature (degC)
            Prec - precipitation rate (mm s-1 = kg m-2 s-1)
            erate - evaporative demand (mm s-1 = kg m-2 s-1)
        Returns:
            self - updated state W
            Trfall - potential infiltration rate to soil profile (mm s-1)
            Interc - interception rate (mm s-1)
            Evap - evaporation from canopy store (mm s-1)
            MBE - mass balance error (mm)
        """
        
        Wmax = self.wmax * self.LAI  # canopy interception capacity mm

        Prec = Prec * dt  # mm

        """ --- initial condition for MBE -"""
        Wo = self.W  # canopy storage

        # interception (mm)
        Interc = (Wmax - self.W) * (1.0 - np.exp(-(1.0 / Wmax) * Prec))
        # new canopy storage, mm        
        self.W = self.W + Interc
        # Throughfall to field layer mm
        Trfall = Prec - Interc

        # evaporate from canopy storage
        if Prec > 0:
            erate = 0.0  # negelect evaporation during precipitation events
        Evap = np.minimum(erate*dt, self.W)  # mm
        self.W = self.W - Evap  # new storage after evaporation

        # canopy storage mass-balance error
        MBE = (self.W  - Wo) - (Prec - Evap - Trfall)

        return Trfall / dt, Interc / dt, Evap / dt,  MBE

    #def run_microclimate(self, Zen, Ib1, Id1, Ib2, Id2, Uo, Ta, D, LWdn0, LWup0, Prec, Rew=1.0, soil_Par_albedo=0.05,
    #                 soil_Nir_albedo=0.35, out_layers=None, plot_figs=False):
    def run_microclimate(self, forcing, out_layers=None, plot_figs=False):
        """
        creates semi-realistic SW, LW and U profiles within canopy and solves throughfall and
        transpiration rate using big-leaf approximation

        Args:
            Zen - solar zenith angle (rad)
            Ib1 - direct Par (Wm-2)
            Id1 - diffuse Par (Wm-2)
            Ib2 - direct Nir (Wm-2)
            Id2 - diffuse Nir (Wm-2)
            Uo - wind speed at canopy top (m/s)
            Ta - air temperature, assumed uniform within canopy (degC)
            D - vapor pressure deficit (Pa)
            LWdn0 - downwelling long-wave radiation at canopy top (Wm-2)
            LWup0 - upwelling (emitted) long-wave radiation at forest floor (Wm-2)
            Rew - relative extractable water, 'bulk root zone' rew = (Theta - Theta_wp) / (Theta_fc - Theta_wp)

            soil_Par_albedo - soil Par albedo (optional, if omitted object property used)
            soil_Nir_albedo - soil Ni albedo (optional, if omitted object property used)
            out_layers [list] - indexes of output layers [0 = forest floor]
            plot_figs - plots profiles'
        Returns:
            U - wind speed (m/s)
            Par, Nir (Wm-2)
            LWdn - downwelling LW (Wm-2)
            trfall - throughfall rate (mm/s)
            evap - evaporation rate (from canopy storage) (mm/s)
            transpi - transpiration rate (mm/s)
            
        """

        """ interception, throughfall and plant transpiration using big-leaf models"""
        Zen = forcing['zenith']
        Ib1 = forcing['direct_par']
        Id1 = forcing['diffuse_par']
        Ib2 = forcing['direct_nir']
        Id2 = forcing['diffuse_nir']
        Uo = forcing['wind_speed']
        Ta = forcing['air_temperature']
        D = forcing['vapor_deficit']
        LWdn0 = forcing['lwdn']
        LWup0 = forcing['lwup']
        Prec = forcing['precipitation']
        Rew = forcing['extractable_water']
        soil_Par_albedo = forcing['ground_par_albedo']
        soil_Nir_albedo = forcing['ground_nir_albedo']

        # --- interception and throughfall rate
        AE = 0.7 * (Ib1 + Id1 + Ib2 + Id2) + eps  # coarse approximation of available energy, Wm-2
        Ga = bulk_aerodynamic_conductance(Uo, self.zm, self.d, self.zom, zos=self.zov)
        Gs_wet = 1e6  # infinite surface conductance when canopy is wet
        erate = penman_monteith(AE, D, Ta, Gs_wet, Ga, units='mm')  # mm/s

        trfall, interc, evap, mbe = self.interception(self.dt, Ta, Prec, erate)  # rates mm/s, mbe in mm

        # --- canopy transpiration rate (mm/s)
        transpi, gs, gc = self.transpiration_rate(Ib1 + Id1, AE, 1e-3*D, Ta, Ga=Ga, f=0.0, Rew=Rew)

        """ microclimate: solve for each canopy layer --- """

        # --- mean horizontal flow velocity
        U = wind_profile(self.z, self.LAI, Uo, self.zm, self.hc, self.d, self.zom, beta=self.beta)

        # --- compute Par and Nir profiles and absorption
        # we need here only SWb & SWd, the incident direct and diffuse radiation (Wm-2(ground))   
        SWb1, SWd1, SWu1, Q_sl1, Q_sh1, q_sl1, q_sh1, q_soil1, f_sl, alb1 = sw_radi(self.LAIz, self.Clump, self.leaf_angle,\
                                                                Zen, Ib1, Id1, self.Par_albedo, soil_Par_albedo, PlotFigs=False)
        SWb2, SWd2, SWu2, Q_sl2, Q_sh2, q_sl2, q_sh2, q_soil2, _, alb2 = sw_radi(self.LAIz, self.Clump, self.leaf_angle,\
                                                                Zen, Ib1, Id1, self.Nir_albedo, soil_Nir_albedo, PlotFigs=False)    

        Par = SWb1 + SWd1
        Nir = SWb2 + SWd2

        # --- compute LW radiation
        T = self.ones * Ta
        _, LWdn, LWup = lw_radi(self.LAIz, self.Clump, self.leaf_angle, T, LWdn0, LWup0, leaf_emi=0.98, soil_emi=0.98, PlotFigs=False)

        if plot_figs:
            plt.figure()
            plt.subplot(221)
            plt.plot(U, self.z/self.hc, 'r-')
            plt.title('U'); plt.ylabel('z/hc'); plt.xlabel('U m/s')

            plt.subplot(222)
            plt.plot(Par, self.z/self.hc, 'g-', Nir, self.z/self.hc, 'r-')
            plt.title('Par(g), Nir(r)'); plt.ylabel('z/hc'); plt.xlabel('Par Wm-2')
            
            plt.subplot(223)
            plt.plot(LWup, self.z/self.hc, 'b-', LWdn, self.z/self.hc, 'c-')
            plt.title('LWup(k), LWdn(c)'); plt.ylabel('z/hc'); plt.xlabel('LW_i Wm-2')

        if out_layers:
            # return values only at output layers
            results = {
                'wind_speed': U[out_layers],
                'par': Par[out_layers],
                'nir': Nir[out_layers],
                'lwdn': LWdn[out_layers],
                'throughfall': trfall[out_layers],
                'interception': interc[out_layers],
                'evaporation': evap[out_layers],
                'transpiration': transpi[out_layers]
                }
        else:
            results = {
                'wind_speed': U,
                'par': Par,
                'nir': Nir,
                'lwdn': LWdn,
                'throughfall': trfall,
                'interception': interc,
                'evaporation': evap,
                'transpiration': transpi
                }

        return results
        # return U, Par, Nir, LWdn, trfall, interc, evap, transpi


def test_model(stime, etime, out_layers=None):
    # out_layers=[0,-1]
    """
    test microclimate simulator
    Args:
        stime, etime = 'yyyy-mm-dd-', period start and end
        out_layers = indices of layers to be returned. 0 = ground, -1 =top (forcing). 
            set =None to get all profiles
    """
    os.chdir('c:\\Projects\\Mxl_Apes\\')
    ffile = os.path.join('c:\\projects\\Mxl_Apes\\Data\\', 'Forcing_2005.csv')

    # --- initializing CanopyModel
    dt = 1800.0
    LAI = 4.0
    loc = {'Lat': 61.51, 'Lon': 24.0, 'Elev': 181.0}
    grid = {'zmax': 15.0, 'Nlayers': 50, 'hc': 15.0}
    
    radi_para = {'Clump': 0.8, 'leaf_angle': 1.0, 'Par_alb': 0.12, 'Nir_alb': 0.55, 
                 'soil_Par_alb': 0.05, 'soil_Nir_alb': 0.35}
    
    # big-leaf model for surface conductance
    phys_para = {'q50': 40.0, 'gsref': 1.9e-3, 'kp': 0.6, 'f': 0.0, 'rw': 0.11, 'rwmin': 0.25} 

    beta = 2.0  # flow parameter
    wmax = 0.20  # mm/LAI    
 
    # forcing   
    Forc, tvec = read_forcing(ffile, loc=loc)
    Forc = Forc[(Forc.index >= stime) & (Forc.index <= etime)]
    Nsteps = len(Forc)
    # -- create model
    cpy = model(dt, loc, LAI, grid, radi_para, phys_para, beta, wmax)
    
    
    # --- create output    
    results = {'U': [], 'Par': [], 'Nir': [], 'LWdn': [], 'Trfall': [], 'Evap': [], 'Transpi': []}
    
    for k in range(0, Nsteps):
        print(' k= ' + str(k))
        Zen = Forc['Zen'].iloc[k]
        Ib1 = Forc['dirPar'].iloc[k]
        Id1 = Forc['diffPar'].iloc[k]
        Ib2 = Forc['dirNir'].iloc[k]
        Id2 = Forc['diffNir'].iloc[k]
        Uo = Forc['U'].iloc[k]
        T = Forc['Ta'].iloc[k]
        D = Forc['VPD'].iloc[k]
        LWdno = Forc['LWin'].iloc[k]
        Prec = Forc['Prec'].iloc[k]
        
        # just to get value here, should be the emitted moss LW from previous dt        
        LWupo = 0.9*LWdno  
        # so should Rew be 'effective root zone relatively extractable water' from previous timestep
        u, par, nir, lwdn, trfall, interc, evap, transpi = cpy.run_microclimate(Zen, Ib1, Id1, Ib2, Id2, Uo,
                                                                        T, D, LWdno, LWupo, Prec, Rew=1.0,
                                                                        soil_Par_albedo=None, 
                                                                        soil_Nir_albedo=None,
                                                                        out_layers=out_layers, plot_figs=False)
        results['U'].append(u)
        results['Par'].append(par)
        results['Nir'].append(nir)
        results['LWdn'].append(lwdn)
        results['Trfall'].append(trfall)
        results['Evap'].append(evap)
        results['Transpi'].append(transpi)

    return results, Forc, cpy
        
def read_forcing(ffile, loc=None):
    """
    reads 30 min forcing data and returns dataframe
    """
    cols = ['year','month', 'day', 'hour', 'minute', 'U', 'ust', 'Ta', 'RH', 'CO2', 'H2O', 'O3',
            'Prec', 'P', 'dirPar', 'diffPar', 'dirNir', 'diffNir', 'Rnet', 'LWin', 'LWout',
            'LWnet', 'Tsh', 'Tsa','Tsc', 'Wh', 'Wa', 'Wc', 'emiatm','cloudfract']

    dat = pd.read_csv(ffile, sep=';', header=None, names=cols)#, parse_dates=[[0,1,2,3,4]])
    tvec = pd.to_datetime(dat[['year', 'month', 'day', 'hour', 'minute']])
    tvec = pd.DatetimeIndex(tvec)
    doy = tvec.dayofyear
    dat['doy'] = doy
    
    # compute solar zenith angle (deg)
    if loc:
        M = len(dat)
        zen = np.zeros(M)
        for k in range(M):
            zen[k], _, _, _, _, _ = solar_dec_azimuth(loc['Lat'], loc['Lon'], dat['doy'][k], dat['hour'][k], dat['minute'][k])

        dat['Zen'] = zen
    # for pre-2009(?) years at hyytiälä, no LWin and LWout so approximate LWin from temperature and atm. emissivity.
    # ea seems to have gaps so lets interpolate them first
    T = dat['Ta']
    ea = dat['emiatm']
    ix = np.where(np.isfinite(ea))[0]
    f = interp1d(ix, ea[ix], kind='nearest', fill_value='extrapolate')
    dat['emiatm'] = f(range(0,len(T)))
    dat['LWin'] = dat['emiatm'] * SIGMA * (T + NT)**4.0
    
    dat['Prec'] = dat['Prec'] / 1800.0  # mm/s
    rh = dat['RH']
    rh[rh > 100] = 100.0
    dat['VPD'] = vpd_from_rh(T, rh)  # Pa
    
    dat.index = tvec
    return dat, tvec
    