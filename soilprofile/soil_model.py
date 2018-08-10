# -*- coding: utf-8 -*-
"""
TEST SCRIPT FOR CREATING SIMPLE SOIL PROFILE OBJECT AND RUNNING SOIL HEAT AND WATER BUDGET MODELS
Uses soil_water and soil_heat modules.

@author: slauniai

LAST EDITS:
    AJ 12.10.2017
    Samuli 25.8.2017.
"""

import numpy as np
import pandas as pd
import os
from datetime import datetime
from matplotlib import pyplot as plt

import soil_heat as hf
import soil_water as wf

eps = np.finfo(float).eps  # machine epsilon

# function calls
# def waterFlow1D(t_final, z, h0, pF, Ksat, Ftop, R, HM=0.0, lbcType='impermeable', lbcValue=None,
#                Wice0=0.0, maxPond=0.0, pond0=0.0, cosalfa=1.0, h_atm=-1000.0, steps=10):
#
# def heatflow_1D(t_final, z, pF, T0, Wliq0, Wice0, ubc, lbc, spara, S=0.0, steps=10):


class SoilModel():
    import soil_water as wf
    import soil_heat as hf
    
    def __init__(self, z, para):
        """
        Initializes SoilModel for 1D solutions of soil water and heat budgets
        Args:
            z (array): vertical grid [m], <0
            para (dict): floats or arrays of len(z)
                * pF (dict): vanGenuchten water retention parameters; scalars or arrays of len(z)
                * Ksat (float/array): saturated vertical hydraulic conductivity [ms-1]
                * Khsat (float/array): saturated horizontal hydraulic conductivity [ms-1]
                * Csv (float/array): dry soil vol. heat capacity [J m-3 (total volume) K-1]
                * vOrg (float(array): organic matter fraction of solid volume [-]
                * vSand (float(array): sand fraction of solid volume [-]
                * vSilt(float(array): silt fraction of solid volume [-]
                * vClay (float(array): clay fraction of solid volume [-]
                * fp (float/array): freezing curve parameter
                * max_pond: [m] maximum pond depth
                * ini_cond (dict): inputs are floats or arrays of len(z)
                    * 'Wtot', vol. water content [-] or
                    * 'h', materix water potential [m] or
                    * 'gwl', ground water level, <0 [m] (assumes steady state in column)
                    * 'T', soil temperature [degC]
                    * 'pond', [m] initial pond depth at surface
                * lbc_heat: lower boundary condition type for heat
                * lbc_water: lower boundary condition type for water (type, value, depth)
                * homogenous (boolean): True assumes vertically homogenous profile and float inputs
                * solve_heat (boolean): True solves heatflow
                * solve_water (boolean): True solves waterflow
            
        Returns:
            object
        """        
        self.para = para

        # --- boundary conditions and flags
        self.lbc_heat = para['lbc_heat']
        self.lbc_water = para['lbc_water']
        
        self.solve_heat = para['solve_heat']
        self.solve_water = para['solve_water']

        N = len(z)
        self.Nlayers = N
        
        # grid
        dz = np.zeros(N)
        dzu = np.zeros(N)
        dzl = np.zeros(N)

        # distances between grid points i-1 and i
        dzu[1:] = z[:-1] - z[1:]
        dzu[0] = -z[0]  # from soil surface to first node, soil surface af z = 0
    
        # compartment thickness (nodes in cell center!! Would be easier to input thicknessess not z)
        dz[0] = 2 * dzu[0]
        for k in range(1, N):
            dz[k] = 2 * dzu[k] - dz[k-1]
    
        # distances between grid points i and i+1
        dzl[:-1] = z[:-1] - z[1:]
        dzl[-1] = dz[-1] / 2.0 #  from last node to bottom surface

        self.grid = {'z': z,
                     'dz': dz,
                     'dzu': dzu,
                     'dzl': dzl}

        # create homogenous soil profile: needs float inputs
        if para['homogenous']:
            # pF -parameters
            self.pF = {key: np.ones(N)*para['pF'][key] for key in para['pF'].keys()}
            # porosity [-]
            self.porosity = np.ones(N)*self.pF['ThetaS']
            # sat. vertical hydr. cond [ms-1]
            self.Ksat = np.ones(N)*para['Ksat']
            # sat. horizontal hydr. cond [ms-1]
            self.Khsat = np.ones(N)*para['Khsat']
            # freezing curve parameter [-]
            fp = np.ones(N)*para['fp']
            # organic matter vol. fraction of solids [-]
            vOrg = np.ones(N)*para['vOrg']
            # sand vol. fraction of solids [-]
            vSand = np.ones(N)*para['vSand']
            # silt vol. fraction of solids [-]
            vSilt = np.ones(N)*para['vSilt']
            # silt vol. fraction of solids [-]
            vClay = np.ones(N)*para['vClay']
            # dry soil vol. heat capacity [J m-3 K-1]           
            Csv = np.ones(N)*para['Csv']

        # create non-homogenous soil profile from inputs
        # this uses information on soil horizons given in 'para', and assigns values for each soil
        # layer above horizon bottom given in para['zh']
        else:
            Nhor = len(para['zh'])
            keys = ['Ksat', 'Khsat', 'Csv', 'fp', 'vOrg', 'vSand', 'vSilt', 'vClay']
            lpara = {key: np.ones(N)*np.NaN for key in keys}
            lpara['pF'] = {key: np.ones(N)*np.NaN for key in para['pF'].keys()}

            ix0 = 0
            for k in range(0, Nhor):
                ix = np.where(z >= para['zh'][k])[0][-1]
                for key in lpara.keys():
                    if key is not 'pF':
                        lpara[key][ix0:ix+1] = para[key][k]
                else:
                    for pp in lpara['pF'].keys():
                        lpara['pF'][pp][ix0:ix+1] = para['pF'][pp][k]
                ix0 = ix + 1

            self.pF = lpara['pF']
            self.porosity = self.pF['ThetaS']

            self.Ksat = lpara['Ksat']
            self.Khsat = lpara['Khsat']

            # freezing curve parameter [-]
            fp = lpara['fp']
            # organic matter vol. fraction of solids [-]
            vOrg = lpara['vOrg']
            # sand vol. fraction of solids [-]
            vSand = lpara['vSand']
            # silt vol. fraction of solids [-]
            vSilt = lpara['vSilt']
            # silt vol. fraction of solids [-]
            vClay = lpara['vClay']
            # dry soil vol. heat capacity [J m-3 K-1]
            Csv = lpara['Csv']

        bedrock = para['Bedrock']
        max_pond = para['max_pond']
        ini_cond = para['ini_cond']
        self.max_pond = max_pond

        # initialize state variables
        T = np.ones(N)*ini_cond['T']
        self.pond = ini_cond['pond']  # pond storage [m]

        if 'h' in ini_cond:
            h = np.ones(N)*ini_cond['h']
            # vol. water content [-]
            Wtot = wf.h_to_cellmoist(self.pF, h, self.grid['dz'])
        elif 'gwl' in ini_cond:
            h = np.ones(N)*(ini_cond['gwl']-z)
            Wtot = wf.h_to_cellmoist(self.pF, h, self.grid['dz'])
        else:
            Wtot = np.ones(N)*ini_cond['Wtot']

        Wliq, Wice, _ = hf.frozen_water(T, Wtot, fp=fp)
        # temperature [degC]
        self.T = T
        # hydraulic head [m]
        self.h = h
        # vol. water content [-]
        self.Wtot = Wtot
        # liq. water content [-]
        self.Wliq = Wliq
        # ice content [-]
        self.Wice = Wice
        # air volume [-]
        self.Wair = self.porosity - self.Wtot

        # ground water level [m]
        gwl = wf.get_gwl(self.h, self.grid['z'])
        self.gwl = gwl
        # print(self.z, self.Wliq, self.T, self.h, self.Wice, self.Ls)

        # vertical & horizontal hydraulic conductivities [ms-1]
        self.Kv = wf.hydraulic_conductivity(self.pF, x=self.h, Ksat=self.Ksat)
        self.Kh = wf.hydraulic_conductivity(self.pF, x=self.h, Ksat=self.Khsat)

        # thermal conductivity [Wm-1K-1]
        self.Lambda = hf.thermal_conductivity(
            self.porosity, self.Wliq, wice=self.Wice,
            vOrg=vOrg, vSand=vSand, vSilt=vSilt, vClay=vClay)

        # In case there is impermeable layer (bedrock bottom), water flow is solved only above this
        self.ix_w = np.where(self.grid['z'] >= self.lbc_water['depth'])
        ix_max = np.max(self.ix_w)
        self.grid_w = {key: self.grid[key][self.ix_w] for key in self.grid.keys()}

        if ix_max < self.Nlayers-1:
            self.h[ix_max + 1:] = np.NaN
            self.Wliq[ix_max + 1:] = 0.0 + eps
            self.Wice[ix_max + 1:] = 0.0 + eps
            self.Wair[ix_max + 1:] = 0.0 + eps
            self.porosity[ix_max + 1:] = 0.0 + eps
            self.Kv[ix_max + 1:] = 0.0 + eps
            self.Kh[ix_max + 1:] = 0.0 + eps
            self.Lambda[ix_max + 1:] = bedrock['Lambda']
            Csv[ix_max + 1:] = bedrock['Cv']
            fp[ix_max + 1:] = 1.0
            vOrg[ix_max + 1:] = 0.0 + eps
            vSand[ix_max + 1:] = 0.0 + eps
            vSilt[ix_max + 1:] = 0.0 + eps
            vClay[ix_max + 1:] = 0.0 + eps

        # input dict for heat flow solver
        self.solids_prop = {'cs': Csv, 'fp': fp,
                            'vOrg': vOrg, 'vSand': vSand, 'vSilt': vSilt, 'vClay': vClay,
                            'bedrockL': bedrock['Lambda']
                            } 

        if self.solve_water:
            # if using equilibrium approach save WstoToGwl
            if 'solve_water_type' in para:
                self.solve_water_type = para['solve_water_type']
                if self.solve_water_type == 'Equilibrium':
                    self.WstoToGwl, self.GwlToWsto = wf.gwl_Wsto(self.grid['dz'], self.pF)
            else:
                self.solve_water_type = None

            # drainage equation
            if 'drainage_equation' in para:
                self.drainage_equation = para['drainage_equation']
            else:
                self.drainage_equation = {'type': None}

            # Keep track of dt used in solving water and heat
            self.dt_water = None
            self.dt_heat = None


    def _run(self, dt, ubc_water, ubc_heat, h_atm=-1000.0, water_sink=None, heat_sink=None, 
                      lbc_water=None, lbc_heat=None):
        """
        Runs soil model for one timestep and updates soil state variables
        
        Returns (dict):
            fluxes
                * 'infiltration' [m]
                * 'evaporation' [m]
                * 'drainage' [m]
                * 'runoff' [m]
                * 'vertical_water_flux'
                * 'vertical_heat_flux'
                * 'water_closure' [m]
            state
                * 'temperature' [degC]
                * 'water_potential' [m]
                * 'volumetric_water_content' [m3m-3]
                * 'ice_content' [m3m-3]
                * 'pond_storage' [m]
                * 'ground_water_level' [m]
                * 'hydraulic_conductivity' [ms-1]
                * 'thermal_conductivity [Wm-1s-1]
        """

        ix = self.ix_w
        pFpara = {'ThetaS': self.pF['ThetaS'][ix], 'ThetaR': self.pF['ThetaR'][ix],
                  'alpha': self.pF['alpha'][ix], 'n': self.pF['n'][ix]}
       
        fluxes = {}

        if water_sink is None:
            water_sink = np.zeros(self.Nlayers)
        if heat_sink is None:
            heat_sink = np.zeros(self.Nlayers)

        # this allows lower bc's be changed during simulation
        if lbc_water:
            self.lbc_water = lbc_water
        if lbc_heat:
            self.lbc_heat = lbc_heat

        if self.solve_water:
            # drainage to ditches
            if self.drainage_equation['type'] == 'Hooghoudt':  # Hooghoudt's drainage
                _, q_drain = wf.drainage_hooghoud(self.grid_w['dz'],
                                              self.Khsat[ix],
                                              self.gwl,
                                              self.drainage_equation['depth'],
                                              self.drainage_equation['spacing'],
                                              self.drainage_equation['width'])
            else:
                q_drain = np.zeros(len(self.grid_w['z']))  # no drainage
            if self.dt_water == None:
                self.dt_water = dt
            if self.dt_heat == None:
                self.dt_heat = dt
            # solve water balance in domain of ix layers
            if self.solve_water_type == 'Equilibrium':  # solving based on equilibrium
                h, Wtot, self.pond, infil, evapo, drainage, trans, roff, fliq, self.gwl, Kv, mbe = \
                    wf.waterStorage1D(t_final=dt,
                                      grid=self.grid_w,
                                      h0=self.h[ix],
                                      pF=pFpara,
                                      Ksat=self.Ksat[ix],
                                      Prec=ubc_water['Prec'],  # flux-based precipitation >0 [m s-1]
                                      Evap=ubc_water['Evap'],  # flux-based evaporation >0 [m s-1]
                                      R=water_sink[ix],  # total water sink (root etc.)
                                      WstoToGwl=self.WstoToGwl,
                                      GwlToWsto=self.GwlToWsto,
                                      HM=q_drain,  # lateral flow 
                                      lbc=self.lbc_water,  # lower boundary condition {'type': xx, 'value': yy}
                                      Wice0=self.Wice[ix],
                                      maxPond=self.max_pond,
                                      pond0=self.pond,
                                      cosalfa=1.0,
                                      h_atm=h_atm)  # atmospheric pressure head [m]; for soil evap. controls
            else:  # Richards equation for solving water flow
                h, Wtot, self.pond, infil, evapo, drainage, trans, roff, fliq, self.gwl, Kv, mbe, self.dt_water = \
                    wf.waterFlow1D(t_final=dt,
                                   grid=self.grid_w,
                                   h0=self.h[ix],
                                   pF=pFpara,
                                   Ksat=self.Ksat[ix],
                                   Prec=ubc_water['Prec'],  # flux-based precipitation >0 [m s-1]
                                   Evap=ubc_water['Evap'],  # flux-based evaporation >0 [m s-1]
                                   R=water_sink[ix],  # total water sink (root etc.)
                                   HM=q_drain,  # lateral flow
                                   lbc=self.lbc_water,  # lower boundary condition {'type': xx, 'value': yy}
                                   Wice0=self.Wice[ix],
                                   maxPond=self.max_pond,
                                   pond0=self.pond,
                                   cosalfa=1.0,
                                   h_atm=h_atm,  # atmospheric pressure head [m]; for soil evap. controls
                                   steps=dt / self.dt_water)
            self.h[ix] = h.copy()
            self.Wtot[ix] = Wtot.copy()
            self.Kv[ix] = Kv.copy()
            self.Wair = self.porosity - self.Wtot

            fluxes.update({'infiltration': infil / dt,
                           'evaporation': evapo / dt,
                           'transpiration': trans /dt,
                           'subsurface_drainage': drainage /dt,
                           'surface_runoff': roff /dt,
                           'total_runoff': (drainage + roff) /dt,
                           'vertical_water_flux': fliq,
                           'MBE': mbe
                          })
        else:
            self.gwl = self.para['ground_water_level']
            self._steady_state_water(pF= self.pF, Ksat=self.Ksat, gwl=self.gwl)

        if self.solve_heat:
            self.T, self.Wliq, self.Wice, fheat, self.Lambda, self.dt_heat, heat_be = \
                hf.heatflow_1D_new(t_final=dt,
                               grid=self.grid,
                               poros=self.porosity,
                               T0=self.T,
                               Wtot=self.Wtot,
                               ubc=ubc_heat,  # upper boundary condition {'type': xx, 'value': yy}
                               lbc=self.lbc_heat,
                               spara=self.solids_prop,  # dict of properties of solid part
                               S=heat_sink,
                               steps= dt / self.dt_heat)

            fluxes.update({'vertical_heat_flux': fheat,
                           'heat_be': heat_be})


        # return state in dictionary
        state = {"water_potential": self.h,
                 "volumetric_water_content": self.Wtot,
                 "column_water_storage": sum(self.Wtot*self.grid['dz']),
                 "ice_content": self.Wice,
                 "pond_storage": self.pond,
                 "ground_water_level": self.gwl,
                 "hydraulic_conductivity": self.Kv,
                 "temperature": self.T,
                 "thermal_conductivity": self.Lambda,
                 }

        return fluxes, state

    def _steady_state_water(self, pF, Ksat, gwl):
        """ Calculates soil water content in steady state with ground water level
        """

        h = np.arange(self.grid['z'][0], gwl, self.grid['z'][1]-self.grid['z'][0])
        h = np.append(np.flip(h, axis=0), np.zeros(self.Nlayers - len(h)))
        self.h = h.copy()

        Wtot = wf.wrc(pF, x=h)
        self.Kv = wf.hydraulic_conductivity(pF, h, Ksat)

        self.Wliq, self.Wice, _ = hf.frozen_water(self.T, Wtot, self.solids_prop['fp'])
        self.Wair = self.porosity - self.Wliq - self.Wice

def update_steady_state(h, T, Wtot, fp, porosity):

    Wliq, Wice = hf.frozen_water(T, Wtot, fp)
    Wair = porosity - Wliq - Wice

    return Wliq, Wice, Wair


#def create_soilmodel():
#    # define grid
#    # [m]
#    z = np.arange(-0.01, -2.0, -0.02)
#    N = len(z)
#
#    # define soil type
#    pF = {'ThetaS': 0.5, 'ThetaR': 0.05, 'alpha': 0.3, 'n': 1.3}
#    poros = pF['ThetaS']
#    vOrg = 0.01
#    vQuartz = 0.1
#    vMineral = 1.0 - poros - vQuartz - vOrg
#    # heat capacity
#    Cvs = hf.volumetric_heat_capacity(poros, orgfract=vOrg)
#    Ksat = 1e-5
#    Khsat = Ksat
#    Ls = hf.thermal_conductivity_solids(
#        poros, vMineral=vMineral, vQuartz=vQuartz, vOrg=vOrg)
#    fp = 2.0
#    ini_cond = {'T': 10.0, 'h': -0.1, 'pond': 0.0}
#
#    # model bcs
#    max_pond = 0.05  # m
#    lbc_heat = {'type': 'temperature', 'value': 10.0}
#    lbc_water = {'type': 'flux', 'value': 0.0}
#
#    model = SoilModel(z, pF, Ksat, Khsat, Cvs, Ls, fp, vOrg, max_pond, ini_cond, lbc_heat, 
#                 lbc_water, homogenous=True, solve_heat=True, solve_water=True)
#
#    return model


class SoilProfile():
    """
    One-dimensional soil profile class definition
    """
    def __init__(self, z, pF, Ksat, Khsat, Cvs, Ls, fp, vOrg, vSand, vSilt, vClay, ini_cond, homogenous=False):
        """
        Initializes soil profile.
        Args:
            z (array): vertical grid [m], <0
            pF (dict): vanGenuchten water retention parameters; scalars or arrays of len(z)
            Ksat (float/array): saturated vertical hydraulic conductivity [ms-1]
            Khsat (float/array): saturated horizontal hydraulic conductivity [ms-1]
            Cvs (float/array): soil solid material vol. heat capacity [J m-3 K-1]
            Ls (float/array): soil solid material thermal conductivity [Wm-1K-1]
            fp (float/array): freezing curve parameter
            ini_cond (dict): inputs are floats or arrays of len(z)
                * 'Wtot', vol. water content [-] or
                * 'h', materix water potential  [m]
                * 'T', soil temperature [degC]
                * 'pond', pond storage [m]
            homogenous (boolean): True assumes vertically homogenous profile and float inputs
        Returns:
            object
        """
        N = len(z)

        # grid
        dz = np.empty(N)
        dzu = np.empty(N)
        dzl = np.empty(N)

        # distances between grid points:
        # dzu is between point i-1 and i, dzl between point i and i+1
        dzu[1:] = z[0:-1] - z[1:N]
        dzu[0] = -z[0]
        dzl[0:-1] = z[0:-1] - z[1:]
        dzl[-1] = (z[-2] - z[-1]) / 2.0

        dz = (dzu + dzl) / 2.0
        dz[0] = dzu[0] + dzl[0] / 2.0

        self.z = z                          # vertical grid [m], <0
        self.dz = dz                        # layer thickness [m]
        del dz, dzu, dzl

        # create homogenous soil profile: needs float inputs
        if homogenous:
            # pF -parameters
            self.pF = {key: np.ones(N)*pF[key] for key in pF.keys()}
            # porosity [-]
            self.porosity = np.ones(N)*self.pF['ThetaS']
            # sat. vertical hydr. cond [ms-1]
            self.Ksat = np.ones(N)*Ksat
            # sat. horizontal hydr. cond [ms-1]
            self.Khsat = np.ones(N)*Khsat
            # soil solid vol. heat capacity [J m-3 K-1]
            self.Cvs = np.ones(N)*Cvs
            # soil solid thermal conductivity [Wm-1K-1]
            self.Ls = np.ones(N)*Ls
            # freezing curve parameter [-]
            self.fp = np.ones(N)*fp
            # organic matter vol. fraction of solids [-]
            self.vOrg = np.ones(N)*vOrg
            # sand vol. fraction of solids [-]
            self.vSand = np.ones(N)*vSand
            # silt vol. fraction of solids [-]
            self.vSilt = np.ones(N)*vSilt
            # silt vol. fraction of solids [-]
            self.vClay = np.ones(N)*vClay

        else:
            self.pF = pF
            self.porosity = self.pF['ThetaS']
            self.Ksat = Ksat
            self.Khsat = Khsat
            self.Cvs = Cvs
            self.Ls = Ls
            self.fp = fp
            self.vOrg = vOrg

        # initialize state variables
        T = np.ones(N)*ini_cond['T']

        if 'h' in ini_cond:
            h = np.ones(N)*ini_cond['h']
            # vol. water content [-]
            Wtot = wf.wrc(self.pF, x=h)
        else:
            Wtot = np.ones(N)*ini_cond['Wtot']

        Wliq, Wice = hf.frozen_water(T, Wtot, fp=self.fp)
        # temperature [degC]
        self.T = T
        # hydraulic head [m]
        self.h = h
        # liq. water content [-]
        self.Wliq = Wliq
        # ice content [-]
        self.Wice = Wice
        # air volume [-]
        self.Wair = self.porosity - self.Wliq - self.Wice
        # ground water level [m]
        self.Gwl = wf.get_gwl(self.h, self.z)

        print(self.z, self.Wliq, self.T, self.h, self.Wice, self.Ls)
        # thermal conductivity [Wm-1K-1]
        self.Lambda = hf.thermal_conductivity_deVries(
            self.porosity, self.Wliq, wice=self.Wice,
            h=self.h, pF=self.pF, T=self.T, ks=self.Ls, vOrg=self.vOrg)
        # hydraulic conductivity [ms-1]
        K = wf.hydraulic_conductivity(self.pF, x=self.Wliq, var='Th')
        K = wf.spatial_average(K, self.z, method='arithmetic')
        self.Khydr = K

        self.pond = ini_cond['pond']        # pond storage [m]
        del K

def create_profile():
    """
    tests creating homogenous soil profile
    """
    # define grid
    # [m]
    z = np.arange(-0.01, -2.0, -0.02)
    N = len(z)

    # define soil type
    pF = {'ThetaS': 0.5, 'ThetaR': 0.05, 'alpha': 0.3, 'n': 1.3}
    poros = pF['ThetaS']
    vOrg = 0.01
    vQuartz = 0.1
    vMineral = 1.0 - vQuartz - vOrg
    # heat capacity
    Cvs = hf.volumetric_heat_capacity(poros, orgfract=vOrg)
    Ksat = 1e-5
    Khsat = Ksat
    Ls = hf.thermal_conductivity_solids(
        poros, vMineral=vMineral, vQuartz=vQuartz, vOrg=vOrg)
    fp = 2.0
    ini_cond = {'T': 10.0, 'h': -0.1, 'pond': 0.0}

    soil0 = SoilProfile(
        z, pF, Ksat, Khsat, Cvs, Ls, fp, vOrg, ini_cond, homogenous=True)

    return soil0

def read_setup(init_file):
    import json
    climoss_path = os.path.join('/projects/Climoss/EnergyBudget/climoss')
    with open(os.path.join(climoss_path, init_file), 'r') as init_file:
        cfg = json.load(init_file)

    general = cfg['general']
    bryo = cfg['bryotype']
    soil = cfg['soil']

    return general, bryo, soil
    
    # dz = np.empty(np.shape(z))
    # dzu = dz.copy(); dzl=dz.copy()
    # dzu[1:] = z[0:-1] - z[1:N]; dzu[0]=-z[0]
    # dzl[0:-1] = z[0:-1] - z[1:]; dzl[-1]=(z[-2] - z[-1])/2.0;

    # dz = (dzu + dzl)/2.0;
    # dz[0] = dzu[0] + dzl[0]/2.0;
    #


# model params

# #%%

# t=np.arange(-10.0,3.0, 0.1)
# w=0.5
# wl,wi=hf.frozen_water(t,w,fp=4.0, To=0.0)

# plt.figure()
# plt.plot(t,wl,'r.-',t,wi,'c.-'); plt.ylabel('water & ice'); plt.xlabel('T')
# #%%

# #w = np.array([0.1,0.2,0.3,0.4,0.5])
# w = np.arange(0.1,0.5,0.05)
# hh = wf.wrc(pF, x=w, var='Th')
# ksolid = hf.thermal_conductivity_solids(poros=0.5, vMineral=0.0, vQuartz=0.0, vOrg=0.5)
#
# ks1=hf.thermal_conductivity_Campbell(
#     pF['ThetaS'], w, wice=0.0, vQuartz=0.05, vClay=0.9, vMineral=0.05)
# ks2=hf.thermal_conductivity_deVries(
#    pF['ThetaS'], w, wice=0.0, h=hh, pF=pF, T=15.0, ks=None, vOrg=0.0, vQuartz=0.05)
# print 'Ks1=' + str(ks1)
# print 'Ks2=' + str(ks2)
#
# #%% heatflow_1D
#
# spara = {'frp': 2.0, 'cs': 2.3, 'ktherm': 2.6}
# T0 = np.ones(np.shape(z))
# W0 = np.zeros(np.shape(z)) + 0.3
# Wi0 = np.zeros(np.shape(z))
#
# # sink term
# S0 = np.zeros(np.shape(z))
# S0[10:20] = 100.0
#
# # heatFlow1D(t_final,z,pF,T0,Wliq0,Wice0, ubc, lbc,spara, S=0.0, steps=10):
# dt0 = 24*60*60 # 1 day
# # dt = 1800
# # Ftop = +1.0e-3/dt0#+0.0e-3#-2e-3/dt0#-2.0e-3/dt0
#
# ubc = {'type': 'flux', 'value': 100.0}
# lbc = {'type': 'temperature', 'value': 1.0}
#
# T_new, W_new, Wi_new, Fsoil, Lambda = hf.heatflow_1D(
#     dt0, z, pF, T0, W0, Wi0, ubc, lbc, spara, S=S0)
#
# plt.close('all')
#
# plt.subplot(121)
# plt.plot(T0,z,'k.-', T_new,z,'r.-'); plt.ylabel('z'); plt.xlabel('T')
# plt.subplot(122)
# plt.plot(W_new,z,'r.-', Wi_new,z,'b.-'); plt.ylabel('z'); plt.xlabel('Wliq (r), Wice(b)')
#
# ##%% test simpler version
# #spara={'frp':2.0, 'Ktherm':0.5, 'cs':1.3e6, 'poros':0.5}

# #T0=np.ones(np.shape(z))*-1.0
# #W0=np.zeros(np.shape(z)) + 0.3

# #S0=np.zeros(np.shape(z))
# ##heatFlow1D(t_final,z,pF,T0,Wliq0,Wice0, ubc, lbc,spara, S=0.0, steps=10):
# #dt0=24*60*60#1 day
# ##dt=1800
# ##Ftop=+1.0e-3/dt0#+0.0e-3#-2e-3/dt0#-2.0e-3/dt0
# #
# #ubc={'type': 'flux', 'value': -1.0}
# #lbc={'type': 'temperature', 'value': 1.0}
# #
# Tnew, W_new, Wi_new, Fsoil, frd, thd = hf.heatFlow1D_Simple(
#     dt0,z,T0,W0, ubc, lbc,spara, S=0.0, steps=10)

# #plt.close('all')

# #plt.plot(T0,z,'k-', Tnew,z,'r.-')
# #%% Test Rankinen et al. soil temperature model

# param={'poros': 0.5, 'Ktherm': 0.7, 'Cs': 1.3e6, 'Cice':10.0e6, 'fs':3.0}

# zlow=-2.0
# z=np.arange(-0.1,zlow,-0.1) #m
# N=len(z)

# To=np.ones(N)*0.0
# Wtot=0.8*param['poros']
# Ds=0.2

# dt0=10*24*3600 #s
# T_sur=-3.0
# T_sur1=T_sur*np.exp(-3.0*Ds)

# Tnew=hf.soil_temperature_Rankinen(dt0, z, To, T_sur, Ds, param, steps=100)
# T_sur=+1.0
# Tnew1=hf.soil_temperature_Rankinen(dt0, z, Tnew, T_sur, Ds, param, steps=100)
