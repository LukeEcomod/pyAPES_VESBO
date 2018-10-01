# -*- coding: utf-8 -*-
"""
Functions and classes for PlantType

@author: L1656
"""
import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from photo import leaf_interface
from phenology import Photo_cycle, LAI_cycle
from rootzone import RootUptake

class PlantType():
    """
    PlantType -class.
    Contains plant-specific properties, state variables and phenology functions
    """

    def __init__(self, z, p, dz_soil, ctr):
        """
        Creates PlantType
        Args:
            z - grid, evenly spaced, np.array
            p - parameters (dict)
            Switch_x - controls
        Returns:
            PlantType instance
        """

        self.Switch_pheno = ctr['pheno_cycle']  # include phenology
        self.Switch_lai = ctr['seasonal_LAI']  # seasonal LAI
        self.Switch_WaterStress = ctr['WaterStress']  # water stress affects stomata

        self.name = p['name']

        # phenology model
        if self.Switch_pheno:
            self.Pheno_Model = Photo_cycle(p['phenop'])  # phenology model instance
            self.pheno_state = self.Pheno_Model.f  # phenology state [0...1]
        else:
            self.pheno_state = 1.0

        # dynamic LAI model
        if self.Switch_lai:
            # seasonality of leaf area
            self.LAI_Model = LAI_cycle(p['laip'])  # LAI model instance
            self.relative_LAI = self.LAI_Model.f  # LAI relative to annual maximum [..1]
        else:
            self.relative_LAI = 1.0

        # physical structure
        self.LAImax = p['LAImax']  # maximum annual 1-sided LAI [m2m-2]
        self.LAI = self.LAImax * self.relative_LAI  # current LAI
        self.lad_normed = p['lad']  # normalized leaf-area density [m-1]
        self.lad = self.LAI * self.lad_normed  # current leaf-area density [m2m-3]

        # root properties
        self.Roots = RootUptake(p['rootp'], dz_soil, self.LAImax)

        self.dz = z[1] - z[0]
        # plant height [m]
        if len(np.where(self.lad_normed > 0)[0]) > 0:
            f = np.where(self.lad_normed > 0)[0][-1]
            self.hc = z[f]
        else:
            self.hc = 0.0
        # leaf gas-exchange parameters
        self.photop0 = p['photop']   # A-gs parameters at pheno_state = 1.0 (dict)
        self.photop = self.photop0.copy()  # current A-gs parameters (dict)
        # leaf properties
        self.leafp = p['leafp']  # leaf properties (dict)
        self.StomaModel = ctr['StomaModel']

        self.Tl_sh = None

    def _update_daily(self, doy, T, PsiL=0.0):
        """
        Updates PlantType pheno_state, gas-exchange parameters, LAI & lad
        Args:
            doy - day of year
            T - daily air temperature [degC]
            Psi_leaf - leaf (or soil) water potential, <0 [MPa]
        NOTE: CALL ONCE PER DAY
        """
        PsiL = np.minimum(-1e-5, PsiL)

        if self.Switch_pheno:
            self.pheno_state = self.Pheno_Model._run(T, out=True)

        if self.Switch_lai:
            self.relative_LAI =self.LAI_Model._run(doy, T, out=True)
            self.LAI = self.relative_LAI * self.LAImax
            self.lad = self.lad_normed * self.LAI

        # scale photosynthetic capacity using vertical N gradient
        f = 1.0
        if 'kn' in self.photop0:
            kn = self.photop0['kn']
            Lc = np.flipud(np.cumsum(np.flipud(self.lad*self.dz)))
            Lc = Lc / np.maximum(Lc[0], eps)
            f = np.exp(-kn*Lc)
        # preserve proportionality of Jmax and Rd to Vcmax
        self.photop['Vcmax'] = f * self.pheno_state * self.photop0['Vcmax']
        self.photop['Jmax'] =  f * self.pheno_state * self.photop0['Jmax']
        self.photop['Rd'] =  f * self.pheno_state * self.photop0['Rd']

        if self.Switch_WaterStress:
            b = self.photop0['drp']
            if 'La' in self.photop0:
                # lambda increases with decreasing Psi as in Manzoni et al., 2011 Funct. Ecol.
                self.photop['La'] = self.photop0['La'] * np.exp(-b*PsiL)
            if 'm' in self.photop0:  # medlyn g1-model, decrease with decreasing Psi  
                self.photop['m'] = self.photop0['m'] * np.maximum(0.05, np.exp(b*PsiL))

    def leaf_gasexchange(self, f_sl, H2O, CO2, T, U, P, Q_sl1, Q_sh1, SWabs_sl, SWabs_sh, LWl, df, Ebal, Tl_ave, gr):
        """
        Compute leaf gas-exchange for PlantType
        """
        
        # initial guess for leaf temperature
        if self.Tl_sh is None or Ebal is False:
            Tl_sh = T.copy()
            Tl_sl = T.copy()
        else:
            Tl_sh = self.Tl_sh.copy()
            Tl_sl = self.Tl_sh.copy()

        # --- sunlit leaves
        sl = leaf_interface(self.photop, self.leafp, H2O, CO2, T, Tl_sl, Q_sl1,
                            SWabs_sl, LWl, U, Tl_ave, gr, P=P, model=self.StomaModel,
                            Ebal=Ebal, dict_output=True)

        # --- shaded leaves
        sh = leaf_interface(self.photop, self.leafp, H2O, CO2, T, Tl_sh, Q_sh1,
                            SWabs_sh, LWl, U, Tl_ave, gr, P=P, model=self.StomaModel,
                            Ebal=Ebal, dict_output=True)

        if Ebal:
            self.Tl_sh= sh['Tl'].copy()
            self.Tl_sl = sl['Tl'].copy()

        # integrate water and C fluxes over all leaves in PlantType, store resuts
        pt_stats = self._integrate(sl, sh, f_sl)

        # --- sink/source terms
        dtsource = df * (f_sl*sl['H'] + (1 - f_sl)*sh['H'])*self.lad  # W m-3
        dqsource = df * (f_sl*sl['E'] + (1.0 - f_sl)*sh['E'])*self.lad  # mol m-3 s-1
        dcsource = - df *(f_sl*sl['An'] + (1.0 - f_sl)*sh['An'])*self.lad  #umol m-3 s-1
        dRstand = np.sum(df * (f_sl*sl['Rd'] + (1.0 - f_sl)*sh['Rd'])*self.lad*self.dz)  # add dark respiration umol m-2 s-1
        Frw = df * (f_sl*sl['Fr'] + (1 - f_sl)*sh['Fr'])*self.lad

        return pt_stats, dtsource, dqsource, dcsource, dRstand, Frw

    def _integrate(self, sl, sh, f_sl):
        """
        integrates layerwise statistics (per unit leaf area) to plant level
        Arg:
            sl, sh - dict of leaf-level outputs for sunlit and shaded leaves:         
            
            x = {'An': An, 'Rd': Rd, 'E': fe, 'H': H, 'Fr': Fr, 'Tl': Tl, 'Ci': Ci,
                 'Cs': Cs, 'gs_v': gsv, 'gs_c': gs_opt, 'gb_v': gb_v}
        Returns:
            y - plant level statistics
        """
        # plant fluxes, weight factors is sunlit and shaded LAI at each layer
        f1 = f_sl*self.lad*self.dz
        f2 = (1.0 - f_sl)*self.lad*self.dz

        keys = ['An', 'Rd', 'E']
        y = {k: np.nansum(sl[k]*f1 + sh[k]*f2) for k in keys}
        del keys

        # effective statistics; layerwise fluxes weighted by An
        g1 = f1 * sl['An'] / np.maximum(eps, np.nansum(f1 * sl['An'] + f2 * sh['An']))
        g2 = f2 * sh['An'] / np.maximum(eps, np.nansum(f1 * sl['An'] + f2 * sh['An']))
        # print sum(g1 + g2)        
        keys = ['Tl', 'Ci', 'Cs', 'gs_v', 'gb_v']

        y.update({k: np.nansum(sl[k]*g1 + sh[k]*g2) for k in keys})
        # print y
        y.update({'Tleaf': f_sl * sl['Tl'] + (1.0 - f_sl) * sh['Tl']})
        y.update({'Tleaf_sl': sl['Tl'] * self.lad})
        y.update({'Tleaf_sh': sh['Tl'] * self.lad})

        return y

