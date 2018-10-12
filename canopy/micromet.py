# -*- coding: utf-8 -*-
"""
.. module: micromet
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes flow and scalar profiles within canopy.
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Oct 02 09:04:05 2018

References:
Launiainen, S., Katul, G.G., Lauren, A. and Kolari, P., 2015. Coupling boreal
forest CO2, H2O and energy flows by a vertically structured forest canopy – 
Soil model with separate bryophyte layer. Ecological modelling, 312, pp.385-405.
"""
import numpy as np
from tools.utilities import central_diff, forward_diff, tridiag, smooth
from constants import *
eps = np.finfo(float).eps  # machine epsilon

class Micromet(object):
    r""" Computes flow and scalar profiles within canopy.
    """
    def __init__(self, z, lad, hc, p):
        r""" Initializes micromet object for computation of flow 
        and scalar profiles within canopy.

        Args:
            z (array): canopy model nodes, height from soil surface (= 0.0) [m]
            lad (array): leaf area density [m2 m-3]
            hc (float): canopy heigth [m]
            p (dict):
                'zos': forest floor roughness length [m]
                'dPdx': horizontal pressure gradient
                'Cd': drag coefficient
                'Utop': ensemble U/ustar [-]
                'Ubot': lower boundary
                'Sc' (dict): {'T','H2O','CO2'} Schmidt numbers
        Returns:
            self (object)
        """

        # parameters
        self.zos = p['zos']  # forest floor roughness length [m]

        self.dPdx = p['dPdx']  # horizontal pressure gradient
        self.Cd = p['Cd']  # drag coefficient
        self.Utop = p['Utop']  # ensemble U/ustar
        self.Ubot = p['Ubot']  # lower boundary
        self.Sc = p['Sc']  # Schmidt numbers
        self.dz = z[1] - z[0]

        # initialize state variables
        _, self.U_n, self.Km_n, _, _, _ = closure_1_model_U(
                z, self.Cd, lad, hc, self.Utop + eps, self.Ubot, dPdx=self.dPdx)

    def normalized_flow_stats(self, z, lad, hc, Utop=None):
        r""" Computes normalized mean velocity profile, shear stress and
        eddy diffusivity within and above horizontally homogenous plant
        canopies using 1st order closure.

        Args:
            z (array): canopy model nodes, height from soil surface (= 0.0) [m]
            lad (array): leaf area density [m2 m-3]
            hc (float): canopy heigth [m]
            Utop (float): U/ustar [-], if None set to self.Utop
        """
        if Utop is None:
            Utop = self.Utop

        _, self.U_n, self.Km_n, _, _, _ = closure_1_model_U(
                z, self.Cd, lad, hc, Utop + eps, self.Ubot, dPdx=self.dPdx)

    def update_state(self, ustaro):
        r""" Updates wind speed profile.
        Args:
            ustaro (float): friction velocity [m s-1]
        """

        U = self.U_n * ustaro
        U[0] = U[1]

        Km = self.Km_n * ustaro + eps
        Km[0] = Km[1]
        self.Km = Km

        return U

    def scalar_profiles(self, gam, H2O, CO2, T, P, source, lbc):
        r""" Solves scalar profiles (H2O, CO2 and T) within canopy.

        Args:
            gam (float): weight for old value in iterations [-]
            H2O (array): water vapor mixing ratio [mol mol-1]
            CO2 (array): carbon dioxide mixing ratio [ppm]
            T (array): ambient air temperature [degC]
            P: ambient pressure [Pa]
            source (dict):
                'H2O' (array): water vapor source [mol m-3 s-1]
                'CO2' (array): carbon dioxide source [umol m-3 s-1]
                'T' (array): heat source [W m-3]
            lbc (dict):
                'H2O' (float): water vapor lower boundary [mol m-2 s-1]
                'CO2' (float): carbon dioxide lower boundary [umol m-2 s-1]
                'T' (float): heat lower boundary [W m-2]

        Returns:
            H2O (array): water vapor mixing ratio [mol mol-1]
            CO2 (array): carbon dioxide mixing ratio [ppm]
            T (array): ambient air temperature [degC]
            err_h2o, err_co2, err_t (floats): maximum error for each scalar
        """
        # previous guess, not values of previous time step!
        H2O_prev = H2O.copy()
        CO2_prev = CO2.copy()
        T_prev = T.copy()

        # --- H2O ---
        H2O = closure_1_model_scalar(dz=self.dz,
                                     Ks=self.Km * self.Sc['H2O'],
                                     Source=source['H2O'],
                                     ubc=H2O[-1],
                                     lbc=lbc['H2O'],
                                     scalar='H2O',
                                     T=T[-1], P=P)
        # relative error
        err_h2o = max(abs((H2O - H2O_prev) / H2O_prev))
        # new H2O
        H2O = gam * H2O_prev + (1 - gam) * H2O

        # --- CO2 ---
        CO2 = closure_1_model_scalar(dz=self.dz,
                                     Ks=self.Km * self.Sc['CO2'],
                                     Source=source['CO2'],
                                     ubc=CO2[-1],
                                     lbc=lbc['CO2'],
                                     scalar='CO2',
                                     T=T[-1], P=P)
        # relative error
        err_co2 = max(abs((CO2 - CO2_prev) / CO2_prev))
        # new CO2
        CO2 = gam * CO2_prev + (1 - gam) * CO2

        # --- T ---
        T = closure_1_model_scalar(dz=self.dz,
                                   Ks=self.Km * self.Sc['T'],
                                   Source=source['T'],
                                   ubc=T[-1],
                                   lbc=lbc['T'],
                                   scalar='T',
                                   T=T[-1], P=P)
        # absolut error
        err_t = max(abs(T - T_prev))
        # new T
        T = gam * T_prev + (1 - gam) * T
# TESTING
#        CO2 = CO2_prev
#        err_co2 = 0.0

        return H2O, CO2, T, err_h2o, err_co2, err_t

def closure_1_model_U(z, Cd, lad, hc, Utop, Ubot, dPdx=0.0, lbc_flux=None):
    """
    Computes normalized mean velocity profile, shear stress and eddy diffusivity
    within and above horizontally homogenous plant canopies using 1st order closure.
    Accounts for horizontal pressure gradient force dPdx, assumes neutral diabatic stability.
    Soleves displacement height as centroid of drag force.
    IN:
       z - height (m), constant increments
       Cd - drag coefficient (0.1-0.3) (-)
       lad - plant area density, 1-sided (m2m-3)
       hc - canopy height (m)
       Utop - U /u* upper boundary
       Uhi - U /u* at ground (0.0 for no-slip)
       dPdx - u* -normalized horizontal pressure gradient
    OUT:
       tau - u* -normalized momentum flux 
       U - u* normalized mean wind speed (-) 
       Km - eddy diffusivity for momentum (m2s-1) 
       l_mix - mixing length (m)
       d - zero-plane displacement height (m)
       zo - roughness lenght for momentum (m)
    CODE:
        Gaby Katul, Samuli Launiainen 2008-2015. Converted to Python 17.5.2017
    """
    lad = 0.5*lad  # frontal plant-area density is half of one-sided
    dz = z[1] - z[2]
    N = len(z)
    U = np.linspace(Ubot, Utop, N)

    nn1 = max(2, np.floor(N/20))  # window for moving average smoothing

    # --- Start iterative solution
    err = 999.9
    eps1 = 0.1
    dPdx_m = 0.0

    while err > 0.01:
        Fd = Cd*lad*U**2  # drag force
        d = sum(z*Fd) / sum(Fd)  # displacement height
        l_mix = mixing_length(z, hc, d)  # m

        # --- dU/dz (m-1)
        y = central_diff(U, dz)
        # y = smooth(y, nn1)

        # --- eddy diffusivity & shear stress
        Km = l_mix**2*abs(y)
        # Km = smooth (Km, nn1)
        tau = -Km * y

        # ------ Set the elements of the Tri-diagonal Matrix
        a1 = -Km
        a2 = central_diff(-Km, dz)
        a3 = Cd*lad*U

        upd = (a1 / (dz*dz) + a2 / (2*dz))  # upper diagonal
        dia = (-a1*2 / (dz*dz) + a3)  # diagonal
        lod = (a1 / (dz*dz) - a2 / (2*dz))  # subdiagonal
        rhs = np.ones(N) * dPdx_m

        # upper BC
        upd[-1] = 0.
        dia[-1] = 1.
        lod[-1] = 0.
        rhs[-1] = Utop

        if not lbc_flux:  # --- lower BC, fixed Ubot
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = Ubot
        else:  # --- lower BC, flux-based
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = 0.  # zero-flux bc
            # how to formulate prescribed flux bc?
            # rhs[0] = lbc_flux

        # --- call tridiagonal solver
        Un = tridiag(lod, dia, upd, rhs)

        err = max(abs(Un - U))

        # --- Use successive relaxations in iterations
        U = eps1*Un + (1.0 - eps1)*U
        dPdx_m = eps1*dPdx + (1.0 - eps1)*dPdx_m

    # ---- return values
    tau = tau / tau[-1]  # normalized shear stress
    zo = (z[-1] - d)*np.exp(-0.4*U[-1])  # roughness length

    y = forward_diff(Un, dz)
    Kmr = l_mix**2 * abs(y)  # eddy diffusivity
    Km = smooth(Kmr, nn1)

    # --- for testing ----
#    plt.figure(101)
#    plt.subplot(221); plt.plot(Un, z, 'r-'); plt.title('U')
#    plt.subplot(222); plt.plot(y, z, 'b-'); plt.title('dUdz')
#    plt.subplot(223); plt.plot(l_mix, z, 'r-'); plt.title('l mix')
#    plt.subplot(224); plt.plot(Km, z, 'r-', Kmr, z, 'b-'); plt.title('Km')    

    return tau, Un, Km, l_mix, d, zo

def closure_1_model_scalar(dz, Ks, Source, ubc, lbc, scalar,
                           T=20.0, P=101300.0, lbc_dirchlet=False):
    """
    Solves stationary scalar profiles in 1-D grid
    INPUT:
        dz - grid size (m), float
        Ks - eddy diffusivity (m2s-1), array
        Source - sink/source term, CO2 (umolm-3s-1), H2O (molm-3s-1), T (Wm-3),
            array
        ubc - upper boundary condition, value of CO2 (ppm), H2O (mol/mol)
            or T (degC) at highest gridpoint
        lbc - lower boundary condition, flux or value:
            flux: CO2 (umol m-2 s-1), H2O (mol m-2 s-1), T (Wm-2).
            value: CO2 (ppm), H2O (mol/mol), T (degC)
        scalar - 'CO2', 'H2O', 'T', string
        T - air temperature (degC), optional
        P - pressure (Pa)
        lbc_dirchlet - True for Dirchlet lower boundary
    OUT:
        Ca - mixing ratio profile (ppm)
    NOTE:
        assumes constant dz and Dirchlet upper boundary condition
    Gaby Katul & Samuli Launiainen, 2009 - 2017
    Kersti: code condensed, forward difference (more stable?)
    """

    dz = float(dz)
    N = len(Ks)
    rho_a = P / (287.05*(T + 273.15))  # air density, kg/m3

#TESTING
#    Source = smooth(Source, 10)

    CF = rho_a / (29e-3)  # mol m-3, molar conc. of air

    if scalar.upper() == 'CO2':  # [umol] -> [mol]
        ubc = 1e-6 * ubc
        Source = 1e-6 * Source
        lbc = 1e-6 * lbc

    if scalar.upper() == 'T':
        CF = 1004.0 * rho_a  # cp*rho_a,  J m-3 K-1

    # --- Set elements of tridiagonal matrix ---
    a = np.zeros(N)  # sub diagonal
    b = np.zeros(N)  # diagonal
    g = np.zeros(N)  # super diag
    f = np.zeros(N)  # rhs

    # intermediate nodes
    a[1:-1] = Ks[1:-1] / (dz**2)
    b[1:-1] = -2 * Ks[1:-1] / (dz**2) - (Ks[2:] - Ks[1:-1]) / (dz**2)
    g[1:-1] = Ks[1:-1] / (dz**2) + (Ks[2:] - Ks[1:-1]) / (dz**2)
    f[1:-1] = -Source[1:-1] / CF

    # uppermost node, Dirchlet boundary
    a[-1] = 0.0
    b[-1] = 1.0
    g[-1] = 0.0
    f[-1] = ubc

    # lowermost node
    if not lbc_dirchlet:  # flux-based
        a[0] = 0.0
        b[0] = 1.
        g[0] = -1.
        f[0] = (lbc / CF)*dz / (Ks[0] + eps)

    else:  #  fixed concentration/temperature
        a[0] = 0.0
        b[0] = 1.
        g[0] = 0.0
        f[0] = lbc

    x = tridiag(a, b, g, f)

    if scalar.upper() == 'CO2':  # [mol] -> [umol]
        x = 1e6*x

    return x

def closure_1_model_scalar2(dz, Ks, Source, ubc, lbc, scalar,
                           T=20.0, P=101300.0, lbc_dirchlet=False):
    """
    Solves stationary scalar profiles in 1-D grid
    INPUT:
        dz - grid size (m), float
        Ks - eddy diffusivity (m2s-1), array
        Source - sink/source term, CO2 (umolm-3s-1), H2O (molm-3s-1), T (Wm-3),
            array
        ubc - upper boundary condition, value of CO2 (ppm), H2O (mol/mol)
            or T (degC) at highest gridpoint
        lbc - lower boundary condition, flux or value:
            flux: CO2 (umol m-2 s-1), H2O (mol m-2 s-1), T (Wm-2).
            value: CO2 (ppm), H2O (mol/mol), T (degC)
        scalar - 'CO2', 'H2O', 'T', string
        T - air temperature (degC), optional
        P - pressure (Pa)
        lbc_dirchlet - True for Dirchlet lower boundary
    OUT:
        Ca - mixing ratio profile (ppm)
    NOTE:
        assumes constant dz and Dirchlet upper boundary condition
    Gaby Katul & Samuli Launiainen, 2009 - 2017
    """
    epsi = eps #1e-5
    dz = float(dz)
    N = len(Ks)
    rho_a = P / (287.05*(T + 273.15))  # air density, kg/m3

    # -----Set elements of tridiagonal matrix
    a1 = Ks
    a2 = np.zeros(N)
    a2[1:] = np.diff(Ks) / dz
    a2[0] = a2[1]
    a3 = np.zeros(N)

    # super-diagonal, diagonal, subdiagonal elements and rhs vector
    upd = a1 / (dz*dz) + a2 / (2*dz)
    dia = -a1*2 / (dz*dz) + a3
    lod = a1 / (dz*dz) - a2 / (2*dz)
    rhs = np.zeros(N)

    # --- uppermost node, Dirchlet boundary
    upd[-1] = 0.
    dia[-1] = 1.
    lod[-1] = 0.

    # CO2
    if scalar.upper() == 'CO2':
        CF = 1e6*rho_a / (29e-3)  # umol m-3, molar conc. of air
        rhs = -Source / CF  # CO2 source/sink profile umol m-3 s-1 --> s-1.'
        rhs[-1] = 1e-6*ubc
        if not lbc_dirchlet:  # --- lower BC, flux-based
            # print 'lbc flux based'
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = (lbc / CF)*dz / (Ks[0] + epsi)
        else:  # --- lower BC, fixed concentration
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = 1e-6*lbc
        # --- call tridiagonal solver
        x = 1e6*tridiag(lod, dia, upd, rhs)  # ppm

    # H2O
    if scalar.upper() == 'H2O':
        CF = rho_a / (29e-3)  # mol m-3, molar conc. of air
        rhs = - Source / CF  # source/sink profile mol m-3 s-1 --> s-1.'
        rhs[-1] = ubc
        if not lbc_dirchlet:  # --- lower BC, flux-based
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = (lbc / CF)*dz / (Ks[0] + epsi)
        else:  # --- lower BC, fixed concentration
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = lbc
        # --- call tridiagonal solver
        x = tridiag(lod, dia, upd, rhs)  # mol/mol

    # Temperature
    if scalar.upper() == 'T':
        CF = 1004.0*rho_a  # cp*rho_a,  J m-3 K-1
        rhs = - Source / CF  # heat source/sink profile W m-3 --> K s-1.'
        rhs[-1] = ubc
        if not lbc_dirchlet:  # --- lower BC, flux-based
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = (lbc / CF)*dz / (Ks[0] + epsi)  # m s-1
            # print rhs
        else:  # --- lower BC, fixed value
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = lbc
        # --- call tridiagonal solver
        x = tridiag(lod, dia, upd, rhs)  # degC
    return x

def mixing_length(z, h, d, l_min=None):
    """
    computes mixing length: linear above the canopy, constant within and
    decreases linearly close the ground (below z< alpha*h/VON_KARMAN)
    IN:
        z - computation grid, m
        h - canopy height, m
        d - displacement height, m
        l_min - intercept at ground, m
    OUT:
        lmix - mixing length, m
    """
    dz = z[1] - z[0]

    if not l_min:
        l_min = dz / 2.0
    
    alpha = (h - d)*VON_KARMAN / h
    I_F = np.sign(z - h) + 1.0
    l_mix = alpha*h*(1 - I_F / 2) + (I_F / 2) * (VON_KARMAN*(z - d))
    
    sc = (alpha*h) / VON_KARMAN
    ix = np.where(z < sc)
    l_mix[ix] = VON_KARMAN*(z[ix] + dz / 2)
    l_mix = l_mix + l_min

    return l_mix

def leaf_boundary_layer_conductance(u, d, Ta, dT, P=101300.):
    """
    Computes 2-sided leaf boundary layer conductance assuming mixed forced and free
    convection form two parallel pathways for transport through leaf boundary layer.
    INPUT: u - mean velocity (m/s)
           d - characteristic dimension of the leaf (m)
           Ta - ambient temperature (degC)
           dT - leaf-air temperature difference (degC)
           P - pressure(Pa)
    OUTPUT: boundary-layer conductances (mol m-2 s-1)
        gb_h - heat (mol m-2 s-1)
        gb_c- CO2 (mol m-2 s-1)
        gb_v - H2O (mol m-2 s-1)
        r - ratio of free/forced convection
    Reference: Campbell, S.C., and J.M. Norman (1998),
    An introduction to Environmental Biophysics, Springer, 2nd edition, Ch. 7
    Gaby Katul & Samuli Launiainen
    Note: the factor of 1.4 is adopted for outdoor environment, see Campbell and Norman, 1998, p. 89, 101.
    """

    u = np.maximum(u, eps)

    # print('U', u, 'd', d, 'Ta', Ta, 'P', P)
    factor1 = 1.4*2  # forced conv. both sides, 1.4 is correction for turbulent flow
    factor2 = 1.5  # free conv.; 0.5 comes from cooler surface up or warmer down

    # -- Adjust diffusivity, viscosity, and air density to pressure/temp.
    t_adj = (101300.0 / P)*((Ta + 273.15) / 293.16)**1.75
    Da_v = MOLECULAR_DIFFUSIVITY_H2O*t_adj
    Da_c = MOLECULAR_DIFFUSIVITY_CO2*t_adj
    Da_T = THERMAL_DIFFUSIVITY_AIR*t_adj
    va = AIR_VISCOSITY*t_adj
    rho_air = 44.6*(P / 101300.0)*(273.15 / (Ta + 273.13))  # [mol/m3]

    # ----- Compute the leaf-level dimensionless groups
    Re = u*d / va  # Reynolds number
    Sc_v = va / Da_v  # Schmid numbers for water
    Sc_c = va / Da_c  # Schmid numbers for CO2
    Pr = va / Da_T  # Prandtl number
    Gr = GRAVITY*(d**3)*abs(dT) / (Ta + 273.15) / (va**2)  # Grashoff number

    # ----- aerodynamic conductance for "forced convection"
    gb_T = (0.664*rho_air*Da_T*Re**0.5*(Pr)**0.33) / d  # [mol/m2/s]
    gb_c=(0.664*rho_air*Da_c*Re**0.5*(Sc_c)**0.33) / d  # [mol/m2/s]
    gb_v=(0.664*rho_air*Da_v*Re**0.5*(Sc_v)**0.33) / d  # [mol/m2/s]

    # ----- Compute the aerodynamic conductance for "free convection"
    gbf_T = (0.54*rho_air*Da_T*(Gr*Pr)**0.25) / d  # [mol/m2/s]
    gbf_c = 0.75*gbf_T  # [mol/m2/s]
    gbf_v = 1.09*gbf_T  # [mol/m2/s]

    # --- aerodynamic conductance: "forced convection"+"free convection"
    gb_h = factor1*gb_T + factor2*gbf_T
    gb_c = factor1*gb_c + factor2*gbf_c
    gb_v = factor1*gb_v + factor2*gbf_v
    # gb_o3=factor1*gb_o3+factor2*gbf_o3

    #r = Gr / (Re**2)  # ratio of free/forced convection

    return gb_h, gb_c, gb_v#, r


def e_sat(T):
    """
    Computes saturation vapor pressure (Pa), slope of vapor pressure curve
    [Pa K-1]  and psychrometric constant [Pa K-1]
    IN:
        T - air temperature (degC)
    OUT:
        esa - saturation vapor pressure in Pa
        s - slope of saturation vapor pressure curve (Pa K-1)
    SOURCE:
        Campbell & Norman, 1998. Introduction to Environmental Biophysics. (p.41)
    """

    esa = 611.0 * np.exp((17.502 * T) / (T + 240.97))  # Pa
    s = 17.502 * 240.97 * esa / ((240.97 + T)**2)

    return esa, s