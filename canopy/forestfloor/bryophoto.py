# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 12:38:05 2018

@author: slauniai
"""
import matplotlib.pyplot as plt
import numpy as np
#from scipy.integrate import odeint
EPS = np.finfo(float).eps  # machine epsilon

# Constants used in the model calculations.
#: [J mol\ :sup:`-1`\ ], latent heat of vaporization at 20\ :math:`^{\circ}`\ C
LATENT_HEAT = 44100.0
# LATENT_HEAT = 2,501e6 #  J kg
#: [kg mol\ :sup:`-1`\ ], molar mass of H\ :sub:`2`\ O
MOLAR_MASS_H2O = 18.015e-3
#: [kg mol\ :sup:`-1`\ ], molar mass of CO\ :sub:`2`\
MOLAR_MASS_CO2 = 44.01e-3
#: [kg mol\ :sup:`-1`\ ], molar mass of C
MOLAR_MASS_C = 12.01e-3
#: [J kg\ :sup:`-1` K\ :sup:`-1`\ ], specific heat of H\ :sub:`2`\ O
SPECIFIC_HEAT_H2O = 4.18e3
#: [J kg\ :sup:`-1` K\ :sup:`-1`\ ], specific heat of organic matter
SPECIFIC_HEAT_ORGANIC_MATTER = 1.92e3
#: [J mol\ :sup:`-1` K\ :sup:`-1`\ ], heat capacity of air at constant pressure
SPECIFIC_HEAT_AIR = 29.3
#: [W m\ :sup:`-2` K\ :sup:`-4`\ ], Stefan-Boltzmann constant
STEFAN_BOLTZMANN = 5.6697e-8
#: [K], zero degrees celsius in Kelvin
DEG_TO_KELVIN = 273.15

#: [K], zero degrees celsius in Kelvin
NORMAL_TEMPERATURE = 273.15
#: [mol m\ :sup:`-3`\ ], density of air at 20\ :math:`^{\circ}`\ C
AIR_DENSITY = 41.6
#: [m\ :sup:`2` s\ :sup:`-1`\ ], kinematic viscosity of air at 20\ :math:`^{\circ}`\ C
AIR_VISCOSITY = 15.1e-6
#: [m\ :sup:`2` s\ :sup:`-1`\ ], thermal diffusivity of air at 20\ :math:`^{\circ}`\ C
THERMAL_DIFFUSIVITY_AIR = 21.4e-6
#: [m\ :sup:`2` s\ :sup:`-1`\ ], molecular diffusvity of CO\ :sub:`2` at 20\ :math:`^{\circ}`\ C
MOLECULAR_DIFFUSIVITY_CO2 = 15.7e-6
#: [m\ :sup:`2` s\ :sup:`-1`\ ], molecular diffusvity of H\ :sub:`2`\ at 20\ :math:`^{\circ}`\ C
MOLECULAR_DIFFUSIVITY_H2O = 24.0e-6
#: [J mol\ :sup:`-1` K\ :sup:``-1], universal gas constant
GAS_CONSTANT = 8.314
#: [kg m\ :sup:`2` s\ :sup:`-1`\ ], standard gravity
GRAVITY = 9.81
#: [kg m\ :sup:`-3`\ ], water density
WATER_DENSITY = 1.0e3


def net_co2_exchange(para, Qp, Ca, T, w, wstar):
    """
    computes net CO2 exchange of moss community
    Args:
        para - parameter dictionary
        Qp - incident PAR umolm-2s-1
        Ca - ambient CO2 (ppm)
        T - moss temperature (degC)
        w - moss water content (g/g)
        wstar - delayed water content (g/g) for desiccation recovery
    Returns:
        An - net CO2 exchange An = -A + Rd, <0 is uptake (umolm-2s-1 or umol g-1 s-1)
        A - photosynthesis rate (umolm-2s-1 or umol g-1 s-1)
        Rd - dark respiration rate (umolm-2s-1 or umol g-1 s-1)
        Cc - internal CO2 (ppm)
        g - total conductance for CO2 (molm-2s-1 or mol g-1 s-1)
    """
    p = para.copy()
    cap, rcap = relative_capacity(p, w, wstar)
    p['Vcmax'] *= cap 
    p['Jmax'] *= cap
    p['alpha'] *= cap
    p['Rd'] *= rcap
        
    # conductance (mol m-2 s-1)
    g = conductance(p, w)
    
    # solve Anet and Cc iteratively until Cc converges
    err = 10^9
    Cc = 0.8*Ca

    while err > 1e-3:
        Cco = Cc
        An, Rd, Av, Aj = photo_farquhar(p, Qp, Cc, T)
        
        Cc = Ca - An / g  # new Cc
        Cc = 0.5*(Cco + Cc)
        err = np.nanmax(abs(Cc - Cco))
    return -An, An - Rd, Rd, Cc, g
    

def conductance(para, w):
    """
    Conductance for CO2 diffusion from bulk air to chloroplast in bryophyte.
    Assumes g = gmax * fw, where gmax is species-specific maximum internal conductance, 
    occurring at w <= wopt, and fw (-) describes decay of conductance due to
    external water (w > wopt).
    
    gmax and wopt are bryophyte traits, while the shape of fw is fixed based on
    data of Williams & Flanagan, 1998. We normalized maximum conductance for Pleurozium and
    Sphagnum to unity, and determined respective wopt's. Then a decreasing exponential 
    function was fitted to data.
    
    Args:
        para - parameters: gmax, wopt, a
        w - water content (g/g)
    Returns:
        g (mol m-2 s-1)
    """
 
    gmax = para['gmax']
    a0 = para['a0']
    a1 = para['a1']
    wopt = para['wopt']
    
    #g = gmax * np.minimum(1.0, a0*np.exp(a1*(w-wopt)) + (1.0 - a0))
    g = gmax * (a0*np.exp(a1*(w-wopt)) + (1.0 - a0)) # this allows g to go above gmax in dry moss
    return g


def relative_capacity(para, w, wstar):
    """
    Relative photosynthetic capacity and dark respiration rate as a function of water content
    Args:
        para - parameter dictionary
        w - current water content (g/g)
        wstar - delayed effective water content for desiccation recovery (g/g)
    Returns:
         rphoto, rrd - relative photosynthetic capacity and dark respiration rate
    """
    
    p = para['CAP_desic']
    
    # r = para['CAP_rewet']
    
    # drying phase is function of w
    cap_dec = np.maximum(0.0, np.minimum(1.0 + p[0]*np.log(w/p[1]), 1.0))
    
    # recovery from desiccation is a function of wstar; now assume reversible
    cap_rw = 1.0
    
    
    rphoto = np.minimum(cap_dec, cap_rw)
    del p #, r
    
    # respiration
    # p = para['Rd_desic']
    # r = para['Rd_rewet]
    
    #rrd = 1.0
    rrd = rphoto # assume to behave as photo
    return rphoto, rrd

def photo_farquhar(photop, Qp, ci, T, co_limi=False):
    """
    Calculates leaf net CO2 exchange and dark respiration rate (umol m-2 s-1).
    INPUT:
        photop - dict with keys:
            Vcmax
            Jmax
            Rd
            qeff
            alpha
            theta
            beta
        Qp - incident Par (umolm-2s-1)
        ci - leaf internal CO2 mixing ratio (ppm)
        T - leaf temperature (degC)
    OUTPUT:
        An - leaf net CO2 exchange (umol m-2 s-1)
        Rd - leaf dark respiration rate (umol m-2 s-1)
        Av - rubisco limited photo (umol m-2 s-1)
        Aj - RuBP -regeneration limited photo (umol m-2 s-1)
        
    """
    Tk = T + 273.15  # K

    # --- constants
    Coa = 2.10e5  # O2 in air (umol/mol)
    TN = 298.15  # reference temperature 298.15 K = 25degC
    R = 8.314427  # gas constant, J mol-1 K-1

    # --- params ----
    Vcmax = photop['Vcmax']
    Jmax = photop['Jmax']
    Rd = photop['Rd']
    alpha = photop['alpha']
    theta = photop['theta']
    beta = photop['beta']  # co-limitation parameter

    # --- CO2 compensation point -------
    Tau_c = 42.75 * np.exp(37830*(Tk - TN) / (TN * R * Tk))

    # ---- Kc & Ko (umol/mol), Rubisco activity for CO2 & O2 ------
    Kc = 404.9 * np.exp(79430.0*(Tk - TN) / (TN * R * Tk))
    Ko = 2.784e5 * np.exp(36380.0*(Tk - TN) / (TN * R * Tk))

    if 'tresp' in photop:  # adjust parameters for temperature
        tresp = photop['tresp']
        Vcmax_T = tresp['Vcmax']
        Jmax_T = tresp['Jmax']
        Rd_T = tresp['Rd']
        Vcmax, Jmax, Rd, Tau_c = photo_temperature_response(Vcmax, Jmax, Rd, Vcmax_T, Jmax_T, Rd_T, Tk)

    Km = Kc*(1.0 + Coa / Ko)
    J = (Jmax + alpha*Qp -((Jmax + alpha*Qp)**2.0 - (4.0*theta*Jmax*alpha*Qp))**0.5) / (2.0*theta)

    # -- rubisco -limited rate
    Av = Vcmax * (ci - Tau_c) / (ci + Km)
    # -- RuBP -regeneration limited rate
    Aj = J/4 * (ci - Tau_c) / (ci + 2.0*Tau_c)

    # An = np.minimum(Av, Aj) - Rd  # single limiting rate
    x = Av + Aj
    y = Av * Aj
    An = (x - (x**2 - 4*beta*y)**0.5) / (2*beta) - Rd  # co-limitation
    return An, Rd, Av, Aj


def photo_temperature_response(Vcmax0, Jmax0, Rd0, Vcmax_T, Jmax_T, Rd_T, T):
    """
    Adjusts Farquhar / co-limitation optimality model parameters for temperature
    INPUT:
        Vcmax0, Jmax0, Rd0 - parameters at ref. temperature 298.15 K
        Vcmax_T, Jmax_T, Rd_T - temperature response parameter lists
        T - leaf temperature (K)
    OUTPUT: Nx1-arrays
        Vcmax, Jmax,Rd (umol m-2(leaf) s-1)
        Gamma_star - CO2 compensation point (ppm)
    CALLED from Farquhar()
    REFERENCES:
        Medlyn et al., 2002.Plant Cell Environ. 25, 1167-1179; based on Bernacchi
        et al. 2001. Plant Cell Environ., 24, 253-260.
    Samuli Launiainen, Luke, 28.3.2017
    """

    # ---constants
    NT = 273.15
    TN = 298.15  # reference temperature 298.15 K = 25degC
    R = 8.314427  # gas constant, J mol-1 K-1

    # --- CO2 compensation point -------
    Gamma_star = 42.75 * np.exp(37830*(T - TN) / (TN * R * T))

    # # ---- Kc & Ko (umol/mol), Rubisco activity for CO2 & O2 ------
    # Kc = 404.9 * np.exp(79430.0*(T - TN) / (TN * R * T))
    # Ko = 2.784e5 * np.exp(36380.0*(T - TN) / (TN * R * T))

    # ------  Vcmax (umol m-2(leaf)s-1) ------------
    Ha = 1e3 * Vcmax_T[0]  # J mol-1, activation energy Vcmax
    Hd = 1e3 * Vcmax_T[1]  # J mol-1, deactivation energy Vcmax
    Sd = Vcmax_T[2]  # entropy factor J mol-1 K-1

    NOM = np.exp(Ha * (T - TN) / (R*NT*T)) * (1.0 + np.exp((NT*Sd - Hd) / (NT*R)))
    DENOM = (1.0 + np.exp((T*Sd - Hd) / (T*R)))
    Vcmax = Vcmax0 * NOM / DENOM

    del Ha, Hd, Sd, DENOM, NOM

    # ----  Jmax (umol m-2(leaf)s-1) ------------
    Ha = 1e3 * Jmax_T[0]  # J mol-1, activation energy Vcmax
    Hd = 1e3 * Jmax_T[1]  # J mol-1, deactivation energy Vcmax
    Sd = Jmax_T[2]  # entropy factor J mol-1 K-1

    NOM = np.exp(Ha * (T - TN) / (R*NT*T)) * (1.0 + np.exp((NT*Sd - Hd) / (NT*R)))
    DENOM = (1.0 + np.exp((T*Sd - Hd) / (T*R)))
    Jmax = Jmax0*NOM / DENOM

    del Ha, Hd, Sd, DENOM, NOM

    # --- Rd (umol m-2(leaf)s-1) -------
    Ha = 1e3 * Rd_T[0]  # J mol-1, activation energy dark respiration
    Rd = Rd0 * np.exp(Ha*(T - TN) / (TN * R * T))

    return Vcmax, Jmax, Rd, Gamma_star


def apparent_photocapacity(b, psi_leaf):
    """
    computes relative photosynthetic capacity as a function of leaf water potential
    Function shape from Kellomäki & Wang, adjustments for Vcmax and Jmax
    IN:
       beta - parameters, 2x1 array
       psi - leaf water potential (MPa)
    OUT:
       f - relative value [0.2 - 1.0]
    """
    psi_leaf = np.array(np.size(psi_leaf), ndmin=1)
    f = (1.0 + np.exp(b[0] * b[1])) / (1.0 + np.exp(b[0] * (b[1] - psi_leaf)))
    f[f < 0.2] = 0.2

    return f


def topt_deltaS_conversion(xin, Ha, Hd, var_in='deltaS'):
    """
    Converts between entropy factor Sv (J mol-1) and temperature optimum
    Topt (K). Medlyn et al. 2002 PCE 25, 1167-1179 eq.19.
    INPUT:
        xin, Ha(J mol-1), Hd(J mol-1)
        input:'deltaS' [Jmol-1] or 'Topt' [K]
    OUT:
        xout - Topt or Sv
    Farquhar parameters temperature sensitivity
    """
    R = 8.314427  # gas constant, J mol-1 K-1
    
    if var_in.lower() == 'deltas':  # Sv --> Topt
        xout = Hd / (xin - R * np.log(Ha / (Hd - Ha)))
    else:  # Topt -->Sv
        c = R * np.log(Ha / (Hd - Ha))
        xout = (Hd + xin * c) / xin
    return xout


def draw_farquhar_curves():
    Vcmax = 30.0
    Jmax = 1.9*Vcmax
    Rd = 0.15*Vcmax
    tresp = {'Vcmax': [69.83, 200.0, 27.56],
             'Jmax': [100.28, 147.92, 19.8],
             'Rd': [33.0], 'include': 'y'}
    
    photop = {'Vcmax': Vcmax, 'Jmax': Jmax, 'Rd': Rd, # umolm-2s-1
              'alpha': 0.24, 'theta': 0.8, 'beta': 0.95, # quantum yield, curvature, co-limitation
              'gmax': 0.06, 'wopt': 7.0, 'a0': 0.7, 'a1': -0.263, 'CAP_desic': [0.53, 7.0],
              'tresp': tresp}

    Qp = np.linspace(0.0, 1600.0, 100)
    cc = 0.8*400
    T = 25.0
    
    An, Rd, Av, Aj = photo_farquhar(photop, Qp, cc, T, co_limi=False)
    photop['beta'] = 1.0
    An1, Rd1, _, _ = photo_farquhar(photop, Qp, cc, T, co_limi=False)

    plt.figure(1)
    # plt.subplot(121)
    plt.plot(Qp, An+Rd, 'k-', Qp, Aj, 'b--', Qp[[0,-1]], [Av,Av], 'r--', Qp, An1+Rd1, 'g--')
    plt.legend(['A', 'A_j', 'A_v', 'A ($\\beta$=1.0)'])
    plt.ylabel('A ($\mu$mol m$^{-2}$ s$^{-1}$)')
    plt.xlabel('absorbed PAR ($\mu$mol m$^{-2}$ s$^{-1}$)')
    plt.title('$C_c=320\, ppm, V_{c,max}=30, J_{max}=1.9\,V_{c,max}, \gamma=0.24, \\theta=0.8, \\beta=0.98, T=25degC$', fontsize=10)


def test():
    
    Vcmax = 15.0
    Jmax = 1.9*Vcmax
    Rd = 0.05*Vcmax
    tresp = {'Vcmax': [69.83, 200.0, 27.56],
             'Jmax': [100.28, 147.92, 19.8],
             'Rd': [33.0], 'include': 'y'}
    
    # pleurozium type
    photop_p = {'Vcmax': Vcmax, 'Jmax': Jmax, 'Rd': Rd, # umolm-2s-1
              'alpha': 0.3, 'theta': 0.8, 'beta': 0.9, # quantum yield, curvature, co-limitation
              'gmax': 0.02, 'wopt': 7.0, 'a0': 0.7, 'a1': -0.263, 'CAP_desic': [0.44, 7.0],
              'tresp': tresp}

    # sphagnum type
    Vcmax = 45.0
    photop_s = {'Vcmax': Vcmax, 'Jmax': 1.9*Vcmax, 'Rd': 0.03*Vcmax, # umolm-2s-1
              'alpha': 0.3, 'theta': 0.8, 'beta': 0.9, # quantum yield, curvature, co-limitation
              'gmax': 0.04, 'wopt': 7.0, 'a0': 0.7, 'a1': -0.263, 'CAP_desic': [0.58, 10.0],
              'tresp': tresp}        
    
    # moisture levels
    w = np.linspace(1, 20 , 50)
    wstar = w.copy()
    
    # make few response plots
    
    # moisture response at 400 and 600 ppm and at light-saturated cond.
    #Qp = np.linspace(0.0, 100, 100)
    Qp = 500.0
    T = 20.0
    fig, (ax) = plt.subplots(ncols=2, nrows=2)
    
    style=['-', '--', '-.']
    k = 0
    for Ca in [400.0, 800.0]:        
        for gmax in [0.03, 0.06]:
            photop_p['gmax'] = gmax
            photop_s['gmax'] = gmax
            A, An, Rd, Cc, g = net_co2_exchange(photop_p, Qp, Ca, T, w, wstar)
            ax[0,0].plot(w, An, linestyle=style[k], label='Ca=' + str(Ca) + ', gmax=' + str(gmax))
            ax[1,0].plot(w, Cc/Ca, linestyle=style[k])
            
            A, An, Rd, Cc, g = net_co2_exchange(photop_s, Qp, Ca, T, w, wstar)
            ax[0,1].plot(w, An, linestyle=style[k], label='Ca=' + str(Ca) + ', gmax=' + str(gmax))
            ax[1,1].plot(w, Cc/Ca, linestyle=style[k])
        k += 1

    ax[0,0].set_title('"Pleurozium"'); 
    ax[0,0].set_ylabel('A umolm-2s-1'); ax[0,0].set_xlabel('w (g/g)')
    ax[1,0].set_ylabel('Cc/Ca'); ax[1,0].set_xlabel('w (g/g)')
    
    ax[0,1].set_title('"Sphagnum"');
    ax[0,1].set_xlabel('w (g/g)')
    ax[1,1].set_xlabel('w (g/g)')
    #ax[0].legend(ncol=1)
    