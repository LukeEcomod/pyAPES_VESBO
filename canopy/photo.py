# -*- coding: utf-8 -*-
"""
Created on Mon May 15 13:43:44 2017

@author: slauniai

photo module contains leaf-scale functions for photosynthesis and stomatal control

"""
import numpy as np
import matplotlib.pyplot as plt
eps = np.finfo(float).eps  # machine epsilon


# ----- CONSTANTS ----------
CP = 29.3  # J/mol/K, heat capacity of air at constant pressure
LMOLAR = 44100.0  # J/mol latent heat of vaporization at 20 deg C
PARCF = 4.56  # conversion from Wm-2 to micromol m-2 s-1
SIGMA = 5.6697e-8  # W m-2 K-4 Stefan-Boltzmann const
a = 1.6  # H2O to CO2 diffusivity ratio [-]
NT = 273.15
R = 8.314427  # universal gas constant, J mol-1 K-1

"""--- scripts for testing functions ---- """


def test_leafscale(method=1):
    Vcmax = 70.75
    Jmax = 78.36
    Rd = 0.023*Vcmax
    tresp = {'Vcmax': [69.83, 200.0, 27.56],
             'Jmax': [100.28, 147.92, 19.8],
             'Rd': [33.0], 'include': 'y'}
    photop = {'Vcmax': Vcmax, 'Jmax': Jmax, 'Rd': Rd,
              'alpha': 0.3, 'theta': 0.8, 'beta': 0.8, 'La': 1600.0, 'm': 2.3,
              'g0': 1e-3, 'kn': 0.0, 'tresp': tresp}

    leafp = {'lt': 0.02, 'emi': 0.98, 'par_alb': 0.12, 'nir_alb': 0.55}

    # env. conditions
    P = 101300.0
    Qp = np.array([1.0, 100., 300., 600., 1200., 1600.])
    #LW = np.array([0.0, -20.0, -30.0, -50.0])
    LW = np.zeros(np.shape(Qp))
    H2O = np.ones(6) * 1.0e3 / P
    CO2 = np.ones(6) * 400.0
    U = 1.0  # np.array([10.0, 1.0, 0.1, 0.01])
    T = 20.0  # np.array([20.0, 19.0, 18.0, 20.0])
    
    SWabs = 0.5*(1-leafp['par_alb'])*Qp / PARCF + 0.5*(1-leafp['nir_alb'])*Qp / PARCF 
    print('SWabs', SWabs)
    if method is not 1:
        x = leaf_interface(photop, leafp, H2O, CO2, T, Qp, SWabs, LW, U, P=101300.0, model=method, Ebal=False, dict_output=True)  
        print x
        plt.figure()
        plt.subplot(221); plt.plot(Qp, x['An'], 'o')
        plt.subplot(222); plt.plot(Qp, x['E'], 'o')
        plt.subplot(223); plt.plot(Qp, x['gs_v'], 'o')
    if method == 1:
        #ci = 300.0
        #Qp = np.arange(10, 1600, 20)
        Qp = 1600.0
        ci = np.arange(200, 700, 10)        
        an, rd, av, aj = photo_farquhar(photop, Qp, ci, T, co_limi=False )
        an1, rd1 = photo_farquhar(photop, Qp, ci, T, co_limi=True)
        plt.figure()
        plt.title('photo.photo_farquhar')
        plt.plot(ci, av, 'k-', ci, aj, 'k--')
        plt.plot(ci, an, 'r.-', ci, an1, 'g.-')
        plt.xlabel('ci (ppm)')
        #plt.ylabel('An (umolm-2s-1)')
        #plt.plot(Qp, an, 'r.-', Qp, an1, 'g.-')
        

def leaf_interface(photop, leafp, H2O, CO2, T, Qp, SWabs, LW, U, P=101300.0, model='CO_OPTI', Ebal=True,
                   dict_output=True):
    """
    Entry-point to coupled leaf gas-exchange and energy balance functions.
    
    CALCULATES leaf photosynthesis (An), respiration (Rd), transpiration (E) and estimates of
    leaf temperature (Tl) and sensible heat fluxes (H) based onleaf energy balance equation coupled with
    leaf-level photosynthesis and stomatal control schemes.
    Energy balance is solved using Taylor's expansion (i.e isothermal net radiation -approximation) which
    eliminates need for iterations with radiation-sceme.
    
    Depending on choise of 'model', photosynthesis is calculated based on biochemical model of Farquhar et
    al. (1980) coupled with various stomatal control schemes (Medlyn, Ball-Woodrow-Berry, Hari, Katul-Vico et al.)
    In all these models, stomatal conductance (gs) is directly linked to An, either by optimal stomatal control principles or
    using semi-empirical models.
    
    INPUT:
        photop - dict. of photoparameters, keys:
            Vcmax - maximum carboxylation velocity (umolm-2s-1)
            Jmax - maximum rate of electron transport (umolm-2s-1)
            Rd - dark respiration rate (umolm-2s-)
            alpha - quantum yield parameter (mol/mol)
            theta - co-limitation parameter of Farquhar-model
            beta - co-limitation parameter of Farquhar-model
            L, m - stomatal parameter (Lambda, m, ...) depending on model
            g0 - residual conductance for CO2 (molm-2s-1)
            
            tresp - dict. of temperature sensitivity parameters, optional.
                omitting neglects temperature adjustments of Vcmax, Jmax, Rd
                Vcmax - [Ha, Hd, Topt]; activation energy (kJmol-1), deactivation energy (kJmol-1), optimum temperature (degC)
                Jmax - [Ha, Hs, Topt];
                Rd - [Ha]; activation energy (kJmol-1)
        leafp - dict of leaf parameters, keys:
            lt - leaf lengthscale (m)
            emi - leaf emissivity (-)
            # Par_alb - leaf Par albedo (-)
            # Nir_alb - leaf Nir albedo (-)
        
        H2O - water vapor mixing ratio (mol/mol)
        CO2 - carbon dioxide mixing ratio (ppm)
        T - ambient air temperature (degC)
        Qp - incident PAR at leaves (umolm-2s-1)
        SWabs - absorbed SW (PAR + NIR) at leaves (Wm-2)
        LW - net isothermal long-wave radiation (Wm-2).
        U - mean wind speed (m/s)
        P - ambient pressure (Pa)

        model - CO_OPTI (Vico et al., 2014)
                MEDLYN (Medlyn et al., 2011 with co-limitation Farquhar)
                BWB (Ball et al., 1987 with co-limitation Farquhar)
                ADD others!!!
        Ebal - True computes leaf temperature by solving energy balance
        dict_output - True returns output as dict, False as separate arrays
      
    OUTPUT:
        x - dict with keys: (or separate arrays in following order)
            An - net CO2 flux (umol m-2 leaf s-1)
            Rd - CO2 respiration (umol m-2 leaf s-1)
            E - H2O flux (transpiration) mol m-2 leaf s-1)
            H - sensible heat flux (W m-2 leaf)
            Fr - non-isothermal radiative flux (W m-2)
            Tl - leaf temperature (degC)
            Ci - leaf internal CO2 mixing ratio (mol/mol)
            Cs - leaf surface CO2 mixing ratio (mol/mol)
            #Cm - mesophyll CO2 mixing ratio (mol/mol)
            gs_v - stomatal conductance for H2O (mol m-2 leaf s-1)
            gs_c - stomatal conductance for CO2 (mol m-2 leaf s-1)
            gb_v - boundary layer conductance for H2O (mol m-2 leaf s-1)
            gbv - leaf boundary layer conductance for H2O (mol m-2 leaf s-1)

    NOTE: Vectorized code can be used in multi-layer sense where inputs are vectors of equal length

    Samuli Launiainen LUKE 3/2011 - 5/2017
    Last edit 16.5.2017
    """
    # -- Parameters -----
    a = 1.6  # H2O to CO2 diffusivity ratio
    lt = leafp['lt']
    ef = leafp['emi']

    # radiative conductance (mol m-2 s-1), Campbell & Norman, 1998
    gr = 4.0*ef*SIGMA*(T+NT)**3 / CP
    gr = 2.0*gr  # LW is emited from both sides

    #
    T = np.array(T)
    Qp = np.array(Qp)
    U = np.array(U)
    H2O = np.array(H2O, ndmin=1)
    CO2 = np.array(CO2, ndmin=1)
    Rabs = np.array(SWabs + LW, ndmin=1)  # isothermal net radiation (Wm-2)
    Tl = T  
    
#    print T, Qp, H2O, CO2, T, Rabs, U, P
    gb_h, gb_v, gb_c, _ = leaf_boundary_layer_conductance(U, lt, T, Tl - T, P)

    err = 999.0
    iterNo = 0
    while err > 0.01 and iterNo < 50:
        iterNo += 1
        #print iterNo
        esat, s = saturation_vapor_pressure(Tl)  
        s = s / P  # mol/mol
        Dleaf = esat / P - H2O  # mol/mol
        Dleaf[Dleaf < 0] = eps

        #print Dleaf
        Told = Tl

        # --- analytical co-limitation model Vico et al. 2013
        if model.upper() == 'CO_OPTI':
            An, Rd, fe, gs_opt, Ci, Cs = photo_c3_analytical(photop, Qp, Tl, Dleaf, CO2, gb_c, gb_v)
            gsv = a*gs_opt
        if model.upper() == 'MEDLYN':
            An, Rd, fe, gs_opt, Ci, Cs = photo_c3_medlyn(photop, Qp, Tl, Dleaf, CO2, gb_c, gb_v, P=P)
            gsv = a*gs_opt
        if model.upper() == 'MEDLYN_FARQUHAR':
            An, Rd, fe, gs_opt, Ci, Cs = photo_c3_medlyn_farquhar(photop, Qp, Tl, Dleaf, CO2, gb_c, gb_v, P=P)
            gsv = a*gs_opt
        if model.upper() == 'BWB':
            rh  = (1 - Dleaf*P / esat)  # rh at leaf (-)
            An, Rd, fe, gs_opt, Ci, Cs = photo_c3_bwb(photop, Qp, Tl, rh, CO2, gb_c, gb_v, P=P)
            gsv = a*gs_opt

        # solve  energy balance
        if Ebal:
            # solve leaf temperature using Taylor's expansion for LW term.
            # For abs(dT)<15 K error is <10Wm-2. Campbell & Norman 1998, p.225
            geff_v = (gb_v*gsv) / (gb_v + gsv)  # molm-2s-1
            Tl = T + (Rabs - LMOLAR*geff_v*Dleaf) / (CP*(gb_h + gr) + LMOLAR*s*geff_v)
            gb_h, gb_v, gb_c, _ = leaf_boundary_layer_conductance(U, lt, T, Tl - T, P)
            err = np.nanmax(abs(Tl - Told))
            # print('err', err)
        else:
            err = 0.0
            H = None
            Fr = None
            E = None

    # outputs
    H = CP*gb_h*(Tl - T)  # Wm-2
    Fr = CP*gr*(Tl - T)  # flux due to radiative conductance (Wm-2)
    E = fe
    
#    if any(np.isnan(An)):
#        print('leafinterface: ', Tl, Qp, Dleaf, CO2, gb_c, gb_v )

    if dict_output:  # return dict
        x = {'An': An, 'Rd': Rd, 'E': fe, 'H': H, 'Fr': Fr, 'Tl': Tl, 'Ci': Ci,
             'Cs': Cs, 'gs_v': gsv, 'gs_c': gs_opt, 'gb_v': gb_v}
        return x
    else:  # return 11 arrays
        return An, Rd, E, H, Fr, Tl, Ci, Cs, gsv, gs_opt, gb_v


def saturation_vapor_pressure(T):
    """
    Computes saturation vapor pressure with respect to free and flat water surface for given temperature T
    INPUT:
        T - temperature (degC), scalar or array
    OUTPUT:
        esat - saturation vapor pressure (Pa)
        delta - slope of saturation vapor pressure curve (Pa degC-1)
    SOURCE:
        Campbell & Norman, 1998. Introduction to Environmental Biophysics.
    """
    
    #constants
    a = 611.0  # Pa
    b = 17.502  # (-)
    c = 240.97  # degC
    
    esat = a*np.exp(b*T / (T+c))  # Pa
    delta = b*c*esat / ((c + T)**2)  # Pa degC-1   
    return esat, delta
    

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
    
    # print('U', u, 'd', d, 'Ta', Ta, 'P', P)
    factor1 = 1.4*2  # forced conv. both sides, 1.4 is correction for turbulent flow
    factor2 = 1.5  # free conv.; 0.5 comes from cooler surface up or warmer down
    
    Da_v = 2.4e-5  # Molecular diffusivity of "water vapor" in air at STP (20C and 11kPa) [m2/s]
    Da_c = 1.57e-5  # Molecular diffusivity of "CO2" in air at STP [m2/s]
    Da_T = 2.14e-5  # Molecular diffusivity of "heat" in air at STP [m2/s]
    va = 1.51e-5  # air viscosity at STP [m2/s]
    g = 9.81  # gravitational constant [m/s2]

    # -- Adjust diffusivity, viscosity, and air density to pressure/temp.
    t_adj = (101300.0 / P)*((Ta + 273.15) / 293.16)**1.75
    Da_v = Da_v*t_adj
    Da_c = Da_c*t_adj
    Da_T = Da_T*t_adj
    va = va*t_adj
    rho_air = 44.6*(P / 101300.0)*(273.15 / (Ta + 273.13))  # [mol/m3]

    # ----- Compute the leaf-level dimensionless groups
    Re = u*d / va  # Reynolds number
    Sc_v = va / Da_v  # Schmid numbers for water
    Sc_c = va / Da_c  # Schmid numbers for CO2
    Pr = va / Da_T  # Prandtl number
    Gr = g*(d**3)*abs(dT) / (Ta + 273.15) / (va**2)  # Grashoff number

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

    r = Gr / (Re**2)  # ratio of free/forced convection

    return gb_h, gb_c, gb_v, r


""" ------- photosynthesis models ------- """

def photo_c3_analytical(photop, Qp, T, VPD, ca, gb_c, gb_v):
    """
    Leaf photosynthesis and gas-exchange by co-limitation optimality model of
    Vico et al. 2013 AFM
    IN:
        photop - parameter dict with keys: Vcmax, Jmax, Rd, alpha, theta, La, tresp
           can be scalars or arrays.
           tresp - dictionary with keys: Vcmax, Jmax, Rd: temperature sensitivity
           parameters. OMIT key if no temperature adjustments for photoparameters.
        Qp - incident PAR at leaves (umolm-2s-1)
        Tleaf - leaf temperature (degC)
        VPD - leaf-air vapor pressure difference(mol/mol)
        ca - ambient CO2 (ppm)
        gb_c - boundary-layer conductance for co2 (mol m-2 s-1)
        gb_v - boundary-layer conductance for h2o (mol m-2 s-1)
    OUT:
        An - net CO2 flux (umolm-2s-1)
        Rd - dark respiration (umolm-2s-1)
        fe - leaf transpiration rate (molm-2s-1)
        gs - stomatal conductance for CO2 (mol/m-2s-1)
        ci - leaf internal CO2 (ppm)
        cs - leaf surface CO2 (ppm)
    """
    Tk = T + 273.15

    # print('T', Tk, 'D', VPD, 'Ca', ca, 'gbc', gb_c, 'qp', Qp)
    # --- const & params---
    MaxIter = 20
    a = 1.6  # H2O to CO2 diffusion rate
    Coa = 2.10e5  # O2 in air (umol/mol)
    TN = 298.15  # reference temperature 298.15 K = 25degC
    R = 8.314427  # gas constant, J mol-1 K-1

    # --- params ----
    Vcmax = photop['Vcmax']
    Jmax = photop['Jmax']
    Rd = photop['Rd']
    alpha = photop['alpha']
    theta = photop['theta']
    La = photop['La']
    g0 = photop['g0']

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

    # --- model parameters k1_c, k2_c [umol/m2/s]
    Km = Kc*(1.0 + Coa / Ko)
    J = (Jmax + alpha*Qp -((Jmax + alpha*Qp)**2.0 - (4*theta*Jmax*alpha*Qp))**(0.5)) / (2.0*theta)

    k1_c = J/4.0
    k2_c = (J/4.0) * Km / Vcmax

    # --- iterative solution for cs
    err = 9999.0
    cnt = 1
    cs = ca  # leaf surface CO2
    while err > 0.01 and cnt < MaxIter:
        NUM1 = -k1_c * (k2_c - (cs - 2*Tau_c))
        DEN1 = (k2_c + cs)**2
        NUM2 = np.sqrt(a*VPD*La*k1_c**2 * (cs - Tau_c)*(k2_c + Tau_c) * ((k2_c + (cs - 2*a*VPD*La))**2)*(k2_c + (cs - a*VPD*La)))
        DEN2 = a*VPD*La*((k2_c + cs)**2) * (k2_c + (cs - a*VPD*La))

        gs_opt = (NUM1 / DEN1) + (NUM2 / DEN2) + eps

        ci = (1. / (2 *gs_opt)) * (-k1_c - k2_c*gs_opt + cs*gs_opt + Rd + np.sqrt((k1_c + k2_c*gs_opt - cs*gs_opt - Rd)**2 \
            - 4*gs_opt*(-k1_c*Tau_c - k2_c*cs*gs_opt - k2_c*Rd)))

        An = gs_opt*(cs - ci)
        An1 = np.maximum(An, 0.0)
        cs0 = cs
        cs = ca - An1 / gb_c

        err = np.nanmax(abs(cs - cs0))
        cnt = cnt + 1
        # print('err', err)


    ix = np.where(An < 0)
    gs_opt[ix] = g0

    if type(ca) is float:
        ci[ix] = ca
        cs[ix] = ca
    else:
        ci[ix] = ca[ix]
        cs[ix] = ca[ix]

    gs_v = a*gs_opt

    geff = (gb_v*gs_v) / (gb_v + gs_v)  # molm-2s-1
    fe = geff*VPD  # leaf transpiration rate

    if len(An) == 1:
        return float(An), float(Rd), float(fe), float(gs_opt), float(ci), float(cs)
    else:
        return An, Rd, fe, gs_opt, ci, cs


def photo_c3_medlyn(photop, Qp, T, VPD, ca, gb_c, gb_v, P=101300.0):
    """
    Leaf gas-exchange by Farquhar-Medlyn model, where co-limitation as in
    Vico et al. 2013 AFM
    IN:
        photop - parameter dict with keys: Vcmax, Jmax, Rd, alpha, theta, La, tresp
           can be scalars or arrays.
           tresp - dictionary with keys: Vcmax, Jmax, Rd: temperature sensitivity
           parameters. OMIT key if no temperature adjustments for photoparameters.
        Qp - incident PAR at leaves (umolm-2s-1)
        Tleaf - leaf temperature (degC)
        VPD - leaf-air vapor pressure difference (mol/mol)
        ca - ambient CO2 (ppm)
        gb_c - boundary-layer conductance for co2 (mol m-2 s-1)
        gb_v - boundary-layer conductance for h2o (mol m-2 s-1)
        P - atm. pressure (Pa)
    OUT:
        An - net CO2 flux (umolm-2s-1)
        Rd - dark respiration (umolm-2s-1)
        fe - leaf transpiration rate (molm-2s-1)
        gs - stomatal conductance for CO2 (mol/m-2s-1)
        ci - leaf internal CO2 (ppm)
        cs - leaf surface CO2 (ppm)
    """
    Tk = T + 273.15
    VPD = 1e-3 * VPD * P  # kPa

    # --- const & params---
    MaxIter = 50
    a = 1.6  # H2O to CO2 diffusion rate
    Coa = 2.10e5  # O2 in air (umol/mol)
    TN = 298.15  # reference temperature 298.15 K = 25degC
    R = 8.314427  # gas constant, J mol-1 K-1

    # --- params ----
    Vcmax = photop['Vcmax']
    Jmax = photop['Jmax']
    Rd = photop['Rd']
    alpha = photop['alpha']
    theta = photop['theta']
    m = photop['m']  # slope parameter
    g0 = photop['g0']

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

    # --- model parameters k1_c, k2_c [umol/m2/s]
    Km = Kc*(1.0 + Coa / Ko)
    J = (Jmax + alpha*Qp -((Jmax + alpha*Qp)**2.0 - (4*theta*Jmax*alpha*Qp))**(0.5)) / (2*theta)
    k1_c = J / 4.0
    k2_c = J / 4.0 * Km / Vcmax

    # --- iterative solution for cs and ci
    err = 9999.0
    cnt = 1
    cs = ca  # leaf surface CO2
    ci = 0.8*ca  # internal CO2
    while err > 0.01 and cnt < MaxIter:
        # CO2 demand (Vico eq. 1) & gs_opt (Medlyn eq. xx)
        An = k1_c * (ci - Tau_c) / (k2_c + ci) - Rd  # umolm-2s-1
        An1 = np.maximum(An, 0.0)
        gs_opt = (1.0 + m / (VPD**0.5)) * An1 / (cs - Tau_c)  # mol m-2s-1
        gs_opt = np.maximum(g0, gs_opt)  # g0 is the lower limit

        # CO2 supply
        cs = np.maximum(ca - An1 / gb_c, 0.5*ca)  # through boundary layer
        ci0 = ci
        ci = np.maximum(cs - An1 / gs_opt, 0.5*ca)  # through stomata

        err = max(abs(ci0 - ci))
        cnt += 1
    # when Rd > photo, assume stomata closed and ci == ca
    ix = np.where(An < 0)
    if type(ca) is float:
        ci[ix] = ca
        cs[ix] = ca
    else:
        ci[ix] = ca[ix]
        cs[ix] = ca[ix]
    gs_opt[ix] = g0
    gs_v = a*gs_opt

    geff = (gb_v*gs_v) / (gb_v + gs_v)  # molm-2s-1
    fe = geff * VPD / (1e-3 * P)  # leaf transpiration rate

    return An, Rd, fe, gs_opt, ci, cs


def photo_c3_medlyn_farquhar(photop, Qp, T, VPD, ca, gb_c, gb_v, P=101300.0):
    """
    Leaf gas-exchange by Farquhar-Medlyn model, where co-limitation as in standard Farquhar-
    model
    IN:
        photop - parameter dict with keys: Vcmax, Jmax, Rd, alpha, theta, La, tresp
           can be scalars or arrays.
           tresp - dictionary with keys: Vcmax, Jmax, Rd: temperature sensitivity
           parameters. OMIT key if no temperature adjustments for photoparameters.
        Qp - incident PAR at leaves (umolm-2s-1)
        Tleaf - leaf temperature (degC)
        VPD - leaf-air vapor pressure difference (mol/mol)
        ca - ambient CO2 (ppm)
        gb_c - boundary-layer conductance for co2 (mol m-2 s-1)
        gb_v - boundary-layer conductance for h2o (mol m-2 s-1)
        P - atm. pressure (Pa)
    OUT:
        An - net CO2 flux (umolm-2s-1)
        Rd - dark respiration (umolm-2s-1)
        fe - leaf transpiration rate (molm-2s-1)
        gs - stomatal conductance for CO2 (mol/m-2s-1)
        ci - leaf internal CO2 (ppm)
        cs - leaf surface CO2 (ppm)
    """
    Tk = T + 273.15
    VPD = 1e-3 * VPD * P  # kPa

    # --- const & params---
    MaxIter = 50
    a = 1.6  # H2O to CO2 diffusion rate
    Coa = 2.10e5  # O2 in air (umol/mol)
    TN = 298.15  # reference temperature 298.15 K = 25degC
    R = 8.314427  # gas constant, J mol-1 K-1

    # --- params ----
    Vcmax = photop['Vcmax']
    Jmax = photop['Jmax']
    Rd = photop['Rd']
    alpha = photop['alpha']
    theta = photop['theta']
    m = photop['m']  # slope parameter
    g0 = photop['g0']
    beta = photop['beta']
    #print beta
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

    # --- model parameters k1_c, k2_c [umol/m2/s]
    Km = Kc*(1.0 + Coa / Ko)
    J = (Jmax + alpha*Qp -((Jmax + alpha*Qp)**2.0 - (4*theta*Jmax*alpha*Qp))**(0.5)) / (2*theta)
    #k1_c = J / 4.0
    #k2_c = J / 4.0 * Km / Vcmax

    # --- iterative solution for cs and ci
    err = 9999.0
    cnt = 1
    cs = ca  # leaf surface CO2
    ci = 0.8*ca  # internal CO2
    while err > 0.01 and cnt < MaxIter:
        # -- rubisco -limited rate
        Av = Vcmax * (ci - Tau_c) / (ci + Km)
        # -- RuBP -regeneration limited rate
        Aj = J/4.0 * (ci - Tau_c) / (ci + 2.0*Tau_c)

        #An = np.minimum(Av, Aj) - Rd  # single limiting rate
        x = Av + Aj
        y = Av * Aj
        An = (x - (x**2.0 - 4.0*beta*y)**0.5) / (2.0*beta) - Rd  # co-limitation

        An1 = np.maximum(An, 0.0)
        #print An1
        # stomatal conductance
        gs_opt = g0 + (1.0 + m / (VPD**0.5)) * An1 / cs
        gs_opt = np.maximum(g0, gs_opt)  # gcut is the lower limit
        #print gs_opt
        # CO2 supply
        cs = np.maximum(ca - An1 / gb_c, 0.5*ca)  # through boundary layer
        ci0 = ci
        ci = np.maximum(cs - An1 / gs_opt, 0.1*ca)  # through stomata

        err = max(abs(ci0 - ci))
        cnt += 1

    # when Rd > photo, assume stomata closed and ci == ca
    ix = np.where(An < 0)
    if type(ca) is float:
        ci[ix] = ca
        cs[ix] = ca
    else:
        ci[ix] = ca[ix]
        cs[ix] = ca[ix]
    gs_opt[ix] = g0
    gs_v = a*gs_opt

    geff = (gb_v*gs_v) / (gb_v + gs_v)  # molm-2s-1
    fe = geff * VPD / (1e-3 * P)  # leaf transpiration rate

    return An, Rd, fe, gs_opt, ci, cs


def photo_c3_bwb(photop, Qp, T, RH, ca, gb_c, gb_v, P=101300.0):
    """
    Leaf gas-exchange by Farquhar-Ball-Woodrow-Berry model, where co-limitation as in
    Vico et al. 2013 AFM
    IN:
        photop - parameter dict with keys: Vcmax, Jmax, Rd, alpha, theta, La, tresp
           can be scalars or arrays.
           tresp - dictionary with keys: Vcmax, Jmax, Rd: temperature sensitivity
           parameters. OMIT key if no temperature adjustments for photoparameters.
        Qp - incident PAR at leaves (umolm-2s-1)
        Tleaf - leaf temperature (degC)
        rh - relative humidity at leaf temperature (-)
        ca - ambient CO2 (ppm)
        gb_c - boundary-layer conductance for co2 (mol m-2 s-1)
        gb_v - boundary-layer conductance for h2o (mol m-2 s-1)
        P - atm. pressure (Pa)
    OUT:
        An - net CO2 flux (umolm-2s-1)
        Rd - dark respiration (umolm-2s-1)
        fe - leaf transpiration rate (molm-2s-1)
        gs - stomatal conductance for CO2 (mol/m-2s-1)
        ci - leaf internal CO2 (ppm)
        cs - leaf surface CO2 (ppm)
    """
    Tk = T + 273.15
    # --- const & params---
    MaxIter = 50
    a = 1.6  # H2O to CO2 diffusion rate
    Coa = 2.10e5  # O2 in air (umol/mol)
    TN = 298.15  # reference temperature 298.15 K = 25degC
    R = 8.314427  # gas constant, J mol-1 K-1

    # --- params ----
    Vcmax = photop['Vcmax']
    Jmax = photop['Jmax']
    Rd = photop['Rd']
    alpha = photop['alpha']
    theta = photop['theta']
    m = photop['m']  # slope parameter
    g0 = photop['g0']

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

    # --- model parameters k1_c, k2_c [umol/m2/s]
    Km = Kc*(1.0 + Coa / Ko)
    J = (Jmax + alpha*Qp -((Jmax + alpha*Qp)**2.0 - (4*theta*Jmax*alpha*Qp))**(0.5)) / (2*theta)
    k1_c = J / 4.0
    k2_c = J / 4.0 * Km / Vcmax

    # --- iterative solution for cs
    err = 9999.0
    cnt = 1
    cs = ca  # leaf surface CO2
    ci = 0.8*ca  # internal CO2
    while err > 0.01 and cnt < MaxIter:
        # CO2 demand (Vico eq. 1) & gs_opt (Medlyn eq. xx)
        An = k1_c * (ci - Tau_c) / (k2_c + ci) - Rd  # umolm-2s-1
        An1 = np.maximum(An, 0.0)
        # bwb -scheme
        gs_opt = g0 + m * An1 / ((cs - Tau_c))*RH
        gs_opt = np.maximum(g0, gs_opt)  # gcut is the lower limit

        # CO2 supply
        cs = np.maximum(ca - An1 / gb_c, 0.5*ca)  # through boundary layer
        ci0 = ci
        ci = np.maximum(cs - An1 / gs_opt, 0.1*ca)  # through stomata

        err = max(abs(ci0 - ci))
        cnt += 1
        
    # when Rd > photo, assume stomata closed and ci == ca
    ix = np.where(An < 0)
    gs_opt[ix] = g0
    ci[ix] = ca[ix]
    cs[ix] = ca[ix]
    gs_v = a*gs_opt

    geff = (gb_v*gs_v) / (gb_v + gs_v)  # molm-2s-1
    esat, _ = saturation_vapor_pressure(T)
    VPD = (1.0 - RH) * esat / P  # mol mol-1
    fe = geff*VPD  # leaf transpiration rate

    return An, Rd, fe, gs_opt, ci, cs


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
        co_limi - True uses co-limitation function of Vico et al., 2014.
    OUTPUT:
        An - leaf net CO2 exchange (umol m-2 leaf s-1)
        Rd - leaf dark respiration rate (umol m-2 leaf s-1)
    NOTE: original and co_limi -versions converge when beta ~ 0.8
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

    if not co_limi:
        # -- rubisco -limited rate
        Av = Vcmax * (ci - Tau_c) / (ci + Km)
        # -- RuBP -regeneration limited rate
        Aj = J/4 * (ci - Tau_c) / (ci + 2.0*Tau_c)

        # An = np.minimum(Av, Aj) - Rd  # single limiting rate
        x = Av + Aj
        y = Av * Aj
        An = (x - (x**2 - 4*beta*y)**0.5) / (2*beta) - Rd  # co-limitation
        return An, Rd, Av, Aj
    else:   # use Vico et al. eq. 1
        k1_c = J / 4.0
        k2_c = (J / 4.0) * Km / Vcmax

        An = k1_c * (ci - Tau_c) / (k2_c + ci) - Rd
        return An, Rd, Tau_c, Kc, Ko, Km, J


def photo_temperature_response(Vcmax0, Jmax0, Rd0, Vcmax_T, Jmax_T, Rd_T, T):
    """
    Adjusts Farquhar / co-limitation optimality model parameters for temperature
    INPUT:
        Vcmax0, Jmax0, Rd0 - parameters at ref. temperature 298.15 K
        Vcmax_T, Jmax_T, Rd_T - temperature response parameter lists
        T - leaf temperature (K)
    OUTPUT: Nx1-arrays
        Vcmax, Jmax,Rd (umol m-2(leaf) s-1)
        Gamma_star - CO2 compensation point
    CALLED from Farquhar(); Opti_C3_Analytical(); Opti_C3_Numerical()
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

    if var_in.lower() == 'deltas':  # Sv --> Topt
        xout = Hd / (xin - R * np.log(Ha / (Hd - Ha)))
    else:  # Topt -->Sv
        c = R * np.log(Ha / (Hd - Ha))
        xout = (Hd + xin * c) / xin
    return xout


def photo_Toptima(T10):
    """
    computes acclimation of temperature optima of Vcmax and Jmax to 10-day mean air temperature
    Args:
        T10 - 10-day mean temperature (degC)
    Returns:
        Tv, Tj - temperature optima of Vcmax, Jmax
        rjv - ratio of Jmax25 / Vcmax25
    Reference: Lombardozzi et al., 2015 GRL, eq. 3 & 4
    """
    # --- parameters
    R = 8.314427  # universal gas constant, J mol-1 K-1
    Hav = 72000.0  # J mol-1
    Haj = 50000.0  # J mol-1
    Hd = 200000.0  # J mol.1

    T10 = np.minimum(40.0, np.maximum(10.0, T10))  # range 10...40 degC
    # vcmax T-optima
    dSv = 668.39 - 1.07*T10  # J mol-1
    Tv = Hd / (dSv - R * np.log(Hav / (Hd - Hav))) - 273.15  # degC
    # jmax T-optima
    dSj = 659.70 - 0.75*T10  # J mol-1
    Tj = Hd / (dSj - R * np.log(Haj / (Hd - Haj))) - 273.15  # degC

    rjv = 2.59 - 0.035*T10  # Jmax25 / Vcmax25

    return Tv, Tj, rjv