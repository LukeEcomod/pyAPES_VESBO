# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 08:16:37 2016


@author: slauniai
"""

import numpy as np
import matplotlib.pyplot as plt
eps = np.finfo(float).eps  # machine epsilon

"""
General functions
"""


def wrc(pF, x=None, var=None):
    """
    vanGenuchten-Mualem soil water retention model (van Genuchten, 1980;
    Schaap and van Genuchten, 2006)

    .. math::
        \\theta(\\psi_s) = \\theta_{res} + \\frac{\\theta_{sat}-\\theta_{res}}
        {(1 + \\lvert \\alpha + \\psi_{s}\\rvert^n)^m}

    where :math:`\\theta_{res}` and :math:`\\theta_{sat}` are residual and saturation
    water contents (m\ :sup:`3` m :sup:`-3`\ ), :math:`\\alpha`\ , *n*, and :math:`m=1-^1/_n`
    are empirical shape parameters.

    Sole input 'pF' draws water retention curve and returns 'None'.
    For drawing give only one pF-parameter set. If several pF-curves are given,
    x can be scalar or len(x)=len(pF). In former case var is pF(x), in latter var[i]=pf[i,x[i]]

    References:
        Schaap and van Genuchten (2005). Vadose Zone 5:27-34
        van Genuchten, (1980). Soil Science Society of America Journal 44:892-898

    Args:
        pF (list/dict):
            0. 'ThetaS' saturated water content [m\ :sup:`3` m :sup:`-3`\ ]
            1. 'ThetaR' residual water content [m\ :sup:`3` m :sup:`-3`\ ]
            2. 'alpha' air entry suction [cm\ :sup:`-1`]
            3. 'n' pore size distribution [-]
        x:
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is volumetric water content
            * [m] if input is water potential
        var: flag for conversion
            * 'Th' for volumetric water content
            * None for water potential
    Returns:
        numpy array float:
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is water potential
            * [m] if input is volumetric water content
    Samuli Launiainen, Luke 2/2016
    """
    if type(pF) is dict:  # dict input
        # Ts, Tr, alfa, n = pF['ThetaS'], pF['ThetaR'], pF['alpha'], pF['n']
        Ts = np.array(pF['ThetaS'])
        Tr = np.array(pF['ThetaR'])
        alfa = np.array(pF['alpha'])
        n = np.array(pF['n'])
        m = 1.0 - np.divide(1.0, n)

    else:  # list input
        pF = np.array(pF, ndmin=1)  # ndmin=1 needed for indexing 0-d arrays
        Ts = pF[:, 0]
        Tr = pF[:, 1]
        alfa = pF[:, 2]
        n = pF[:, 3]
        m = 1.0 - np.divide(1.0, n)

    def theta_psi(x):
        # converts water content (m3m-3) to potential (m)
        x = np.minimum(x, Ts)
        x = np.maximum(x, Tr)  # checks limits
        s = (Ts - Tr) / ((x - Tr) + eps)
        Psi = -1e-2 / alfa*(s**(1.0 / m) - 1.0)**(1.0 / n)  # m
        Psi[np.isnan(Psi)] = 0.0
        return Psi

    def psi_theta(x):
        # converts water potential (m) to water content (m3m-3)
        x = 100*np.minimum(x, 0)  # cm
        Th = Tr + (Ts - Tr) / (1 + abs(alfa*x)**n)**m
        return Th

    # This does all the work
    if x is None and np.size(Ts) == 1:  # draws pf-curve
        xx = -np.logspace(-4, 5, 100)  # cm
        yy = psi_theta(xx)
        #  field capacity and wilting point
        fc = psi_theta(-1.0)
        wp = psi_theta(-150.0)

        fig = plt.figure()
        fig.suptitle('vanGenuchten-Mualem WRC', fontsize=16)
        ttext = r'$\theta_s=$' + str(Ts) + r', $\theta_r=$' + str(Tr) +\
                r', $\alpha=$' + str(alfa) + ',n=' + str(n)

        plt.title(ttext, fontsize=14)
        plt.semilogx(-xx, yy, 'g-')
        plt.semilogx(1, fc, 'ro', 150, wp, 'ro')  # fc, wp
        plt.text(1, 1.1*fc, 'FC'), plt.text(150, 1.2*wp, 'WP')
        plt.ylabel(r'$\theta$  $(m^3m^{-3})$', fontsize=14)
        plt.xlabel('$\psi$ $(m)$', fontsize=14)
        plt.ylim(0.8*Tr, min(1, 1.1*Ts))

        del xx, yy
        return None

    elif x is None:
        print 'soil_cores.wrc: To draw curve give only one pF -parameter set'
        return None

    if var is 'Th':
        y = theta_psi(x)  # 'Theta-->Psi'
    else:
        y = psi_theta(x)  # 'Psi-->Theta'

    return y

def h_to_cellmoist(pF,h,dzu,dzl):
    
    dzu = dzu / 2
    dzl = dzl / 2
    dzu[0] = dzu[0] * 2
    dzl[-1] = dzl[-1] * 2
    
    # Ts, Tr, alfa, n = pF['ThetaS'], pF['ThetaR'], pF['alpha'], pF['n']
    Ts = np.array(pF['ThetaS'])
    Tr = np.array(pF['ThetaR'])
    alfa = np.array(pF['alpha'])
    n = np.array(pF['n'])
    m = 1.0 - np.divide(1.0, n)
    
    # Theta based on cell center h
    x = np.minimum(h, 0)
    theta = Tr + (Ts - Tr) / (1 + abs(alfa*100*x)**n)**m    
    
    # Correct moisture of partly unsatrurated cells
    ix1 = np.where(h < dzu )
    ix2 = np.where(h > -dzl)
    ix = np.intersect1d(ix1, ix2)
    del ix1, ix2
    # moisture of unsaturated part
    x[ix] = (dzu[ix] - h[ix]) / 2
    theta[ix] = Tr[ix] + (Ts[ix] - Tr[ix]) / (1 + abs(alfa[ix]*100*x[ix])**n[ix])**m[ix]
    # total moisture
    theta[ix] = (theta[ix] * (dzu[ix] - h[ix]) + Ts[ix] * (dzl[ix] + h[ix]))/(dzu[ix] + dzl[ix])
    
    return theta

def diff_wcapa(pF,h,dzu,dzl):

    dh = 1e-5
    theta_plus = h_to_cellmoist(pF,h + dh,dzu,dzl)
    theta_minus = h_to_cellmoist(pF,h - dh,dzu,dzl)

    dwcapa = (theta_plus - theta_minus) / (2 * dh)
    
    return dwcapa

def effSat(pF, x, var=None):
    """
    Effective saturation of a soil layer
    Args:
        pF (list/dict):
            0. 'ThetaS' saturated water content [m\ :sup:`3` m :sup:`-3`\ ]
            1. 'ThetaR' residual water content [m\ :sup:`3` m :sup:`-3`\ ]
            2. 'alpha' air entry suction [cm\ :sup:`-1`]
            3. 'n' pore size distribution [-]
        x:
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is volumetric water content
            * [m] if input is water potential
        var: flag for conversion
            * 'Th' for volumetric water content
            * None for water potential
    Returns:
        np.array:
            * effective saturation [-] s = (x - ThetaR)/(ThetaS - ThetaR)
    """
    if type(pF) is dict:  # dict input
        Ts = np.array(pF['ThetaS'])
        Tr = np.array(pF['ThetaR'])
        alfa = np.array(pF['alpha'])
        n = np.array(pF['n'])
    else:  # list input
        pF = np.array(pF, ndmin=1)  # ndmin=1 for indexing of 0-dim arrays
        Ts = pF[0]
        Tr = pF[1]
        alfa = pF[2]
        n = pF[3]

    if var is not None or var is 'Th':  # x=Th
        s = np.minimum((x - Tr) / (Ts - Tr + eps), 1.0)
    else:
        #x = 100*np.minimum(x, 0)  # cm                                         #Turha toistaa?
        #x = Tr + (Ts - Tr) / (1.0 + abs(alfa*x)**n)**(1.0 - 1.0 / n)           #Turha toistaa?
        s = np.minimum((wrc(pF, x=x) - Tr) / (Ts - Tr + eps), 1.0)

    return s


def hydraulic_conductivity(pF, x=None, var=None, Ksat=1):
    """
    Hydraulic conductivity following vanGenuchten-Mualem -model. Sole input 'pF'
    draws relative conductivity curve.
    Args:
        pF (dict/list):
            0. 'ThetaS' saturated water content [m\ :sup:`3` m :sup:`-3`\ ]
            1. 'ThetaR' residual water content [m\ :sup:`3` m :sup:`-3`\ ]
            2. 'alpha' air entry suction [cm\ :sup:`-1`]
            3. 'n' pore size distribution [-]
        x (float or np.array):
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is vol.water content
            * [m] if input is water potential
        var (str): flag for conversion
            * 'Th' for volumetric water content
            * None for water potential
        Ksat (float or np.array): saturated hydraulic conductivity [units]
    Returns:
        Kh - hydraulic conductivity (if Ksat ~=1 then in [units], else relative [-])
    """

    if type(pF) is dict:  # dict input
        alfa = np.array(pF['alpha'])
        n = np.array(pF['n'])
        m = 1.0 - np.divide(1.0, n)
    else:  # list input
        pF = np.array(pF, ndmin=1)
        alfa = pF[:, 2]
        n = pF[:, 3]
        m = 1.0 - np.divide(1.0, n)

    def relcond(x):
        nm = (1.0 - abs(alfa*x)**(n - 1.0) * (1 + abs(alfa*x)**n)**(-m))**2
        dn = (1.0 + abs(alfa*x)**n)**(m / 2.0)
        r = nm / (dn + eps)
        return r

    if x is None and np.size(alfa) == 1:  # draws curve
        xx = -np.logspace(-4, 5, 100)  # cm
        yy = relcond(xx)

        fig = plt.figure()
        fig.suptitle('Hydr. cond. (vanGenuchten-Mualem)', fontsize=16)
        ttext = r'$K_{sat}=$' + str(Ksat) + r', $\alpha=$' + str(alfa) \
                + ', n=' + str(n)

        plt.title(ttext, fontsize=14)
        plt.semilogx(-xx, yy, 'g-')
        plt.ylabel(r'K_{h} / K_{sat}', fontsize=14)
        plt.xlabel('$\psi$ $(cm)$', fontsize=14)

        del xx, yy
        return None

    elif x is None:
        print 'hydrCond: To draw curve give only one pF -parameter set'
        return None

    # this computes and returns
    x = np.array(x)
    if x is not None and var is 'Th':
        x = wrc(pF, x=x, var='Th')

    Kh = Ksat*relcond(100.0*np.minimum(x, 0))                                   # x includes heads > 0

    return Kh


def rew(x, fc, wp):
    """
    Relative plant extractable water at the prevailing vol. water content
    Args:
        float, list or np.array:
            x: vol. water content [vol/vol]
            fc: field capacity [vol/vol]
            wp: wilting point [vol/vol]
    Returns:
        rew (np.array) [-]
    """
    rew = np.minimum((x - wp) / (fc - wp + eps), 1.0)
    rew = np.maximum(rew, 0.0)
    return rew


def moisture_convert(x, bd, var=None):
    """
    Converts between volumetric [m3m-3] and gravimetric [kg kg-1] water content.
    Args :
        x (float, np.array):
            * volumetric [vol/vol] or
            * gravimetric water cont. on dry-mass basis [H2O/dryweight]
        bd (float, np.array): soil bulk density (dry density)
            [kg/m3 = 0.001 g/cm3]. Thus 1 g/cm3=1000kg/m3
        var (str): 'Th' if x is vol. water content
    Returns:
        (np.array):
        converted water content [kg kg-1] or [m3m-3]
    """
    dw = 1000.0  # water density [kg/m3]
    if var is 'Th':
        y = np.array(x)*np.divide(dw, bd)  # return grav. water content
    else:
        y = np.array(x)*bd / dw  # return vol. water content
    return y

"""
Water flow in soils in 1D

"""


def waterFlow1D(t_final, z, h0, pF, Ksat, Prec, Evap, R, HM=0.0, lbc={'type': 'impermeable', 'value': None},
                Wice0=0.0, maxPond=0.0, pond0=0.0, cosalfa=1.0, h_atm=-1000.0, steps=10):
    """
    Solves soil water flow in 1-D using implicit, backward finite difference solution of Richard's
    equation.

    Args:
        t_final (float): solution timestep [s]
        z (array): grid,<0 (soil surface = 0), monotonically decreasing [m] 
                    (all the nodes, also top and bottom node, are in the centre of the soil compartments)
        h0 (array): initial hydraulic head [m]
        pF (dict): dict of vanGenuchten soil pF-parameters; pF.keys()=['ThetaR', 'ThetaS', 'n', 'alpha']
        Ksat (array): saturated hydr. cond. [ms-1]
        Prec (float): potential top boundary flux, >0 [m s-1]                   # Eritelty sateeksi ja haihdunnaksi koska vain toinen näistä potentiaalista
        Evap (float): potential evaporation from surface, >0 [m s-1]
        R (array): local sink/source array due root uptake & well [s-1], >0 for sink
        HM (array): net lateral flux array, e.g. to ditches [s-1], >0 for net outflow
        lbc: lower bc {'type': 'impermeable', 'flux', 'free_drain', 'head', and 'value'; give for 'head' and 'flux'}, head [m], flux [m s-1] <0 for outflow
        Wice0 (array): vol. ice content [m3m-3] - not needed now; could be used to scale hydr.conductivity
        maxPond (float): maximum allowed pond depth at surface [m]
        pond0 (float): initial pond depth [m]
        cosalfa - 1 for vertical water flow, 0 for horizontal transport
        h_atm - hydraulic head [m] in equilibrium with air humidity - used to compute soil evaporation supply  # Ei vakio?
        steps (int): initial number of subtimesteps used to proceed to 't_final'
    Returns:
        h (array): new hydraulic head [m]
        W (array): new total water content [m3m-3]
        h_pond (float): new ponding depth [m]
        C_inf (float): total infiltration [m], <=0
        C_eva (float): total evaporation [m],>=0
        C_dra (float): total drainage (to ditches and through bottom) [m], >=0 from profile
        C_trans (float): total root uptake [m], >=0 from profile
        C_roff (float): total surface runoff [m]
        Fliq (array): vertical water fluxes [ms-1] at t_final
        gwl (float): ground water level [m]; if not in computational layer then assumed hydrostatic equilibrium with node -1
        Kv (array): hydraulic conductivity [ms-1]
        mbe (float): total mass balance error [m]
    REFERENCES:
        vanDam & Feddes (2000): Numerical simulation of infiltration, evaporation and shallow
        groundwater levels with the Richards equation, J.Hydrol 233, 72-85.
    CODE:
        Samuli Launiainen, Luke 8.4.2016. Converted from Matlab (APES SoilProfile.WaterFlow)
        Kersti Haahti, 29.12.2017: 
            - Upper bc switching between head and flux as in vanDam & Feddes (2000) 
            - small fixes (dz, dzu, dzl, KLh (spatial_average, hydraulic_conductivity: saturated cells = Ksat), get_gwl, output, inputs)
            - adjustments to selecting dt: minimally 30s, iterNo can reach 20, still neads checking..
            - lbc options checked, adjustments of heads in get_gwl removed
    NOTE:
        (8.4.2016): upper bc restriction checks needs to be tested
        (   -"-  ): include macropore adjustment as in APES-code?
        (29.12.2017): What is Fliq used for? Calculation migth need checking..
    """

    Conv_crit = 1.0e-5
    S = R + HM  # net imposed sink/source term (root uptake + lateral flow, e.g drainage)

    # -------------Get computation grid -------------
    N = len(z)  # nr of nodal points, 0 is top

    dz = np.empty(N)
    dzu = np.empty(N)
    dzl = np.empty(N)

    # distances between grid points: dzu is between point i-1 and i, dzl between point i and i+1
    dzu[1:] = z[:-1] - z[1:]
    dzu[0] = -z[0]  # from soil surface to first node, soil surface af z = 0
    dzl[:-1] = z[:-1] - z[1:]
    dzl[-1] = (z[-2] - z[-1]) / 2.0 #  from last node to bottom surface
    # compartment thickness
    dz = (dzu + dzl) / 2.0
    dz[0] = dzu[0] + dzl[0] / 2.0
    dz[-1] = dzu[-1] / 2.0 + dzl[-1] 
    
    # ----soil variables and save intial conditions--
    if type(Ksat) is float:
        Ksat = np.zeros(N) + Ksat

    poros = pF['ThetaS']
    W_ini = h_to_cellmoist(pF,h0,dzu,dzl) #wrc(pF, x=h0)  # m3m-3
    # h_ini = h0
    pond_ini = pond0

    # these change during solution
    W = W_ini.copy()
    h = h0.copy()
    h_pond = pond0

    # -------- find solution at t_final --------

    dto = t_final / steps  # initial time step [s]
    t = 0.0  # running time
    dt = dto
    
    # cumulative boundary fluxes
    C_inf = 0.0
    C_eva = 0.0
    C_dra = 0.0
    C_trans = 0.0
    C_roff = 0.0
    # Rsink = sum(S*dz)

    while t < t_final:  # loop until solution timestep
        # these will stay constant during iteration over time step "dt"
        h_old = h.copy()
        W_old = W.copy()

        # these change during iteration
        h_iter = h.copy()
        W_iter = W.copy()
        afp_iter = poros - W_iter  # air filled prosity [m3m-3]
        
        KLh = hydraulic_conductivity(pF, x=h_iter, Ksat=Ksat)  # [m s-1]
        KLh = spatial_average(KLh, method='arithmetic')                     # Returns K at i-1/2

        err1 = 999.0
        err2 = 999.0
        iterNo = 0

        Conv_crit2 = 1.0

        # start iterative solution of Richards equation
        pass_flag = True
        while (err1 > Conv_crit or err2 > Conv_crit2):  # and pass_flag is False:
            # print pass_flag
            iterNo += 1

            # ------ lower bc check -------
            if lbc['type'] == 'free_drain':
                q_bot = -KLh[-1]*cosalfa
                # if h_iter(N)<-1, q_bot=0; end % allow free drainage until field capacity
            elif lbc['type'] == 'impermeable':
                q_bot = 0.0
            elif lbc['type'] == 'flux':
                q_bot = max(lbc['value'], -KLh[-1]*cosalfa)

            elif lbc['type'] == 'head':
                h_bot = lbc['value']
                q_bot = -KLh[-1]*(h_iter[-1] - h_bot) / dzl[-1] - KLh[-1]*cosalfa # note, this is approximation for output

            # ------- upper bc check ---------
            # using swiching between flux and head bc as in Dam and Feddes (2000)
            # note: q0<0 infiltration, q0>0 evaporation
            q0 = Evap - Prec - h_pond / dt  # potential rate m/s
            MaxInf = max(-KLh[0]*(h_pond - h_iter[0] - z[0]) / dzu[0], -Ksat[0])  # max infiltration rate to top node
            MaxEva = -KLh[0]*(h_atm - h_iter[0] - z[0]) / dzu[0]  # limited by soil supply [m/s]
            Qin = (q_bot - sum(S*dz) - q0) * dt  # net flow to column
            
            Airvol = max(0.0, sum((poros - W_old)*dz))  # maximum inflow (m) tolerated by soil column
            if q0 < 0:  # case infiltration
                if Airvol <= eps:  # initially saturated profile
                    if Qin >= 0: #and Airvol / sum(poros*dz) < 1.0e-3:  # inflow exceeds outflow, matrix stays saturated
                        print 'saturated soil column, ponding water, h_sur = ' + str(min(Qin, maxPond)) + ' h = ' + str(h_iter[0])  + ' W = ' + str(W_old[0])
                        h_sur = min(Qin, maxPond)
                        ubc_flag = 'head'

                    else:  # outflow exceeds inflow
                        print 'outflow exceeds inflow' + ' q0 = ' + str(q0) + ' Qin = ' + str(Qin) + ' h_pond = ' + str(h_pond) 
                        q_sur = q0
                        ubc_flag = 'flux'
                        if iterNo ==1:
                            h_iter -= dz
                            W_iter = h_to_cellmoist(pF,h_iter,dzu,dzl)

                else:  # initially unsaturated profile
                    if Qin >= Airvol:  # only part fits into profile
                        print 'only part fits into profile, h_sur = ' + str(min(Qin - Airvol, maxPond)) + ' h = ' + str(h_iter[0])  + ' W = ' + str(W_old[0])
                        h_sur = min(Qin - Airvol, maxPond)
                        ubc_flag = 'head'
                        
                    else:  # all fits into profile
                        print 'all fits into profile, q0 = ' + str(q0) + ' Airvol = ' + str(Airvol)  + ' h_pond = ' + str(h_pond)+ ' MaxInf = ' + str(MaxInf)
                        if iterNo ==1 and Airvol < 1e-3:
                            h_iter -= dz
                            W_iter = h_to_cellmoist(pF,h_iter,dzu,dzl)
                        if q0 < MaxInf:
                            h_sur = h_pond
                            ubc_flag = 'head'
                        else:
                            q_sur = q0
                            ubc_flag = 'flux'

            else:  # case evaporation
                if iterNo ==1 and Airvol < 1e-3:
                    h_iter -= dz
                    W_iter = h_to_cellmoist(pF,h_iter,dzu,dzl)
                if q0 > MaxEva:
                    print 'case evaporation, limited by atm demand, q0 = ' + str(q0) + ' MaxEva = ' + str(MaxEva)  + ' h_pond = ' + str(h_pond)
                    h_sur = h_atm
                    ubc_flag = 'head'
                else:
                    print 'case evaporation, no limit ' + 'h_pond = ' + str(h_pond)+ ' Qin = ' + str(Qin)+ ' Airvol= ' + str(Airvol) + ' q0 = ' + str(q0) + ' dt= ' + str(dt)
                    q_sur = q0
                    ubc_flag = 'flux'

                    
                    
            if Qin >= Airvol:
                Conv_crit2 = 1e-8
            else:
                Conv_crit2 = 1e-6

            # ----------------------------------------------
            C = diff_wcapa(pF,h_iter,dzu,dzl)  #diff_capa(pF, h_iter)  # differential water capacity

            # set up tridiagonal matrix
            a = np.zeros(N)  # sub diagonal
            b = np.zeros(N)  # diagonal
            g = np.zeros(N)  # super diag
            f = np.zeros(N)  # rhs

            # ------ intermediate nodes for j in range(1,N-1):
            b[1:-1] = C[1:-1] + dt / dz[1:-1] * (KLh[1:N-1] / dzu[1:-1] + KLh[2:N] / dzl[1:-1])
            a[1:-1] = - dt / (dz[1:-1] * dzu[1:-1]) * KLh[1:N-1]
            g[1:-1] = - dt / (dz[1:-1] * dzl[1:-1]) * KLh[2:N]

            f[1:-1] = C[1:-1] * h_iter[1:-1] - (W_iter[1:-1] - W_old[1:-1]) + dt / dz[1:-1]\
                      * (KLh[1:N-1] - KLh[2:N]) * cosalfa - S[1:-1] * dt

            # ---------- top node (j=0) ---------
            if ubc_flag != 'head':  # flux bc
                b[0] = C[0] + dt / (dz[0] * dzl[0]) * KLh[1]
                a[0] = 0.0
                g[0] = -dt / (dz[0] * dzl[0]) * KLh[1]

                f[0] = C[0] * h_iter[0] - (W_iter[0] - W_old[0]) + dt / dz[0] * (-q_sur - KLh[1] * cosalfa) - S[0] * dt
            
            else:  # head boundary
                b[0] = C[0] + dt / dz[0] *(KLh[0] / dzu[0] + KLh[1] / dzl[0])
                a[0] = 0;
                g[0] = -dt / (dz[0] * dzl[0]) * KLh[1];
                f[0] = C[0] * h_iter[0] - (W_iter[0] - W_old[0])  + dt / dz[0] * ((KLh[0] - KLh[1]) * cosalfa + KLh[0] / dzu[0] * h_sur) - S[0]*dt;

            # ------ bottom node (j=N) ---------
            if lbc['type'] != 'head':  # flux bc
                b[-1] = C[-1] + dt / (dz[-1] * dzu[-1]) * KLh[N-1]
                a[-1] = -dt / (dz[-1] * dzu[-1]) * KLh[N-1]  # -dt/(dz[-1]*dzu[-1])*KLh[-2]
                g[-1] = 0.0

                f[-1] = C[-1] * h_iter[-1] - (W_iter[-1] - W_old[-1])\
                        + dt / dz[-1] * (KLh[N-1] * cosalfa + q_bot) - S[-1] * dt

            else:  # head boundary, fixed head "at node N+1"
                b[-1] = C[-1] + dt / dz[-1] *(KLh[N-1] / dzu[-1] +  KLh[N] / dzl[-1])
                a[-1] = - dt / (dz[-1] * dzu[-1]) * KLh[N-1]
                g[-1] = 0.0

                f[-1] =  C[-1] * h_iter[-1] - (W_iter[-1] - W_old[-1]) + dt / dz[-1] * ((KLh[N-1] - KLh[N]) * cosalfa + KLh[N] / dzl[-1] * h_bot) - S[-1] * dt

            # ---------save old and update iteration values
            h_iterold = h_iter.copy()
            W_iterold = W_iter.copy()

            h_iter = thomas(a,b,g,f)
            # h_iter[h_iter>0]=0.0
            # print h_iter
            # find GWL
            _, gwl = get_gwl(h_iter, z) 
            W_iter = h_to_cellmoist(pF,h_iter,dzu,dzl) #wrc(pF, x=h_iter)
            afp_iter = poros - W_iter

            if any(abs(h_iter - h_iterold) > 1.0):
                print 'break'
                print (a)
                print (b)
                print (g)
                print (f)
                print (C)
                print (h)
                print (h_iter)
                break
            if iterNo == 20:
                dt = dt / 3.0
                if dt < 10:
                    dt = max(dt,30)
                    iterNo = 0
                    print 'More than 20 iterations, new dt = ' + str(dt)
                    continue  # re-try with smaller timestep
                else:
                    break
                    print 'Problem with solution'
            elif pass_flag is False or any(np.isnan(h_iter)):
                #print pass_flag,  ', dt: ' + str(dt) + ' qsur = ' + str(q_sur) + ' qbot = ' + str(q_bot)
                if any(np.isnan(h_iter)):
                    print 'nan found'
                    break  # break while loop

            err1 = max(abs(W_iter - W_iterold))
            err2 = max(abs(h_iter - h_iterold))
            print 'err1 = ' + str(err1) + ' err2 = ' + str(err2)

        # ------- ending iteration loop -------
        # print 'Pass: t = ' + str(t) + ' dt = ' + str(dt) + ' iterNo = ' + str(iterNo)
        # new state at t
        h = h_iter.copy()
        W = h_to_cellmoist(pF,h_iter,dzu,dzl)  #wrc(pF, x=h_iter)
        # calculate q_sur and q_bot in case of head boundaries
        if ubc_flag == 'head':
            q_sur = -KLh[0] * (h_sur - h[0]) / dzu[0] - KLh[0]
        if lbc['type'] == 'head':
            q_bot = -KLh[-1]*(h[-1] - h_bot) / dzl[-1] - KLh[-1]*cosalfa
            
        # ---- update cumulative Infiltration, evaporation, drainage and h_pond [m] ----
        if q_sur <= eps:        
            h_pond -= (-Prec + Evap - q_sur)*dt
            #h_pond = max(0, h_pond)                                                # Voiko jättää???
            rr = max(0, h_pond - maxPond)  # create runoff if h_pond>maxpond
            h_pond -= rr
            C_roff += rr        
            del rr
            C_inf += (q_sur - Evap)*dt
            C_eva += Evap*dt
        else:
            C_eva += (q_sur + Prec)*dt + h_pond
            h_pond = 0.0
        
#        if abs(h_pond) < eps:
#            h_pond = 0.0

        C_dra += -q_bot*dt + sum(HM*dz)*dt  # drainage + net lateral flow
        C_trans += sum(R*dz)*dt

        t += dt  # solution time & new initial timestep
        
        print 't = ' + str(t) + ' dt = ' + str(dt) + ' iterNo = ' + str(iterNo)
        
        dt_old = dt
        
        if iterNo <= 2:
            dt = dt * 1.25
        elif iterNo > 4:
            dt = dt / 1.25
        
        dt = max(dt,360)#30)

        if dt_old == t_final or t_final > t:
            dto = min(dt, t_final)
        
        dt = min(dt, t_final-t)
        

    # ----- at t = t_final: ending while loop & compute outputs
                
    # vertical fluxes
    KLh = hydraulic_conductivity(pF, x=h, Ksat=Ksat)
    #KLh = spatial_average(KLh, method='arithmetic')

    Fliq = nodal_fluxes(z, h, KLh)  # ms-1

    # mass balance error [m]
 #   mbe = Prec*t_final + (h_to_colmoist(pF,h0,dzu,dzl) - h_to_colmoist(pF,h,dzu,dzl)) + (pond_ini - h_pond) - C_dra - C_trans - C_roff - C_eva
    mbe = Prec*t_final + (sum(W_ini*dz) - sum(W*dz)) + (pond_ini - h_pond) - C_dra - C_trans - C_roff - C_eva
    # print 'mbe:' +str(mbe)

    return h, W, h_pond, C_inf, C_eva, C_dra, C_trans, C_roff, Fliq, gwl, KLh, mbe, dto


""" Utility functions """

def get_gwl(head, x):
    """
    finds ground water level and adjusts head in saturated nodes by adding hydrostatic pressure
    Args:
        head (array): heads in nodes [m]
        x (array): grid,<0, monotonically decreasing [m]                                      
    Returns:
        head (array): adjusted heads [m]
        gwl (float): ground water level in column [m]
    """
    sid = find_index(head, lambda x: x >= 0)  # returns indices of heads >=0
    if len(sid) > 0:
        if sid[0] > 0:  # gwl below first node
            gwl = x[sid[0]-1] + head[sid[0]-1]
        else:
            gwl = x[sid[0]] + head[sid[0]]  # gwl in first node, limited to ground surface (at z = 0)

        #head[sid] = gwl - x[sid]  # m, >0                                      # Messes up calculation when using other tah impermeable lower boundary 
    else:
        gwl = head[-1] + x[-1]  # if gwl not in profile, assume hydr.
                                # equilibrium between last node and gwl

    return head, gwl


def thomas(a, b, C, D):
    """
    Tridiagonal matrix algorithm of Thomas
    a=subdiag, b=diag, C=superdiag, D=rhs
    """
    n = len(a)
    V = np.zeros(n)
    G = np.zeros(n)
    U = np.zeros(n)
    x = np.zeros(n)

    V[0] = b[0].copy()
    G[0] = C[0] / V[0]
    U[0] = D[0] / V[0]

    for i in range(1, n):  # nr of nodes
        V[i] = b[i] - a[i]*G[i - 1]
        U[i] = (D[i] - a[i]*U[i - 1]) / V[i]
        G[i] = C[i] / V[i]

    x[-1] = U[-1]
    inn = n - 2
    for i in range(inn, -1, -1):
        x[i] = U[i] - G[i] * x[i + 1]
    return x

def diff_capa(pF, head):
    """
    Derivative of vGenuchten soil water retention curve [m-1]
    Args: 
        pF (dict): dict of vanGenuchten soil pF-parameters; pF.keys()=['ThetaR', 'ThetaS', 'n', 'alpha']
        head (array): head [m]
    Returns: 
        x (array): dW/dhead, derivative of vGenuchten soil water retention curve for given heads [m-1]
    """

    h = -100.0*head  # cm
    ts = pF['ThetaS']
    tr = pF['ThetaR']
    n = pF['n']
    m = 1.0 - np.divide(1, n)
    alfa = pF['alpha']

    # print ts, tr, n, m, alfa                                                  # for x <= 0 "RuntimeWarning: invalid value encountered in power"
    x = 100.0*(ts - tr)*(n - 1.0)*alfa**n*h**(n - 1.0) / ( (1.0 + (alfa*h)**n)**(m + 1.0))
    x[h <= 0.0] = 0.0
    return x


def find_index(a, func):
    """
    finds indexes or array elements that fill the condition
    call as find_index(a, lambda x: criteria)
    """
    return [i for (i, val) in enumerate(a) if func(val)]

def spatial_average(y, x=None, method='arithmetic'):
    """
    Calculates spatial average of quantity y, from node points to soil compartment edges
    Args: 
        y (array): quantity to average
        x (array): grid,<0, monotonically decreasing [m]
        method (str): flag for method 'arithmetic', 'geometric','dist_weighted'
    Returns: 
        f (array): averaged y, note len(f) = len(y) + 1
    """

    N = len(y)
    f = np.empty(N+1)  # Between all nodes and at surface and bottom
    if method is 'arithmetic':
        f[1:-1] = 0.5*(y[:-1] + y[1:])
        f[0] = y[0]
        f[-1] = y[-1]

    elif method is 'geometric':
        f[1:-1] = np.sqrt(y[:-1] * y[1:])
        f[0] = y[0]
        f[-1] = y[-1]

    elif method is 'dist_weighted':                                             # En ymmärrä, ei taida olla käyttössä
        a = (x[0:-2] - x[2:])*y[:-2]*y[1:-1]
        b = y[1:-1]*(x[:-2] - x[1:-1]) + y[:-2]*(x[1:-1] - x[2:])

        f[1:-1] = a / b
        f[0] = y[0]
        f[-1] = y[-1]

    return f


def nodal_fluxes(x, y, K):
    """
    Calculates fluxes between nodal points in 1-D grid, f = -K*dh/dx            # now at nodal points?
    Args: 
        x (array): grid, monotonic!                                             # Elevation of nodes [m]?
        y (array): heads in nodes [m]
        K (array): hydraulic conductivity at nodes [ms-1]
    Returns: 
        f (array): flux [ms-1]
    """
    f = np.empty(np.shape(y))
    yprim = central_difference(x, y)                                            # yprim = (y_left-y_rigth)/2dx

    f = -K*yprim  # flux
    return f


def central_difference(x, y):
    """
    Derivative by central difference
    Args: 
        x (array): grid,  monotonic?
        y (array): values
    Returns: 
        yprim (array): derivative
    """
    yprim = np.empty(np.shape(y))

    yprim[1:-1] = (y[2:] - y[0:-2]) / (x[2:] - x[0:-2])  # central difference
    yprim[0] = (y[1] - y[0]) / (x[1] - x[0])  # forward difference at left boundary
    yprim[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])  # backward difference at right boundary

    return yprim


""" 1D transient water flow solved by Crank-Nicholson scheme, following FEMMA -code (Koivusalo, Lauren et al.) """
""" TÄTÄ PITÄÄ ITEROIDA. MASSA EI SÄILY, ALARAJAN REUNAEHTO EI TOIMI EIKÄ NUMERIIKKA PELAA."""

def waterFlow1D_CN(dt0, z, h0, pF, Ksat, Ftop, R, F=0.0, lbcType='impermeable',lbcValue=None, Wice=0.0, maxPond=0.0, pond0=0.0, h_atm=-1000.0, Implic=0.5, steps=10, MaxIter=100, IterLim=5.0e-5):
    """
    Solves soil water flow in 1-D using Crank-Nicholson (predictor-corrector) finite difference solution of Richard's equation.
    Reference: Koivusalo (2009) FEMMA-document:
    IN:
        t_final - solution timestep [s]
        z - grid,<0, monotonically decreasing
        h0 - initial hydraulic head [m]
        pF - dict of vanGenuchten soil pF-parameters; pF.keys()=['ThetaR', 'ThetaS', 'n', 'alpha']
        Ksat - saturated hydr. cond. [ms-1]
        Ftop - potential top boundary flux (<0 for infiltration)
        R - local sink/source array due root uptake & well [s-1], <0 for sink
        HM - net lateral flux array [s-1], <0 for net outflow
        lbcType - lower bc type: 'impermeable', 'flux', 'free_drain', 'head'
        lbcValue - lower bc value; give for 'head' and 'flux'
        Wice0 - vol. ice content [m3m-3] - not needed now; could be used to scale hydr.conductivity
        maxPond - maximum allowed pond depth at surface [m]
        pond0 - initial pond depth [m]
        cosalfa - 1 for vertical water flow, 0 for horizontal transport
        h_atm - hydraulic head [m] in equilibrium with air humidity - used to compute soil evaporation supply
        steps - initial subtimesteps used to proceed to 't_final'
    OUT:
        h - new hydraulic head [m]
        W - new total water content [m3m-3]
        h_pond - new ponding depth [m]
        C_inf - total infiltration [m], <=0
        C_eva - total evaporation [m],>=0
        C_dra - total drainage from profile [m], <0 from profile
        C_roff - total surface runoff [m]
        Fliq - vertical water fluxes [ms-1] at t_final
        gwl - ground water level [m]; if not in computational layer then assumed hydrostatic equilibrium with node -1
        mbe - total mass balance error [m]
    CODE:
        Samuli Launiainen, Luke 8.4.2016. Converted from Matlab (APES SoilProfile.WaterFlow)
    NOTE:
        (8.4.2016): upper bc restriction checks needs to be tested
        (   -"-  ): include macropore adjustment as in APES-code?

    prepare sink term S -> (Eact(lyr) + (HorizLatFlowOut(lyr)) + HorizDeepLatFlow(lyr) _
        + HorizDrFlow(lyr)) / Dz(1)
    prepare Qinf = Precip + SurfSto / Dt

    """
    from soil_core import hydrCond
    from soil_core import wrc

    #lowerBC's
    if lbcType is 'impermeable': Qbot=0.0
    if lbcType is 'flux': Qbot=lbcValue
    if lbcType is 'free_drain': Qbot=None #value setlater
    if lbcType is 'head': h_bot=lbcValue

    S=R+F #net imposed sink/source term (root uptake + lateral flow)

    #-------------Get computation grid -------------
    N=len(z) #nr of nodal points, 0 is top
    dz=np.empty(N); dzu=np.empty(N); dzl=np.empty(N)
    #if any(z)>0: z=-z; #must be monotonic negative values

    #distances between grid points: dzu is between point i-1 and i, dzl between point i and i+1
    dzu[1:]=z[0:-1] - z[1:N]; dzu[0]=-z[0]
    dzl[0:-1]=z[0:-1] - z[1:]; dzl[-1]=(z[-2] - z[-1])/2.0;

    dz=(dzu + dzl)/2.0;
    dz[0]=dzu[0] + dzl[0]/2.0;
    print 'dz = ', dz
    print 'dzu = ', dzu
    print 'dzl = ', dzl

    #----soil variables and save intial conditions--
    if type(Ksat) is float: Ksat=np.zeros(N)+Ksat
    Ws=np.array(pF['ThetaS'])
    S=np.array(R)+np.array(F); del R, F #lumped sink/source term
    W_ini=wrc(pF,x=h0); #m3m-3

    #these change during process
    h_new=h0.copy()
    h_pond=pond0;


    #cumulative boundary fluxes
    C_inf=0.0;
    C_eva=0.0;
    C_dra=0.0;
    C_roff=0.0;

    dt=dt0/steps # internal timestep [s]
    for idt in range(1,steps): #internal time loop

        print 'idt: ' +str(idt)
        h_old = h_new.copy()
        h_iter = h_new.copy()
        iterNo=0

        for ic in range(MaxIter+1): #iteration loop
            iterNo +=1;

            KLh = hydrCond(pF, x = h_iter, Ksat=Ksat) #hydr conductivity m/s
            KLh=spatialAverage(KLh,method='arithmetic')
            C=diffCapa(pF, h_iter)    #differential water capacity

            #check upper boundary condition
            #print airVol
            if Ftop==0.0: Qtop=0.0;
            if Ftop<0: #case infiltration
                airVol = sum((Ws - wrc(pF, x=h_old))*dz)  #air filled porosity
                MaxInf=-Ksat[0]/dz[0]*(Implic *h_iter + (1.0-Implic)*h_old)

                Qtop =max(-airVol / dt, Ftop, MaxInf)
                #print 'infiltr'
            if Ftop>=0: #case evaporation
                MaxEva=-KLh[0]/dz[0]*(h_atm - h_old[0] - 1) #maximum evaporation by soil supply
                Qtop=max(Ftop,MaxEva)
                #print 'evap'
            #print Qtop


            #set up tridiagonal matrix
            a= np.zeros(N); b = np.zeros(N); g = np.zeros(N); f = np.zeros(N)
            #mid layers
            a[1:-1] = -Implic * KLh[0:-2] / ( dz[1:-1] * dzu[1:-1]) #subdiag
            b[1:-1] = C[1:-1] / dt + Implic * ( KLh[0:-2] / ( dz[1:-1] * dzu[1:-1]) + KLh[1:-1] / ( dz[1:-1] * dzl[1:-1]) )   #diag
            g[1:-1] = - Implic *KLh[1:-1] / (dz[1:-1]*dzl[1:-1]) #superdiag

            f[1:-1] = C[1:-1]*h_old[1:-1] / dt - (1-Implic) *( KLh[0:-2]/dz[1:-1]* ( (h_old[1:-1]-h_old[0:-2]) /dzu[1:-1] - 1.0) +\
                    KLh[1:-1]/dz[1:-1]* ( (h_old[2:]-h_old[1:-1]) /dzl[1:-1] - 1.0) ) - S[1:-1] #RHS

            #bottom bc
            if lbcType=='free_drain': Qbot=-KLh[-1];

            #flux boundary (impermeable: Qbot=0, prescribed flux: Qbot=lbcValue, free drainage: Qbot=k[-1])
            if lbcType is not 'head':
                #print dzz
                a[-1] = -Implic * KLh[-1]/ (dz[-1] *dzu[-1])
                b[-1] = C[-1]/dt + Implic*KLh[-1] / (dz[-1] *dzu[-1])
                g[-1] = 0.0
                f[-1] = C[-1]*h_old[-1] /dt - (1- Implic) * KLh[-1]/dz[-1]* ( (h_old[-1]-h_old[-2])/dzu[-1] - 1.0 ) + Qbot/dz[-1] - S[-1]

            else:   #fixed head lbcValue
                a[-1] = -Implic * KLh[-1] / ( dz[-2] * dzu[-1]) #subdiag
                b[-1] = C[-1] / dt + Implic * ( KLh[-2] / ( dz[-1] * dzu[-1]) + KLh[-1] / ( dz[-1] * dzl[-1]) )   #diag
                g[-1] = - Implic *KLh[-1] / (dz[-1]*dzl[-1]) #superdiag

                f[-1] = C[-1]*h_old[-1] / dt - (1-Implic) *( KLh[-2]/dz[-1]* ( (h_old[-1]-h_old[-2]) /dzu[-1] - 1.0) +\
                        KLh[-1]/dz[-1]* ( (h_bot-h_old[-1]) /dzl[-1] - 1.0) ) - S[-1] #RHS

            # top bc is flux-based
            a[0] = 0.0
            b[0] = C[0]/dt + Implic*KLh[1]/(dz[0]*dzl[0])
            g[0] = -Implic*KLh[1]/(dz[0]*dzl[0])

            f[0] = C[0]*h_old[0]/dt + (1 - Implic)* ( KLh[1]/dz[0]*( ( h_old[1] - h_old[0])/dzl[0] -1))-S[0] -Qtop/dz[0]

            #call tridiagonal solver
            h_iterold = h_iter.copy()
            h_iter=thomas(a,b,g,f);
            h_iter, gwl=getGwl(h_iter, z)

            err=max(abs(h_iter - h_iterold))
            #print err
            if err < IterLim:
                print iterNo
                print h_iter
                break
            #reset matrix coefficients
            #a.fill(np.NaN); C.fill(np.NaN); b.fill(np.NaN); D.fill(np.NaN)
    #update state variable and cumulative fluxes
    #print 'ic :' + str(ic)

        h_new=h_iter.copy()
        W_new=wrc(pF, x=h_new)
        C_inf += min(0.0, Qtop)*dt
        C_eva += max(0.0, Qtop)*dt
        C_dra += Qbot*dt

        if Qtop<=0:
            h_pond=max(0, h_pond - (Ftop - Qtop)*dt)
            rr=max(0, h_pond - maxPond) #create runoff if h_pond>maxpond
            h_pond=h_pond - rr;
            C_roff += rr; del rr

       # if ic==MaxIter: psiiNew=psiiOld0
    mbe=None

    return h_new, W_new, h_pond, C_inf, C_eva, C_dra, C_roff, mbe

def solveRichardsSteady(z, gwl, pF, Ksv, Ksh=None, h_root=-50.0, figSwitch=False):
    """
    Computes steady-state solution of Richards equation between ground water level and bottom of root zone.
    IN:
        z - grid [m], <0
        gwl - depth of ground water level [m]
        pF - dict of vanGenuchten pF parameters or list in order [ThetaS, ThetaR, alpha, n]
        Ksv - vertical sat. hydr. conductity [ms-1]
        Ksh - horizontal sat. hydr. conductity [ms-1]. None if Ksh = Ksv
        h_root - suction at upper boundary (root zone bottom) [m]
        figSwitch - True plots figures
    OUT:
        X - hydraulic head profile [m]
        UpFlux - capillary upflux to root zone [ms-1]
    """

    Omega = 0.5
    Conv_crit = 0.00001
    maxIter = 500

    # -------- Get computation grid -------------
    N = len(z)  # nr of nodal points, 0 is top
    dz = np.empty(N)
    dzu = np.empty(N)
    dzl = np.empty(N)

    # distances between grid points: dzu is between point i-1 and i, dzl between point i and i+1
    dzu[1:] = z[0:-1] - z[1:N]
    dzu[0] = -z[0]
    dzl[0:-1] = z[0:-1] - z[1:]
    dzl[-1] = dzl[-2]

    dz = (dzu + dzl) / 2.0
    dz[0] = dzu[0] + dzl[0] / 2.0
    #    print 'z = ', z
    #    print 'dz = ', dz
    #    print 'dzu = ', dzu
    #    print 'dzl = ', dzl

    if type(Ksv) is float:
        Ksv = np.ones(N)*Ksv
    if type(Ksh) is float:
        Ksh = np.ones(N)*Ksh
    elif Ksh is None:
        Ksh = Ksv

    # tridiagonal elements
    A = np.zeros(N)
    B = np.zeros(N)
    C = np.zeros(N)
    D = np.zeros(N)

    X = np.zeros(N)
    X_old = np.zeros(N)

    # initial condition
    X[z <= gwl] = gwl - z[-1]
    X[z > gwl] = 0.5*h_root

#    dd = [z[0], gwl,z[-1]]
#    yy = [h_root,  0.0, 0.0]
#    iPsii= interp1d(dd, yy); del yy, dd
#    X=iPsii(z).copy()
    # print 'X=', X

    # X_ini=X.copy() #save for possible smaller Omega

    # bc at root zone bottom
    A[0] = 0.0
    B[0] = 1.0
    C[0] = 0.0
    D[0] = h_root

    # bc at gwl
    A[-1] = 0.0
    B[-1] = 1.0
    C[-1] = 0.0
    D[-1] = max(0.0, gwl - z[-1])

    err = 9999.99
    iterNo = 0
    while err > Conv_crit and iterNo < maxIter:
        iterNo += 1
        Con = hydraulic_conductivity(pF, x=X, Ksat=Ksv)  # hydr conductivity m/s in this iteration
        Con = spatial_average(Con, 'arithmetic')
        X_old = X.copy()

        # tridiag elements, middle nodes
        A[1:-1] = Con[1:-1] / dzu[1]
        C[1:-1] = Con[2:] / dzl[1]
        B[1:-1] = - C[1:-1] - A[1:-1]
        D[1:-1] = - (Con[1:-1] - Con[2:])

        Xi = thomas(A, B, C, D)

        X = Omega*Xi + (1.0 - Omega)*X_old

        err = max(abs(X - X_old))
        # print err
        if iterNo == maxIter:
            Omega = 0.1
            iterNo = 0
            # X=X_ini.copy() #reset to initial condition

    # ------ convergence, compute upflux profile
    # print 'iterNo=' + str(iterNo) +',Omega=' + str(Omega)

    flx = np.zeros(N)
    for k in range(0, N-1):
        flx[k] = - Con[k]*((X[k] - X[k+1]) / dz[k])
    flx[X >= 0.0] = 0.0
    UpFlux = flx[1]  # assume that flux at node i=1 equals capillary rise to root zone

#    flux using central difference
#    xx=X.copy(); xx[xx>0.0]=0.0
#    flx=nodalFluxes(z,xx,Con) #compute nodal fluxes
#    flx[X>0]=0.0
#    UpFlux=flx[1] #assume that flux at node i=1 equals capillary rise to

    if figSwitch is True:
        plt.figure()
        plt.subplot(121); plt.plot(X,z,'r.-'); plt.ylabel('z'); plt.xlabel('h [m]')
        plt.subplot(122); plt.plot(flx,z,'r.-'); plt.ylabel('z'); plt.xlabel('Flux [ms-1]')
        print 'upflux = ' + str(UpFlux)

    return X, UpFlux

"""
************ Drainage equations ********************

"""


def drainage_linear(zs, Ksat, GWL, DitchDepth, DitchSpacing):
    """"
    Calculates drainage from soil profile to ditch using simple linear equation,
    i.e. accounts only drainage from layers where GWL<DitchDepth
    INPUT:
       zs - depth of soil node (m), array, zs<0
       Ksat - saturated hydraulic conductivity (m/s),array
       GWL - ground water level below surface (m), GWL<0
       DitchDepth (m), depth of drainage ditch bottom
       DitchSpacing (m), horizontal spacing of drainage ditches
    OUTPUT:
       Q - total drainage from soil profile (m/s), >0 is outflow
       Qz_drain - drainage from each soil layer (m m-1s-1), i.e. sink term to Richard's eq.
    """

    dx = 1.0  # unit length of horizontal element (m)
    N = len(zs)

    Qz_drain = np.zeros(N)
    dz = np.zeros(N)

    dz[1:] = (zs[0:-1] - zs[1:])  # m
    dz[0] = -2*zs[0]

    Keff = Ksat * dz / dx  # transmissivity m s-1 in each layer

    # (-), positive gradient means flow towards ditches, return flow neglected
    hgrad = np.max([(GWL + DitchDepth) / (0.5*DitchSpacing), 0.0])

    ix1 = np.where((zs - GWL) < 0)
    ix2 = np.where(zs > -DitchDepth)  # layers above ditch bottom where drainage is possible
    ix = np.intersect1d(ix1, ix2)
    del ix1, ix2

    Qz_drain[ix] = Keff[ix]*hgrad  # layerwise drainage ms-1
    Q = sum(Qz_drain)
    Qz_drain = Qz_drain / dz  # s-1, sink term

    return Q, Qz_drain


def drainage_hooghoud(zs, Ksat, GWL, DitchDepth, DitchSpacing, DitchWidth, Zbot=None):
    """
    Calculates drainage to ditch using Hooghoud's drainage equation,
    i.e. accounts only drainage from saturated layers above and below ditch bottom
    Args:
       zs (array): grid,<0, monotonically decreasing [m]
       Ksat (array): saturated hydr. cond. [ms-1]
       GWL (float): ground water level below surface, <0 [m]
       DitchDepth (float): depth of drainage ditch bottom, >0 [m]
       DitchSpacing (float): horizontal spacing of drainage ditches [m]
       DitchWidth (float): ditch bottom width [m]
       Zbot (float): distance to impermeable layer, >0 [m]
    Returns:
       Q (float): total drainage from soil profile, >0 is outflow [ms-1]
       Qz_drain (array): drainage from each soil layer, i.e. sink term to Richard's eq. [m m-1s-1]
    Reference:
       Follows Koivusalo, Lauren et al. FEMMA -document. Ref: El-Sadek et al., 2001.
       J. Irrig.& Drainage Engineering.

    Samuli Launiainen, Metla 3.11.2014.; converted to Python 14.9.2016
    Kersti Haahti, 29.12.2017. Code checked, small corrections
    """
    
    N = len(zs)
    Qz_drain = np.zeros(N)
    Qa = 0.0
    Qb = 0.0

    dz = np.empty(N)  # Depth of soil compartment
    dz[1:-1] = (zs[0:-2] - zs[1:-1]) / 2 + (zs[1:-1] - zs[2:]) / 2 # m
    dz[0] = -2*zs[0]
    dz[-1] = zs[-2] - zs[-1]
    
    if Zbot is None or Zbot > sum(dz):  # Can't be lower than soil profile bottom
        Zbot = sum(dz)

    Hdr = np.maximum(0, GWL + DitchDepth)  # depth of saturated layer above ditch bottom
    Hdr = min(Hdr, DitchDepth)
    
    # print Hdr
    if Hdr > 0:
        Trans = Ksat[:]*dz  # transmissivity of layer, m2s-1

        # -------- drainage from saturated layers above ditch base
        
        # layers above ditch bottom where drainage is possible
        ix = np.intersect1d(np.where((zs - GWL) < 0), np.where(zs > -DitchDepth))  

        if ix.size > 0:
            Ka = sum(Trans[ix]) / sum(dz[ix])  # effective hydraulic conductivity ms-1
            Qa = 4 * Ka * Hdr**2 / (DitchSpacing**2)  # m s-1, total drainage above ditches
            # sink term s-1, rel=Trans(dra)/sum(Trans(dra)) partitions Qa by relative transmissivity of layer
            Qz_drain[ix] = Qa*Trans[ix] / sum(Trans[ix]) / dz[ix]
            del ix

        # -----drainage from saturated layers below ditch base

        # layers below ditch bottom where drainage is possible
        ix = np.intersect1d(np.where(zs <= -DitchDepth), np.where((zs - GWL) < 0))

        Kb = sum(Trans[ix]) / sum(dz[ix])  # effective hydraulic conductivity ms-1

        # compute equivalent depth Deq
        Dbt = Zbot - DitchDepth  # distance from impermeable layer to ditch bottom
        A = 3.55 - 1.6 * Dbt / DitchSpacing + 2 * (2.0 / DitchSpacing)**2.0
        Reff = DitchWidth / 2.0  # effective radius of ditch
        

        if Dbt / DitchSpacing <= 0.3:
            Deq = Dbt / (1.0 + Dbt / DitchSpacing * (8 / np.pi * np.log(Dbt / Reff) - A))  # m
        else:
            Deq = np.pi * DitchSpacing / (8 * (np.log(DitchSpacing / Reff) - 1.15))  # m

        Qb = 8 * Kb * Deq * Hdr / DitchSpacing**2  # m s-1, total drainage below ditches
        Qz_drain[ix] = Qb * Trans[ix] / sum(Trans[ix]) / dz[ix]  # sink term s-1
        del ix

    Q = Qa + Qb  # total drainage m s-1, positive is outflow to ditch
    # print Hdr, Qa, Qb, Q
        # plt.figure(1)
    # plt.subplot(121); plt.plot([0 1],[GWL GWL],'r-',[0 1],[-DitchDepth -DitchDepth],'k-')
    # plt.subplot(122); plt.plot(Qz_drain,zs,'r-'); ylabel('zs'); xlabel('Drainage profile s-1')
    return Q, Qz_drain


""" ************** Functions for testing stuff in this module **************"""


def test_drainage(DitchDepth, DitchSpacing):
    """ tests drainage equations"""

    DitchWidth = 0.8
    Zbot = 1.0

    zs = -np.arange(0.1, 2.0, 0.01)  # grid
    Ksat = 1e-5  # m/s

    gwl = np.arange(-DitchDepth, max(zs), 0.01)
    N = len(gwl)
    Drain = np.zeros([N, 1])
    DrainHoog = np.zeros([N, 1])
    Qprofile = []
    for k in range(0, N):
        Q, Qz = drainage_linear(zs, Ksat, gwl[k], DitchDepth, DitchSpacing)
        Drain[k] = Q
        Qprofile.append(Qz)
        del Q, Qz
        Q, Qz = drainage_hooghoud(zs, Ksat, gwl[k], DitchDepth, DitchSpacing, DitchWidth, Zbot)
        DrainHoog[k] = Q

    plt.figure()
    plt.subplot(121); plt.plot(gwl,Drain,'r.-', gwl, DrainHoog, 'g.-'); plt.title('Q ms-1')
    plt.subplot(122); plt.plot(gwl,DrainHoog,'c.-')

    return Drain, Qprofile
