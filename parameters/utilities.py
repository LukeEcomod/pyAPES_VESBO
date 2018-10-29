# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:05:05 2018

@author: L1656
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
eps = np.finfo(float).eps  # machine epsilon

import seaborn as sns

def fit_pF(head, watcont, fig=False, labels=None, percentage=False, kPa=False):
    """

    """

    if fig:
        plt.figure()
        colors = sns.color_palette("hls", len(watcont))
        c = 0
    head = np.array(head)
    if kPa:
        head = head * 10  # kPa -> cm
    vg_ini = (0.88, 0.09, 0.03, 1.3)
    bounds = ((eps, eps, eps, 1.0), (1.0, 1.0, 5.0, 5.0))
    van_g = lambda h, *p:   p[1] + (p[0] - p[1]) / (1. + (p[2] * h) **p[3]) **(1. - 1. / p[3])
    vgen_all = []

    for k in range(0, len(watcont)):
        Wcont = np.array(watcont[k])
        ix = np.where(Wcont >= 0)
        if percentage:
            Wcont[ix] = Wcont[ix] / 100  # % -> fraction
        if labels is None:
            label = ''
        else:
            label = labels[k] + ': '
        try:
            vgen, _ = curve_fit(van_g, head[ix], Wcont[ix], p0=vg_ini, bounds=bounds)
            label+='Ts=%5.3f, Tr=%5.3f, alfa=%5.3f, n=%5.3f' % tuple(vgen)
        except RuntimeError:
            vgen = [-1, -1, -1, -1]
            label+='No fit!'
        vgen_all.append(vgen)

        if fig:

            plt.semilogy(Wcont[ix], head[ix], '.',color = colors[c])
            xx = np.logspace(-1, 5.0, 100)
            plt.semilogy(van_g(xx, *vgen), xx, '-',color = colors[c],
                         label=label)
            c += 1

    if fig:
        plt.xlabel(r'$\theta$  $(m^3m^{-3})$', fontsize=14)
        plt.ylabel('$-h$ $(cm)$', fontsize=14)
        plt.ylim(xx[0], xx[-1])
        plt.xlim(0.0, 1.0)
        plt.legend()

    return vgen_all

def peat_hydrol_properties(x, unit='g/cm3', var='bd', ptype='A', fig=False, labels=None):
    """
    Peat water retention and saturated hydraulic conductivity as a function of bulk density
    Päivänen 1973. Hydraulic conductivity and water retention in peat soils. Acta forestalia fennica 129.
    see bulk density: page 48, fig 19; degree of humification: page 51 fig 21
    Hydraulic conductivity (cm/s) as a function of bulk density(g/cm3), page 18, as a function of degree of humification see page 51 
    input:
        - x peat inputvariable in: db, bulk density or dgree of humification (von Post)  as array \n
        - bulk density unit 'g/cm3' or 'kg/m3' \n
        - var 'db' if input variable is as bulk density, 'H' if as degree of humification (von Post) \n
        - ptype peat type: 'A': all, 'S': sphagnum, 'C': Carex, 'L': wood, list with length of x 
    output: (ThetaS and ThetaR in m3 m-3)
        van Genuchten water retention parameters as array [ThetaS, ThetaR, alpha, n] \n
        hydraulic conductivity (m/s)
    Arin koodi
    """
    #paras is dict variable, parameter estimates are stored in tuples, the model is water content = a0 + a1x + a2x2, where x is
    para={}  #'bd':bulk density in g/ cm3; 'H': von Post degree of humification
    para['bd'] ={'pF0':(97.95, -79.72, 0.0), 'pF1.5':(20.83, 759.69, -2484.3),
            'pF2': (3.81, 705.13, -2036.2), 'pF3':(9.37, 241.69, -364.6),
            'pF4':(-0.06, 249.8, -519.9), 'pF4.2':(0.0, 174.48, -348.9)}
    para['H'] ={'pF0':(95.17, -1.26, 0.0), 'pF1.5':(46.20, 8.32, -0.54),
            'pF2': (27.03, 8.14, -0.43), 'pF3':(17.59, 3.22, -0.07),
            'pF4':(8.81, 3.03, -0.10), 'pF4.2':(5.8, 2.27, -0.08)}
    
    intp_pF1={}  # interpolation functions for pF1        
    intp_pF1['bd'] = interp1d([0.04,0.08,0.1,0.2],[63.,84.,86.,80.],fill_value='extrapolate')
    intp_pF1['H'] = interp1d([1.,4.,6.,10.],[75.,84.,86.,80.],fill_value='extrapolate')
    
    #Saturatated hydraulic conductivity parameters
    Kpara ={'bd':{'A':(-2.271, -9.80), 'S':(-2.321, -13.22), 'C':(-1.921, -10.702), 'L':(-1.921, -10.702)}, 
            'H':{'A':(-2.261, -0.205), 'S':(-2.471, -0.253), 'C':(-1.850, -0.278), 'L':(-2.399, -0.124)}}

    x = np.array(x)
    prs = para[var]
    pF1=intp_pF1[var]
    if unit=='kg/m3'and var=='db': x=x/1000.
    if  np.shape(x)[0] >1 and len(ptype)==1:
        ptype=np.repeat(ptype, np.shape(x)[0])

    wcont = lambda x, a0, a1, a2: a0 + a1*x + a2*x**2.
    potentials =np.array([0.001, 1.,3.2, 10.,100.,1000.,1500.])
    wc = (np.array([wcont(x,*prs['pF0']), pF1(x), wcont(x,*prs['pF1.5']), wcont(x,*prs['pF2']),
               wcont(x,*prs['pF3']), wcont(x,*prs['pF4']),wcont(x,*prs['pF4.2'])]))
    wc = np.transpose(wc)
    pF_para = fit_pF(potentials, wc, fig=fig, labels=labels, percentage=True, kPa=True)

    Ksat = np.zeros((np.size(x)))
    K = lambda x, a0, a1: 10.**(a0 + a1*x) / 100.   # to m/s
    for i, a, pt in zip(range(len(x)), x, ptype):
        Ksat[i] = K(a, *Kpara[var][pt])  # hydraulic conductivity (cm/s -> m/s)

    return pF_para, Ksat

""" Functions for computing lad profiles """

def lad_profiles(grid, dbhfile, quantiles, hs, plot=False):
    """
    Leaf area density (lad) profiles for tree species and understory shrubs
    Args:
        grid (dict):
            'zmax': heigth of grid from ground surface [m]
            'Nlayers': number of layers in grid [-]
        dbhfile (str): file path to dbhfile
        quantiles (list): cumulative frequency limits for grouping trees
        normed (boolean): True returns sum(lad*dz) normalized to unity
    Returns:
        lad_p, lad_s, lad_d, lad_g (arrays): leaf area density profiles
            for model treegroups and understory shrubs (m2/m3)
    """
    # grid [m]
    z = np.linspace(0, grid['zmax'], grid['Nlayers'])
    # tress species lad profiles
    lad_p, lad_s, lad_d, _, _, _, lai_p, lai_s, lai_d = model_trees(z, quantiles, normed=True, dbhfile=dbhfile, plot=plot)
    # understory shrubs
    lad_g = np.zeros([len(z), 1])
    lad_g[z < hs] = 1.0
    lad_g[0] = 0.0
    if sum(lad_g) < 1.0:
        lad_g[1] = 1.0
    lad_g = lad_g / np.maximum(sum(lad_g * z[1]), eps)

    return lad_p, lad_s, lad_d, lad_g, lai_p, lai_s, lai_d

def model_trees(z, quantiles, normed=False,
                dbhfile='c:\\projects\\MLM_Hyde\\Data\\hyde_runkolukusarjat.txt',
                plot=False):
    """
    reads runkolukusarjat from Hyde and creates lad-profiles for pine, spruce and decid.
    Args:
        z - grid (m)
        quantiles - cumulative frequency limits for grouping trees
        normed - True returns sum(lad*dz) normalized to unity
    Returns:
        lad_p, lad_s, lad_d - leaf-area density profiles for model treegroups (m2/m3)
        n_p, n_s, n_d - trees / ha in model treegroups
    """
    dat = np.loadtxt(dbhfile, skiprows=1)
    dz = z[1]-z[0]

    M = len(quantiles)
    # year 2008 data
    pine = dat[:, [0, 1]]
    spruce = dat[:, [0, 2]]
    decid = dat[:, [0, 3]]

    # pines
    h, hb, mleaf, L, a = profiles_hyde(pine, 'pine', z)
    n = pine[:, 1]
    c = np.cumsum(n) / np.maximum(sum(n), eps)  # relative frequency
    m = 0.0
    lad_p = np.zeros([len(z), M])
    n_p = np.zeros(M)
    lai_p = np.zeros(M)
    for k in range(M):
        f = np.where((c > m) & (c <= quantiles[k]))[0]
        lad_p[:, k] = np.sum(a[:, f], axis=1)
        n_p[k] = np.sum(n[f])
        lai_p[k] = sum(dz*lad_p[:,k])
        m = quantiles[k]
        if normed:
            lad_p[:, k] = lad_p[:, k] / np.maximum(np.sum(lad_p[:, k] * dz), eps)

    # spruces
    h, hb, mleaf, L, a = profiles_hyde(spruce, 'spruce', z)
    n = spruce[:, 1]
    c = np.cumsum(n) / np.maximum(sum(n), eps)  # relative frequency
    m = 0.0
    lad_s = np.zeros([len(z), M])
    n_s = np.zeros(M)
    lai_s = np.zeros(M)
    for k in range(M):
        f = np.where((c > m) & (c <= quantiles[k]))[0]
        lad_s[:, k] = np.sum(a[:, f], axis=1)
        n_s[k] = np.sum(n[f])
        lai_s[k] = sum(dz*lad_s[:,k])
        m = quantiles[k]
        if normed:
            lad_s[:, k] = lad_s[:, k] / np.maximum(np.sum(lad_s[:, k] * dz), eps)

    # decid
    h, hb, mleaf, L, a = profiles_hyde(decid, 'birch', z)
    n = decid[:, 1]
    c = np.cumsum(n) / np.maximum(sum(n), eps)  # relative frequency
    m = 0.0
    lad_d = np.zeros([len(z), M])
    n_d = np.zeros(M)
    lai_d = np.zeros(M)
    for k in range(M):
        f = np.where((c > m) & (c <= quantiles[k]))[0]
        lad_d[:, k] = np.sum(a[:, f], axis=1)
        n_d[k] = np.sum(n[f])
        lai_d[k] = sum(dz*lad_d[:,k])
        m = quantiles[k]
        if normed:
            lad_d[:, k] = lad_d[:, k] / np.maximum(np.sum(lad_d[:, k] * dz), eps)

    if plot:
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        plt.figure(figsize=(2.5,3.5))
        for k in range(M):
            plt.plot(lad_p[:, k],z,color=colors[0], label='pine, %.2f m$^2$m$^{-2}$' % lai_p[k])#,lad_g,z)
            plt.plot(lad_s[:, k],z,color=colors[1], label='spruce, %.2f m$^2$m$^{-2}$' % lai_s[k])
            plt.plot(lad_d[:, k],z,color=colors[2], label='decid, %.2f m$^2$m$^{-2}$' % lai_d[k])
        plt.title("  ")#dbhfile.split("/")[-1])
        plt.ylabel('height [m]')
        if normed:
            plt.xlabel('normalized lad [-]')
        else:
            plt.xlabel('lad [m$^2$m$^{-3}$]')
        plt.tight_layout()

    return lad_p, lad_s, lad_d, n_p, n_s, n_d, lai_p, lai_s, lai_d

def profiles_hyde(data, species, z):
    # height h model based on Näslund equation parameterized using Hyytiälä stand inventory from 2008
    # trunk base ht model based on Tahvanainen & Forss, 2008 For. Ecol. Manag. Fig 4:
    d = data[:,0]  # dbh, cm
    N = data[:,1]  # trees ha-1

    # compute 
    if species == 'pine':
        h = 1.3 + d**2.0 / (1.108 + 0.197*d)**2.0
        ht = -3.0 + 0.76*h
        ht = np.maximum(0.0, ht)
        
        y, L = marklund(d, species)
        mleaf = y[3]*N  # kg/ha
        L = L*1e-4*N  # leaf area m2/m2
        
    if species == 'birch':
        h = 1.3 + d**2.0 / (0.674 + 0.201*d)**2.0
        ht = -2.34 + 0.58*h
        ht = np.maximum(0.0, ht)

        y, L = marklund(d, species)
        mleaf = y[3]*N  # kg/ha
        L = L*1e-4*N  # leaf area m2/m2
        
    if species == 'spruce':
        h = 1.3 + d**3.0 / (1.826 + 0.303*d)**3.0
        ht = -2.34 + 0.58*h
        ht = np.maximum(0.0, ht)        

        y, L = marklund(d, species)
        mleaf = y[3]*N  # kg/ha
        L = L*1e-4*N  # leaf area m2/m2

    # get crown leaf mass profiles
    lad = np.zeros([len(z),len(d)])
    for k in range(len(d)):
        a, _ = crown_biomass_distr(species,z,htop=h[k],hbase=ht[k])
        lad[:,k] =a*L[k]

    return h, ht, mleaf, L, lad

def marklund(d, species):
    """
    marklund(d, species):
    Computes tree biomasses for Scots pine, Norway spruce or birch based on Marklund et al. 1988.
    Returns biomasses of each tree compartment (in kg), based on diameter at breast height (d, in cm)\
    INPUT:
        d - DBH, diameter at breast height (cm)
        species - 'pine', 'spruce', 'birch'.
    OUTPUT:
        array y:
        y[0] - stem wood (kg)
        y[1] - stem bark
        y[2] - living branches incl. needles / leaves
        y[3] - needles /leaves
        y[4] - dead branches
        y[5] - stumps
        y[6] - roots, >=5cm
        y[7] - roots, <5cm
        L - one-sided leaf area (m2)\n
    SOURCE:
        Marklund L.G.(1988): Biomassafunktioner for tall, gran och björk i Sverige. SLU Rapport 45, 73p
        Kärkkäinen L., 2005 Appendix 4.
        Kellomäki et al. (2001). Atm. Env.\n
    AUTHOR:
        Samuli Launiainen, 25.4.2014
    """
    SLA = {'pine': 6.8, 'spruce': 4.7, 'decid': 14.0}  # Härkönen et al. 2015 BER 20, 181-195
    #d=np.array(float(d))

    if species.lower()=='pine':
        y1=np.exp(11.4219*(d/(d+14))-2.2184)    #stem wood
        y2=np.exp(8.8489*(d/(d+16))-2.9748)     #stem bark
        y3=np.exp(9.1015*(d/(d+10))-2.8604)     #living branches incl. needles
        y4=np.exp(7.7681*(d/(d+7))-3.7983)      #needles
        y5=np.exp(9.5938*(d/(d+10))-5.3338)     #dead branches
        y6=np.exp(11.0481*(d/(d+15))-3.9657)    #stumps
        y7=np.exp(13.2902*(d/(d+9))-6.3413)     #roots, >=5cm
        y8=np.exp(8.8795*(d/(d+10))-3.8375)     #roots, <5cm

        y=[y1, y2, y3, y4, y5, y6, y7, y8] # array of biomasses (kg)
        L=y4*SLA['pine'] # leaf area (m2) based on specific foliage area (SLA)
        
        return y, L
    
    elif species.lower()=='spruce':
        y1=np.exp(11.4873*(d/(d+14))-2.2471)   #stem wood
        y2=np.exp(9.8364*(d/(d+15))-3.3912)    #stem bark
        y3=np.exp(8.5242*(d/(d+13))-1.2804)    #living branches incl. needles
        y4=np.exp(7.8171*(d/(d+12))-1.9602)    #needles
        y5=np.exp(9.9550*(d/(d+18))-4.3308)    #dead branches
        y6=np.exp(10.6686*(d/(d+17))-3.3645)   # stumps
        y7=np.exp(13.3703*(d/(d+8))-6.3851)    #roots, >=5cm
        y8=np.exp(7.6283*(d/(d+12))-2.5706)    #roots, <5cm
        
        y=[y1, y2, y3, y4, y5, y6, y7, y8] # array of biomasses (kg)
        L=y4*SLA['spruce'] # leaf area (m2) based on specific foliage area (SLA)   
        
        return y, L
    
    elif species.lower()=='birch':
        #silver and downy birch
        y1=np.exp(10.8109*(d/(d+11))-2.3327)    #stem wood
        y2=np.exp(10.3876*(d/(d+14))-3.2518)    #stem bark
        y3=np.exp(10.2806*(d/(d+10))-3.3633)    #living brances excluding leaves
        y4=np.exp(8.0580*(d/(d+8))-3.9823)      #Kellomäki et al. 2001 Atm. Env.
        y5=np.exp(7.9266*(d/(d+5))-5.9507)      #dead branches

        y=[y1, y2, y3, y4, y5] # array of biomasses (kg)
        
        L=y4*SLA['decid'] # leaf area (m2) based on specific foliage area (SLA) 
        
        return y, L
        
    else:
        print('Vegetation.marklund: asked species (pine, spruce, birch) not found')
        return None

def crown_biomass_distr(species,z,htop=1,hbase=0,PlotFigs="False"):
    """
    crown_biomass_distr(species,z,htop=1,hbase=0,PlotFigs="False"):    
    Computes vertical crown and needle biomass profiles for pine, spruce or birch using models 
    of Fors & Tahvanainen (2008).
    INPUT:
        species - 'pine','spruce','birch'\n
        z - height array (m), constant increments\n
        htop - tree height (m),optional (default = 1, when z has to be relative height within crown [0...1] \n
        hbase - crown base height (m), optional \n
        PlotFigs - "True" plots density profiles\n
    OUTPUT:
        Ldens - relative leaf biomass density profile \n
        Cdens - relative crown biomass density profile \n
    SOURCE:
        Tahvanainen & Fors, 2008. For.Ecol.Manag. 255, 455-467\n
    AUTHOR:
        Samuli Launiainen, Metla, 26.6.2014\n
    NOTES:
        Leaf biomass profiles are available only for pine and spruce, use 
        crown biomass profile as surrogate for birch\n
    """
    species=species.lower()
    z=np.array(z)
    dz=z[1]-z[0] #m
    htop=float(htop) #m
    hbase=float(hbase) #m       
    N=np.size(z) #nodes
    
    #hrel= np.maximum(0,(z - hbase)/(htop-hbase)) #relative height within crown
    hrel=(z - hbase)/(htop-hbase)
    hrel[hrel>1]=1 
    hrel[hrel<0]=0

    # dictionarys of parameter values, key is species string
    # table8, uses crown biomass as surrogate for birch leaf biomass density
    NeedleBiomasspara = {'pine':[0.04276, 1.1365,3.43024,4.98748,0.99285],
                         'spruce':[0.04166, 1.31352, 2.62216, 3.96188, 0.98435],
                         'birch':[0.11902, 0.96061, 3.90684, 3.65942, 0.98929]}
    # table7
    CrownBiomasspara = {'pine':[0.07881, 1.03037,3.59557,3.66652,0.99071],
                        'spruce':[0.05892, 1.18495, 2.58915, 2.76521,0.98635],
                        'birch':[0.11902, 0.96061, 3.90684, 3.65942, 0.98929]}
                      
    #compute cumulative biomasses and vertical profiles
    cnd=np.zeros(N)
    Ldens=np.zeros(N)   
    ccd=np.zeros(N)
    Cdens=np.zeros(N)
  
    #needle/leaf density
    beta=NeedleBiomasspara[species]
        
    cnd=beta[0] + beta[1]*(1-np.exp(-beta[2]*hrel))**beta[3]*beta[4]
    cnd[hrel==0]=0
    cnd=cnd/max(cnd)    #cumulative needle biomass [0...1]
    
    #needle density profile, integral(Ldens*dz) over z gives unity
    Ldens[1:N]=np.diff(cnd,1)/dz
    Ldens[Ldens<0]=0

    #eliminate kink at crown base
    ind = (Ldens > 0.0).nonzero()[0][0]
    Ldens[ind]=np.mean(Ldens[ind-1:ind])       
    Ldens=Ldens/sum(Ldens*dz) #normalize integral to unity    

    del ind, beta
    
    #print sum(Ldens*dz)

    #crown biomass density                          
    beta=CrownBiomasspara[species]
        
    ccd=beta[0] + beta[1]*(1-np.exp(-beta[2]*hrel))**beta[3]*beta[4]
    ccd[hrel==0]=0
    ccd=ccd/max(cnd)    #cumulative crown biomass [0...1]
    
    #crown density profile, integral(Ldens*dz) over z gives unity
    Cdens[1:N]=np.diff(ccd,1)/dz
    Cdens[Cdens<0]=0
    
    #eliminate kink at crown base
    ind = (Cdens > 0.0).nonzero()[0][0]
    Cdens[ind]=np.mean(Ldens[ind-1:ind])       
    Cdens=Cdens/sum(Cdens*dz) #normalize integral to unity    
    #print sum(Cdens*dz)
    
    del ind, beta

    if PlotFigs=="True":
        
        plt.figure(99)            
        plt.plot(Ldens,z,'b.-',Cdens,z,'r.-')
        plt.title("Vegetation.crown_biomass_distr")
        plt.ylabel("z (m)")
        plt.xlabel("density")    
        plt.legend((species +' leafmass',species +' crownmass'),loc='best')      
    
    return Ldens,Cdens 