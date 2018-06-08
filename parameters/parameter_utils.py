# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:05:05 2018

@author: L1656
"""
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
eps = np.finfo(float).eps  # machine epsilon

def fit_pF(head, watcont, fig=False):
    """

    """

    if fig:
        plt.figure()
        colors = ['r', 'b', 'g', 'm', 'c']
        c = 0
    head = np.array(head)
    head = head * 10  # kPa -> cm
    vg_ini=(0.88,	 0.09, 0.03, 1.3)
    van_g = lambda h, *p:   p[1] + (p[0] - p[1]) / (1. + (p[2] * h) **p[3]) **(1. - 1. / p[3])
    vgen_all = []

    for k in range(0, len(watcont)):
        Wcont = np.array(watcont[k])
        ix = np.where(Wcont >= 0)
        Wcont[ix] = Wcont[ix] / 100  # % -> fraction
        try:
            vgen, _ = curve_fit(van_g, head[ix], Wcont[ix], p0=vg_ini)
            label='pF: Ts=%5.3f, Tr=%5.3f, alfa=%5.3f, n=%5.3f' % tuple(vgen)
        except RuntimeError:
            vgen = [-1, -1, -1, -1]
            label='No fit!'
        vgen_all.append(vgen)

        if fig:

            plt.semilogy(Wcont[ix], head[ix], '.',color = colors[c])
            xx = np.logspace(-1, 4.2, 100)
            plt.semilogy(van_g(xx, *vgen), xx, '-',color = colors[c],
                         label=label)
            c += 1
            if c > 4:
                c = 0

    if fig:
        plt.xlabel(r'$\theta$  $(m^3m^{-3})$', fontsize=14)
        plt.ylabel('$-h$ $(cm)$', fontsize=14)
        plt.ylim(xx[0], xx[-1])
        plt.xlim(0.0, 1.0)
        plt.legend()

    return vgen_all

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
    lad_g = np.ones([len(z), 1])
    lad_g[z > hs] = 0.0
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
        plt.figure(figsize=(3,4))
        for k in range(M):
            plt.plot(lad_p[:, k],z,color=colors[0])#,lad_g,z)
            plt.plot(lad_s[:, k],z,color=colors[1])
            plt.plot(lad_d[:, k],z,color=colors[2])
            plt.legend(['pine','spruce','decid'])
        plt.title("  ")#dbhfile.split("/")[-1])
        plt.ylabel('height [m]')
        if normed:
            plt.xlabel('normalized lad [-]')
        else:
            plt.xlabel('lad [m$^2$m$^{-3}$]')
            plt.legend(['pine, %.2f m$^2$m$^{-2}$' % lai_p[0],
                        'spruce, %.2f m$^2$m$^{-2}$' % lai_s[0],
                        'decid, %.2f m$^2$m$^{-2}$' % lai_d[0]],
                        frameon=False, labelspacing=0.1, borderpad=0.0)
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
        print 'Vegetation.marklund: asked species (pine, spruce, birch) not found'
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