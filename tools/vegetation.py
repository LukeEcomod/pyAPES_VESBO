# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 12:07:47 2014

@author: slauniai
"""
    
import numpy as np
import matplotlib.pylab as plt
         

"""
Package Vegetation:
    Contains general methods related to vegetation biomass amount and its variation among sites.\n
    Call all methods as e.g. "cVege.fUndestoryBiomassModel(args)", where "cVege" is instance of cVegetation -class. \n
    See documentation of interface of each method for inputs, outputs and datatypes.\n
    Uses numpy for array manipulation & calculation.\n
FUNCTIONS:
    marklund- tree allometry for Scots pine, Norway spruce and birch (Marklund et al. 1988).\n
    crown_biomass_distr- relative leaf and crown (leaf+branch) biomass distribution (Tahvanainen & Forss, 2008). \n
    understory_biomass - understory biomass per plant group in different stands (Muukkonen & Mäkipää, 2006).\n
AUTHOR:
    Samuli Launiainen, METLA 6/2014 \n
VERSION:
    16.11.2016 (SL): small edits, sent to Aleksi \n     
"""

def model_trees(dbhfile, z, quantiles=[0.85, 1.0], normed=False):
    """
    reads runkolukusarjat from Hyde and creates lad-profiles for pine, spruce and decid.
    Uses 2008 data
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
    
    # use year 2008 data
    pine = dat[:, [0, 2]]
    spruce = dat[:, [0, 4]]
    decid = dat[:, [0, 6]]

    # pines
    h, hb, mleaf, L, a = profiles_hyde(pine, 'pine', z)
    n = pine[:, 1]
    c = np.cumsum(n) / sum(n)  # relative frequency
    m = 0.0
    lad_p = np.zeros([len(z), M])
    n_p = np.zeros(M)
    for k in range(M):
        f = np.where((c > m) & (c <= quantiles[k]))[0]
        lad_p[:, k] = np.sum(a[:, f], axis=1)
        n_p[k] = np.sum(n[f])
        m = quantiles[k]
        if normed:
            lad_p[:, k] = lad_p[:, k] / np.sum(lad_p[:,k] * dz)

    # spruces
    h, hb, mleaf, L, a = profiles_hyde(spruce, 'spruce', z)
    n = spruce[:, 1]
    c = np.cumsum(n) / sum(n)  # relative frequency
    m = 0.0
    lad_s = np.zeros([len(z), M])
    n_s = np.zeros(M)
    for k in range(M):
        f = np.where((c > m) & (c <= quantiles[k]))[0]
        lad_s[:, k] = np.sum(a[:, f], axis=1)
        n_s[k] = np.sum(n[f])
        m = quantiles[k]
        if normed:
            lad_s[:, k] = lad_s[:, k] / np.sum(lad_s[:, k] * dz)

    # decid
    h, hb, mleaf, L, a = profiles_hyde(decid, 'birch', z)
    n = decid[:, 1]
    c = np.cumsum(n) / sum(n)  # relative frequency
    m = 0.0
    lad_d = np.zeros([len(z), M])
    n_d = np.zeros(M)
    for k in range(M):
        f = np.where((c > m) & (c <= quantiles[k]))[0]
        lad_d[:, k] = np.sum(a[:, f], axis=1)
        n_d[k] = np.sum(n[f])
        m = quantiles[k]
        if normed:
            lad_d[:, k] = lad_d[:, k] / np.sum(lad_d[:, k] * dz)

    return lad_p, lad_s, lad_d, n_p, n_s, n_d


def naslund_heightcurve(d, species, p=None):
    """
    Computes tree height from dbh using Näslund -equation parameterized by
    Siipilehto & Kangas, 2015 Metsätalouden Aikakausikirja
    Args:
        d - DBH, diameter at breast height (cm)
        species - 'pine', 'spruce', 'birch'
    Returns:
        z - scalar or array of tree heights (m)
    """
    para = {'pine': [2.0, 1.181, 0.247], 'spruce': [3.0, 1.703, 0.327], 'birch': [2.0, 1.014, 0.238]}
    if p is None:    
        p = para[species]
    
    z = 1.3 + d**p[0] / ( p[1] + p[2]*d)**p[0]

    return z    


def crownheight(h, species):
    """
    Computes crown height from relationships derived from Tahvanainen & Forss, 2008 For. Ecol. Manag.
    Fig 4
    Args:
        h - height (m)
    Returns:
        ht - height to crown base
    """
    
    if species == 'pine':
        ht = -3.0 + 0.76*h
        ht = np.maximum(0.0, ht)
    if species == 'spruce':
        ht = -1.99 + 0.51*h - 7.7e-3*h**2.0
        ht = np.maximum(0.0, ht)
    if species == 'birch':
        ht = -2.23 + 0.58*h
        ht = np.maximum(0.0, ht)

    return ht

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
    SLA = {'pine': 6.2, 'spruce': 4.9, 'decid': 12.0}  # Majasalmi et al. 2013 ForEcol, Luoma 1997
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

def crown_biomass_distr(species,z,htop=1,hbase=0,PlotFigs=False):
    """
    crown_biomass_distr(species,z,htop=1,hbase=0,PlotFigs=False):    
    Computes vertical crown and needle biomass profiles for pine, spruce or birch using models 
    of Fors & Tahvanainen (2008).
    INPUT:
        species - 'pine','spruce','birch'\n
        z - height array (m), constant increments\n
        htop - tree height (m),optional (default = 1, when z has to be relative height within crown [0...1] \n
        hbase - crown base height (m), optional \n
        PlotFigs - True plots density profiles\n
    OUTPUT:
        Ldens - relative leaf biomass density profile \n
        Cdens - relative crown biomass density profile \n
    SOURCE:
        Tahvanainen & Fors, 2008. For.Ecol.Manag. 255, 455-467\n
    AUTHOR:
        Samuli Launiainen, Metla, 26.6.2014\n
    NOTES:
        Leaf biomass profiles are available only for pine and spruce, use crown biomass profile as surrogate for birch\n
        
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

    plt.plot(Ldens,z,'b.-',Cdens,z,'r.-')

    if PlotFigs:
        
        plt.figure(99)            
        plt.plot(Ldens,z,'b.-',Cdens,z,'r.-')
        plt.title("Vegetation.crown_biomass_distr")
        plt.ylabel("z (m)")
        plt.xlabel("density")    
        plt.legend((species +' leafmass',species +' crownmass'),loc='best')      
    
    return Ldens,Cdens 


def lad_weibul(z, LAI, h, b=None, c=None, species=None):
    """
    Generates leaf-area density profile from Weibull-distribution
    INPUT:
        z: height array (m), monotonic and constant steps
        LAI: leaf-area index (m2m-2)
        h: canopy height (m), scalar
        b: Weibull shape parameter 1, scalar
        c: Weibull shape parameter 2, scalar
        species: 'pine', 'spruce', 'birch' to use table values
    OUTPUT:
        LAD: leaf-area density (m2m-3), array \n
    SOURCE:
        Teske, M.E., and H.W. Thistle, 2004, A library of forest canopy structure for use in interception modeling.
        Forest Ecology and Management, 198, 341-350. 
        Note: their formula is missing brackets for the scale param. 
    AUTHOR:
        Gabriel Katul, 2009. Coverted to Python 16.4.2014 / Samuli Launiainen
    """
    
    para = {'pine': [0.906, 2.145], 'spruce': [2.375, 1.289], 'birch': [0.557, 1.914]} 

    if b is None or c is None:
        b, c = para[species]
    
    z = np.array(z)
    dz = abs(z[1]-z[0])
    N = np.size(z)
    LAD = np.zeros(N)

    # dummy variables
    a = np.zeros(N)
    x = z[z <= h] / h  # normalized heigth
    n = np.size(x)

    # weibul-distribution
    a[0:n] = -(c/b)*(((1.0-x)/b)**(c-1.0))*(np.exp(-((1.0-x) / b)**c)) / (1.0-np.exp(-(1.0 / b)**c))
    # a[0] = 0.0  # no leaves at ground
    a = a / sum(a*dz)

    LAD = LAI*a

    # plt.figure(1)
    # plt.plot(LAD,z,'r-')      
    return LAD
    

def understory_biomass(site, age, x=[]):
    """
    understory_biomass(site, age, x=[]):    
    Computes understory biomasses using models of Muukkonen & Makipaa, 2006 Bor. Env. Res.\n
    INPUT:
        site - site type string:
            'pine_upland'
            'spruce_upland'
            'broadleaved_upland'
            'spruce_mire'
            'pine_mire'
        age - stand age (yrs), scalar \n
        x - array of independent variables (optional, if not provided age-based model is used):
            x[0]=lat (degN)
            x[1]=lon (degN) 
            x[2]=elev (m)
            x[3]=temperature sum (degC)
            x[4]=site nutrient level (-) 
            x[5]=stem vol. (m3 ha-1)
            x[6]=stem nr (ha-1)
            x[7]=basal area (m2 ha-1)
            x[8]=site drainage status,integer
    OUTPUT:
        y - dry biomasses (kg ha-1) of different groups\n
        gname - array of group names corresponding to y columns. Bottom layer includes Mosses + Lichens, Field layer includes Herbs&Grasses, Shrubs \n
    SOURCE:
        Muukkonen & Makipaa, 2006. Bor.Env.Res. 11, 355-369.\n
    AUTHOR:
        Samuli Launiainen 18.06.2014 \n
    NOTE:
         For age-only models: Total biomasses and group-level biomasses may be highly different for spruce_upland and broadleaved_upland! \n
         Multi-regression models not yet tested!
         In model equations independent variables named differently to M&M (2006): here x[0] = z1, x[1]=z2, ... x[7]=z8 and x[8]=z10\n
         \n
         Site nutrient level x[4] at upland sites:
             1: herb-rich forest 
             2: herb-rich heat f. 
             3: mesic heath f. 
             4: sub-xeric heath f.
             5: xeric heath f. 
             6: barren heath f.
             7: rock,cliff or sand f. 
         Site nutrient level x[4] at mires:\n
             1: herb-rich hw-spruce swamps, pine mires, fens, 
             2: V.myrtillus / tall sedge spruce swamps, tall sedge pine fens, tall sedge fens,
             3: Carex clobularis / V.vitis-idaea swamps, Carex globularis pine swamps, low sedge (oligotrophic) fens,
             4: Low sedge, dwarf-shrub & cottongrass pine bogs, ombo-oligotrophic bogs,
             5: S.fuscum pine bogs, ombotrophic and S.fuscum low sedge bogs.
         Drainage status x[8] at mires (Paavilainen & Paivanen, 1995):
             1: undrained
             2: Recently draines, slight effect on understory veg., no effect on stand
             3: Transforming drained mires, clear effect on understory veg and stand
             4: Transforming drained mires, veget. resembles upland forest site type, tree-stand forest-like.
  
    """
        
    
    age=np.array(age)
    x=np.array(x)
    
    if np.size(x)==0:
        model = 'StandAge'
    else:
        model = 'MultiRegression'    
    
    print('Vegetation.understory_biomass: computing sitetype ' +site +' using ' +model +' -based model')

    if site.lower()=='pine_upland':
        
        gname=['Dwarf shrubs','Herbs&Grasses','Mosses','Lichens','Field layer total','Bottom layer total','Total']
        y=np.zeros(7)
        
        if np.size(x)==0:  #only age as explaning variable
            # dependent variable is sqrt(biomass -0.5)
            y[0] = 16.68 + 0.129*age - 4e-4*age**2  #Dwarf shrubs
            y[1] = 11.725 -0.098*age +2e-4*age**2   #Herbs&grasses
            y[2] = 27.329 + 0.138*age -5e-4*age**2  #mosses
            y[3] = 7.975 -2e-4*age**2               #lichens
            y[4] = 22.521 + 0.069*age -2e-4*age**2  #Field layer total
            y[5] = 32.952 + 0.085*age -6e-7*age**3  #Bottom layer total
            y[6] = 42.641 + 0.094*age -8e-7*age**3  #Total
            
            cf=[126.91, 55.02, 361.44, 260.56, 96.72, 355.13, 231.56] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1
            
            return Y, gname
        else:
            # dependent variable is sqrt(biomass -0.5)
            y[0] = -10.328 + 0.005*x[1]*age -4e-4*age**2 -3e-5*x[3] +9e-4*x[0]*x[3]  #Dwarf shrubs
            y[1] = 15.223 + 0.023*x[4]*age -0.062*x[4]*x[0] +9e-6*x[3]**2 -1e-4*x[3]*age -0.028*x[4]*x[7]  #Herbs&grasses
            y[2] = 68.365 +3e-4*x[3]*age -3e-4*x[5]**2 -0.07*x[4]*age +2e-4*x[5]*x[2] -0.2*x[4]*x[0] -4e-5*x[3]**2 +0.014*x[4]*x[3] +0.01*x[7]**2  #mosses
            y[3] = -53.196 +0.378*x[4]*x[0] -0.014*x[4]*x[3] + 2e-5*x[3]**2 #lichens
            y[4] = 13.865 +0.013*x[0]*x[1] -2.969*x[4] +3e-5*x[3]*age #Field layer total
            y[5] = 8.623 +0.09*x[4]*x[0] +0.004*x[1]*age +3e-5*x[2]*x[3] -3e-4*x[5]**2 -5e-4*age**2 + 8e-4*x[5]*age #Bottom layer total
            y[6] = 22.523 + 0.084*x[4]*x[0] + 0.01*x[1]*age -0.031*x[4]*age -7e-4*age**2 -3e-4*x[5]**2 +6e-4*x[5]*age  #Total
            
            cf=[126.91, 55.02, 361.44, 260.56, 96.72, 355.13, 231.56] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1                
            return Y, gname
            
    if site.lower()=='spruce_upland':
        
        gname=['Dwarf shrubs','Herbs&Grasses','Mosses','Field layer total','Total']
        y=np.zeros(5)         
        
        if np.size(x)==0:  #only age as explaning variable
            # dependent variable is sqrt(biomass -0.5)
            y[0] = 10.375 -0.033*age + 0.001*age**2 -6e-6*age**3  #Dwarf shrubs
            y[1] = 15.058 -0.113*age + 3e-4*age**2  #Herbs&grasses
            y[2] = 19.282 + 0.164*age -1e-6*age**3  #mosses
            y[3] = 15.399 +0.036*age  #Field layer total
            y[4] = 27.349 + 0.157 -1e-6*age**3  #Total
            
            cf=[87.14, 44.60, 264.82, 67.15, 206.67] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1
            return Y, gname 
        else:
            y[0] = 10.903 +0.027*x[4]*age + 0.045*x[4]*x[0] -6e-6*x[3]**2  #Dwarf shrubs
            y[1] = 21.49 -0.05*x[4]*x[0] -0.006*x[4]*x[5] -0.008*x[4]*age  #Herbs&grasses
            y[2] = 9.672 +0.029*x[4]*age +0.078*x[4]*x[0] +0.186*x[7] #mosses
            y[3] = -42.593 + 0.981*x[0] -0.008*x[7]**2 +0.002*x[7]*age  #Field layer total
            y[4] = 22.522 +0.026*x[4]*age +0.11*x[4]*x[0] -0.003*x[4]*x[3]  #Total
            
            cf=[87.14, 44.60, 264.82, 67.15, 206.67] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1
            return Y, gname    

    if site.lower()=='broadleaved_upland':
        
        gname=['Dwarf shrubs','Herbs&Grasses','Mosses','Field layer total','Total']
        y=np.zeros(5)         
        
        if np.size(x)==0:  #only age as explaning variable
            # dependent variable is sqrt(biomass -0.5)
            y[0] = 7.102 +4e-4*age**2 #Dwarf shrubs
            y[1] = 20.58 - 0.423*age + 3e-4*age**2 -2e-5*age**3  #Herbs&grasses
            y[2] = 13.555 -0.056*age #mosses
            y[3] = 18.831 + 2e-4*age**2  #Field layer total
            y[4] = 25.645 + 0.037*age #Total
            
            cf=[77.67, 55.55, 236.6, 55.40, 156.51] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1 
            return Y, gname                                   

        else:
            # dependent variable is sqrt(biomass -0.5)
            y[0] =  3.217 +0.034*x[4]*x[2] -3e-4*x[2]**2 +5e-4*x[2]*age -0.077*x[4]*x[7] #Dwarf shrubs
            y[1] =  -192.32 -3.451*x[4] + 0.117*x[0]*x[1] -1e-5*x[2]*x[3] -0.065*x[1]**2 +0.002*x[0]*x[3] -0.003*x[1]*x[3] #Herbs&grasses
            y[2] = 20.931 +0.096*x[4]*x[0] -6e-4*x[1]*x[3] #mosses
            y[3] = -95.393 +0.094*x[0]*x[1] -1e-6*x[6]*x[3] -0.106*x[1]**2 +5e-4*x[0]*x[3] #Field layer total
            y[4] = 19.8 +0.691*x[4]*x[0] -38.578*x[4]#Total
            
            cf=[77.67, 55.55, 236.6, 55.40, 156.51] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1
            return Y, gname   
               
    if site.lower()=='spruce_mire': #hardwood-spruce mires & paludified forests
        
        gname=['Bottom layer total','Field layer total','Total']
        y=np.zeros(3)         
        if np.size(x)==0:  #only age as explaning variable
            # dependent variable is sqrt(biomass -0.5)
            y[0] = 25.923 + 1e-6*age**3 #Bottom layer total
            y[1] = 19.334 + 0.094*age #Field layer total
            y[2] = 36.655 + 4e-4*age**2  #Total
            
            cf=[98.10, 162.58, 116.54] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1
            return Y, gname  
        else:
            y[0] =  -3.182 + 0.022*x[0]*x[1] +2e-4*x[2]*age -0.077*x[4]*x[1] -0.003*x[1]*x[5] + 2e-4*x[5]**2  #Bottom layer total
            y[1] =  23.24 -1.163*x[9]**2 +1.515*x[4]*x[9] -2e-5*x[5]*x[6] +8e-5*x[3]*age +1e-5*x[6]*x[2] #Field layer total
            y[2] =  35.52 +0.001*x[1]*x[2] -1.1*x[8]**2 -2e-5*x[5]*x[6] +4e-5*x[6]*age +0.139*x[1]*x[8] #Total
            
            cf=[98.10, 162.58, 116.54] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1               
            return Y, gname             
    
    if site.lower()=='pine_mire': #pine mires
        
        gname=['Bottom layer total','Field layer total','Total']
        y=np.zeros(3)         
        if np.size(x)==0:  #only age as explaning variable
            # dependent variable is sqrt(biomass -0.5)
            y[0] = 36.039 +3e-6*age**3 #Bottom layer total
            y[1] = 35.861 + 0.026*age #Field layer total
            y[2] = 51.877 + 0.042*age  #Total
            
            cf=[222.22, 133.26, 167.40] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1
            return Y, gname                 
        else: 
            y[0] =  31.809 +0.008*x[1]*x[2] -3e-4*x[6]*x[7] +6e-5*x[6]*age -0.188*x[2] #Bottom layer total
            y[1] =  48.12 -1e-5*x[3]**2 +0.013*x[4]*age -0.04*x[5]*age +0.026*x[4]*x[5] #Field layer total
            y[2] =  50.098 +0.005*x[1]*x[2] -1e-5*x[5]*x[6] +0.026*x[4]*age -1e-4*x[2]*x[3] -0.014*x[5]*x[8]#Total           

            cf=[222.22, 133.26, 167.40] #correction factor
            Y = y**2 -0.5 +cf #kg ha-1                
            return Y, gname  
