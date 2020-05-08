# -*- coding: utf-8 -*-
"""
General utility functions

Note:
    migrated to python3
    - nothing changed

@author: Kersti Haahti
"""

import numpy as np

def forward_diff(y, dx):
    """
    computes gradient dy/dx using forward difference
    assumes dx is constatn
    """
    N = len(y)
    dy = np.ones(N) * np.NaN
    dy[0:-1] = np.diff(y)
    dy[-1] = dy[-2]
    return dy / dx


def central_diff(y, dx):
    """
    computes gradient dy/dx with central difference method
    assumes dx is constant
    """
    N = len(y)
    dydx = np.ones(N) * np.NaN
    # -- use central difference for estimating derivatives
    dydx[1:-1] = (y[2:] - y[0:-2]) / (2 * dx)
    # -- use forward difference at lower boundary
    dydx[0] = (y[1] - y[0]) / dx
    # -- use backward difference at upper boundary
    dydx[-1] = (y[-1] - y[-2]) / dx

    return dydx

def tridiag(a, b, C, D):
    """
    tridiagonal matrix algorithm
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
        V[i] = b[i] - a[i] * G[i - 1]
        U[i] = (D[i] - a[i] * U[i - 1]) / V[i]
        G[i] = C[i] / V[i]

    x[-1] = U[-1]
    inn = n - 2
    for i in range(inn, -1, -1):
        x[i] = U[i] - G[i] * x[i + 1]
    return x

def smooth(a, WSZ):
    """
    smooth a by taking WSZ point moving average.
    NOTE: even WSZ is converted to next odd number.
    """
    WSZ = int(np.ceil(WSZ) // 2 * 2 + 1)
    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid') / WSZ
    r = np.arange(1, WSZ-1, 2)
    start = np.cumsum(a[:WSZ-1])[::2] / r
    stop = (np.cumsum(a[:-WSZ:-1])[::2] / r)[::-1]
    x = np.concatenate((start, out0, stop))
    return x

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

def lad_weibul(z, LAI, h, hb=0.0, b=None, c=None, species=None):
    """
    Generates leaf-area density profile from Weibull-distribution
    Args:
        z: height array (m), monotonic and constant steps
        LAI: leaf-area index (m2m-2)
        h: canopy height (m), scalar
        hb: crown base height (m), scalar
        b: Weibull shape parameter 1, scalar
        c: Weibull shape parameter 2, scalar
        species: 'pine', 'spruce', 'birch' to use table values
    Returns:
        LAD: leaf-area density (m2m-3), array \n
    SOURCE:
        Teske, M.E., and H.W. Thistle, 2004, A library of forest canopy structure for 
        use in interception modeling. Forest Ecology and Management, 198, 341-350. 
        Note: their formula is missing brackets for the scale param.
        Here their profiles are used between hb and h
    AUTHOR:
        Gabriel Katul, 2009. Coverted to Python 16.4.2014 / Samuli Launiainen
    """
    
    para = {'pine': [0.906, 2.145], 'spruce': [2.375, 1.289], 'birch': [0.557, 1.914]} 
    
    if (max(z) <= h) | (h <= hb):
        raise ValueError("h must be lower than uppermost gridpoint")
        
    if b is None or c is None:
        b, c = para[species]
    
    z = np.array(z)
    dz = abs(z[1]-z[0])
    N = np.size(z)
    LAD = np.zeros(N)

    a = np.zeros(N)

    # dummy variables
    ix = np.where( (z > hb) & (z <= h)) [0]
    x = np.linspace(0, 1, len(ix)) # normalized within-crown height

    # weibul-distribution within crown
    cc = -(c / b)*(((1.0 - x) / b)**(c - 1.0))*(np.exp(-((1.0 - x) / b)**c)) \
            / (1.0 - np.exp(-(1.0 / b)**c))

    a[ix] = cc
    a = np.abs(a / sum(a*dz))    

    LAD = LAI * a

    # plt.figure(1)
    # plt.plot(LAD,z,'r-')      
    return LAD

def lad_constant(z, LAI, h, hb=0.0):
    """
    creates constant leaf-area density distribution from ground to h.

    INPUT:
        z: height array (m), monotonic and constant steps
        LAI: leaf-area index (m2m-2)
        h: canopy height (m), scalar
        hb: crown base height (m), scalar
     OUTPUT:
        LAD: leaf-area density (m2m-3), array
    Note: LAD must cover at least node 1
    """
    if max(z) <= h:
        raise ValueError("h must be lower than uppermost gridpoint")

    z = np.array(z)
    dz = abs(z[1]-z[0])
    N = np.size(z)
    
#    # dummy variables
#    a = np.zeros(N)
#    x = z[z <= h] / h  # normalized heigth
#    n = np.size(x)
#
#    if n == 1: n = 2
#    a[1:n] = 1.0
    
    # dummy variables
    a = np.zeros(N)
    ix = np.where( (z > hb) & (z <= h)) [0]
    if ix.size == 0:
        ix = [1]

    a[ix] = 1.0
    a = a / sum(a*dz)
    LAD = LAI * a
    return LAD