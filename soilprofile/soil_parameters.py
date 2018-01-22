# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 14:05:24 2017

@author: L1656
"""
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

def fit_pF(head, watcont, fig=False):
    """

    """

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
        plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left")

    return vgen_all

thickness = [0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0]
nodes = [10, 5, 5, 2, 2, 5, 5]
watcont = [[94.69, 49.42, 29.61, 21.56, 20.05, 17.83, 16.54],
           [91.41, 66.26, 56.98, 45.58, 41.44, 39.32, 37.89],
           [89.12, -999, 72.83, 63.97, 54.40, 50.15, 48.80],
           [89.46, -999, 82.46, 76.79, 66.93, 63.61, 62.53],
           [89.46, -999, 82.46, 76.79, 66.93, 63.61, 62.53],
           [89.46, -999, 82.46, 76.79, 66.93, 63.61, 62.53],
           [89.46, -999, 82.46, 76.79, 66.93, 63.61, 62.53]]
head = [0.0001, 0.3, 0.981, 4.905, 9.81, 33.0, 98.1]
Kvsat = [2e-4, 2e-5, 5e-6, 3e-6, 5e-6, 3e-6, 1e-6]
Khmult = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

pF_para = fit_pF(head, watcont, fig=False)

z = np.array([])
ThetaS = np.array([])
ThetaR = np.array([])
alpha = np.array([])
n = np.array([])

Ksat = np.array(Kvsat)
Khsat = np.array(Kvsat)*np.array(Khmult)

zh = np.array([0.0])

for k in range(0, len(thickness)):
    z = np.append(z, -np.arange(zh[k] + thickness[k]/(2 * nodes[k]), zh[k] + thickness[k], thickness[k] / nodes[k]))
    ThetaS = np.append(ThetaS, pF_para[k][0])
    ThetaR = np.append(ThetaR, pF_para[k][1])
    alpha = np.append(alpha, pF_para[k][2])
    n = np.append(n,  pF_para[k][3])
    zh = np.append(zh, zh[k] + thickness[k])

zh = -zh[1:]
"""
para = {
        'z': z,
        'zh':zh,
        'pF': {
                'ThetaS': ThetaS, 
                'ThetaR': ThetaR, 
                'alpha': alpha,     # 1/cm
                'n': n
                },                  # (dict): vanGenuchten water retention parameters; scalars or arrays of len(z)
        'Ksat': Ksat,               # (float/array): saturated vertical hydraulic conductivity [ms-1]
        'Khsat': Khsat,              # (float/array): saturated horizontal hydraulic conductivity [ms-1] - used in drainage equation
        'Csv': -np.ones(len(zh)),                # (float/array): dry soil vol. heat capacity [J m-3 (total volume) K-1]
        'vOrg': -np.ones(len(zh)),               # (float/array): organic matter fraction of solid volume [-]
        'vSand': -np.ones(len(zh)),              # (float/array): sand fraction of solid volume [-]
        'vSilt': -np.ones(len(zh)),              # (float(array): silt fraction of solid volume [-]
        'vClay': -np.ones(len(zh)),              # (float(array): clay fraction of solid volume [-]
        'fp': -np.ones(len(zh)),                 # (float/array): freezing curve parameter
        'max_pond': 0.01,           # (float) maximum pond depth [m]
        'ini_cond': {               # (dict): inputs are floats or arrays of len(z)
                'gwl': -0.2,        # (float) [m] or (float/array) Wtot', vol. water content [-] or 'h', matrix water potential [m]
                'T': -1.0,          # soil temperature [degC]
                'pond': 0.0
                },                  # [m] initial pond depth at surface
        'lbc_heat': -1.0,           # lower boundary condition type for heat
        'lbc_water': {
                'type': 'impermeable',
                'value': None,
                'depth': -2.0,
                },                  # lower boundary condition type for water (type, value, depth)
        'homogenous': False,         # (boolean): True assumes vertically homogenous profile and float inputs
        'solve_heat': False,        # (boolean): True solves heatflow
        'solve_water': True,        # (boolean): True solves waterflow
        'solve_water': True,        # (boolean): True solves waterflow
        'solve_water_type':'Equilibrium',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
        'Bedrock': {
                'Cv': 2160000.0,
                'Lambda': 3.0
                },
        'drainage_equation': {
                'type': 'Hooghoudt',
                'depth': 1.0,
                'spacing': 40.0,
                'width': 1.0,
                }                   # (dict): Drainage equation and drainage parameters
        }

"""
para = {
        'z': z, # np.arange(-0.005, -2.0, -0.01),
        'pF': {
                'ThetaS': 0.88, 
                'ThetaR': 0.093, 
                'alpha': 0.029,     # 1/cm
                'n': 1.34
                },                  # (dict): vanGenuchten water retention parameters; scalars or arrays of len(z)
        'Ksat': 1e-5,               # (float/array): saturated vertical hydraulic conductivity [ms-1]
        'Khsat': 1e-5,              # (float/array): saturated horizontal hydraulic conductivity [ms-1] - used in drainage equation
        'Csv': -1.0,                # (float/array): dry soil vol. heat capacity [J m-3 (total volume) K-1]
        'vOrg': -1.0,               # (float/array): organic matter fraction of solid volume [-]
        'vSand': -1.0,              # (float/array): sand fraction of solid volume [-]
        'vSilt': -1.0,              # (float(array): silt fraction of solid volume [-]
        'vClay': -1.0,              # (float(array): clay fraction of solid volume [-]
        'fp': -1.0,                 # (float/array): freezing curve parameter
        'max_pond': 0.01,           # (float) maximum pond depth [m]
        'ini_cond': {               # (dict): inputs are floats or arrays of len(z)
                'gwl': -0.2,        # (float) [m] or (float/array) Wtot', vol. water content [-] or 'h', matrix water potential [m]
                'T': -1.0,          # soil temperature [degC]
                'pond': 0.0
                },                  # [m] initial pond depth at surface
        'lbc_heat': -1.0,           # lower boundary condition type for heat
        'lbc_water': {
                'type': 'impermeable',
                'value': None,
                'depth': -2.0,
                },                  # lower boundary condition type for water (type, value, depth)
        'homogenous': True,         # (boolean): True assumes vertically homogenous profile and float inputs
        'solve_heat': False,        # (boolean): True solves heatflow
        'solve_water': True,        # (boolean): True solves waterflow
        'solve_water': True,        # (boolean): True solves waterflow
        'solve_water_type':'Equilibrium',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
        'Bedrock': {
                'Cv': 2160000.0,
                'Lambda': 3.0
                },
        'drainage_equation': {
                'type': 'Hooghoudt',
                'depth': 1.0,
                'spacing': 40.0,
                'width': 1.0,
                }                   # (dict): Drainage equation and drainage parameters
        }

        