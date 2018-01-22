# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:05:05 2018

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

    plt.xlabel(r'$\theta$  $(m^3m^{-3})$', fontsize=14)
    plt.ylabel('$-h$ $(cm)$', fontsize=14)
    plt.ylim(xx[0], xx[-1])
    plt.xlim(0.0, 1.0)
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left")

    return vgen_all
