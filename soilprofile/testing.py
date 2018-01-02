# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:15:20 2017

@author: L1656
"""
t_max=24
t_final=3600;
z = np.arange(-0.01, -2.0, -0.02);
h0 = np.arange(-0.5, 1.5, 0.02);
h_pond0 = 0.0
_, gwl = get_gwl(h0,z)
plt.plot(h0,z,'k-')
plt.hold(True)
pF = {'ThetaS': 0.88, 'ThetaR': 0.093, 'alpha': 0.029, 'n': 1.34};
Ksat = np.ones(len(z))*1e-4;
Prec=0.000002
Evap=0.0 #000001
R=np.zeros(len(z))
R[0:9]=0.0000001
_, Qz_drain=drainage_hooghoud(z, Ksat, gwl, 1.0, 40.0, 1.0);
#Qz_drain=Qz_drain*0.0;
MBE=np.zeros(t_max)
GWL=np.zeros(t_max)
H_pond=np.zeros(t_max)
Theta=np.zeros(t_max)
for k in range(0, t_max):
    h, W, h_pond, C_inf, C_eva, C_dra, C_trans, C_roff, Fliq, gwl, KLh, mbe = waterFlow1D(t_final, z, h0, pF, Ksat, Prec, Evap, R,HM = Qz_drain, lbc={'type': 'impermeable', 'value': None}, maxPond=0.01, pond0=h_pond0, steps=10)
    h0=h.copy()
    h_pond0=h_pond
    plt.plot(h,z,'r-')
    _, Qz_drain=drainage_hooghoud(z, Ksat, gwl, 1.0, 40.0, 1.1, 2.0);
    MBE[k]=mbe
    GWL[k]=gwl
    H_pond[k]=h_pond
    Theta[k]=sum(W*0.02)


#plt.plot(MBE)