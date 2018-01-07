# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:15:20 2017

@author: L1656
"""
t_max=100
t_final=3600;
z = np.arange(-0.01, -2.0, -0.02)

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

h0 = -0.0 - z
h_pond0 = 0.0
_, gwl = get_gwl(h0,z)
plt.plot(h0,z,'k-')
plt.hold(True)
pF = {'ThetaS': np.ones(len(z))*0.88, 'ThetaR': np.ones(len(z))*0.093, 'alpha': np.ones(len(z))*0.029, 'n': np.ones(len(z))*1.34};
Ksat = np.ones(len(z))*1e-5;
Prec=0.000000
Evap=0.0000001
R=np.zeros(len(z))
#R[0:9]=0.000001
_, Qz_drain=drainage_hooghoud(z, Ksat, gwl, 1.0, 40.0, 1.0);
#Qz_drain=Qz_drain*0.0;
Dra=np.zeros(t_max)
MBE=np.zeros(t_max)
GWL=np.zeros(t_max)
H_pond=np.zeros(t_max)
Theta=np.zeros(t_max)
for k in range(0, t_max):
    print 'k = ' + str(k)
    h, W, h_pond, C_inf, C_eva, C_dra, C_trans, C_roff, Fliq, gwl, KLh, mbe, dto = waterFlow1D(t_final, z, h0, pF, Ksat, Prec, Evap, R,HM = Qz_drain, lbc={'type': 'flux', 'value': 0.0}, maxPond=0.01, pond0=h_pond0, steps=2)
    h0=h.copy()
    h_pond0=h_pond
    plt.plot(h,z,'r-')
    _, Qz_drain=drainage_hooghoud(z, Ksat, gwl, 1.0, 40.0, 1.0);
    Dra[k]=C_dra
    MBE[k]=mbe
    GWL[k]=gwl
    H_pond[k]=h_pond
    Theta[k]=sum(W*dz)
    if k == 30:
        Prec=0.000001


#plt.plot(MBE)