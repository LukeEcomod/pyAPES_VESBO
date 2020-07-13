# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:49:30 2020

@author: 03110850
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyAPES_utilities.parameter_utilities import aleksis_combination, crown_biomass_distr
eps = np.finfo(float).eps  # machine epsilon

data = pd.read_csv(r'O:\Projects\Antoine_SLU\Stand data\TreeInventory_Selection_Footprint80.csv',
                   sep=',', header='infer', encoding = 'ISO-8859-1')

data['Height'] = data['Height'] / 10.0

data['basal_area'] = np.pi * (data['Diameter'] / 100.0 / 2)**2

data['species_short'] = 'birch'
data.loc[(data['Species'] == 'PinusSylvestris'),'species_short'] = 'pine'
data.loc[(data['Species'] == 'PiceaAbies'),'species_short'] = 'spruce'

# plt.figure()
# for species in set(data['species_short']):
#     plt.plot(data[data['species_short'] == species]['Diameter'],data[data['species_short'] == species]['Height'],'o')

# Crown base height (m): Tahvanainen & Forss, 2008 For. Ecol. Manag. Fig 4
data['crown_base_height'] = -1.9 + 0.70 * data['Height']  # pine
data.loc[(data['species_short'] == 'spruce'),'crown_base_height'] = -0.4 + 0.28 * data['Height']  # spruce
data.loc[(data['species_short'] == 'birch'),'crown_base_height'] = -0.8 + 0.48 * data['Height']  # birch
data['crown_base_height'] = np.maximum(0.0, data['crown_base_height'])

# plt.figure()
# for species in set(data['species_short']):
#     plt.plot(data[data['species_short'] == species]['Height'],data[data['species_short'] == species]['crown_base_height'],'o')

# Leaf area index (m2)
data['LAI'] = np.nan
data['foliage_biomass'] = np.nan
for species in set(data['species_short']):
    ix=(data['species_short'] == species)
    y, data.loc[ix,'LAI'] = aleksis_combination(d=data['Diameter'][ix].values,
                                                h=data['Height'][ix].values,
                                                ch=data['crown_base_height'][ix].values,
                                                species=species)
    data.loc[ix,'foliage_biomass']=y[3]

# plt.figure()
# for species in set(data['species_short']):
#     plt.plot(data[data['species_short'] == species]['Diameter'],data[data['species_short'] == species]['LAI'],'o')

# plt.figure()
# for species in set(data['species_short']):
#     plt.plot(data[data['species_short'] == species]['foliage_biomass'],data[data['species_short'] == species]['BiomassFoliage'],'o')

data.loc[len(data)] = 0
data.loc[len(data)-1,'StandId'] = 1042
data.loc[len(data)-1,'species_short'] = 'birch'
data.loc[len(data)] = 0
data.loc[len(data)-1,'StandId'] = 244
data.loc[len(data)-1,'species_short'] = 'birch'

z = np.linspace(0, 32, 101)

# Crown leaf mass profiles
data['lad'] = data['LAI'] * data.apply(lambda row: crown_biomass_distr(species=row['species_short'], z=z,
                                                         htop=row['Height'],
                                                         hbase=row['crown_base_height'])[0],
                         axis=1)

data_gr = data[['StandId','species_short','LAI','basal_area']].groupby(['StandId','species_short']).sum()

# plot area
data_gr['LAI'] = data_gr['LAI'] / (np.pi*10**2)
data_gr['basal_area'] = data_gr['basal_area'] / (np.pi*10**2 / 10000)

data_gr['lad'] = len(data_gr['LAI'])*[np.ones(101)]

plt.figure(figsize=(20,12))
i=1
for plot in set(data['StandId']):
    if i == 1:
        ax=plt.subplot(2,5,i)
    else:
        plt.subplot(2,5,i,sharex=ax)
    plt.title(plot)
    for species in set(data[data['StandId']==plot]['species_short']):
        data_gr.loc[(plot, species),'lad'] = np.sum(
            np.stack(data[(data['StandId']==plot) & (data['species_short']==species)]['lad'].to_numpy()),axis=0) / (np.pi*10**2)  # plot area
        plt.plot(data_gr.loc[(plot, species),'lad'], z, label='%s, %.2f m$^2$m$^{-2}$, %.2f m$^2$ha$^{-1}$' % (species, data_gr.loc[(plot, species),'LAI'], data_gr.loc[(plot, species),'basal_area']))
    lad_total = np.sum(
            np.stack(data[(data['StandId']==plot)]['lad'].to_numpy()),axis=0) / (np.pi*10**2)
    plt.plot(lad_total, z, ':k', label='total, %.2f m$^2$m$^{-2}$, %.2f m$^2$ha$^{-1}$' % (data_gr.loc[(plot),'LAI'].sum(),data_gr.loc[(plot),'basal_area'].sum()))
    plt.legend()
    if i == 1 or i == 6:
        plt.ylabel('Height (m)')
    if i >= 6:
        plt.xlabel('Leaf area density (m$^2$m$^{-3}$)')
    i += 1

plt.tight_layout()

data_gr['lad_normed'] = data_gr['lad'] / (data_gr['LAI'] + eps)

data_gr = data_gr.unstack(level=0)
