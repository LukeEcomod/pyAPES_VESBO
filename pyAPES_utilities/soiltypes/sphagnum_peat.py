"""
Created on Fri Mar 13

@author: Antti-Jussi Kieloaho

Parameters for Degero peat soil
"""

# depth of layer bottom [m], soil surface at 0.0
zh = [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1., -1.5, -2.0]
N = len(zh)

# TEST
porosity = [0.945, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918]
residual_water_content = [0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098]
pf_alpha = [0.338, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072]
pf_n = [1.402, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371]

Kvsat = [30*8.99E-05, 20*2.98E-05, 10*9.86E-06, 3.27E-06, 1.08E-06, 3.58E-07, 1.19E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07]
Khmult = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat

Khsat = [Kvsat[i] * Khmult[i] for i in range(N)]

soil_properties = {
    'pF': {  # vanGenuchten water retention parameters
        'ThetaS': porosity,
        'ThetaR': residual_water_content,
        'alpha': pf_alpha,  # [cm-1]
        'n': pf_n
    },
    'saturated_conductivity_vertical': Kvsat,  # saturated vertical hydraulic conductivity [m s-1]
    'saturated_conductivity_horizontal': Khsat,  # saturated horizontal hydraulic conductivity [m s-1]
    'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
    'solid_composition': {  # fractions of solid volume [-]
        'organic': [1.0 for i in range(N)],
        'sand': [0.0 for i in range(N)],
        'silt': [0.0 for i in range(N)],
        'clay': [0.0 for i in range(N)]
    },
    'freezing_curve': [0.5 for i in range(N)],  # freezing curve parameter
    'bedrock': {
        'solid_heat_capacity': 2.16e6,  # [J m-3 (solid) K-1]
        'thermal_conductivity': 3.0  # thermal conductivity of non-porous bedrock [W m-1 K-1]
    }
}

# original based on parametrization of sphagnum peat in SpaFHy_peat: Kersti Haahti
sphagnum  = {
    'soil_id': 1.0,
    'z': [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1., -1.5, -2.0],
    'pF': {  # vanGenuchten water retention parameters
        'ThetaS': [0.945, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918],
        'ThetaR': [0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098],
        'alpha': [0.338, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072],
        'n': [1.402, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371]},
    'saturated_conductivity': [30*8.99E-05, 20*2.98E-05, 10*9.86E-06, 3.27E-06, 1.08E-06, 3.58E-07, 1.19E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07],
} 