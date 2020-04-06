"""
Created on Tue Mar 24.3.2020

This module contains tools to perform sensititivity analysis, e.g. Morris

Note:
Sensitivity tools are build based on 3rd party library SALib
(https://salib.readthedocs.io/en/latest/)

Structure of this module:
1. sensitivity sampling and analysis
2. data handling tools for parameter dictionaries and results
3. plotting
4. input/output

Workflow:
1. run simulations
    - first check right amount of optimal trajectories
    - use more than 2 levels in final simulations
2. read results
    - check that all the results are meaningfull and do not contain nans
    - aggregate results (later use roll to aggragate)
3. read samples
4. analyse results using samples
    - use aggregated results
    - visualize


TODO:
- redo sensitivity sampling as there is now fixed soil type

@author: Antti-Jussi Kieloaho
"""

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import xarray as xr

import json
from copy import deepcopy
import time

from SALib.analyze.morris import analyze as mrr_analyze
from SALib.sample.morris import sample as mrr_sample
from SALib.sample import saltelli

from parameters.parameter_tools import get_parameter_list

#
# --- 1. sensitivity sampling and analysis ---
#

def sensitivity_sampling(
    parameterset,
    moss_ranges,
    trajectories=600,
    optimal_trajectories=10,
    num_levels=4,
    grid_jump=1,
    save_samples=None):
    """ Creates Morris analysis samples for sensitivity analysis from pyAPES parametersets

    1. Parametersets as an argument
    2. ranges for as an argument
    3. soil type as a part of parametersets
    4. check how to convert from SALib works to pyAPES

    Note this is used in parallelAPES driver to get parameters space
    
    If problems with SALib, check Consice API Reference for updates 
    (https://salib.readthedocs.io/en/latest/)
    
    In newer version of SALib, Args of Morris sampling:
        problem (dict): SALib parameters (the problem definition)
        N (int): the number of trajectories to generate
        num_levels (int, default=4): the number of optimal trajectories to sample (from 2 to N)
        local_optimization (bool, default=True): Flag whether to use local optimization

    Args:
        parameterset (tuple(str, dict)): parameterset name and pyAPES parametersets
        moss_ranges (dict): ranges of parameters used in the analysis
        trajectories (int): number of random trajectories
        optimal_trajectories (int): number of optimal trajectories (between 2 and N)
        num_levels (int): number of levels (should be even)
        grid_jump (int): (in older version of SALib)
        save_samples (bool): flag for saving parameter space
    """
    
    scenario = parameterset['scenario']
    moss_type = moss_ranges['moss_type']
    del moss_ranges['moss_type']
    
    ranges = convert_to_salib(moss_ranges)

    # First try newer version of SALib
    try:
        samples = mrr_sample(
            problem=ranges,
            N=trajectories,
            num_levels=num_levels,
            local_optimization=True,
            optimal_trajectories=optimal_trajectories
        )
    # Then try older version
    except TypeError:
        samples = mrr_sample(
            problem=ranges,
            N=trajectories,
            num_levels=num_levels,
            grid_jump=grid_jump,
            local_optimization=True,
            optimal_trajectories=optimal_trajectories
        )
    # If there is problem in something else
    except:
        raise

    if save_samples:
        file_name = scenario + '_{}_trajectories_{}_levels'.format(
        optimal_trajectories,
        num_levels)
        
        jsonify(samples, file_name + 'params')
        jsonify(ranges, file_name + '_ranges')

    # how here onward should be modified?
    new_parameterset = convert_to_pyapes(samples, ranges, moss_type)
    
    # merges new_parameterset to parameterset given as arguments
    new_parameterset = merge(new_parameterset, parameterset, override=True)
    
    # iterate parameterset with default parameters
    parameter_list = get_parameter_list(new_parameterset)
    
    return parameter_list


def analyse_output(output, param_space, ranges, window_size=3, num_levels=4):
    """ Morris analysis of one output variable
    
    Args:
        output (DataArray): output variable
        param_space (list): parameters space containing samples
        ranges (list): definition of input variables
        num_levels (int): number of levels in Morris samples
    """
    
    samples = np.array(param_space)
    
    # resampling to daily values
    output = output.resample(date='24H').mean(dim='date')
    
    # sliding window with DataArray.rolling
    rolling_mean = output.rolling(date=window_size, center=True).mean().dropna('date')
    
    # Morris analysis
    cases = []
    mu_star = []
    mu_star_conf = []
    dates = []
    
    for day in rolling_mean:
        case = mrr_analyze(
            problem=ranges,
            X=samples,
            Y=day.values,
            num_levels=num_levels
        )
        
        cases.append(case)
        mu_star.append(case['mu_star'])
        mu_star_conf.append(case['mu_star_conf'])
        
        date_list = [None] * len(ranges['names'])
        for x in range(len(date_list)):
            date_list[x] = day.date.values
    
        dates.append(np.array(date_list))
    
    results = {
        'mu_star': mu_star,
        'mu_star_conf': mu_star_conf,
        'dates': dates,
        'cases': cases
    }
    
    return results

#
# --- 2. data handling tools for parameter dictionaries and results ---
#

def merge_parametersets(overriding, default, override=True):
    """ Merges overriding parameterset a to parameterset b.
    
    In case of clashing keys parameters are taken from overriding parameterset
    
    Args:
        overriding (dict): parameterset that overrides default
        default (dict): default parameterset
    """

    if not isinstance(overriding, dict) or not isinstance(default, dict):
        raise ValueError('Parametersets have to be dicts!')
        
    new_parameterset = merge(overriding, default, override=override)
    
    return new_parameterset
        
def merge(a, b, path=None, override=False):
    """ Merges dictionary b into a.
    """
    
    if path is None: path = []

    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)], override=override)

            elif a[key] == b[key]:
                pass  # same leaf value

            else:
                if override:
                    b[key] = a[key]

                else:
                    raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))

        else:
            a[key] = b[key]
    
    return a


def convert_to_pyapes(samples, ranges, moss_type,):
    """ Converts SALib samples to pyAPES parameterset.

    Args:
        samples (dict): SALib simulation values
        ranges (dict): SALib simulation ranges
        moss_type (str): 'forest_moss' or 'sphagnum'

    """

    parameterset = {
        'count': len(samples),
        'canopy': {
            'forestfloor': {
                'bottom_layer_types': {
                    moss_type: {
                        #'coverage': 1.0,
                    },
                },
            },
        }
    }

    for sample in samples:
        for idx in range(ranges['num_vars']):
            name = ranges['names'][idx]

            if name in ('albedo_PAR', 'albedo_NIR'):

                if 'optical_properties' not in parameterset['canopy']['forestfloor']['bottom_layer_types'][moss_type]:
                    parameterset['canopy']['forestfloor']['bottom_layer_types'][moss_type]['optical_properties'] = {}
            
                if name not in parameterset['canopy']['forestfloor']['bottom_layer_types'][moss_type]['optical_properties']:
                    parameterset['canopy']['forestfloor']['bottom_layer_types'][moss_type]['optical_properties'][name] = tuple()
                
                parameterset['canopy']['forestfloor']['bottom_layer_types'][moss_type]['optical_properties'][name] += tuple([sample[idx]])
                
            else:

                if name not in parameterset['canopy']['forestfloor']['bottom_layer_types'][moss_type]:
                    parameterset['canopy']['forestfloor']['bottom_layer_types'][moss_type][name] = tuple()
               
                parameterset['canopy']['forestfloor']['bottom_layer_types'][moss_type][name] += tuple([sample[idx]])

    return parameterset


def convert_to_salib(params):
    """ Converts pyAPES parameters to SALib form.

    Parameters needs to be flat (not a nested dictionary)

    Args:
        params (dict): parameter ranges [min, max]
    """

    param_dict = {
        'num_vars': len(params),
        'names': [],
        'bounds': []
    }

    for key, item in params.items():
        param_dict['names'].append(key)
        param_dict['bounds'].append(item)

    return param_dict


def sort_Si(Si, key, sortby='mu_star'):
    """ Sorts the analysis results, i.e. for plotting
    """
    return np.array([Si[key][x] for x in np.argsort(Si[sortby])])


def arrange_measures(results):
    """ Arrange Morris results by output variable
    """
    measures = {}

    for var in variables:
        names = results[var]['names']

        mus = []
        mu_stars = []
        sigmas = []

        for idx in range(len(names)):
            mus.append((names[idx], results[var]['mu'][idx]))
            mu_stars.append((names[idx], results[var]['mu_star'][idx]))
            sigmas.append((names[idx], results[var]['sigma'][idx]))

        measures[var] = {
                'sigma': sigmas,
                'mu': mus,
                'mu_star': mu_stars}

    return measures


def rank_parameters(results, names):
    """ Ranks parameters 
    1. Sort Morris results for each variable
    2. Use index of sorted list as a rank for each parameter
    3. Sum ranks of each parameters together (NOTE: the largest rank is the best one)
    4. Sort ranks and find lowest score
    """

    mustar_ranks = {name: 0 for name in names}
    sigma_ranks = {name: 0 for name in names}

    for var in results.keys():

        sorted_mustar = sort_Si(results[var], 'names', 'mu_star')
        mustar_values = sort_Si(results[var], 'mu_star', 'mu_star')

        sorted_sigma = sort_Si(results[var], 'names', 'sigma')
        sigma_values = sort_Si(results[var], 'sigma', 'sigma')

        for idx in range(len(names)):

            mustar_name = sorted_mustar[idx]
            max_value = max(mustar_values)

            if mustar_values[idx] > 0.0:
                mustar_ranks[mustar_name] += mustar_values[idx] / max_value
            else:
                mustar_ranks[mustar_name] += 0.0

            sigma_name = sorted_sigma[idx]
            max_value = max(sigma_values)

            if sigma_values[idx] > 0.0:
                sigma_ranks[sigma_name] += sigma_values[idx] / max_value
            else:
                sigma_ranks[sigma_name] += 0.0

    return {'mu_star': mustar_ranks, 'sigma': sigma_ranks}

#
# --- 3. plotting ---
#

def horizontal_bar_plot(ax, Si, sortby='mu_star', txt=None, x_label=True):
    """
    """

    names_sorted = sort_Si(Si, 'names', sortby)
    names = []
    for name in names_sorted:
        names.append(param_names[name])

    mu_star_sorted = sort_Si(Si, 'mu_star', sortby)
    mu_star_conf_sorted = sort_Si(Si, 'mu_star_conf', sortby)

    y_pos = np.arange(len(mu_star_sorted))

    ax.barh(y_pos,
            mu_star_sorted,
            xerr=mu_star_conf_sorted,
            align='center',
            ecolor='black')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names)
    if x_label:
        ax.set_xlabel(r'$\mu^\star$') # add unit if available

    ax.set_ylim(min(y_pos)-1, max(y_pos)+1)

    if txt != None:
        ax.text(0.92, 0.05, txt, transform=ax.transAxes)

def covariance_plot(ax, Si, x='mu_star', y='sigma', txt=None, x_label=True):
    """
    """

    out = ax.scatter(Si[x], Si[y], c=u'k', marker=u'o')
    ax.set_ylabel(y)

    if y == 'sigma':
        ax.set_ylim(0.0, )

    if x == 'mu_star':
        ax.set_xlim(0.0, )

    x_axis_bounds = np.array(ax.get_xlim())

    line1, = ax.plot(x_axis_bounds, x_axis_bounds, 'k-')
    line2, = ax.plot(x_axis_bounds, 0.5 * x_axis_bounds, 'k--')
    line3, = ax.plot(x_axis_bounds, 0.1 * x_axis_bounds, 'k-.')

    ax.legend((line1, line2, line3),
              (r'$\sigma / \mu^{\star} = 1.0$',
               r'$\sigma / \mu^{\star} = 0.5$',
               r'$\sigma / \mu^{\star} = 0.1$'),
              loc='best')

    if x_label:
        ax.set_xlabel(x) # add unit if available

    ax.set_ylim(0.0-(0.01 * np.array(ax.get_ylim()[1])), )

    for i, t in enumerate(Si['names']):
        a = abs(Si[y][i])
        b =abs(Si[x][i])
        if b != 0.0:
            ratio = a / b
        else:
            ratio = 0.0

#        if ratio >= .5 and Si[x][i] > 15.:
        ax.annotate(t, (Si[x][i], Si[y][i]))

    if txt != None:
        ax.text(0.05, 0.95, txt, transform=ax.transAxes)

    return out


def wedge_plot(ax, Si, traj):
    """ wedge plot for visual interpretation of the analysis results
    """

    ax.scatter(Si['mu'], Si['sigma'])

    x_axis_bounds = np.array(ax.get_xlim())
    y_axis_bounds = np.array(ax.get_ylim())
    ax.set_ylim(0.0, )

    line1, = ax.plot(2.0*(y_axis_bounds/np.sqrt(traj)), y_axis_bounds, 'k--')
    line2, = ax.plot(-2.0*(y_axis_bounds/np.sqrt(traj)), y_axis_bounds, 'k--')

    ax.set_xlabel(r'$\mu$ ')  # add unit if available
    ax.set_ylim(0.0-(0.01 * np.array(ax.get_ylim()[1])), )

    for i, txt in enumerate(Si['names']):
        ax.text(Si['mu'][i], Si['sigma'][i], txt)

        
#  -----------------
# | 4. input/output |
#  -----------------

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def jsonify(to_be_dumped, file_name, file_path=None):
    """ To save dictionaries in json-format

    Args:
        to_be_dumped (dict): dictionary to be saved in file
        file_name (str): file name without .json
        file_path (str): file path where file will be saved
    """
    from os import path, makedirs, getcwd

    if not file_name or not isinstance(file_name, str):
        print(type(file_name), file_name)
        raise ValueError('Give file name without .json')

    if not file_path:
        file_path = 'results/'
        pyAPES_folder = getcwd()
        file_path = path.join(pyAPES_folder, file_path)

    if not path.exists(file_path):
        makedirs(file_path)
    
    file_name = path.join(file_path, file_name + '.json')
    
    with open(file_name, 'w') as fp:
        json.dump(to_be_dumped, fp, cls=NumpyEncoder)
