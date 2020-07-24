""" Created on Thu Mar 26.3.2020

This module contains tools to parametrize pyAPES model
Module contains tools:
    to construct default parameters
    to construct parametersets for scenarios


Reasoning:
    Parameter handling should be separated from driver
    Driver has sole responsible to run the pyAPES model and handle IO of the model

Logging configuration remains inside driver for now

@author: Antti-Jussi Kieloaho
"""
from copy import deepcopy as copy

from .parametersets import get_parameters
from tools.iotools import read_forcing

def get_parameter_list(default_params, modifications):
    """ returns list of modified default parameters and their forcing.

    If using name of scenarios, check available scenarios from parametersets module.
    If using own modifications, follow dictionary structure found in parametersets module.

    Args:
        scenario (str/dict/):
            - name of predefined scenarios or
            - modifications as a dict or
            - list of modification dicts
    Returns:
        param_list (list): list of modified parameters
    """

    if isinstance(modifications, str):
        modifications = get_parameters(modifications)

    if not isinstance(modifications, dict):
        raise ValueError('Give scenario name or modification dict (parameterset)')

    Nsim = modifications['count']  # number of simulations

    # change default parameters according to modifications defined in parametersets
    param_list = [iterate_parameters(modifications, copy(default_params), count) for count in range(Nsim)]

    if len(set([param_list[i]['general']['start_time'] for i in range(Nsim)])) > 1:
        raise ValueError('starting time cannot be modified')

    if len(set([param_list[i]['general']['end_time'] for i in range(Nsim)])) > 1:
        raise ValueError('end time cannot be modified')

    if len(set([param_list[i]['general']['dt'] for i in range(Nsim)])) > 1:
        raise ValueError('dt cannot be modified')

    start_time = param_list[0]['general']['start_time']
    end_time = param_list[0]['general']['end_time']
    dt = param_list[0]['general']['dt']

    for i in range(Nsim):

        forc_filename = param_list[i]['general']['forc_filename']

        forcing = read_forcing(
            forc_filename=forc_filename,
            start_time=start_time,
            end_time=end_time,
            dt=dt
        )

        param_list[i]['forcing'] = forcing
        param_list[i]['nsim'] = i

    return param_list

def iterate_parameters(parameters, default, count):
    """ Goes through recursively senstivity nested parameter dictionary.
    Args:
        paramters (dict): nested dictionary
        default (dict): default parameter dictionary
        count (int): number of scenarios
    """

    for key, value in parameters.items():
        if key == 'count' or key == 'scenario':
            continue
        elif isinstance(value, dict):
            if key in default:
                default[key] = iterate_parameters(value, default[key], count)
            else:
                print(key + ' not in parameters')

        else:
            if isinstance(value, tuple):
                default[key] = value[count]
            else:
                default[key] = value

    return default
