#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 11:07:09 2018
TODO:
    - dump parameter space to file
    - check if filter/adapeter can be used for configure loggers Formatter
    (There is need for add nsim, process id can be added somwhere directily)
@author: ajkieloaho
"""

import os
import sys
import multiprocessing as mp
from threading import Thread
from multiprocessing import Process, Queue, Pool  # , cpu_count
#from psutil import cpu_count
from copy import deepcopy

from tools.iotools import initialize_netcdf, write_ncf
from pyAPES import Model

import time

import logging
import logging.handlers
import logging.config


def _result_writer(ncf):
    """
    Args:
        ncf: NetCDF4 file handle
    """

    logger = logging.getLogger()
    logger.info("Writer is ready!")

    while True:
        # results is tuple (Nsim, data)
        results = writing_queue.get()

        if results is None:
            ncf.close()
            logger.info("NetCDF4 file is closed. and Writer closes.")
            break

        logger.info("Writing results of simulation {}".format(results[0]))
        write_ncf(nsim=results[0], results=results[1], ncf=ncf)
# logging to a single file from multiple processes
# https://docs.python.org/dev/howto/logging-cookbook.html#logging-to-a-single-file-from-multiple-processes

def _logger_listener():
    """
    Args:
        queue (Queue): logging queue
    """

    while True:
        record = logging_queue.get()

        if record is None:
            # print('logger done')
            break

        logger = logging.getLogger(record.name)
        logger.handle(record)


def _worker():
    """
    Args:
        task_queue (Queue): queue of task initializing parameters
        result_queue (Queue): queue of model calculation results
        logging_queue (Queue): queue for model loggers
    """

    # --- LOGGING ---
    qh = logging.handlers.QueueHandler(logging_queue)
    root = logging.getLogger()

    # !!! root level set should be in configuration dictionary!!!
    root.handlers = []
    root.setLevel(logging.INFO)
    root.addHandler(qh)

    # --- TASK QUEUE LISTENER ---
    while True:
        task = task_queue.get()

        if task is None:
            root.info('Worker done')
            break

        root.info("Creating simulation {}".format(task['nsim']))

        try:
            model = Model(
                dt=task['general']['dt'],
                canopy_para=task['canopy'],
                soil_para=task['soil'],
                forcing=task['forcing'],
                outputs=output_variables['variables'],
                nsim=task['nsim'],
            )

            result = model.run()
            writing_queue.put((task['nsim'], result))

        except:
            message = 'FAILED: simulation {}'.format(task['nsim'])
            root.info(message + '_' + sys.exc_info()[0])
        # can return something if everything went right


def driver(ncf_params,
           logging_configuration,
           N_workers):
    """
    Args:
        ncf_params (dict): netCDF4 parameters
        logging_configuration (dict): parallel logging configuration
        N_workers (int): number of worker processes
    """

    # --- PROCESSES ---
    running_time = time.time()

    workers = []
    for k in range(N_workers):
        workers.append(
            Process(
                target=_worker,
            )
        )

        task_queue.put(None)
        workers[k].start()

    # --- NETCDF4 ---
    ncf, _ = initialize_netcdf(
        variables=ncf_params['variables'],
        sim=ncf_params['Nsim'],
        soil_nodes=ncf_params['Nsoil_nodes'],
        canopy_nodes=ncf_params['Ncanopy_nodes'],
        planttypes=ncf_params['Nplant_types'],
        groundtypes=ncf_params['Nground_types'],
        time_index=ncf_params['time_index'],
        filepath=ncf_params['filepath'],
        filename=ncf_params['filename'])

    writing_thread = Thread(
        target=_result_writer,
        args=(ncf,)
    )

    writing_thread.start()

    # --- LOGGING ---
    logging.config.dictConfig(logging_configuration)

    logging_thread = Thread(
        target=_logger_listener,
    )

    logging_thread.start()

    # --- USER INFO ---

    logger = logging.getLogger()
    logger.info('Number of worker processes is {}, number of simulations: {}'.format(N_workers, ncf_params['Nsim']))

    # --- CLOSE ---

    # join worker processes
    for w in workers:
        w.join()

    logger.info('Worker processes have joined.')
    logger.info('Running time %.2f seconds' % (time.time() - running_time))

    # end logging queue and join
    logging_queue.put_nowait(None)
    logging_thread.join()

    # end writing queue and join
    writing_queue.put_nowait(None)
    writing_thread.join()

    logger.info('Results are in path: ' + ncf_params['filepath'])

    return ncf_params['filepath']

if __name__ == '__main__':
    """
    SL changes for hydetrends
    """
    import argparse

    from parameters.outputs_S1 import parallel_logging_configuration, output_variables
    from parameters.SmearII import gpara, cpara, spara
    from parameters.parametersets_S1 import get_parameter_list_S1
    from parameters.parametersets_S2 import get_parameter_list_S2
    from parameters.parametersets_S3 import get_parameter_list_S3

    
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpu', help='number of cpus to be used', type=int)
    parser.add_argument('--scenario', help='scenario name', type=str)
    parser.add_argument('--year', help='year as forcing', type=int)
    
    args = parser.parse_args()

    # --- Queues ---
    manager = mp.Manager()
    logging_queue = Queue()
    writing_queue = Queue()
    task_queue = Queue()

    # --- TASKS ---
    scen = args.scenario
    year = args.year
    
    # list of parameters
    parameters = {
        'general': gpara,
        'canopy': cpara,
        'soil': spara
        }

    #tasks = get_parameter_list(parameters, scen)
    if scen.upper()== 'S1':
        tasks = get_parameter_list_S1(year)
    elif scen.upper()== 'S2':
        tasks = get_parameter_list_S2('S2', years=[year, year])
    elif scen.upper()== 'S3':
        tasks = get_parameter_list_S3(year)
        
    # ncf parameters
    ncf_params = {
        'variables': output_variables['variables'],
        'Nsim': len(tasks),
        'Nsoil_nodes': len(tasks[0]['soil']['grid']['dz']),
        'Ncanopy_nodes': tasks[0]['canopy']['grid']['Nlayers'],
        'Nplant_types': len(tasks[0]['canopy']['planttypes']),
        'Nground_types': 1,  # 
        'time_index': tasks[0]['forcing'].index,
        #'filename': time.strftime('%Y%m%d%H%M_') + scen + '_pyAPES_results.nc',
        'filename': scen + '_%d.nc' %year,
        'filepath': tasks[0]['general']['results_directory'],
    }

    for para in tasks:
        task_queue.put(deepcopy(para))

    # --- Number of workers ---
    Ncpu = args.cpu

    N_workers = Ncpu - 1

    parallel_logging_configuration['handlers']['parallelAPES_file']['filename'] = time.strftime('%Y%m%d%H%M_') + scen + '.log'

    # --- DRIVER CALL ---
    outputfile = driver(
        ncf_params=ncf_params,
        logging_configuration=parallel_logging_configuration,
        N_workers=N_workers)
    
#    if scen.upper() == 'S1':
#        ncfiles = [f for f in glob.glob(r'results/Hyytiala/S1*.nc')]
#        ds = xarray.open_dataset(r'results/Hyytiala/S1.nc')
    print(outputfile)