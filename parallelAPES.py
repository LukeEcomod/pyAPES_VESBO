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
from tools.iotools import jsonify
from pyAPES import Model

import time

import logging
import logging.handlers
import logging.config


def _result_writer(ncf):
    """
    Args:
        queue (Queue): results queue
        ncf_param: parameters to initialize NetCDF4 file
    """

    while True:
        # results is tuple (Nsim, data)
        results = writing_queue.get()

        if results is None:
            ncf.close()
            print('results done')
            break

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
            print('logger done')
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

#    from pyAPES_utilities.spinup import soil_spinup

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
            print('worker done')
            break

        root.info("Creating simulation {}".format(task['nsim']))

        try:
#            soil_temperature = soil_spinup(
#                task['general'],
#                task['canopy'],
#                task['soil'],
#                task['forcing'],
#                moss_type,
#                0.01
#            )
#
#            task['soil']['heat_model']['initial_condition']['temperature'] = soil_temperature

            model = Model(
                task['general'],
                task['canopy'],
                task['soil'],
                task['forcing'],
                task['nsim']
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
        task_queue (queue): simulation parameter queue
        logging_queue (queue): queue for logging events
        writing_queue (queue): queue for writing netcdf results
        N_workers (int): number of worker processes
    """

    # --- PROCESSES ---
    running_time = time.time()

    workers = []
    for k in range(N_workers):
        workers.append(
            Process(
                target=_worker,
        #       args=(
        #           task_queue,
        #           writing_queue,
        #           logging_queue,
        #        )
            )
        )

        task_queue.put(None)
        workers[k].start()
    #pool = Pool(N_workers, _worker,)
    #_ = pool.apply_async(_worker, ())

    # --- NETCDF4 ---
    ncf, _ = initialize_netcdf(
        variables=ncf_params['variables'],
        sim=ncf_params['Nsim'],
        soil_nodes=ncf_params['Nsoil_nodes'],
        canopy_nodes=ncf_params['Ncanopy_nodes'],
        plant_nodes=ncf_params['Nplant_types'],
        forcing=ncf_params['forcing'],
        filename=ncf_params['file_name'])

    writing_thread = Thread(
        target=_result_writer,
        args=(ncf,)
        #args=(writing_queue, ncf,)
    )

    writing_thread.start()

    # --- LOGGING ---
    logging.config.dictConfig(logging_configuration)

    logging_thread = Thread(
        target=_logger_listener,
        #args=(logging_queue,)
    )

    logging_thread.start()

    # --- USER INFO ---

    logger = logging.getLogger()
    logger.info('Number of worker processes is {}, number of simulations: {}'.format(N_workers, Nsim))

    # --- CLOSE ---

    # join worker processes
    for w in workers:
        print("working")
        w.join()

    logger.info('Running time %.2f seconds' % (time.time() - running_time))

    # end logging queue and join
    logging_queue.put_nowait(None)
    logging_thread.join()

    # end writing queue and join
    writing_queue.put_nowait(None)
    writing_thread.join()

    logger.info('Results are in path: ' + ncf_params['output_path'])

    return ncf_params['output_path']


def get_tasks(scenario='all'):
    """ Creates parameters space for tasks
    """
#    from parameters.ebal import sensitivity_sampling
#    from parameters.parametersets import parameters, iterate_parameters
    from parameters.parametersets import get_parameters, iterate_parameters
    from copy import deepcopy as copy
    parametersets = get_parameters(scenario)

    from parameters.general import gpara
    from parameters.canopy import cpara
    from parameters.soil import spara

    from tools.iotools import read_forcing
#
#    if predefined is None:
#        save_samples = True
#    else:
#        save_samples = False

    default_params = {
            'general': gpara,
            'canopy': cpara,
            'soil': spara
            }

    Nsim = parametersets['count']
    para_space = [iterate_parameters(parametersets, copy(default_params), count) for count in range(Nsim)]

#    para_space, filename = sensitivity_sampling(
#        moss_type=moss_type,
#        soil_type=soil_type,
#        optimal_trajectories=optimal_trajectories,
#        num_levels=num_levels,
#        save_samples=save_samples,
#        use_predefined=predefined
#    )

    forcing = read_forcing(
        gpara['forc_filename'],
        gpara['start_time'],
        gpara['end_time'],
        dt=gpara['dt']
    )

    # save parameters into json file
    #jsonify(para_space, file_name='results/'+timestr+'_pyAPES_parameters.json')

    task_list = []
    for k, para in enumerate(para_space):
        para.update({
            'nsim': k,
            'forcing': forcing
        })

        task_list.append(deepcopy(para))

    filename = time.strftime('%Y%m%d%H%M') + '_pyAPES_results.nc'

    pyAPES_folder = os.getcwd()
    filepath = os.path.join(pyAPES_folder, "results", filename)

    ncf = {
        'variables': gpara['variables'],
        'Nsim': len(task_list),
        'Nsoil_nodes': len(spara['grid']['dz']),
        'Ncanopy_nodes': cpara['grid']['Nlayers'],
        'Nplant_types': len(cpara['planttypes']),
        'forcing': forcing,
        'file_name': filename,
        'output_path': filepath,
    }

    return task_list, ncf


if __name__ == '__main__':
    import argparse
    from parameters.general import parallel_logging_configuration
#    from parameters.parametersets import parameters
    #mp.set_start_method('spawn')

    parser = argparse.ArgumentParser()
    parser.add_argument('--cpu', help='number of cpus to be used', type=int)
    parser.add_argument('--scenario', help='scenario name (all, control, partial, clearcut)', type=str)
#    parser.add_argument('--moss_type', help='hylocomium or sphagnum', type=str)
#    parser.add_argument('--soil_type', help='organic or mineral', type=str)
#    parser.add_argument('--trajectories', help='number of optimal trajectories to be used', type=int)
#    parser.add_argument('--levels', help='number of levels to be used', type=int)
#    parser.add_argument('--samples', help='predefined sample space', type=str)
#
    args = parser.parse_args()
#
#    timestr = time.strftime('%Y%m%d%H%M')

    # --- Queues ---
    manager = mp.Manager()
    logging_queue = Queue()
    writing_queue = Queue()
    task_queue = Queue()

    # --- TASKS ---
#    moss_type = args.moss_type
#    soil_type = args.soil_type
#    optimal_trajectories = args.trajectories
#    num_levels = args.levels
#    predefined = args.samples

    tasks, ncf_params = get_tasks(args.scenario)

    Nsim = len(tasks)

    for para in tasks:
        task_queue.put(deepcopy(para))

    # --- Number of workers ---
    Ncpu = args.cpu


    if Ncpu is None:
        #Ncpu = cpu_count(logical=False)
        Ncpu = 1

#   if Nsim > (Ncpu - 1):
#       N_workers = Ncpu - 1
#   else:
    N_workers = Ncpu - 1

#    parallel_logging_configuration['handlers']['parallelAPES_file']['filename'] = 'sensitivity_'+moss_type+'_'+soil_type+'.log'

    # --- DRIVER CALL ---
    outputfile = driver(
        ncf_params=ncf_params,
        logging_configuration=parallel_logging_configuration,
#        task_queue=task_queue,
#        logging_queue=logging_queue,
#        writing_queue=writing_queue,
        N_workers=N_workers)

    print(outputfile)