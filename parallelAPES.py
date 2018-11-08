#!/usr/bin/env python2
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
from multiprocessing import Process, Queue  # , cpu_count
from psutil import cpu_count
from copy import deepcopy

from tools.iotools import initialize_netcdf, write_ncf
from tools.iotools import jsonify
from pyAPES import Model

import time

import logging
import logging.handlers
import logging.config
import threading


def _result_writer(queue, ncf_param):
    """
    Args:
        queue (Queue): results queue
        ncf_param: parameters to initialize NetCDF4 file
    """
    ncf, _ = initialize_netcdf(
        variables=ncf_param['variables'],
        sim=ncf_param['Nsim'],
        soil_nodes=ncf_param['Nsoil_nodes'],
        canopy_nodes=ncf_param['Ncanopy_nodes'],
        plant_nodes=ncf_param['Nplant_types'],
        forcing=ncf_param['forcing'],
        filename=ncf_param['file_name'])

    while True:
        # results is tuple (Nsim, data)
        results = queue.get()

        if results is None:
            ncf.close()
            break

        write_ncf(nsim=results[0], results=results[1], ncf=ncf)


# logging to a single file from multiple processes
# https://docs.python.org/dev/howto/logging-cookbook.html#logging-to-a-single-file-from-multiple-processes


def _logger_listener(queue):
    """
    Args:
        queue (Queue): logging queue
    """

    while True:
        record = queue.get()

        if record is None:
            break

        logger = logging.getLogger(record.name)
        logger.handle(record)


def _worker(task_queue, result_queue, logging_queue):
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
    root.setLevel(logging.DEBUG)
    root.addHandler(qh)
    # --- TASK QUEUE LISTENER ---
    while True:
        task = task_queue.get()

        if task is None:
            break

        model = Model(task['general'],
                      task['canopy'],
                      task['soil'],
                      task['forcing'],
                      task['nsim'])

        result = model.run()

        result_queue.put((task['nsim'], result))

        # can return something if everything went right


def driver(create_ncf=True,
           dbhfile='letto2014.txt',
           logging_configuration=None,
           Ncpu=None):
    """
    Args:
        create_ncf (boolean)
        dbhfile (str)
        logging_configuration (dict)
        ncpu (int)
    """
    from parameters.general import gpara
    from parameters.canopy import get_cpara
    from parameters.soil import get_spara

    from parameters.sensitivity import get_parameters, iterate_parameters

    from tools.iotools import read_forcing

    # from logging.config import dictConfig
    # --- LOGGING ---

    logging_queue = Queue()

    logging.config.dictConfig(logging_configuration)
    logging_thread = threading.Thread(target=_logger_listener,
                                      args=(logging_queue,))
    logging_thread.start()

    logger = logging.getLogger()

    # --- PARAMETERS ---
    # Refactor: parameter handling in separate function
    parameters = get_parameters(name=None)
    Nsim = parameters['count']

    cpara = get_cpara(dbhfile=dbhfile)
    spara = get_spara('organic')

    default_para = {
            'canopy': cpara,
            'soil': spara,
            'general': gpara
            }

    para_space = [iterate_parameters(
        parameters,
        deepcopy(default_para),
        count) for count in range(Nsim)]

    timestr = time.strftime('%Y%m%d%H%M')

    # save parameters into json file
    jsonify(para_space, file_name='results/'+timestr+'_pyAPES_parameters.json')

    # --- FORCING ---
    forcing = read_forcing(gpara['forc_filename'],
                           gpara['start_time'],
                           gpara['end_time'],
                           dt=gpara['dt'])

    # --- TASKS ---
    tasks = []
    for k, para in enumerate(para_space):
        para.update({'nsim': k,
                     'forcing': forcing})
        tasks.append(deepcopy(para))

    # --- NETCDF4 ---
    writing_queue = Queue()

    filename = timestr + '_pyAPES_results.nc'

    pyAPES_folder = os.getcwd()
    filepath = os.path.join(pyAPES_folder, "results", filename)

    ncf_param = {'variables': gpara['variables'],
                 'Nsim': Nsim,
                 'Nsoil_nodes': len(spara['grid']['dz']),
                 'Ncanopy_nodes': cpara['grid']['Nlayers'],
                 'Nplant_types': len(cpara['planttypes']),
                 'forcing': forcing,
                 'file_name': filename,
                 'output_path': filepath,
                 }

    writing_thread = threading.Thread(target=_result_writer,
                                      args=(writing_queue, ncf_param,))
    writing_thread.start()

    # --- PROCESSES ---
    task_queue = Queue()
    for para in para_space:
        task_queue.put(para)

    if Ncpu is None:
        Ncpu = cpu_count(logical=False)

    if Nsim > (Ncpu - 1):
        N_workers = Ncpu - 1
    else:
        N_workers = Nsim

    logger.info('Number of worker processes is {}, number of simulations: {}'.format(N_workers, Nsim))

    running_time = time.time()

    workers = []
    for k in range(N_workers):
        workers.append(Process(target=_worker,
                               args=(task_queue,
                                     writing_queue,
                                     logging_queue,
                                     )))
        task_queue.put(None)
        workers[k].start()

    # --- CLOSE ---

    # join worker processes
    for w in workers:
        w.join()

    # end logging queue and join
    logging_queue.put(None)
    logging_thread.join()

    # end writing queue and join
    writing_queue.put(None)
    writing_thread.join()

    logger.info('Running time %.2f seconds' % (time.time() - running_time))
    logger.info('Results are in path: ' + ncf_param['output_path'])

    return ncf_param['output_path']


if __name__ == '__main__':
    import argparse
    from parameters.general import parallel_logging_configuration

    parser = argparse.ArgumentParser()
    parser.add_argument('--cpu', help='number of cpus to be used', type=int)

    args = parser.parse_args()

    if args.cpu:
        outputfile = driver(create_ncf=True,
                            dbhfile='letto2014.txt',
                            logging_configuration=parallel_logging_configuration,
                            Ncpu=args.cpu)
    else:
        outputfile = driver(create_ncf=True,
                            dbhfile="letto2014.txt",
                            logging_configuration=parallel_logging_configuration)

    print(outputfile)
