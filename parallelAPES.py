#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 11:07:09 2018

@author: ajkieloaho
"""

import os
from multiprocessing import Process, JoinableQueue, Queue, cpu_count
import multiprocessing
from copy import deepcopy
from pyAPES import initialize_netcdf, _write_ncf
from pyAPES import Model

import time

import logging
import logging.handlers
import Threading


def _result_writer(queue, ncf_param):
    """
    """

    ncf, _ = initialize_netcdf(
        variables=ncf_param['variables'],
        sim=ncf_param['Nsim'],
        soil_nodes=ncf_param['Nsoil_nodes'],
        canopy_nodes=ncf_param['Ncanopy_nodes'],
        plant_nodes=ncf_param['Nplant_types'],
        forcing=ncf_param['forcing'],
        filename=ncf_param['file_name'])

    print(" Queue consumer is listening")

    flag = True
    while flag:

        results = queue.get()

        if results['Nsim'] == -999:
            print('Writing of results is DONE!')
            ncf.close()
            queue.task_done()
            flag = False
        else:
            print(' Processing simulation number: {}'.format(results['Nsim']))
            _write_ncf(nsim=results['Nsim'], results=results['data'], ncf=ncf)
            queue.task_done()

# logging to a single file from multiple processes
# https://docs.python.org/dev/howto/logging-cookbook.html#logging-to-a-single-file-from-multiple-processes


def logger_thread(queue):
    while True:
        record = queue.get()
        if record is None:
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)


def _worker(task_queue, result_queue, logging_queue):
    """
    """
    # setting up logging
    # might be worth of doing inside while loop
    qh = logging.handlers.QueueHandler(logging_queue)
    root = logging.getLogger()

    # !!! root level set should be in configuration dictionary!!!

    root.setLevel(logging.DEBUG)
    root.addHandler(qh)

    # task queue sentinel
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


def driver(create_ncf=True, dbhfile='letto2014.txt', logging_configuration=None):
    """
    """
    from parameters.canopy import get_cpara
    from parameters.sensitivity import parameters, iterate_parameters
    from parameters.general import gpara
    from parameters.soil import spara
    from tools.iotools import read_forcing

    # --- LOGGING ---

    logging_queue = Queue()

    logging.config.dictConfig(logging_configuration)
    logging_thread = Threading.Thread(target=logger_thread, args=(logging_queue))
    logging_thread.start()

    # --- TASKS ---
    # Refactor: parameter handling in separate function

    Nsim = parameters['count']

    cpara = get_cpara(dbhfile=dbhfile)

    default_para = {
            'canopy': cpara,
            'soil': spara,
            'general': gpara
            }

    para_space = [iterate_parameters(parameters, deepcopy(default_para), count) for count in range(Nsim)]

    forcing = read_forcing(gpara['forc_filename'],
                           gpara['start_time'],
                           gpara['end_time'],
                           dt=gpara['dt'])


    model_input = {'forcing': forcing}

    tasks = []
    for k in range(Nsim):
        model_input.update({'nsim': k})
        model_input.update(para_space[k])
        tasks.append(deepcopy(model_input))

    # --- NETCDF4 ---
    results_queue = JoinableQueue()

    timestr = time.strftime('%Y%m%d%H%M')
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


    writing_process = Process(target=_result_writer, args=(results_queue, ncf_param))

    writing_process.deamon = True
    writing_process.start()

    # --- PROCESSES ---
    task_queue = Queue()
    for para in para_space:
        task_queue.put(para)

    N_workers = cpu_count() - 2

    workers = []
    for k in range(N_workers):
        workers.append(Process(target=_worker, args=(task_queue, results_queue, logging_queue)))
        task_queue.put(None)
        workers[k].start()

#    #--- RUNS IN MAIN PROCESS and WRITE IN ANOTHER ---
#    for task in tasks:
#
#        model = Model(task['general'],
#                      task['canopy'],
#                      task['soil'],
#                      task['forcing'],
#                      task['nsim'])
#        result = model.run()
#        results_queue.put({'Nsim': task['nsim'], 'data': deepcopy(result)})
#
#        del result

    # --- CLOSE ---

    # end worker processes
    for w in workers:
        w.join()

    # end writing queue
    results_queue.put({'Nsim': -999, 'data': 'DONE'})
    results_queue.join()

    # end logging queue
    logging_queue.put(None)
    logging_thread.join()


    return ncf_param['output_path']


if __name__ == '__main__':
    import logging.config

    logging_configuration = {
        'version': 1,
#        'disable_existing_loggers': False,
        'formatters': {
                'default': {'format': '%(asctime)s %(levelname)s %(name)s %(message)s'},
                'model': {'format': '%(levelname)s %(name)s %(funcName)s %(message)s'},
                },
        'handlers': {
                'console': {
                        'class' : 'logging.StreamHandler',
                        'formatter': 'model',
                        'level': 'INFO'  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        },
                'file': {
                        'class': 'logging.FileHandler',
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'formatter': 'model',
                        'filename': 'pyAPES.log',
                        'mode': 'w',  # a == append, w == overwrite
                        },
                },
        'loggers': {
                'pyAPES': {
                        'handlers': ['file', 'console'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
#                        'propagate': True,
                        },
                'canopy':{
                        'handlers': ['file', 'console'],
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
#                        'propagate': True,
                        },
                'soil':{
                        'handlers': ['file', 'console'],
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
#                        'propagate': True,
                        },
                },
        'root': {
                'level': 'DEBUG',
                'handlers': ['file', 'console']
                }
        }

    outputfile = driver(create_ncf=True,
                        dbhfile="letto2014.txt",
                        logging_configuration=logging_configuration)

    print(outputfile)

