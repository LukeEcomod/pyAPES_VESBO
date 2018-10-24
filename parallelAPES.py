#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 11:07:09 2018

@author: ajkieloaho
"""

import os
from multiprocessing import Process, JoinableQueue, cpu_count, Manager
import multiprocessing
from copy import deepcopy
from pyAPES import initialize_netcdf, _write_ncf
from pyAPES import Model

import time

import logging

results_queue = JoinableQueue()

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

def _poolworker(task):
    """
    Args:
        task (dict)
    """

    print("about to create a model...")

    model = Model(task['general'],
                  task['canopy'],
                  task['soil'],
                  task['forcing'],
                  task['nsim'])
    print(model.Nsim)
    print('about to run the model!')
    result = model.run()

    print("work with simulation: {} is done".format(task['nsim']))

    return (task['nsim'], result)


def _result_handler(result):
    """
    Args:
        result (dict)
    """
    print("result_handler")
    results_queue.put({'Nsim': result[0], 'data': deepcopy(result[1])})


def do_work(in_queue, out_queue):
    """ https://stackoverflow.com/questions/11996632...
    .../multiprocessing-in-python-while-limiting-the-number-of-running-processes
    """

    raise NotImplementedError

# logging to a single file from multiple processes
# https://docs.python.org/dev/howto/logging-cookbook.html#logging-to-a-single-file-from-multiple-processes

import logging
import logging.handlers
import Threading

class QueueHandler(logging.Handler):
    """
    This handler sends events to a queue. Typically, it would be used together
    with a multiprocessing Queue to centralise logging to file in one process
    (in a multi-process application), so as to avoid file write contention
    between processes.
    This code is new in Python 3.2, but this class can be copy pasted into
    user code for use with earlier Python versions.
    """

    def __init__(self, queue):
        """
        Initialise an instance, using the passed queue.
        """
        logging.Handler.__init__(self)
        self.queue = queue

    def enqueue(self, record):
        """
        Enqueue a record.
        The base implementation uses put_nowait. You may want to override
        this method if you want to use blocking, timeouts or custom queue
        implementations.
        """
        self.queue.put_nowait(record)

    def prepare(self, record):
        """
        Prepares a record for queueing. The object returned by this
        method is enqueued.
        
        The base implementation formats the record to merge the message
        and arguments, and removes unpickleable items from the record
        in-place.
        
        You might want to override this method if you want to convert
        the record to a dict or JSON string, or send a modified copy
        of the record while leaving the original intact.
        """
        # The format operation gets traceback text into record.exc_text
        # (if there's exception data), and also puts the message into
        # record.message. We can then use this to replace the original
        # msg + args, as these might be unpickleable. We also zap the
        # exc_info attribute, as it's no longer needed and, if not None,
        # will typically not be pickleable.
        self.format(record)
        record.msg = record.message
        record.args = None
        record.exc_info = None
        return record

    def emit(self, record):
        """
        Emit a record.
        Writes the LogRecord to the queue, preparing it first.
        """
        try:
            self.enqueue(self.prepare(record))
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


def logger_thread(queue):
    while True:
        record = queue.get()
        if record is None:
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)


# !!! THIS WILL BE NEW WORKER PROCESS !!!
def worker_process(task, queue):
    qh = logging.handlers.QueueHandler(queue)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.addHandler(qh)

    print("about to create a model...")

    model = Model(task['general'],
                  task['canopy'],
                  task['soil'],
                  task['forcing'],
                  task['nsim'])

    print('about to run the model!')

    result = model.run()

    print("work with simulation: {} is done".format(task['nsim']))

    return (task['nsim'], result)



def driver(create_ncf=True, dbhfile='letto2014.txt', logging_configuration=None):
    """
    """
    from parameters.canopy import get_cpara
    from parameters.sensitivity import parameters, iterate_parameters
    from parameters.general import gpara
    from parameters.soil import spara
    from tools.iotools import read_forcing

    cpara = get_cpara(dbhfile=dbhfile)

    logging_queue = JoinableQueue()

    logging.config.dictConfig(logging_configuration)
    logging_thread = Threading.Thread(target=logger_thread, args=(logging_queue))
    logging_thread.start()

    Nsim = parameters['count']

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

    # --- POOL OF WORKERS ---
    # does not work because logger is not pickable
    Ncpu = cpu_count()
    pool = multiprocessing.Pool(4)
    for task in tasks:
      pool.apply_async(worker_process, (task, logging_queue), callback=_result_handler)

    pool.close()
    pool.join()

    # --- FREE PROCESSES ---
    # there is risk of too many running processes (fork bomb)
    # needs 'pool' -producer-consumer -pattern

#    processes = []
#    for task in tasks:
#        process = pool.Process(target=_worker, args=(task,))
#        process.start()
#        processes.append(process)
#
#    for process in processes:
#        process.join()


    # end writing process
    results_queue.put({'Nsim': -999, 'data': 'DONE'})
    results_queue.join()

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

