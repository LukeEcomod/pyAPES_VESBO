#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 11:07:09 2018

@author: ajkieloaho
"""

from multiprocessing import Process, JoinableQueue
from copy import deepcopy
from pyAPES import driver, initialize_netcdf, _write_ncf


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


def _worker(queue, task):
    """
    Args:
        queue (JoinableQueue): for results to write netCDF4 file
        task (pyAPES.Model): pyAPES model instance
    """
    results = task.run()
    print("work with simulation: {} is done".format(task.Nsim))
    queue.put({'Nsim': task.Nsim, 'data': deepcopy(results)})
    del results


def drive_parallel(Nsim, tasks, ncf_param):
    result_queue = JoinableQueue()
    writing_process = Process(target=_result_writer, args=(result_queue, ncf_param))

    writing_process.deamon = True
    writing_process.start()

#    for task in tasks:
#
#        results = task.run()
#        result_queue.put({'Nsim': task.Nsim, 'data': deepcopy(results)})
#
#        del results

    processes = []
    for task in tasks:
        process = Process(target=_worker, args={result_queue, task})
        process.start()
        processes.append(process)

    for process in processes:
        process.join()

    result_queue.put({'Nsim': -999, 'data': 'DONE'})
    result_queue.join()

    return 'DONE'


class ParallelAPES(object):
    """
    """

    def __init__(self,
                 Nsim,
                 tasks):
        self.Nsim = Nsim
        self.tasks = tasks

    def run(self, ncf_param):
        """
        """
        result_queue = JoinableQueue()
        writing_process = Process(target=_result_writer, args=(result_queue, ncf_param))

        writing_process.deamon = True
        writing_process.start()

        for task in self.tasks:

            results = task.run()
            result_queue.put({'Nsim': task.Nsim, 'data': deepcopy(results)})

            del results

#        processes = []
#        for task in self.tasks:
#            process = Process(target=_worker, args={task, result_queue})
#            process.start()
#            processes.append(process)
#
#        for process in processes:
#            process.join()

        result_queue.put({'Nsim': -999, 'data': 'DONE'})
        result_queue.join()

        print('DONE')

        return ncf_param['output_path']


if __name__ == '__main__':
    outputfile = driver(create_ncf=True, parallel=True, dbhfile="letto2014.txt")
    print(outputfile)