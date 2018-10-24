#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:23:50 2018

Note:
    migrated to python3
    - no changes

@author: ajkieloaho
"""

import numpy as np

def solver_array(func, init, time, method):
    """ Solves system of ODEs, with m equations and m unknows

    examples of func:

        def problem1(u, t):
            r = [u[1], -u[0]]
            return np.asarray(r)

        def problem3(u, t):
            dudt[0] = u[1]
            dudt[1] =-u[0]
            return dudt

        NOTE:
            Every time the function problem1 and problem2 is called
            (and that happens twice at each time level!),
            a new array must be made from a list.
            This can be avoided by implementing a class.
            In the class, a numpy array is allocated
            for the right-hand side and reused this in subsequent calls:

        class Problem3:
            def __init__(self):
                # Allocate an array for dudt for efficiency
                self.dudt = np.zeros(2)

            def __call__(self, u, t):
                self.dudt[0] = u[1]
                self.dudt[1] = -u[0]
                return self.dudt

    Args:
        func (obj): function that defines problem to be solved
        init (float/list): initial value
        method (obj): numerical method to advance u one time step.
        time (array): an array of time points.

    Returns:
        u (array): integration results (du/dt)
        t (array): array of time points

    """

    time = np.asarray(time)
    N = len(time) - 1

    if isinstance(init, (float, int)):
        init = [init]  # wrap in list, which then will be array

    init = np.asarray(init)

    if not isinstance(func(init, 0), np.ndarray):
        raise TypeError('func (%s) must return numpy array' % func.__name__)

    u = np.zeros((N + 1, len(init)))
    u[0] = init[:]

    for idx in range(N):
        u[idx + 1] = method(u, idx, time, func)

    return u, time


def RK2_array(u, idx, time, func):
    """ 2nd order Runge-Kutta

    Under construction

    ATTENTION:
        This does not work at the moment

        !!!check MEANING OF dt!!!
        SHOULD USE TIME STEP BUT CONSTRUCTED FOR TIME

    Args:
        u (array): integration results (du/dt) from previous steps
        idx (int): time point in the u array
        time (array): an array of time points.
        func (obj): function that defines problem to be solved
    """

    raise NotImplemented

    dt = time[idx + 1] - time[idx]

    K1 = func(u[idx], 0.5 * dt)
    K2 = func(u[idx] + dt * K1, dt)

    new_u = u[idx] + dt * (K1 + K2) / 2.0

#    K1 = dt * func(u[idx], time[idx])
#    K2 = dt * func(u[idx] + 0.5 * K1, 0.5 * time[idx])
#
#    new_u = u[idx] + K2

#    K1 = dt * func(u[idx], dt)
#    K2 = dt * func(u[idx] + 0.5 * K1, 0.5 * dt)

    # h is time interval
    # as we proceed one time step per time then dt is t
    # Heun's method with single corrector
    # Middle point
#    K1 = func(u[idx], 0.5 * dt)
#    K2 = func(u[idx] + 0.5 * dt * K1, dt)
#
#    new_u = u[idx] + 0.5 * (K1 + K2) * dt
#
#    K1 = func(u[idx], 0.5 * dt)
#    K2 = func(u[idx] + dt * K1, dt + 0.5 * dt)
#    new_u = u[idx] + 0.5 * (K1 + K2) * dt

    return new_u

def RK4_array(u, idx, time, func):
    """ 4nd order Runge-Kutta

    Under construction

    ATTENTION:
        This does not work at the moment

        !!!check MEANING OF dt!!!
        SHOULD USE TIME STEP BUT CONSTRUCTED FOR TIME

    Args:
        u (array): integration results (du/dt) from previous steps
        idx (int): time point in the u array
        time (array): an array of time points.
        func (obj): function that defines problem to be solved
    """

    raise NotImplemented

    dt = time[idx + 1] - time[idx]
    dt2 = dt * 0.5

#    K1 = dt * func(u[idx], time[idx])
#    K2 = dt * func(u[idx] + 0.5 * K1, dt2)
#    K3 = dt * func(u[idx] + 0.5 * K2, dt2)
#    K4 = dt * func(u[idx] + K3, time[idx] + dt)

    K1 = dt * func(u[idx], dt)
    K2 = dt * func(u[idx] + 0.5 * K1, dt2)
    K3 = dt * func(u[idx] + 0.5 * K2, dt2)
    K4 = dt * func(u[idx] + K3, dt)

    new_u = u[idx] + (1./6.) * (K1 + 2. * K2 + 2. * K3 + K4)
    return new_u


def ForwardEuler_array(u, idx, time, func):
    """ Forward Euler method

    It works fine if 20 timesteps are used in bryotype heat and water.

    Only problems occurs if something sudden phenomena in forcing occurs.

    Args:
        u (array): integration results (du/dt) from previous steps
        idx (int): time point in the u array
        time (array): an array of time points.
        func (obj): function that defines problem to be solved
    """
    dt = time[idx + 1] - time[idx]
    new_u = u[idx] + dt * func(u[idx], dt)

    return new_u