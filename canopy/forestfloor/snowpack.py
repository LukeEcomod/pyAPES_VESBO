#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 08:54:17 2018

@author: ajkieloaho
"""

import numpy as np

EPS = np.finfo(float).eps  # machine epsilon

class SnowpackModel(object):
    """ Represents snow cover-soil-atmosphere interactions

    """

    def __init__(self, properties):
        self.properties = properties

    def run(self, dt, forcing):
        """ Calculates one timestep and updates states of SnowpackModel instance
        """

        raise NotImplemented

# EOF