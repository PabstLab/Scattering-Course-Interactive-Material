#!/usr/bin/python

"""
Scattering amplitude and intensity functions for cylinders

"""

import math
import numpy as np
from scipy.special import spherical_jn, jv

########################################################################################
"""Scattering amplitude of a cylinder with scalar t"""
def cylinder_amp_scalar (q,R,H,t):
    with np.errstate(divide='ignore', invalid='ignore'):
        Beta = np.where( t == 1, 1, 2 * jv(1, q*R*np.sqrt(1-t**2)) / (q*R*np.sqrt(1-t**2)) )
    return (np.pi*R**2*H) * spherical_jn(0, q*t*H/2.) * Beta
    
########################################################################################
"""Scattering amplitude of a cylinder"""
def cylinder_amp (q,R,H,t):
    with np.errstate(divide='ignore', invalid='ignore'):
        Beta = np.where( t == 1, 1, 2 * jv(1, q[:,None]*R*np.sqrt(1-t[None,:]**2)) / (q[:,None]*R*np.sqrt(1-t[None,:]**2)) )
    return (np.pi*R**2*H) * spherical_jn(0, q[:,None]*t[None,:]*H/2.) * Beta

########################################################################################
"""Scattering intensity of a cylinder"""
def cylinder_int (q,Dr,R,H):

    ### t-array for cylinder orientation
    Nt = 90
    t = np.linspace(0, 1, Nt)
    I_partial = np.empty([q.size,Nt],dtype=float)
    I_partial = cylinder_amp(q,R,H,t)**2

    ### Calculating orientation average
    I = np.trapz(I_partial, t, 1/(Nt-1), axis=1)
    
    return Dr**2 * I
