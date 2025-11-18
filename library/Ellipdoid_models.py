#!/usr/bin/python

"""
Scattering amplitude and intensity functions for ellipdoids

"""

import math
import numpy as np
from scipy.special import spherical_jn, jv

"""Size of orientation array"""
Nt = 90
"""Size of PDF array"""
Np = 30

########################################################################################
"""Normal (Gaussian) probability density function (PDF)"""
def PDF_norm(x, mu, sig) :
    return np.exp( -(x-mu)**2 / (2*sig**2) ) /( sig*np.sqrt(2*np.pi) )

########################################################################################
"""Scattering amplitude of a ellipsoid"""
def ellipsoid_amp (q,R,e,t):

    u = np.sqrt( (R*t)**2 + (R*e*np.sqrt(1-t**2))**2 )
    A = np.empty([q.size, Nt], dtype=float)
    A = (4/3*np.pi*R**3*e) * 3 * ( np.sin(q[:,None]*u[None,:]) - q[:,None]*u[None,:]*np.cos(q[:,None]*u[None,:]) ) / (q[:,None]*u[None,:])**3
    return A

########################################################################################
"""Scattering intensity of a ellipsoid"""
def ellipsoid_int (n, q, Dr, R, e, con):

    ### t-array for cylinder orientation
    t = np.linspace(0, 1, Nt, endpoint=False)
    I_partial = np.empty([q.size, Nt], dtype=float)
    I_partial = ellipsoid_amp(q, R, e,t)**2

    ### Calculating orientation average
    I = np.trapz(I_partial, t, 1/(Nt-1), axis=1)
    
    return n * Dr**2 * I + con

########################################################################################
"""Scattering intensity of a polydisperse system of ellipsoid"""
"""Normal distribution"""
def ellipsoid_int_normal (n, q, Dr, R, rel, e, con):

    A = np.empty([q.size, Np, Nt], dtype=float)

    R_array = np.linspace(R-3*R*rel, R+3*R*rel, Np, endpoint=False)
    t = np.linspace(0, 1, Nt, endpoint=False)

    u = np.empty([Np, Nt], dtype=float)
    u = np.sqrt( (R_array[:,None]*t[None,:])**2 + (R_array[:,None]*e*np.sqrt(1-t[None,:]**2))**2 )

    A = (4./3.*np.pi*R_array[None,:,None]**3*e) * 3 * ( np.sin(q[:,None,None]*u[None,:,:]) - q[:,None,None]*u[None,:,:]*np.cos(q[:,None,None]*u[None,:,:]) ) / (q[:,None,None]*u[None,:,:])**3
    
    I_full = A**2 * PDF_norm(R_array[None,:,None], R, R*rel)
    I_poly = np.trapz(I_full, t, 1/(Nt-1), axis=2)
    I = np.trapz(I_poly, R_array, 6*R*rel/Np, axis=1)

    return n * Dr**2 * I + con

########################################################################################
"""Scattering intensity of a polydisperse system of core-shell ellipsoid"""
"""Normal distribution on R"""
def ellipsoid_coreshell_int_normal (n, q, r0, r1, rs, R0, DR, rel, e, con):

    A = np.empty([q.size, Np, Nt], dtype=float)
    A0 = np.empty_like(A)
    A1 = np.empty_like(A)

    R_array = np.linspace(R0-3*R0*rel, R0+3*R0*rel, Np, endpoint=False)
    t = np.linspace(0, 1, Nt, endpoint=False)

    u0 = np.empty([Np, Nt], dtype=float)
    u0 = np.sqrt( (R_array[:,None]*t[None,:])**2 + (R_array[:,None]*e*np.sqrt(1-t[None,:]**2))**2 )
    u1 = np.empty([Np, Nt], dtype=float)
    u1 = np.sqrt( ((R_array[:,None]+DR)*t[None,:])**2 + ((R_array[:,None]+DR)*e*np.sqrt(1-t[None,:]**2))**2 )

    A0 = (4./3.*np.pi*R_array[None,:,None]**3*e) * 3 * ( np.sin(q[:,None,None]*u0[None,:,:]) - q[:,None,None]*u0[None,:,:]*np.cos(q[:,None,None]*u0[None,:,:]) ) / (q[:,None,None]*u0[None,:,:])**3
    A1 = (4./3.*np.pi*R_array[None,:,None]**3*e) * 3 * ( np.sin(q[:,None,None]*u1[None,:,:]) - q[:,None,None]*u1[None,:,:]*np.cos(q[:,None,None]*u1[None,:,:]) ) / (q[:,None,None]*u1[None,:,:])**3

    A = (r1-rs)*A1 + (r0-r1)*A0

    I_full = A**2 * PDF_norm(R_array[None,:,None], R0, R0*rel)
    I_poly = np.trapz(I_full, t, 1/(Nt-1), axis=2)
    I = np.trapz(I_poly, R_array, 6*R0*rel/Np, axis=1)

    return n * I + con
