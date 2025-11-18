#!/usr/bin/python

"""
Scattering amplitude and intensity functions for monodisperse and polydisperse spheres

"""

import math
import numpy as np

########################################################################################
"""Scattering amplitude of a sphere"""
def sphere_amp (x,R):
    y = x*R
    #A = 0
    #if y != 0 :   A = 3 * ( np.sin(y) - y*np.cos(y) ) / y**3
    #elif y == 0 : A = 1
    return np.where ( y==0, 1, 3 * ( np.sin(y) - y*np.cos(y) ) / y**3 )

########################################################################################
"""Scaled scattering amplitude of a sphere"""
def sphere_amp_scaled (y):
    return 3 * ( np.sin(y) - y*np.cos(y) ) / y**3

########################################################################################
"""Scattering intensity of a sphere"""
def sphere_int (x,R,Dr):
    return Dr**2 * (4/3*np.pi*R**3)**2 * sphere_amp (x,R)**2

########################################################################################
"""Normal (Gaussian) probability density function (PDF)"""
def PDF_norm(x, mu, sig) :
    return np.exp( -(x-mu)**2 / (2*sig**2) ) /( sig*np.sqrt(2*np.pi) )

########################################################################################
"""Lognormal probability density function (PDF)"""
def PDF_lognorm(x, mu, sig) :
    return np.exp(-(np.log(x)-mu)**2 / (2*sig**2) ) /( x*sig*np.sqrt(2*np.pi))

########################################################################################
"""Scattering intensity of a polydisperse system of spheres"""
"""Normal distribution"""
def sphere_int_normal (x,n,Rm,sR,Dr):

    Ipoly = np.empty([x.size])
    if sR == 0 :
        Ipoly = sphere_int(x,Rm,Dr)
    else:    
        N = 30
        R = np.linspace(Rm-3*sR, Rm+3*sR, N)
        pdf = PDF_norm(R, Rm, sR)
    
        I = sphere_int(x[:,None],R[None,:],Dr)
    
        I[:,] *= pdf
        Ipoly = np.trapz(I, R, 6*sR/N, axis=1)

    return n * Ipoly
########################################################################################
"""Scattering intensity of a polydisperse system of spheres"""
"""Normal distribution"""
"""for 2D plots"""
def sphere_int_normal_2D (q,n,Rm,sR,Dr):

    Ipoly = np.empty_like(q)
    if sR == 0 :
        Ipoly = sphere_int(q,Rm,Dr)
    else:    
        N = 30
        R = np.linspace(Rm-3*sR, Rm+3*sR, N)
        pdf = PDF_norm(R, Rm, sR)
      
        I = sphere_int(q[:,:,None],R[None,:],Dr)
       
        I[:,:,] *= pdf
        Ipoly = np.trapz(I, R, 6*sR/N, axis=2)

    return n * Ipoly
    
########################################################################################
"""Scattering intensity of a polydisperse system of spheres"""
"""Log-Normal distribution"""
def sphere_int_lognorm (x,n,m,s,Dr):

    Rm = math.exp(m+s**2/2)
    sR = math.sqrt( ( math.exp(s**2) -1 ) * math.exp(2*m+s**2) )
   
    N = 60
    R = np.linspace(1e-12, (Rm+8*sR), N)
    pdf = PDF_lognorm(R, m, s)

    I = sphere_int(x[:,None],R[None,:],Dr)

    Ipoly = np.empty([x.size])
    I[:,] *= pdf
    Ipoly = np.trapz(I, R, (Rm+8*sR)/N, axis=1)

    return n * Ipoly

########################################################################################
"""Scattering intensity of a spherical core-shell"""
def sphere_core_shell_int (x,R0,DR,rs,r0,r1):

    P0 = (r0-r1) * (4/3*np.pi*R0**3) * sphere_amp(x,R0)
    P1 = (r1-rs) * (4/3*np.pi*(R0+DR)**3) * sphere_amp(x,(R0+DR))
    
    return (P0+P1)**2

########################################################################################
"""Scattering intensity of a spherical core-shell"""
"""Normal distribution"""
def sphere_core_shell_int_normal (x,n,R0,sigma,DR,rs,r0,r1):

    P0 = (r0-r1) * (4/3*np.pi*R0**3) * sphere_amp(x,R0)
    P1 = (r1-rs) * (4/3*np.pi*(R0+DR)**3) * sphere_amp(x,(R0+DR))

    Ipoly = np.empty([x.size])
    if sigma != 0 :

        N = 40
        R = np.linspace(R0-3*sigma, R0+3*sigma, N)
        pdf = PDF_norm(R, R0, sigma)
    
        I = sphere_core_shell_int(x[:,None],R[None,:],DR,rs,r0,r1)
   
        I[:,] *= pdf
        Ipoly = np.trapz(I, R, 6*sigma/N, axis=1)

    elif sigma == 0:
        Ipoly = sphere_core_shell_int(x,R0,DR,rs,r0,r1)

    return n * Ipoly
    
########################################################################################

########################################################################################
"""Scattering amplitude of a spherical 3xcore-shell"""
def sphere_core_3xshell_amp (x,R0,d_h,d_t,r_h,r_t,r_s):

    P0 = (r_s-r_h) * (4/3*np.pi*R0**3) * sphere_amp(x,R0)
    P1 = (r_h-r_t) * (4/3*np.pi*(R0+d_h)**3) * sphere_amp(x,(R0+d_h))
    P2 = (r_t-r_h) * (4/3*np.pi*(R0+d_h+2*d_t)**3) * sphere_amp(x,(R0+d_h+2*d_t))
    P3 = (r_h-r_s) * (4/3*np.pi*(R0+d_h+2*d_t+d_h)**3) * sphere_amp(x,(R0+d_h+2*d_t+d_h))
    
    return (P0+P1+P2+P3)
    
########################################################################################
"""Scattering intensity of a spherical 3xcore-shell"""
def sphere_core_3xshell_int (x,R0,d_h,d_t,r_h,r_t,r_s):
    
    return sphere_core_3xshell_amp(x,R0,d_h,d_t,r_h,r_t,r_s)**2

########################################################################################
"""Scattering intensity of a 3xspherical core-shell"""
"""Normal distribution"""
def sphere_core_3xshell_int_normal (x,n,R0,sigma,d_h,d_t,r_h,r_t,r_s):

    Ipoly = np.empty([x.size])
    if sigma != 0 :

        N = 40
        R = np.linspace(R0-3*sigma, R0+3*sigma, N)
        pdf = PDF_norm(R, R0, sigma)
    
        I = sphere_core_3xshell_int(x[:,None],R[None,:],d_h,d_t,r_h,r_t,r_s)
   
        I[:,] *= pdf
        Ipoly = np.trapz(I, R, 6*sigma/N, axis=1)

    elif sigma == 0:
        Ipoly = sphere_core_3xshell_int(x,R0,d_h,d_t,r_h,r_t,r_s)

    return n * Ipoly
    
########################################################################################
"""Scattering amplitude of a spherical 5xcore-shell"""
def sphere_core_5xshell_amp (x,R0,d_h,d_t,d_m,r_h,r_t,r_m,r_s):

    P0 = (r_s-r_h) * (4/3*np.pi*R0**3) * sphere_amp(x,R0)
    P1 = (r_h-r_t) * (4/3*np.pi*(R0+d_h)**3) * sphere_amp(x,(R0+d_h))
    P2 = (r_t-r_m) * (4/3*np.pi*(R0+d_h+d_t)**3) * sphere_amp(x,(R0+d_h+d_t))
    P3 = (r_m-r_t) * (4/3*np.pi*(R0+d_h+d_t+d_m)**3) * sphere_amp(x,(R0+d_h+d_t+d_m))
    P4 = (r_t-r_h) * (4/3*np.pi*(R0+d_h+d_t+d_m+d_t)**3) * sphere_amp(x,(R0+d_h+d_t+d_m+d_t))
    P5 = (r_h-r_s) * (4/3*np.pi*(R0+d_h+d_t+d_m+d_t+d_h)**3) * sphere_amp(x,(R0+d_h+d_t+d_m+d_t+d_h))
    
    return (P0+P1+P2+P3+P4+P5)
    
########################################################################################
"""Scattering intensity of a spherical 3xcore-shell"""
def sphere_core_5xshell_int (x,R0,d_h,d_t,d_m,r_h,r_t,r_m,r_s):
    
    return sphere_core_5xshell_amp(x,R0,d_h,d_t,d_m,r_h,r_t,r_m,r_s)**2

########################################################################################
"""Scattering intensity of a 5xspherical core-shell"""
"""Normal distribution"""
def sphere_core_5xshell_int_normal (x,n,R0,sigma,d_h,d_t,d_m,r_h,r_t,r_m,r_s):

    Ipoly = np.empty([x.size])
    if sigma != 0 :

        N = 40
        R = np.linspace(R0-3*sigma, R0+3*sigma, N)
        pdf = PDF_norm(R, R0, sigma)
    
        I = sphere_core_5xshell_int(x[:,None],R[None,:],d_h,d_t,d_m,r_h,r_t,r_m,r_s)
   
        I[:,] *= pdf
        Ipoly = np.trapz(I, R, 6*sigma/N, axis=1)

    elif sigma == 0:
        Ipoly = sphere_core_5xshell_int(x,R0,d_h,d_t,d_m,r_h,r_t,r_m,r_s)

    return n * Ipoly
    
########################################################################################



