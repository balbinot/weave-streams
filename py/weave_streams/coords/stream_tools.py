#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import matplotlib.pyplot as p
from scipy import linalg
from numpy import cos, sin

RAD = np.pi/180.

Ag = np.array([[-0.0548755604, +0.4941094279, -0.8676661490],
               [-0.8734370902, -0.4448296300, -0.1980763734],
               [-0.4838350155, +0.7469822445, +0.4559837762]])
k = 4.74047

def Mrot(alpha_pole, delta_pole, phi1_0):

    """
    Rotation matrix given a pole in Equatorial
    """

    alpha_pole *= RAD
    delta_pole = (90.-delta_pole)*RAD
    phi1_0 *= RAD

    M1 = np.array([[np.cos(alpha_pole),np.sin(alpha_pole),0.],
                   [-np.sin(alpha_pole),np.cos(alpha_pole),0.],
                   [0.,0.,1.]])

    M2 = np.array([[np.cos(delta_pole),0.,-np.sin(delta_pole)],
                   [0.,1.,0.],
                   [np.sin(delta_pole),0.,np.cos(delta_pole)]])

    M3 = np.array([[np.cos(phi1_0),np.sin(phi1_0),0.],
                   [-np.sin(phi1_0),np.cos(phi1_0),0.],
                   [0.,0.,1.]])

    return np.dot(M3,np.dot(M2,M1))

def Mrot3(a,b,c):
    """
        Three angle version (e.g. LMC, Orphan, Sag). Is it equivalent to Mrot?
    """
    a *= RAD
    b *= RAD
    c *= RAD
    M00 = cos(c)*cos(a) - cos(b)*sin(a)*sin(c)
    M01 = cos(c)*sin(a) + cos(b)*sin(a)*sin(c)
    M02 = sin(c)*sin(b)
    M10 = -sin(c)*cos(a) - cos(b)*sin(a)*cos(c)
    M11 = -sin(c)*sin(a) + cos(b)*cos(a)*cos(c)
    M12 = cos(c)*sin(b)
    M20 = sin(b)*sin(a)
    M21 = -sin(b)*cos(a)
    M22 = cos(b)

    M = np.array([[M00, M01, M02], [M10, M11, M12], [M20, M21, M22] ])

    return M



def phi12(alpha, delta, alpha_pole, delta_pole, phi1_0, inv=False):

    vec_radec = np.array([np.cos(alpha*RAD)*np.cos(delta*RAD),
                          np.sin(alpha*RAD)*np.cos(delta*RAD),
                          np.sin(delta*RAD)])

    M = Mrot(alpha_pole,delta_pole,phi1_0)
    if inv:
        M = np.linalg.inv(M)
    vec_phi12 = np.dot(M, vec_radec).T

    phi1 = np.arctan2(vec_phi12[:,1],vec_phi12[:,0])*(1./RAD)
    phi2 = np.arcsin(vec_phi12[:,2])*(1./RAD)

    return [phi1,phi2]

def MM(lon, lat, inv=False):
    cl = np.cos(lon*RAD)
    cb = np.cos(lat*RAD)
    sl = np.sin(lon*RAD)
    sb = np.sin(lat*RAD)
    if inv==True:
        return np.array([[cl*cb, -sl, -cl*sb],
                         [sl*cb,  cl, -sl*sb],
                         [sb   ,   0,  cb]]).T
    else:
        return np.array([[cl*cb, -sl, -cl*sb],
                         [sl*cb,  cl, -sl*sb],
                         [sb   ,   0,  cb]])

def mulb2UVW(l, b, r, mul, mub, vr):

    if vr==0:
        vr = np.zeros_like(mura)
        r = np.ones_like(mura)

    v_vec = np.array([vr, k*r*mul, k*r*mub])
    V_vec = np.dot(v_vec, np.dot(Ag.T, MM(l,b))).T
    return V_vec

def muradec2UVW(ra, dec, r, mura, mudec, vr):

    if vr==0:
        vr = np.zeros_like(mura)
        r = np.ones_like(mura)

    v_vec = np.array([vr, k*r*mura, k*r*mudec])
    V_vec = np.dot(v_vec, np.matmul(Ag.T, MM(ra, dec))).T

    return V_vec

def UVW2phi12(U, V, W, phi1, phi2, r, alpha_pole, delta_pole, phi1_0):

    if r==1:
        r = np.ones_like(phi1)

    cra = np.cos(phi1*RAD)
    cdec = np.cos(phi2*RAD)
    sra = np.sin(phi1*RAD)
    sdec = np.sin(phi2*RAD)

    R = Mrot(alpha_pole, delta_pole, phi1_0)

    MUVWphi12 = np.array([[cra*cdec, -sra, -cra*sdec],
                          [sra*cdec,  cra, -sra*sdec],
                          [sdec   ,   0,    cdec]])

    UVWvec = np.array([U,V,W])
    M = np.dot(Ag, np.dot(MUVWphi12.T,R))
    vec = np.dot(M, UVWvec).T

    vr  = vec[:,0]
    muphi1 = vec[:,1]/(k*r)
    muphi2 = vec[:,2]/(k*r)

    return (vr, muphi, muphi2)

def ellR(ra, dec, ra0, dec0, PA, e, sky=False):
    cPA = np.cos(RAD*PA)
    sPA = np.sin(RAD*PA)
    if sky:
        X = (ra - ra0)*np.cos(np.deg2rad(dec0))
    else:
        X = (ra - ra0)
    Y = (dec - dec0)
    R = np.sqrt(((X*cPA - Y*sPA)/(1-e))**2 + (X*sPA+ Y*cPA)**2)
    return R
