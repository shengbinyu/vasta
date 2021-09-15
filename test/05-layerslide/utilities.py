import numpy as np
import math

def lineEquation(x, T1, T2):
    if (abs(T2[0] - T1[0])< 1e-6):
        if (x > T1[0]):
            y = True
        else:
            y = False
    else:
        k = (T2[1] - T1[1]) / (T2[0] - T1[0])
        y = k * (x - T1[0]) + T1[1]
    return y

def lineBoundary(x,y, T1, T2):
    if (abs(T2[0] - T1[0])< 1e-6):
        if ((x - T1[0]) > 1e-3):
            outofline = True
        else:
            outofline = False
    else:
        k = (T2[1] - T1[1]) / (T2[0] - T1[0])
        y0 = k * (x - T1[0]) + T1[1]
        if ((y - y0) > 1e-3):
            outofline = True
        else:
            outofline = False
    return outofline

def bToc(mesh_b,b1,b2):
    number_mesh = len(mesh_b)
    mesh_c = np.zeros(mesh_b.shape)
    for ii in range(number_mesh):
        k1 = mesh_b[ii,0]
        k2 = mesh_b[ii,1]
        c1 = k1 * b1[0] + k2 * b2[0]
        c2 = k1 * b1[1] + k2 * b2[1]
        mesh_c[ii,0]=c1
        mesh_c[ii,1]=c2
    return mesh_c

def listSupercell(super):
    if (super % 2 == 1):
        super = int((super + 1) / 2)
        super_ss = np.arange(-1 * super + 1, super)
    else:
        super = int((super) / 2)
        super_ss = np.arange(-1 * super, super)
    return super_ss

def get_twist_angle(index):
    n = index[0]
    m = index[1]
    cstheta = (m**2+n**2+4*m*n)/(m**2+n**2+m*n) * 0.5
    theta=math.acos(cstheta)*180/math.pi
    return theta
