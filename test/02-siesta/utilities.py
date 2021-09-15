import numpy as np
import math

def lineEquation(x, T1, T2):
    k = (T2[1] - T1[1]) / (T2[0] - T1[0])
    y = k * (x - T1[0]) + T1[1]
    return y

def listSupercell(super):
    if (super % 2 == 1):
        super = int((super + 1) / 2)
        super_ss = np.arange(-1 * super + 1, super)
    else:
        super = int((super) / 2)
        super_ss = np.arange(-1 * super, super)
    return super_ss
