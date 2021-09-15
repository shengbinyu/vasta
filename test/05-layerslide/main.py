from readfile import *
from writefile import *
from system import *
from lattice import *
from utilities import *
import numpy as np
import copy

def getpwd():
    pwd=os.popen('pwd').read().rstrip()
    return pwd

def intlist(listin):
    listout = []
    for x in listin:
        listout.append(int(x))
    return listout

def floatlist(listin):
    listout = []
    for x in listin:
        listout.append(float(x))
    return listout


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    latt = read_POSCAR('./')
    # print(latt.strprint())
    dfx=np.array([3/6,3/12],dtype='float')
    atomslist = ['Mo', 'S', 'Se']
    indexlist = intlist([2, 2, 2])
    scale_positions_new = latt.layer_slide_change(atomslist, indexlist, dfx)
    latt_new = copy.deepcopy(latt)
    latt_new.update_scale_positions(scale_positions_new)
    write_POSCAR('./', latt_new)
    latt_new.strprint()
