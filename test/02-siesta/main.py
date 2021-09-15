from readfile import *
from writefile import *
from system import *
from lattice import *
import numpy as np

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    latt = read_POSCAR('./')
    print(latt.strprint())
    write_fdf(latt)
