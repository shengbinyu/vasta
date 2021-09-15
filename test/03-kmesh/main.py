from readfile import *
from writefile import *
from system import *
from lattice import *
from utilities import *

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    latt = read_POSCAR('./')
    print(latt.strprint())
    write_kmesh('./',latt,1e-3)
