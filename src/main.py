from readfile import *
from writefile import *
from system import *
from lattice import *
import numpy as np


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    latt = read_POSCAR('./')
    print(latt.strprint())
    twistJanus = System(lattice=latt,supercell=[8,8,1])
    twist_positions = twistJanus.twist_layer(21.787,15.0)
    # print(twistJanus.strprint())
    # twistJanus.plot_2D_atoms(15.0)
    twistlatt=twistJanus.get_Twist_lattice([1,2])
    twistlatt.strprint()
    write_POSCAR('./',twistlatt)