from readfile import *
from writefile import *
from system import *
from lattice import *
from utilities import *

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    latt = read_POSCAR('./')
    # print(latt.strprint())
    index=[29,30]
    twist_angle=get_twist_angle(index)
    print('twist angle:', twist_angle)
    twistJanus = System(lattice=latt,supercell=[150,150,1])
    twist_positions = twistJanus.twist_layer(twist_angle,15.0)
    # twistJanus.plot_2D_atoms(15.0)
    twistlatt=twistJanus.get_Twist_lattice(index)
    twistlatt.strprint()
    write_POSCAR('./',twistlatt)
    write_fdf('./',twistlatt)