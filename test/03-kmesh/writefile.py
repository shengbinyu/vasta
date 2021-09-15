import numpy as np
import math
from lattice import *
from vector import *
from utilities import *
import  matplotlib.pyplot as plt

"""
    @author:usb
    @date:2021/5/18
    @version:1.0
    write structure file ,band data ,and dos data,format as :
    POSCAR.vasp  write_POSCAR()
    band.dat     write_band()
    siesta.fdf   write_fdf()
    Kmesh.dat    write_kmesh()
"""

def write_POSCAR(pwd,atoms):
    symbols=atoms.symbols
    filedir=pwd+'/'+symbols+'.vasp'
    file = open(filedir,'w+')
    file.writelines(symbols)
    file.writelines('\n')          #The first line
    file.writelines(' 1.000\n')        #The second line

    #The third line
    celllist = atoms.cell.tolist()
    for j in range(3):
        cellstr = [str(i) for i in celllist[j]]
        for k in range(3):
            file.writelines('  ')
            file.writelines(cellstr[k].ljust(20,"0"))
            file.writelines('  ')
        file.writelines('\n')

    #The forth line 
    elements = atoms.elements
    for i in range(len(elements)):
        file.writelines('   ')
        file.writelines(elements[i].ljust(5," "))
    file.writelines('\n')

    #The fifth line
    elements_numbers = atoms.numbers.tolist()
    numbers_str = [str(i) for i in elements_numbers]
    for i in range(len(elements)):
        file.writelines('   ')
        file.writelines(numbers_str[i].ljust(5," "))
    file.writelines('\n')

    #The sixth to last line
    atoms_numbers = atoms.count_atoms()
    file.write('Direct\n')
    scale_positions = atoms.scale_positions.tolist()
    for j in range(atoms_numbers):
        scale_positions_str = [str(i) for i in scale_positions[j]]
        for k in range(3):
            file.writelines('  ')
            file.writelines(scale_positions_str[k].ljust(20,"0"))
            file.writelines('  ')
        file.writelines('\n')
    file.close()


def write_band(pwd,band):
    """
        This is used for writing 'band.dat' file
    """
    filedir=pwd+'/band.dat'
    file = open(filedir,'w+')
    file.write('#  Kstep\tE(eV)\n')    #The first line

    Kstep   = band.Kstep.tolist()
    Kenergy = band.Kenergy.tolist()
    Knumber = band.Knumber
    Nband   = band.Nband
    ispin   = band.ispin
    #band.strprint()
    for i in range(Nband):
        file.writelines('#  band '+str(i+1)+'\n')
        for j in range(Knumber):
            file.writelines('  ')
            file.writelines(str(round(Kstep[j],5)).ljust(7,"0"))
            file.writelines('\t')
            for k in range(ispin):
                file.writelines(str(Kenergy[j][i][k]).ljust(10,"0"))
                file.writelines('\t')
            file.writelines('\n')
    file.close()


def write_fdf(pwd,atoms):
    """
        can write some format of file:
        siesta input file:{label}.fdf <- my.fdf is used
        VASP structure file POSCAR
    """

    elements_list = ['H','C','N','O','F','S','Cl','Cr','Se','Br','Mo','Te','I','W']
    elements_number = ['1','6','7','8','9','16','17','24','34','35','42','52','53','74']
    species = atoms.elements
    species_numbers = []
    for i in range(len(species)):
        list_index = elements_list.index(species[i])
        species_numbers.append(elements_number[list_index])

    numbers = atoms.numbers
    positions = atoms.positions
    lattice_constant = atoms.cell_norm()
    lattice_angle    = atoms.vecangle
    lattice_cell = atoms.cell
    lattice = np.hstack((lattice_constant,lattice_angle))
    total_numbers = numbers.sum()
    species_counts = len(numbers)
    filename=pwd+'/'+atoms.symbols+'.fdf'

    file = open(filename,'w')
    file.write('#------------------ General system descriptors --------------#\n')
    file.write('SystemName        '+atoms.symbols+'\n')
    file.write('SystemLabel       '+atoms.symbols+'\n')
    file.write('\n')
    file.write('NumberOfSpecies   '+str(species_counts)+'\n')
    file.write('NumberOfAtoms     '+str(total_numbers)+'\n')
    file.write('\n')
    file.write('%block ChemicalSpeciesLabel\n')
    for i in range(len(species)):
        file.write('\t'+str(i+1)+'\t'+str(species_numbers[i])+'\t'+species[i]+'\n')
    file.write('%endblock ChemicalSpeciesLabel\n')

    file.write('\n')
    file.write('LatticeConstant        1.0 Ang\n')
    # file.write('%block LatticeParameters\n')
    # for i in range(6):
    #     file.write('\t'+str(round(lattice[i],4)))
    # file.write('\n')
    # file.write('%endblock LatticeParameters\n')
    file.write('%block LatticeVectors\n')
    for i in range(3):
        for j in range(3):
            file.write('\t'+str(round(lattice_cell[i,j],4)).ljust(8,"0"))
        file.write('\n')
    file.write('%endblock LatticeVectors\n')
    file.write('\n')
    file.write('AtomicCoordinatesFormat  Ang\n')
    file.write('\n')

    file.write('%block AtomicCoordinatesAndAtomicSpecies\n')
    for i in range(len(numbers)):
        if (i == 0):
              index = 0
        else:
            index=index+numbers[i-1]  #have problem!
        for j in range(numbers[i]):
            numindex = index + j
            for k in range(3):
                file.write('\t'+str(round(positions[numindex,k],5)).ljust(8,"0"))
            file.write('\t'+str(i+1)+'\n')
    file.write('%endblock AtomicCoordinatesAndAtomicSpecies\n')
    file.write('#------------------------------------------------------------#\n')
    file.write('\n')
    file.close()

    # os.popen('cat my.fdf>>temp.fdf')
    # open both files, copy content from my.fdf file to temp.fdf file
    with open('my.fdf', 'r') as firstfile, open(filename, 'a') as secondfile:
        # read content from first file
        for line in firstfile:
            # append content to second file
            secondfile.write(line)
    secondfile.close()


def write_kmesh(pwd, atoms,pho):
    """
          This method is used to write Kmesh.dat file in the first Brillouin zone:
                KPOINTS  format
          Only applied to two dimensional hexagonal lattice!
      """

    recipcell = atoms.recipcell
    cell_norm = atoms.cell_norm()
    a_norm = cell_norm[0]
    # b_norm = np.linalg.norm(recipcell[0,:])
    alpha = atoms.vecangle[2]
    b1 = recipcell[0, :]
    b2 = recipcell[1, :]

    # six points in concer of BZ
    if (abs(alpha-120) < 1e-3):
        K1 = (1/3)*(b1+b2)
        K1_vec = vector(array=K1)
        K2 = K1_vec.vec_roate(-60).transarr()
        K3 = K1_vec.vec_roate(-120).transarr()
        K4 = K1_vec.vec_roate(-180).transarr()
        K5 = K1_vec.vec_roate(-240).transarr()
        K6 = K1_vec.vec_roate(-300).transarr()
    elif (abs(alpha-60) < 1e-3):
        K1 = (1 / 3) * (b1 + 2* b2)
        K1_vec = vector(array=K1)
        K2 = K1_vec.vec_roate(-60).transarr()
        K3 = K1_vec.vec_roate(-120).transarr()
        K4 = K1_vec.vec_roate(-180).transarr()
        K5 = K1_vec.vec_roate(-240).transarr()
        K6 = K1_vec.vec_roate(-300).transarr()

    # generate mesh in the bais of b1 and b2
    k_number = int(math.sqrt(4 / pho))
    k_line = np.linspace(-1, 1, k_number)
    mesh_b = np.array([], dtype='float')
    for ii in k_line:
        for jj in k_line:
            mesh_b = np.append(mesh_b, [ii, jj,0.0])

    mesh_b = np.reshape(mesh_b, (-1, 3))
    mesh_c = bToc(mesh_b, b1, b2)

    #  The points of mesh in the first of BZ
    mesh_bz = np.array([],dtype='float')
    for ii in range(len(mesh_c)):
        x = mesh_c[ii,0]
        y = mesh_c[ii,1]
        #boudary of BZ
        flag12 = lineBoundary(x,y, K1, K2)
        flag23 = lineBoundary(x,y, K2, K3)
        flag34 = lineBoundary(x,y, K3, K4)
        flag45 = lineBoundary(x,y, K4, K5)
        flag56 = lineBoundary(x,y, K5, K6)
        flag61 = lineBoundary(x,y, K6, K1)
        if (flag12*flag45 == -1) and (flag23*flag56 == -1) and (flag34*flag61 == -1):
            mesh_bz = np.append(mesh_bz,mesh_c[ii,:])

    mesh_bz = np.reshape(mesh_bz,(-1,3))
    mesh_number = len(mesh_bz)
    mesh_bz_w = mesh_bz/(2*math.pi)

    #write mesh to file
    filename = pwd + '/'+'Kmesh.dat'
    file = open(filename, 'w')
    file.write('First Brillouin zone generated by vasta \n')
    file.write(str(mesh_number) +'\n')
    file.write('Cartesian\n')
    for i in range(mesh_number):
        for j in range(3):
            file.write('   '+str(round(mesh_bz_w[i,j],6)).ljust(10,"0"))
        file.write('    1.000000\n')
    file.close()

    # plot to check
    print('Number of points in Kmesh:', mesh_number)
    print('Plot Kmesh in first BZ to check...')
    plt.figure()
    plt.scatter(mesh_bz[:, 0], mesh_bz[:, 1],s=4, c='r')
    plt.plot([K1[0], K1[0]], [K2[1], K2[1]], c='k')
    plt.plot([K2[0], K2[0]], [K3[1], K3[1]], c='k')
    plt.plot([K3[0], K3[0]], [K4[1], K4[1]], c='k')
    plt.plot([K4[0], K4[0]], [K5[1], K5[1]], c='k')
    plt.plot([K5[0], K5[0]], [K6[1], K6[1]], c='k')
    plt.plot([K6[0], K6[0]], [K1[1], K1[1]], c='k')
    plt.axis('equal')
    plt.show()