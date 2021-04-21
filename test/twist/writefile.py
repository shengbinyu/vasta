import numpy as np
import os
from lattice import *

"""
    @author:usb
    @date:2018/4/18
    @version:2.0
    write structure file ,band data ,and dos data,format as :
    POSCAR.vasp  write_POSCAR()
    band.dat     write_band()
    siesta.fdf   write_fdf()
    !dos.dat     write_dos() ->next
"""

def write_POSCAR(pwd,atoms):
    filedir=pwd+'/POSCAR.vasp'
    file = open(filedir,'w+')
    file.writelines(atoms.symbols)
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
    flag = input("Direct corridnate or Cartesian corrdinate('D''d'[default] or 'C''c' ):\n")
    if (flag == 'd' or flag == 'D' or flag == '\n'):
        file.write('Direct\n')
        scale_positions = atoms.scale_positions.tolist()
        for j in range(atoms_numbers):
            scale_positions_str = [str(i) for i in scale_positions[j]]
            for k in range(3):
                file.writelines('  ')
                file.writelines(scale_positions_str[k].ljust(20,"0"))
                file.writelines('  ')
            file.writelines('\n')
    else:
        file.write('Cartesian\n')
        positions = atoms.positions.tolist()
        for j in range(atoms_numbers):
            positions_str = [str(round(i,15)) for i in positions[j]]
            for k in range(3):
                file.writelines('  ')
                file.writelines(positions_str[k].ljust(18,"0"))
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


def write_fdf(atoms):
    """
        can write some format of file:
        siesta input file:{label}.fdf <- my.fdf is used
        VASP structure file POSCAR
    """

    elements_list = ['H','C','N','O','F','S','Se','Mo']
    elements_number = ['1','6','7','8','9','16','34','42']
    species = atoms.elements
    species_numbers = []
    for i in range(len(species)):
        list_index = elements_list.index(species[i])
        species_numbers.append(elements_number[list_index])

    numbers = atoms.numbers
    positions = atoms.positions
    lattice_constant = atoms.cell_norm()
    lattice_angle    = atoms.cell_angle()
    lattice = np.hstack((lattice_constant,lattice_angle))
    total_numbers = numbers.sum()
    species_counts = len(numbers)

    file = open('temp.fdf','w')
    file.write('#------------------ General system descriptors --------------#\n')
    file.write('SystemName        '+atoms.symbols[0]+'\n')
    file.write('SystemLabel       '+atoms.symbols[0]+'\n')
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
    file.write('%block LatticeParameters\n')
    for i in range(6):
        file.write('\t'+str(round(lattice[i],3)))
    file.write('\n')
    file.write('%endblock LatticeParameters\n')
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
                file.write('\t'+str(round(positions[numindex,k],3)))
            file.write('\t'+str(i+1)+'\n')
    file.write('%endblock AtomicCoordinatesAndAtomicSpecies\n')
    file.write('#------------------------------------------------------------#\n')
    file.write('\n')
    os.popen('cat my.fdf>>temp.fdf')
    file.close()