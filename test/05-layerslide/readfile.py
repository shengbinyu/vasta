import numpy as np
from  lattice import *
""" 
    read file with vasp package,such as :
    POSCAR     read_POSCAR()
    
"""

def get_line(file,split):
    strdata = file.readline() 
    strdata = strdata.strip('\n')
    strdata = strdata.split(split)
    strlist = []
    for x in strdata:
        if (x == ''):
            continue
        strlist.append(x.strip())
    return strlist
       
def read_POSCAR(pwd):
    filedir=pwd+'/POSCAR'
    file = open(filedir,'r')
    strdata = get_line(file,' ')
    symbols = strdata[0]
    strdata = get_line(file,' ')
    scale = float(strdata[0])

    param = []
    for i in range(3):
        strlist = get_line(file,'  ')
        param.append(strlist)
    param = np.array(param,dtype = 'float')
    param.shape = -1,3
    strlist = get_line(file,' ')
    Element_species = strlist
    strlist = get_line(file,' ')
    Element_number  = strlist
    Element_number  = np.array(Element_number,dtype='int32')

    Element_positions = []
    strlist = get_line(file,' ')
    if ('S' in strlist[0] or 's' in strlist[0]):
        cordlabel=get_line(file, ' ')
    else:
        cordlabel = strlist
    for i in range(len(Element_number)):
        for j in range(Element_number[i]):
            strdata = get_line(file,' ')
            Element_positions.append(strdata[0:3])

    Element_positions = np.array(Element_positions,dtype='float')
    file.close()

    param = scale*param
    if ('D' in cordlabel[0] or 'd' in cordlabel[0]):
        scale_poisitions = Element_positions
        atoms = Lattice(symbols=symbols,elements=Element_species,numbers=Element_number
                ,cell=param,scale_positions=Element_positions)
    else:
        poisitions = Element_positions
        atoms = Lattice(symbols=symbols,elements=Element_species,numbers=Element_number,
                cell=param,positions=Element_positions)

    return atoms

