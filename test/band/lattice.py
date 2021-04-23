import numpy as np
import copy
from vector import *

class Lattice(object):
    def __init__(self,symbols=None,elements=None,numbers=None,totatoms=None,atomindex=None,cell=None,recipcell=None,positions=None,
            scale_positions=None,vecangle=None):
        if symbols is not None:
            self.symbols = symbols
        if elements is not None:
            self.elements = elements
        if numbers is not None:
            self.numbers = numbers
        if cell is not None:
            self.cell = cell
            if recipcell is None:
                self.recipcell = self.trans_recipcell()
        if recipcell is not None:
            self.recipcell = recipcell
            if cell is None:
                self.cell = self.trans_cell()
        if positions is not None:
            self.positions = positions
            if scale_positions is None:
                self.scale_positions = self.trans_scale_positions()
        if scale_positions is not None:
            self.scale_positions = scale_positions
            if positions is None:
                self.positions = self.trans_positions()
        if totatoms is not None:
            self.totatoms = totatoms
        elif atomindex is None:
            self.totatoms = self.count_atoms()
        if atomindex is not None:
            self.atomindex = atomindex
        elif atomindex is None:
                self.atomindex = self.get_atomindex()
        if vecangle is not None:
            self.vecangle = vecangle
        elif vecangle is None:
                self.vecangle = self.get_vecangle()


    def get_symbols(self):
        return self.symbols

    def get_elements(self):
        return self.elements

    def get_numbers(self):
        return self.cell

    def get_positions(self):
        return self.positions

    def get_scale_positions(self):
        return self.scale_positions

    def set_symbols(self,symbols):
        self.symbols = symbols

    def set_elements(self,elements):
        self.elements = elements

    def set_numbers(self,numbers):
        self.numbers = numbers

    def set_positions(self,positions):
        self.positions = positions

    def set_scale_positions(self,scale_positions):
        self.scale_positions = scale_positions

    def update_positions(self,position_new):
        self.positions = position_new
        self.scale_positions = self.trans_scale_positions()

    def update_scale_positions(self,scale_positions_new):
        self.scale_positions = scale_positions_new
        self.positions = self.trans_positions()

    def update_cell(self,cell_new):
        self.cell = cell_new
        self.recipcell = self.trans_recipcell()

    def update_cell(self,recipcell_new):
        self.recipcell = recipcell_new
        self.cell = self.trans_cell()

    def update_elements_numbers(self,elements_new,numbers_new):
        self.elements = elements_new
        self.numbers= numbers_new

    def count_atoms(self):
        numbers = self.numbers
        atoms_count = 0
        for i in range(len(numbers)):
            for j in range(numbers[i]):
                atoms_count += 1
        return atoms_count

    def get_atomindex(self):
        numbers = self.numbers
        totatoms = self.totatoms
        atomindex = np.array([],dtype='int32')
        for ii in range(len(numbers)):
            index_tmp = np.ones(numbers[ii],dtype='int32')*ii
            atomindex = np.append(atomindex,index_tmp)
        return atomindex

    def cell_norm(self):
        cell = self.cell
        v = []
        for i in range(3):
            v.append(vector(x=cell[i,0],y=cell[i,1],z=cell[i,2]).norm())
        v_norm = np.array(v,dtype = 'float')
        return v_norm

    def get_vecangle(self):
        cell = np.array(self.cell,dtype='float')
        a1 = vector(cell[0, 0], cell[0, 1], cell[0, 2])
        a2 = vector(cell[1, 0], cell[1, 1], cell[1, 2])
        a3 = vector(cell[2, 0], cell[2, 1], cell[2, 2])
        alpha = vector.vec_angle(a2,a3)
        beta  = vector.vec_angle(a1,a3)
        gamma = vector.vec_angle(a1,a2)
        angle = np.array([alpha,beta,gamma],dtype='float')
        return angle

    def get_volume(self):
        cell = self.cell
        a1   = vector(cell[0,0],cell[0,1],cell[0,2])
        a2   = vector(cell[1,0],cell[1,1],cell[1,2])
        a3   = vector(cell[2,0],cell[2,1],cell[2,2])
        volume = a1.dot(a2.cross(a3))
        return volume

    def trans_scale_positions(self):
        positions = self.positions
        atoms_count = self.count_atoms()
        cell_t = np.array(self.cell,dtype='float')
        cell_t = np.linalg.inv(cell_t)
        scale_positions = np.zeros((atoms_count,3),dtype='float')
        for i in range(atoms_count):
            scale_positions[i,:] = np.dot(positions[i,:],cell_t)
        return scale_positions

    def trans_positions(self):
        scale_positions = self.scale_positions
        atoms_count = self.count_atoms()
        cell = self.cell
        positions = np.zeros((atoms_count,3),dtype='float')
        for i in range(atoms_count):
            for j in range(3):
                positions[i,:] += scale_positions[i,j]*cell[j,:]
        return positions

    def trans_recipcell(self):
        volume = self.get_volume()
        cell = self.cell
        a1   = vector(cell[0,0],cell[0,1],cell[0,2])
        a2   = vector(cell[1,0],cell[1,1],cell[1,2])
        a3   = vector(cell[2,0],cell[2,1],cell[2,2])
        b1   = a2.cross(a3)*2*pi/volume
        b2   = a3.cross(a1)*2*pi/volume
        b3   = a1.cross(a2)*2*pi/volume
        recipcell = np.zeros((3,3),dtype='float')
        recipcell[0,:] = b1.transarr()
        recipcell[1,:] = b2.transarr()
        recipcell[2,:] = b3.transarr()
        #recipcell[0,:] = [b1.x,b1.y,b1.z]
        #recipcell[1,:] = [b2.x,b2.y,b2.z]
        #recipcell[2,:] = [b3.x,b3,y,b3.z]
        return  recipcell

    def find_position(self,atomslist,indexlist):
        species = self.elements
        species_number = self.numbers
        positions = self.positions
        atomslistlen = len(atomslist)
        indexlistlen = len(indexlist)
        if (atomslistlen != indexlistlen):
            print('Error:length of ',atomslist,'is not same as ',indexlist)
            sys.exit(0)
        listlen = atomslistlen
        indexlist = np.array(indexlist,dtype='int')
        indexlist = indexlist - 1
        index_set = []
        posit = np.zeros((listlen,3),dtype='float')
        for i in range(listlen):
            for j in range(len(species)):
                if (atomslist[i] == species[j]):
                    index = 0
                    for k in range(j):
                        index += species_number[k]
                    index = index + indexlist[i]
                    index_set.append(index)
                    posit[i,:]= positions[index,:]
        return index_set,posit
        
    def atoms_distance(self,atomslist,indexlist):
        index_set,posit = self.find_position(atomslist,indexlist)
        atoms1 = vector(array=posit[0,:])
        atoms1.strprint()
        atoms2 = vector(array=posit[1,:])
        atoms2.strprint()
        distance = atoms1.get_distance(atoms2)
        return distance

    def vert_distance(self,atomslist,indexlist):
        index_set,posit = self.find_position(atomslist,indexlist)
        heigh = abs(posit[0,2] - posit[1,2])
        return heigh

    def atoms_angle(self,atomslist,indexlist):
        index_set,posit = self.find_position(atomslist,indexlist)
        atom1 = vector(array=posit[0,:])
        atom1.strprint()
        atom2 = vector(array=posit[1,:])
        atom2.strprint()
        atom3 = vector(array=posit[2,:])
        atom3.strprint()
        v21 = atom1.sub(atom2)
        v23 = atom3.sub(atom2)
        angle = v21.get_angle(v23)
        return angle

    def layer_distance_change(self,atomslist,indexlist,d_change):
        positions_new = copy.deepcopy(self.positions)
        index_set,posit = self.find_position(atomslist,indexlist)
        listlen = len(atomslist)
        posit_new = np.zeros((listlen,3),dtype='float')
        for i in range(listlen):
            posit_new[i,:] = posit[i,:]
            posit_new[i,2] = posit[i,2]+d_change
        for i in range(listlen):
            positions_new[index_set[i],:] = posit_new[i,:]
        return positions_new

    def delete_atoms(self,atomslist,indexlist):
        positions = copy.deepcopy(self.positions)
        elements  = copy.deepcopy(self.elements)
        numbers   = copy.deepcopy(self.numbers)
        index_set,posit = self.find_position(atomslist,indexlist)
        positions_new = np.delete(positions,index_set,axis=0)
        index_zero = []
        for i in range(len(elements)):
            count_i   = atomslist.count(elements[i])
            numbers[i] = numbers[i] - count_i
            if numbers[i] == 0 :
                index_zero.append(i)
        if len(index_zero) != 0:
            for j in range(len(index_zero)):
                del elements[index_zero[j]]
            numbers = np.delete(numbers,index_zero)
        return positions_new,elements,numbers

    def strprint(self):
        volume = self.get_volume()
        print('Sysmbol:\n',self.symbols)
        print('cell:\n',self.cell)
        print('recipcell:\n',self.recipcell)
        print('sepcies:\n',self.elements)
        print('sepcies_number:\n',self.numbers)
        print('scale_positions:\n',self.scale_positions)
        print('positions:\n',self.positions)
        print('cell_norm:\n',self.cell_norm())
        print('volume:\n',volume)
        print('vecangle:\n', self.vecangle)