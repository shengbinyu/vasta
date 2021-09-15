import numpy as np
import copy
from math import *
from vector import *
from utilities import *
from lattice import *
import matplotlib.pyplot as plt

class System(object):
    def __init__(self,lattice,symbols=None,supercell=None,elements=None,numbers=None,totatoms=None,atomindex=None,positions=None):
        if lattice is not None:
            self.lattice = lattice
        elif lattice is None:
            print('System have none lattice!')

        if symbols is not None:
            self.symbols = symbols
        elif symbols is None:
            self.symbols = self.lattice.symbols

        if supercell is not None:
            self.supercell = supercell
        elif supercell is None:
            self.supercell = np.array([1,1,1])

        if elements is not None:
            self.elements = elements
        elif elements is None:
            self.elements = self.get_Elements()

        if numbers is not None:
            self.numbers = numbers
        elif elements is None:
            self.numbers = self.get_Numbers()

        if totatoms is not None:
            self.totatoms = totatoms
        elif elements is None:
            self.totatoms = self.get_TotalAtoms()

        if atomindex is not None:
            self.atomindex = atomindex

        if positions is not None:
            self.positions = positions
        elif positions is None:
            self.positions,self.atomindex = self.get_Positions()


    def get_Elements(self):
        return self.lattice.elements

    def get_Numbers(self):
        numbers = self.lattice.numbers
        supercell = self.supercell
        count = supercell[0]*supercell[1]*supercell[2]
        numbers = numbers*count
        return numbers

    def get_TotalAtoms(self):
        numbers = self.numbers
        totatoms= int(np.sum(numbers))
        return totatoms

    def get_Positions(self):
        positions = self.lattice.positions
        super=np.array(self.supercell,dtype='int32')
        super_x = listSupercell(super[0])
        super_y = listSupercell(super[1])
        super_z = listSupercell(super[2])
        cell = self.lattice.cell
        numbers = np.array(self.lattice.numbers)
        numbers_atom = np.sum(numbers)
        atomindex_cell = self.lattice.atomindex
        super_index = np.array([],dtype='int32')
        for x in super_x:
            for y in super_y:
                for z in super_z:
                    tmp = [x,y,z]
                    super_index = np.append(super_index,tmp)
        super_index=super_index.reshape(-1,3)
        super_positions = np.array([],dtype='float')
        atom_index_super = np.array([],dtype='int')
        for atom_index in range(int(numbers_atom)):
            for index in super_index:
                v_tmp = np.array(cell[0,:]*index[0]+cell[1,:]*index[1]+cell[2,:]*index[2])
                positions_tmp=positions[atom_index]+v_tmp
                super_positions = np.append(super_positions,positions_tmp)
                atom_index_super = np.append(atom_index_super,atomindex_cell[atom_index])
        super_positions = super_positions.reshape(-1,3)
        return super_positions,atom_index_super

    def twist_layer(self,twist_angle,layer_index):
        twist_angle = twist_angle * pi /180
        atoms_tot_number  = self.totatoms
        positions = self.positions
        rx = positions[:,0]
        ry = positions[:,1]
        rz = positions[:,2]
        twist_atom_index_TF = rz > layer_index
        twist_atom_index = np.arange(atoms_tot_number,dtype='int32')[twist_atom_index_TF]
        for index in twist_atom_index:
            positions_x = rx[index] * cos(twist_angle) - ry[index] * sin(twist_angle)
            positions_y = rx[index] * sin(twist_angle) + ry[index] * cos(twist_angle)
            positions[index,0]=positions_x
            positions[index,1]=positions_y
        twist_positions = positions
        self.positions = positions
        return twist_positions

    def strprint(self):
        print('Supercell:\n', self.supercell)
        print('Elements:\n',  self.elements)
        print('Numbers:\n' ,  self.numbers)
        print('Total Atoms:\n', self.totatoms)
        # print('AtomIndex:\n', self.atomindex)
        print('Position:\n',  self.positions)

    def plot_2D_atoms(self,layer_index):
        positions = self.positions
        rz = positions[:, 2]
        atom_index_TF_top = rz > layer_index
        atom_index_TF_bottom = rz < layer_index
        top_layer = positions[atom_index_TF_top,:]
        bottom_layer = positions[atom_index_TF_bottom,:]
        plt.figure()
        plt.scatter(bottom_layer[:, 0], bottom_layer[:, 1], c='b')
        plt.scatter(top_layer[:, 0], top_layer[:, 1],c='r')
        plt.scatter(0,0,c='k',marker='*')
        plt.axis('equal')
        plt.show()

    def get_Twist_lattice(self,vec_index):
        positions = self.positions
        elements = self.elements
        totatmos  = self.totatoms
        atom_index = self.atomindex
        # vec_angle = self.lattice.vecangle[0]
        vec_org = self.lattice.cell
        m,n = vec_index
        symbols = 'twist'+str(m)+str(n)
        T0 = np.zeros(3,dtype='float')
        T1 = m*vec_org[0,:]+n*vec_org[1,:]
        vec1 = vector(array=T1)
        vec2 = vec1.vec_roate(120)
        T2 = vec2.transarr()
        Td = T1+T2
        atoms_positions = np.array([],dtype='float')
        atoms_twist_index = np.array([], dtype='int32')
        for ii in range(totatmos):
            x = positions[ii,0]
            y = positions[ii,1]
            y_tmp1 = lineEquation(x, T0, T1)
            y_tmp2 = lineEquation(x, T0, T2)
            y_tmp3 = lineEquation(x, T1, Td)
            y_tmp4 = lineEquation(x, T2, Td)

            prec = 1e-3
            comp1eq = fabs(y - y_tmp1) < prec
            comp2eq = fabs(y - y_tmp2) < prec
            comp1gq = (y - y_tmp1) > prec
            comp2gq = (y - y_tmp2) > prec
            comp3lq = (y_tmp3 - y) > prec
            comp4lq = (y_tmp4 - y) > prec
            comp1 = comp1gq or comp1eq
            comp2 = comp2gq or comp2eq
            comp3 = comp3lq
            comp4 = comp4lq

            if (comp1 and comp2 and comp3 and comp4):
                atoms_positions  =  np.append(atoms_positions,positions[ii,:])
                atoms_twist_index = np.append(atoms_twist_index,atom_index[ii])

        atoms_positions = atoms_positions.reshape(-1,3)
        twist_cell = np.vstack((T1,T2,vec_org[2,:]))
        numbers = np.zeros(len(elements),dtype='int32')
        for jj in range(len(elements)):
            numbers[jj]=len(atoms_positions[atoms_twist_index == jj])
        twistlattice = Lattice(symbols=symbols,elements=elements,numbers=numbers,cell=twist_cell,positions=atoms_positions)
        # plot to check
        plt.figure()
        plt.scatter(0, 0, c='k', marker='*')
        plt.scatter(positions[:,0],positions[:,1],c='r')
        plt.scatter(atoms_positions[:, 0], atoms_positions[:, 1], c='b')
        plt.plot([T1[0],T0[0]],[T1[1],T0[1]],c='k')
        plt.plot([T2[0], T0[0]], [T2[1], T0[1]], c='k')
        plt.plot([Td[0], T1[0]], [Td[1], T1[1]], c='k')
        plt.plot([Td[0], T2[0]], [Td[1], T2[1]], c='k')
        plt.axis('equal')
        plt.show()
        
        self.twistlattice = twistlattice
        return twistlattice