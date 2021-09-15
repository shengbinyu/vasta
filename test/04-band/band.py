import numpy as np
import sys
import matplotlib.pyplot as plt
from math import *

class Band(object):
    """Band Object
    Author   : yushengbin
    Date     : 2017/4/29
    Purpose  : get bandgap;get band structure
    Parmeters: fermi/Kcord/Kpath/Knumber/Kpoints/Kenergy/Kstep
               /Kstepcord/ispin/Nband/bandgap/vbm_locate/cbm_locate
    """
    def __init__(self,fermi=None,ispin=None,Nband=None,Kcord=None,
            Kcordstep=None,Kstep=None,Kpath=None,Knumber=None,
            Kpoints=None,Kenergy=None,bandgap=None):
        if fermi is not None:
            self.fermi = fermi
        if Kcord is not None:
            self.Kcord = Kcord
        if Kpath is not None:
            self.Kpath = Kpath
        if Knumber is not None:
            self.Knumber = Knumber
        if Kpoints is not None:
            self.Kpoints = Kpoints
        if Kenergy is not None:
            self.Kenergy = Kenergy
        if Kstep is not None:
            self.Kstep = Kstep
        else:
            self.Kstep = self.trans_Kstep()
        if Kcordstep is not None:
            self.Kcordstep = Kcordstep
        else:
            self.Kcordstep = self.trans_Kcordstep()
        if ispin is not None:
            self.ispin = ispin 
        else:
            self.ispin = self.get_ispin()
        if Nband is not None:
            self.Nband = Nband
        else:
            self.Nband = self.get_Nband()
        if bandgap is not None:
            self.bandgap = bandgap
        else:
            bandgap,vbm_locate,cbm_locate = self.get_bandgap()
            self.set_bandgap(bandgap,vbm_locate,cbm_locate)



    def set_fermi(self,fermi):
        self.fermi = fermi

    def set_ispin(self,ispin):
        self.ispin = ispin

    def set_Kcord(self,Kcord):
        self.Kcord = Kcord

    def set_Kcordstep(self,Kcordstep):
        self.Kcordstep = Kcordstep

    def set_Kpath(self,Kpath):
        self.Kpath = Kpath

    def set_Kstep(self,Kstep):
        self.Kstep = Kstep

    def set_Knumber(self,Knumber):
        self.Knumber = Knumber

    def set_Kpoints(self,Kpoints):
        self.Kpoints = Kpoints
    
    def set_Kenergy(self,Kenergy):
        self.Kenergy = Kenergy

    def set_bandgap(self,bandgap,vbm_locate,cbm_locate):
        self.bandgap = bandgap
        self.vbm_locate = vbm_locate
        self.cbm_locate = cbm_locate

    def strprint(self):
        print('fermi:\n',self.fermi)
        print('ispin:\n',self.ispin)
        print('Nband:\n',self.Nband)
        print('Kcord: \n',self.Kcord)
        print('Kcordstep: \n',self.Kcordstep)
        print('Knumber: \n',self.Knumber)
        print('Kpath: \n',self.Kpath)
        print('Kpoinits:\n',self.Kpoints)
        print('Kenergy:\n ',self.Kenergy)
        print('Kstep:\n ',self.Kstep)
        print('bandgap:\n',self.bandgap)
        print('vbm_locate:\n',self.vbm_locate)
        print('cbm_locate:\n',self.cbm_locate)
        #print('len(step)= ',self.Kstep.shape,'\nshape(Kenergy)= '
                #,self.Kenergy.shape)

    def get_ispin(self):
        Knumber,Nband,ispin = self.Kenergy.shape
        return ispin

    def get_Nband(self):
        Knumber,Nband,ispin = self.Kenergy.shape
        return Nband

    def trans_Kstep(self):
        from readfile import read_POSCAR
        from main import getpwd
        Kpoints = self.Kpoints 
        Knumber = self.Knumber
        pwd = str(getpwd())
        atoms = read_POSCAR(pwd)   #!!need atoms information
        Kstep = np.zeros((Knumber),dtype = 'float')
        for i in range(Knumber-1):
            Kspc = Kpoints[i+1,:] - Kpoints[i,:]
            Kspc = np.dot(Kspc,atoms.recipcell)
            Kstep[i+1] =Kstep[i]+sqrt(Kspc[0]**2+Kspc[1]**2+Kspc[2]**2)
        return Kstep

    def trans_Kcordstep(self):
        Kstep = self.trans_Kstep()
        Kcord = self.Kcord
        Knumber  = self.Knumber
        Kinteral = int(Knumber/(len(Kcord)-1))
        Kcordstep= np.zeros((len(Kcord)),dtype = 'float')
        for i in range(len(Kcord)):
            if (i == 0):
                Kcordstep[i] = 0.00
                continue
            index = i * Kinteral - 1
            Kcordstep[i] = Kstep[index]
        return Kcordstep

    def get_bandgap(self):
        fermi = self.fermi
        Kenergy = self.Kenergy
        Knumber = self.Knumber
        Nband   = self.Nband
        ispin   = self.ispin
        vbmtemp = -10000
        cbmtemp = 10000
        if (ispin == 1):
            for i in range(Knumber):
                index = -1
                for j in range(Nband):
                    if (Kenergy[i,j,0] <= fermi):
                        index += 1
                if ((Kenergy[i,index,0] <= fermi) and (Kenergy[i,index,0] >= vbmtemp)):
                        vbmtemp = Kenergy[i,index,0]
                        vbm_locate = [i,index]
                if ((Kenergy[i,index+1,0] >= fermi ) and (Kenergy[i,index+1,0] <= cbmtemp)):
                        cbmtemp = Kenergy[i,index+1,0]
                        cbm_locate = [i,index+1]
            vbm = np.array(vbmtemp,dtype = 'float')
            cbm = np.array(cbmtemp,dtype = 'float' )
            bandgap = cbm - vbm
        elif (ispin == 2):
            vbm = np.zeros((2),dtype = 'float')
            cbm = np.zeros((2),dtype = 'float')
            vbm_locate = np.zeros((2,2),dtype = 'int')
            cbm_locate = np.zeros((2,2),dtype = 'int')
            for k in range(2):
                vbmtemp = -10000
                cbmtemp = 100000
                for i in range(Knumber):
                    index = -1
                    for j in range(Nband):
                        if (Kenergy[i,j,k] <= fermi):
                            index +=1
                    if ((Kenergy[i,index,k] <= fermi) and (Kenergy[i,index,k] >= vbmtemp)):
                        vbmtemp = Kenergy[i,index,0]
                        vbm_locate[k,:] = np.array([i,index])
                    if ((Kenergy[i,index+1,k] >= fermi ) and (Kenergy[i,index+1,k] <= cbmtemp)):
                        cbmtemp = Kenergy[i,index+1,0]
                        cbm_locate[k,:] = np.array([i,index+1])
                vbm[k] = vbmtemp
                cbm[k] = cbmtemp
            bandgap = cbm - vbm
        return bandgap,vbm_locate,cbm_locate

    def pltband(self,Elim):
        """ plot band structure and save figure  """
        """set fermi level is zero"""
        Kcordstep = self.Kcordstep
        Knumber = self.Knumber
        Kpath   = self.Kpath
        Kstep   = self.Kstep
        Kenergy = self.Kenergy - self.fermi
        ispin   = self.ispin
        Nband   = self.Nband
        bandgap = self.bandgap
        vbm_locate = self.vbm_locate
        cbm_locate = self.cbm_locate
        Emin   = Elim[0]
        Emax   = Elim[1]
        if (Emax < Emin):
            print('Error: Emax(',Emax,') is smaller Emin(',Emin,')')
            sys.exit(0)
        outText = []
        for x in Kpath:
            if (x == '\Gamma'):
                outText.append("$"+x+"$")
                continue
            outText.append(x)
        if (ispin == 1):
            for i in range(Nband):
                plt.plot(Kstep,Kenergy[:,i,:],'b')
        elif(ispin == 2):
            for i in range(Nband):
                plt.plot(Kstep,Kenergy[:,i,0],'g')
                plt.plot(Kstep,Kenergy[:,i,1],'b')
        plt.ylabel('Energy/eV')
        plt.xticks(Kcordstep,outText)
        plt.xlim(0,Kstep[-1])
        plt.ylim(Emin,Emax)

        E_fermi = np.zeros((len(Kstep)))
        plt.plot(Kstep[:],E_fermi,'r--')
        YK=np.arange(Emin,Emax,0.1)
        Kcordstep.shape = -1,1
        K =Kcordstep.repeat(len(YK),axis=1)
        for i in range(1,len(Kcordstep)-1):
            plt.plot(K[i,:],YK,'k')

        #tx1 = Kstep[vbm_locate[0,0]]
        #ty1 = Kenergy[vbm_locate[0,0],vbm_locate[0,1],0]
        #x0  = Kstep[cbm_locate[0,0]]
        #y0  = Kenergy[cbm_locate[0,0],cbm_locate[0,1],0]
        #print('xy=(',x0,',',y0,')\nxytext=(',tx1,',',ty1,')')
        #plt.annotate('',xy=(x0,y0,),xytext=(tx1,ty1),arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))
        plt.savefig('band.png',dpi=150)


