
class Element(object):
    def __init__(self, symbol, number, mass=None,en_name=None,chs_name=None):
        self.symbol=symbol
        self.number=number
        self.mass = mass
        self.en_name=en_name
        self.chs_name = chs_name

    def strprint(self):
        print('Element number:', self.number)
        print('Element symbol:', self.symbol)
        print('Element atom mass:', self.mass)
        print('Element English name:', self.en_name)
        print('Element Chinese name:', self.chs_name)

# Periodic Table
#First line
H=Element('H',1,mass=1.008,en_name='Hydrogen',chs_name='氢')
He=Element('He',2,mass=4.0026,en_name='Helium',chs_name='氦')
#Second line
Li=Element('Li',3,mass=6.94,en_name='Lithium',chs_name='锂')
Be=Element('Be',4,mass=9.0122,en_name='Beryllium',chs_name='铍')
B=Element('B',5,mass=10.81,en_name='Boron',chs_name='硼')
C=Element('C',6,mass=12.011,en_name='Carbon',chs_name='碳')
N=Element('N',7,mass=14.007,en_name='Nitrogen',chs_name='氮')
O=Element('O',8,mass=15.999,en_name='Oxygen',chs_name='氧')
F=Element('F',9,mass=18.998,en_name='Fluorine',chs_name='氟')
Ne=Element('Ne',10,mass=20.180,en_name='Neon',chs_name='氖')
#Third line

