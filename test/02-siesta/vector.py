import numpy as np
from math import *

class vector(object):
    def __init__(self,x=0,y=0,z=0,array=None):
        if array is None:
            self.x = x
            self.y = y
            self.z = z
        else:
            self.x = array[0]
            self.y = array[1]
            self.z = array[2]

    def __add__(self,num):
        return vector(self.x+num,self.y+num,self.z+num)

    def __sub__(self,num):
        return vector(self.x-num,self.y-num,self.z-num)

    def __mul__(self,num):
        return vector(self.x*num,self.y*num,self.z*num)

    def __truediv__(self,num):
        return vector(self.x/num,self.y/num,self.z/num)

    def __pow__(self,num):
        return vector(self.x**num,self.y**num,self.z**num)

    def add(self,obj):
        return vector(self.x+obj.x,self.y+obj.y,
                self.z+obj.z)

    def sub(self,obj):
        return vector(self.x-obj.x,self.y-obj.y,
                self.z-obj.z)

    def mul(self,obj):
        return vector(self.x*obj.x,self.y*obj.y,
                self.z*obj.z)

    def dot(self,obj):
        return self.x*obj.x+self.y*obj.y+ self.z*obj.z

    def cross(self,obj):
        result = vector()
        result.x = self.y*obj.z - self.z*obj.y
        result.y = self.z*obj.x - self.x*obj.z
        result.z = self.x*obj.y - self.y*obj.x
        return result

    def norm(self):
        return  sqrt(self.x**2+self.y**2+self.z**2)

    def transarr(self):
        return np.array([self.x,self.y,self.z],dtype='float')

    def get_distance(self,obj):
        dx = self.x - obj.x
        dy = self.y - obj.y
        dz = self.z - obj.z
        distance = sqrt(abs(dx**2+dy**2+dz**2))
        return distance

    def get_angle(self,obj1,obj2):
        v1 = obj1.sub(self)
        v2 = obj2.sub(self)
        angle = acos((v1.dot(v2))/(v1.norm()*v2.norm()))
        angle = angle*180/ pi
        return angle

    def vec_angle(self,v2):
        v1 = self
        angle = acos((v1.dot(v2)) / (v1.norm() * v2.norm()))
        angle = angle * 180 / pi
        return angle

    def vec_roate(self,angle_z):
        theta = angle_z *pi / 180
        roate_matrix = np.array([[cos(theta),-sin(theta),0],[sin(theta),cos(theta),0],[0,0,1]],dtype='float')
        vec = self.transarr()
        vec = vec.T
        vecnew = np.dot(roate_matrix,vec).T
        vecroate = vector(array=vecnew)
        return vecroate

    def strprint(self):
        print('('+str(self.x)+', '+str(self.y)+', '+str(self.z)+')')

