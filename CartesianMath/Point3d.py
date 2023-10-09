# -*- coding: utf-8 -*-
"""
CIS PA #1

Kiana Bronder, kbronde1
"""
import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

# 1a. 3D points, rotations, and frame transformations

#############################################
# Copyright (c) 2018 Fabricio JC Montenegro #
# Version 1.0                               #
# https://github.com/SplinterFM/Point3d     #
#############################################

# I took this code off this Github and added a frame variable

class Point3d(object):
    def __init__(self, frame, x=0, y=0, z=0):
        if type(x) == list:
            # if receives a list
            try:
                assert len(x) == 3
                self.coords = np.array(x+[1],dtype=np.float64)
            except:
                raise ValueError('To create a Point3d with a list it must have length 3')
        elif type(x) == tuple:
            # if receives a tuple
            try:
                assert len(x) == 3
                self.coords = np.array(list(x)+[1],dtype=np.float64)
            except:
                raise ValueError('To create a Point3d with a tuple it must have length 3')
        elif type(x) == np.ndarray:
            try:
                assert x.ndim == 1 and x.size == 4 and x[3] == 1
                self.coords = np.array(x, dtype=np.float64)
            except:
                print ("ndim", x.ndim)
                print ("size", x.size)
                print ("[3]", x[3])
                raise ValueError('To create a Point3d with an np.array it must have ndim 1, size 4, and the last element must be 1')
        elif type(x) == Point3d:
            self.coords = np.array(x.coords, dtype=np.float64)
        else:
            self.coords = np.array([x,y,z,1], dtype=np.float64)

        self.iter = 0
        self.frame = frame

    def copy(self):
        return Point3d(self.frame, self.coords)


    def __getattr__(self, name):
        """ When the user asks for attributes x, y, or z, we return
            coords[0], [1], and [2] """
        if name == 'x':
            return self.coords[0]
        elif name == 'y':
            return self.coords[1]
        elif name == 'z':
            return self.coords[2]
        elif name == 'frame':
            return self.frame

    '''
    def __getitem__(self, key):
        return self.coords[key]

    def __setattr__(self, name, value):
        """ For x, y, and z sets coords[0], [1], and [2].
            For everything else does the normal set. """
        if name == 'x':
            self.coords[0] = value
        elif name == 'y':
            self.coords[1] = value
        elif name == 'z':
            self.coords[2] = value
        elif name == 'frame':
            self.frame = value
        else:
            self.__dict__[name] = value

    def __setitem__(self, key, value):
        self.coords[key] = value


    def __repr__(self):
        return "Point3d({0}, {1}, {2})".format(self.x, self.y, self.z)

    def __str__(self):
        return "({0}, {1}, {2})".format(self.x, self.y, self.z)
    '''
# UPDATE TO INCLUDE FRAME CONSIDERATION
    def __add__(self, other):
        try:
            assert isinstance(other, Point3d)
            assert self.frame == other.__getattr__('frame')
            x = self.x + other.x
            y = self.y + other.y
            z = self.z + other.z
            return Point3d(x,y,z)
        except:
            raise TypeError("Arithmetic requires 2 Points in same frame")

    def __sub__(self, other):
        try:
            assert isinstance(other, Point3d)
            assert self.frame == other.__getattr__('frame')
            x = self.x - other.x
            y = self.y - other.y
            z = self.z - other.z
            return Point3d(x,y,z)
        except:
            raise TypeError("Arithmetic requires 2 Points in same frame")

    def __mul__(self, other):
        try:
            assert isinstance(other, Point3d)
            assert self.frame == other.__getattr__('frame')
            x = self.x * other
            y = self.y * other
            z = self.z * other
            return Point3d(x,y,z)
        except:
            raise TypeError("Arithmetic requires 2 Points in same frame")

    def dot(self, other):
        return self.coords.dot(other.coords)

    def cross(self, other):
        a = np.array([self.x, self.y, self.z])
        b = np.array([other.x, other.y, other.z])
        c = np.cross(a,b)
        return Point3d(list(c))

    def __pos__(self):
        return self

    def __neg__(self):
        return self * -1

    def __iter__(self):
        return self

    def next(self):
        if self.iter >= self.coords.size-1:
            self.iter = 0
            raise StopIteration
        else:
            self.iter += 1
            return self.coords[self.iter-1]

    def length(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def translate(self, vector):
        self.x += vector.x
        self.y += vector.y
        self.z += vector.z

    def rotateX(self, radians):
        c = np.cos(radians)
        s = np.sin(radians)
        xmat = np.array([[1, 0, 0, 0],
                        [0, c,-s, 0],
                        [0, s, c, 0],
                        [0, 0, 0, 1]])
        self.coords = xmat.dot(self.coords)

    def rotateY(self, radians):
        c = np.cos(radians)
        s = np.sin(radians)
        ymat = np.array([[c, 0, s, 0],
                        [ 0, 1, 0, 0],
                        [-s, 0, c, 0],
                        [ 0, 0, 0, 1]])
        self.coords = ymat.dot(self.coords)

    def rotateZ(self, radians):
        c = np.cos(radians)
        s = np.sin(radians)
        zmat = np.array([[c, -s, 0, 0],
                        [s, c, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])
        self.coords = zmat.dot(self.coords)
        
    def visualize(self, other):
        x = np.array([self.x, other.__getattr__('x')])
        y = np.array([self.y, other.__getattr__('y')])
        z = np.array([self.z, other.__getattr__('z')])
        return x,y,z