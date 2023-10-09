import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from Point3d import Point3d

class Frame(object):
    
    def __init__(self, name, R, p, neighbors=[]):
        """
            
        Parameters
        ----------
        name : string or character
            Name of the frame (e.g., 'A').
        R : numpy ndarray
            Frame's rotation matrix.
        p : numpy array
            Frame's origin point.
        neighbors : list, optional
            List containing all adjacent Frames. The default is [].
    
        Returns
        -------
        None.
    
        """
        self.name = name
        self.R = R
        self.p = p
        self.neighbors = neighbors
    
    def inverse(self):
        invR = np.transpose(self.R) # R^-1 = R^T
        invP = (-1) * np.transpose(self.R) * self.p
        out = Frame(self.name, invR, invP, self.neighbors)
        return out
    
    def __mul__(self, other): # self * other
        if isinstance(other, Frame): # F_self * F_other
            R = self.R * other.R
            p = self.R * other.p + self.p
            name = self.name+" to "+other.name
            F_transformation = Frame(name,R,p) # set no neighbors?
            return F_transformation
        elif isinstance(other, Point3d): # F_self * p_other
            return (self.R * other + self.p)
        else:
            raise TypeError('Cannot multiply Frame by non-Point object')
    
    # TO-DO: finish this function (currently only adds target to path)
    def __findFrame__(self,point,path):
        for frame in self.neighbors: # iterate through all neighbors
            if frame.name == point.frame: # if target frame is a neighbor
                path.append(frame)
                return path
        # otherwise search neighbors' children
        for frame in self.neighbors: # separate loop so it's bfs not dfs
            frame.__findFrame__(point, path)
        return path
    
    def transformPoint(self,point):
        # check frame
        if point.__getattr__('frame') == self.name:
            return point
        
        # transform point to current frame
        path = self.__findFrame__(point, [])
        
        ## TO DO:
        # loop through path
        # multiply frame by frame to get transformation frame matrix F
        new_point = F*point
        return new_point
    
    def visualizePoint(self, point, frame):
        p = self.transformPoint(p,frame)
        x = point.__getattr__('x')
        y = point.__getattr__('y')
        z = point.__getattr__('z')
        return x,y,z