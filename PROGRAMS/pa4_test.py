#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 09:12:49 2023

@author: kianabronder
"""
import numpy as np
from FileIO import read_calbody, read_calreadings, read_empivot, read_optpivot, read_output1, read_ctfiducials 
from FileIO import read_emfiducials, read_emnav, read_output2, read_mesh, read_output3, read_samplereadings, read_probbody
from Registration import registrationArunMethod
from EMPivotCalibration import pivotCalibration, GtoEM
from OpticalPivotCalibration import opticalCalibration
from Point3d import Point3d
from DistortionCorrection import bernstein, calcDistortionCorrection, correctDistortion
from ClosestPointOnTriangle import findClosestPointOnTriangle
from Mesh import Mesh, insert, search
from Frame import Frame

def calcDistance(x, y):
    return np.sqrt((x[0] - y[0])**2 + (x[1] - y[1])**2 + (x[2] - y[2])**2)

# i do not think this is right but the paper did not specify how to check for convergence
def checkConvergence(w, converged): 
    if sum(w)/w.shape[0] >= 0.5:
        converged = True
    return converged

## Test file
dataset = "PA345 Student Data/PA4-A-Debug"
print("PA4 Output Errors w KdTree Search:")

# test I/O functions
A, A_tip, N_A = read_probbody("PA345 Student Data/Problem4-BodyA.txt")
B, B_tip, N_B = read_probbody("PA345 Student Data/Problem4-BodyB.txt")
V, ind, n = read_mesh("PA345 Student Data/Problem4MeshFile.sur")
mesh = Mesh(V, ind, n)
a, b = read_samplereadings(dataset+"-SampleReadingsTest.txt", N_A, N_B)
s_exp, c_exp, mag = read_output3(dataset+"-Output.txt") # same file format as PA3
s_error = 0
c_error = 0
mag_error = 0

d_k = np.empty((a.shape[0], 3))
A_tip = A_tip.transpose()[..., np.newaxis]
for i in range(a.shape[0]): # N_samples:
    F_Ai = registrationArunMethod(A, a[i], "A")
    F_Bi = registrationArunMethod(B, b[i], "B")
    F_BA = F_Bi.inverse() * F_Ai
    d_k[i] = (F_BA.R @ A_tip)[:,0] + F_BA.p.coords
    # d_error += calcDistance(d_exp[i], d_k[i])

# for PA4, iteratively find F_reg
s_k = np.empty((a.shape[0], 3))

# kdtree search to find the closest points c_k to s_k
root = None
for i in range(ind.shape[0]): # build tree
    verts = mesh.getVerticesOfTriangle(i)
    root = insert(root, mesh.calcCentroid(verts), verts, 0, i)

converged = False
F_reg = Frame("reg", np.identity(3), Point3d("reg", 0, 0, 0)) # Assume F_reg = I for initial guess
d_max = 0.5
while not converged:
    s_k = F_reg.R.dot(d_k) + F_reg.p.coords
    c_k = np.zeros((a.shape[0], 3))
    w = np.zeros(d_k.shape[0])
    for i in range(d_k.shape[0]): # search tree
        nearest_node = search(root, s_k[i])
        c_k[i] = findClosestPointOnTriangle(s_k[i], nearest_node.triangle)
        # if i < (d_k.shape[0] - 1):
        #     s_k[i+1] = F_reg.R.dot(d_k[i+1]) + F_reg.p.coords
        ## I'm definitely creating F_reg wrong
        # F_reg = registrationArunMethod(d_k, c_k, "reg")
        # print(nearest_node.idx)
        if np.abs(c_k[i] - (F_reg.R.dot(s_k[i]) + F_reg.p.coords)) <= d_max:
            w_i = 1
    F_reg = ## TODO: argmin function
    converged = checkConvergence(w, converged)
    
for i in range(d_k.shape[0]): # calculate errors
    c_error += calcDistance(c_exp[i], c_k[i])
    s_error += calcDistance(s_exp[i], s_k[i])
    shortest_dist = calcDistance(s_k[i], c_k[i])
    mag_error += (np.abs(mag[i]-shortest_dist))

# Slight difference in these values
print("s_k error: ", s_error/s_k.shape[0])
print("c_k error: ", c_error/c_k.shape[0])
print("mag error: ", mag_error/c_k.shape[0])