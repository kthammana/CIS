"""
CIS PA #4

Kiana Bronder, kbronde1
Keerthana Thammana, lthamma1
"""
import numpy as np
from FileIO import read_probbody, read_mesh, read_samplereadings, read_output3
from Registration import registrationArunMethod
from Point3d import Point3d
from ClosestPointOnTriangle import findClosestPointOnTriangle
from Mesh import Mesh, insert, search
from testing import calcDistance
from Frame import Frame

import os
cwd = os.getcwd()
print(cwd)

# change input and output filenames
dataset = "PA345 Student Data/PA4-J-Unknown"
output_filename = "PA4-J-Unknown-Output.txt"

f = open("../OUTPUT/"+output_filename, "w")

A, A_tip, N_A = read_probbody("PA345 Student Data/Problem4-BodyA.txt")
B, B_tip, N_B = read_probbody("PA345 Student Data/Problem4-BodyB.txt")
a, b = read_samplereadings(dataset+"-SampleReadingsTest.txt", N_A, N_B)

# Calculate d_k
A_tip = A_tip.transpose()[..., np.newaxis]
d_k = np.empty((a.shape[0], 3))
for i in range(a.shape[0]): # N_samples:
    F_Ai = registrationArunMethod(A, a[i], "A")
    F_Bi = registrationArunMethod(B, b[i], "B")
    F_BA = F_Bi.inverse() * F_Ai
    d_k[i] = (F_BA.R @ A_tip)[:,0] + F_BA.p.coords
    
f.write(str(a.shape[0]) + " " + output_filename+ "\n") # first line of file

# Construct KD tree
V, ind, n = read_mesh("PA345 Student Data/Problem4MeshFile.sur")
mesh = Mesh(V, ind, n)
root = None
for i in range(ind.shape[0]):
    verts = mesh.getVerticesOfTriangle(i)
    root = insert(root, mesh.calcCentroid(verts), verts, 0, i)
    
# Find c_k and calculate s_k
s_k = np.empty((a.shape[0], 3))
c_k = np.zeros((a.shape[0], 3))
converged = False
F_reg = Frame("reg", np.identity(3), Point3d("reg", 0, 0, 0)) # Assume F_reg = I for initial guess
iterations = 0
d_max = 0.3
while not converged and iterations < 100: # max num of iterations is 100
    for i in range(d_k.shape[0]): # search tree
        s_k[i] = F_reg.R.dot(d_k[i]) + F_reg.p.coords
        nearest_node = search(root, s_k[i])
        c_k[i] = findClosestPointOnTriangle(s_k[i], nearest_node.triangle)
        
    F_reg = registrationArunMethod(d_k, c_k, "reg")

    # Update s_k and calculate || s_k - c_k ||
    mag = np.empty(d_k.shape[0])
    for i in range(s_k.shape[0]):
        s_k[i] = F_reg.R.dot(d_k[i]) + F_reg.p.coords
        mag[i] = calcDistance(s_k[i], c_k[i])

    if np.mean(mag) < d_max:
        converged = True 

    iterations += 1

# Write s_k, c_k, and || s_k - c_k || to file
for i in range(d_k.shape[0]):
    point_s = Point3d("",s_k[i])
    point_c = Point3d("",c_k[i])
    string = point_s.__str__() + "     " + point_c.__str__()
    string += "     {:.3f}".format(round(mag[i], 3))
    f.write(string + "\n")
    
print('Wrote', output_filename)