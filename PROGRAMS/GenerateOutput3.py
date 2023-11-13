"""
CIS PA #3

Kiana Bronder, kbronde1
Keerthana Thammana, lthamma1
"""
import numpy as np
from FileIO import read_probbody, read_mesh, read_samplereadings, read_output3
from Registration import registrationArunMethod
from Point3d import Point3d
from ClosestPointOnTriangle import findClosestPointOnTriangle

import os
cwd = os.getcwd()
print(cwd)

A, A_tip, N_A = read_probbody("PA345 Student Data/Problem3-BodyA.txt")
B, B_tip, N_B = read_probbody("PA345 Student Data/Problem3-BodyB.txt")

a, b = read_samplereadings("PA345 Student Data/PA3-A-Debug-SampleReadingsTest.txt", N_A, N_B)

d_k = np.empty((a.shape[0], 3))
for i in range(a.shape[0]): # N_samples:
    F_A = registrationArunMethod(a[i], A, "A")
    F_B = registrationArunMethod(b[i], B, "B")
    F_BA = F_B.inverse() * F_A
    point = F_BA.R * A_tip + F_BA.p.coords
    print(point)
    d_k[i] = F_BA.R * A_tip + F_BA.p.coords