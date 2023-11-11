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

A, A_tip = read_probbody("PA345 Student Data/Problem5-BodyA.txt")
print(A_tip)
