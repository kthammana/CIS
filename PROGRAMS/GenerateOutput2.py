"""
CIS PA #2

Kiana Bronder, kbronde1
Keerthana Thammana, lthamma1
"""
import numpy as np
from Point3d import Point3d
from FileIO import read_calbody, read_calreadings, read_empivot, read_optpivot
from FileIO import read_ctfiducials, read_emfiducials, read_emnav
from Registration import registrationArunMethod
from EMPivotCalibration import pivotCalibration
from OpticalPivotCalibration import opticalCalibration

import os
cwd = os.getcwd()
print(cwd)

# change input and output filenames
filename = "pa1-unknown-h-output1.txt"
d, a, c = read_calbody("PA1 Student Data/pa1-unknown-h-calbody.txt")
G = read_empivot("PA1 Student Data/pa1-unknown-h-empivot.txt")
D,H = read_optpivot("PA1 Student Data/pa1-unknown-h-optpivot.txt")
D, A, C = read_calreadings("PA1 Student Data/pa1-unknown-h-calreadings.txt")

# output file:
    # N_C , N_frames, NAME-OUTPUT1.TXT
    # P_EM
    # P_opt
    # C_1 through C_N_frames

# output file name and directory
f = open("../OUTPUT/"+filename, "w")

# Calculate P_EM
P_EM = pivotCalibration(G)

# Calculate P_opt
P_opt = opticalCalibration(d,D,H)

# Calculate C_exp and write all variables to output file
N_frames = D.shape[0]
N_C = c.shape[0]
f.write(str(N_C) + ", " + str(N_frames) + ", " + filename + "\n")
f.write(P_EM.__str__() + "\n")
f.write(P_opt.__str__() + "\n")
for i in range(N_frames): # N_frames
    F_D = registrationArunMethod(d, D[i], "D")
    F_A = registrationArunMethod(a, A[i], "A")
    F_DA = F_D.inverse() * F_A
    for j in range(N_C): # N_C
        C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
        point = Point3d("EM",C0_exp[0][0], C0_exp[1][0], C0_exp[2][0])
        f.write(point.__str__() + "\n")
f.close()