"""
CIS PA #2

Kiana Bronder, kbronde1
Keerthana Thammana, lthamma1
"""
import numpy as np
from FileIO import read_calbody, read_calreadings, read_empivot, read_output1, read_ctfiducials, read_emfiducials, read_emnav, read_output2
from Registration import registrationArunMethod
from EMPivotCalibration import pivotCalibration, GtoEM
from Point3d import Point3d
from DistortionCorrection import calcDistortionCorrection, correctDistortion


import os
cwd = os.getcwd()
print(cwd)

# change input and output filenames
dataset = "PA2 Student Data/pa2-unknown-j"
filename = "pa2-unknown-j-output2.txt"
print(filename)
d, a, c = read_calbody(dataset+"-calbody.txt")
G = read_empivot(dataset+"-empivot.txt")
D, A, C = read_calreadings(dataset+"-calreadings.txt")
# C_exp,P_em,P_opt = read_output1(dataset+"-output1.txt")

# output file:
    # N_frames, NAME-OUTPUT2.TXT
    # v_x,i , v_y,i , v_z,i

# output file name and directory
f = open("../OUTPUT_PA2/"+filename, "w")

# calculating expected Cs
C_expected = np.zeros(C.shape)
for i in range(D.shape[0]):
    F_D = registrationArunMethod(d, D[i], "D")
    F_A = registrationArunMethod(a, A[i], "A")
    F_DA = F_D.inverse() * F_A
    for j in range(c.shape[0]):
        C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
        C_expected[i][j] = [C0_exp[0][0], C0_exp[1][0], C0_exp[2][0]]

coef, q_min, q_max = calcDistortionCorrection(np.vstack(C_expected), np.vstack(C), 5)
# C_errors = np.zeros(C.shape[0:2]) # stores error of each frame
# for i in range(C.shape[0]): # N_frames
#     for j in range(C.shape[1]): # N_C
#         P_Cexp = Point3d("C", C_exp[i][j])
#         C_errors[i][j] = P_Cexp.error(C_expected[i][j])
# print('Average calibration error per point:', np.mean(C_errors), 'mm')

# EM pivot calibration
G = read_empivot(dataset+"-empivot.txt")
G_corr = np.empty(G.shape)
for i in range(G_corr.shape[0]):
    G_corr[i] = correctDistortion(G[i], coef, q_min, q_max, 5)
P_em_exp, P_tip, g = pivotCalibration(G_corr)

# Fiducial point registration
b = read_ctfiducials(dataset+"-ct-fiducials.txt")
G_EM = read_emfiducials(dataset+"-em-fiducialss.txt")
G_EM_corr = np.empty(G_EM.shape)
for i in range(G_EM_corr.shape[0]):
    G_EM_corr[i] = correctDistortion(G_EM[i], coef, q_min, q_max, 5)
B = GtoEM(G_EM_corr, P_tip, g)

# Calculate F_reg
F_reg = registrationArunMethod(B, b, "EM")

# Compute tip location w.r.t. CT image
G_nav = read_emnav(dataset+"-EM-nav.txt")
G_nav_corr = np.empty(G_nav.shape)
for i in range(G_nav_corr.shape[0]):
    G_nav_corr[i] = correctDistortion(G_nav[i], coef, q_min, q_max, 5)
V = GtoEM(G_nav_corr, P_tip, g)
v_exp = np.empty([G_nav.shape[0],3])
for i in range(V.shape[0]):
    v_exp[i] = F_reg.R.dot(V[i]) + F_reg.p.coords

# Calculate v_exp and write to output file
N_frames = G_nav_corr.shape[0]
f.write(str(N_frames) + ", " + filename + "\n")
# f.write('N_frames:' + str(N_frames) + ", " + filename + "\n")
for i in range(N_frames):
    point = Point3d("EM",v_exp[i])
    f.write(point.__str__() + "\n")
    # f.write(point.__str__() + "\\\\" + "\n")
    
# Calculate average v error
# v = read_output2(dataset+"-output2.txt")
# errors = np.empty(v.shape[0])
# for i in range(v.shape[0]):
#     p_v = Point3d("CT",v[i])
#     errors[i] = p_v.error(v_exp[i])
# print('Avg v error:', np.mean(errors))
f.close()