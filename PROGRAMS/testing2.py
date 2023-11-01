import numpy as np
from FileIO import read_calbody, read_calreadings, read_empivot, read_optpivot, read_output1
from Registration import registrationArunMethod
from EMPivotCalibration import pivotCalibration
from OpticalPivotCalibration import opticalCalibration
from Point3d import Point3d
from DistortionCorrection import calcDistortionCorrection, correctDistortion

C_exp,P_em,P_opt = read_output1("PA1 Student Data/pa1-debug-g-output1.txt")

# calculating expected Cs
d, a, c = read_calbody("PA1 Student Data/pa1-debug-g-calbody.txt")
D, A, C = read_calreadings("PA1 Student Data/pa1-debug-g-calreadings.txt")
C_expected = np.zeros(C.shape)
for i in range(D.shape[0]):
    F_D = registrationArunMethod(d, D[i], "D")
    F_A = registrationArunMethod(a, A[i], "A")
    F_DA = F_D.inverse() * F_A
    for j in range(c.shape[0]):
        C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
        C_expected[i][j] = [C0_exp[0][0], C0_exp[1][0], C0_exp[2][0]]

coef, q_min, q_max = calcDistortionCorrection(np.vstack(C_expected), np.vstack(C), 3)

# EM pivot calibration
G = read_empivot("PA1 Student Data/pa1-debug-g-empivot.txt")
G_corr = np.empty(G.shape)
for i in range(G_corr.shape[0]):
    G_corr[i] = correctDistortion(G[i], coef, q_min, q_max, 3)
print(G)
P_em_exp = pivotCalibration(G_corr)
print('EM Calculated output:', P_em_exp.__str__())
print('EM Expected output:',P_em)
print('EM Pivot Error:',P_em_exp.error(P_em),'mm')