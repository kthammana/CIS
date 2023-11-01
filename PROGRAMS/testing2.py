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
C_errors = np.zeros(C.shape[0:2]) # stores error of each frame
for i in range(C.shape[0]): # N_frames
    for j in range(C.shape[1]): # N_C
        # P_C = Point3d(C[i][j])
        P_Cexp = Point3d("C", C_exp[i][j])
        C_errors[i][j] = P_Cexp.error(C[i][j])
print('Average calibration error per point:', np.mean(C_errors), 'mm')

### DEBUGGING REGISTRATION
# use given d and artificial D
artiD = np.array([[0.3404, 0.5853, 0.2238],[-100.2885, -80.1432, 76.7532], 
              [111.0513, -62.0671, 79.7082],[10.4224, -142.7955, 156.2376],
              [11.1531, -109.2218, -101.3906],[-89.4758, -189.9503, -24.8612],
              [121.8640, -171.8742, -21.9062],[21.2351, -252.6026, 54.6232]])
artiF_D = registrationArunMethod(d, artiD, "D")
# artificial R and p calculated on MATLAB
artiR = np.array([[0.0721, 0.7381, -0.6709],[-0.7320, -0.4177, -0.5382],[-0.6774, 0.5299, 0.5102]])
print('Calculated R:', artiF_D.R)
print('Expected R:',artiR)
print('Calculated p:',artiF_D.p.__str__())
print('Expected p: 0.3404, 0.5853, 0.2238')
reg_error = 0
for i in range(3):
    temp_point = Point3d("D",F_D.R[i])
    reg_error += temp_point.error(artiR[i])
    # print(i,':',reg_error)
reg_error += artiF_D.p.error([0.3404, 0.5853, 0.2238])
print('Registration error:',reg_error)

# EM pivot calibration
G = read_empivot("PA1 Student Data/pa1-debug-g-empivot.txt")
G_corr = np.empty(G.shape)
for i in range(G_corr.shape[0]):
    G_corr[i] = correctDistortion(G[i], coef, q_min, q_max, 3)
# print(G)
P_em_exp = pivotCalibration(G_corr)
print('EM Calculated output:', P_em_exp.__str__())
print('EM Expected output:',P_em)
print('EM Pivot Error:',P_em_exp.error(P_em),'mm')