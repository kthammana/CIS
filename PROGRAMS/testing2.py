import numpy as np
from FileIO import read_calbody, read_calreadings, read_empivot, read_optpivot, read_output1, read_ctfiducials, read_emfiducials, read_emnav, read_output2
from Registration import registrationArunMethod
from EMPivotCalibration import pivotCalibration, GtoEM
from OpticalPivotCalibration import opticalCalibration
from Point3d import Point3d
from DistortionCorrection import calcDistortionCorrection, correctDistortion

# ### DEBUGGING REGISTRATION
# # use given d (from PA1 data set g) and artificial D
# artiD = np.array([[0.3404, 0.5853, 0.2238],[-100.2885, -80.1432, 76.7532], 
#               [111.0513, -62.0671, 79.7082],[10.4224, -142.7955, 156.2376],
#               [11.1531, -109.2218, -101.3906],[-89.4758, -189.9503, -24.8612],
#               [121.8640, -171.8742, -21.9062],[21.2351, -252.6026, 54.6232]])
# artiF_D = registrationArunMethod(d, artiD, "D")
# # artificial R and p calculated on MATLAB
# artiR = np.array([[0.0721, 0.7381, -0.6709],[-0.7320, -0.4177, -0.5382],[-0.6774, 0.5299, 0.5102]])
# print('Calculated R:', artiF_D.R)
# print('Expected R:',artiR)
# print('Calculated p:',artiF_D.p.__str__())
# print('Expected p: 0.3404, 0.5853, 0.2238')
# # reg_error = 0
# # for i in range(3):
# #     temp_point = Point3d("D",F_D.R[i])
# #     reg_error += temp_point.error(artiR[i])
# #     # print(i,':',reg_error)
# # reg_error += artiF_D.p.error([0.3404, 0.5853, 0.2238])
# # print('Registration error:',reg_error)


def printPA2OutputErrors(dataset):

<<<<<<< HEAD
coef, q_min, q_max = calcDistortionCorrection(np.vstack(C_expected), np.vstack(C), 3)
C_errors = np.zeros(C.shape[0:2]) # stores error of each frame
for i in range(C.shape[0]): # N_frames
    for j in range(C.shape[1]): # N_C
        P_Cexp = Point3d("C", C_exp[i][j])
        C_errors[i][j] = P_Cexp.error(C_expected[i][j])
print('Average calibration error per point:', np.mean(C_errors), 'mm')
=======
    C_exp,P_em,P_opt = read_output1(dataset+"-output1.txt")
>>>>>>> 1a6b0c8585bdb54b697993c2af8a66a4a7de2a3e

    # calculating expected Cs
    d, a, c = read_calbody(dataset+"-calbody.txt")
    D, A, C = read_calreadings(dataset+"-calreadings.txt")
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
            P_Cexp = Point3d("C", C_exp[i][j])
            C_errors[i][j] = P_Cexp.error(C_expected[i][j])
    print('Average calibration error per point:', np.mean(C_errors), 'mm')

    # EM pivot calibration
    G = read_empivot(dataset+"-empivot.txt")
    G_corr = np.empty(G.shape)
    for i in range(G_corr.shape[0]):
        G_corr[i] = correctDistortion(G[i], coef, q_min, q_max, 3)
    P_em_exp, P_tip, g = pivotCalibration(G_corr)
    print('EM Pivot Error:',P_em_exp.error(P_em),'mm')

    # Fiducial point registration
    b = read_ctfiducials(dataset+"-ct-fiducials.txt")
    G_EM = read_emfiducials(dataset+"-em-fiducialss.txt")
    G_EM_corr = np.empty(G_EM.shape)
    for i in range(G_EM_corr.shape[0]):
        G_EM_corr[i] = correctDistortion(G_EM[i], coef, q_min, q_max, 3)
    B = GtoEM(G_EM_corr, P_tip, g)

    # Calculate F_reg
    F_reg = registrationArunMethod(B, b, "EM")

    # Compute tip location w.r.t. CT image
    G_nav = read_emnav(dataset+"-EM-nav.txt")
    G_nav_corr = np.empty(G_nav.shape)
    for i in range(G_nav_corr.shape[0]):
        G_nav_corr[i] = correctDistortion(G_nav[i], coef, q_min, q_max, 3)
    V = GtoEM(G_nav_corr, P_tip, g)
    v_exp = np.empty([G_nav.shape[0],3])
    for i in range(V.shape[0]):
        v_exp[i] = F_reg.R.dot(V[i]) + F_reg.p.coords

    # Calculate average v error
    v = read_output2(dataset+"-output2.txt")
    errors = np.empty(v.shape[0])
    for i in range(v.shape[0]):
        p_v = Point3d("CT",v[i])
        errors[i] = p_v.error(v_exp[i])
    print('Avg v error:', np.mean(errors))

printPA2OutputErrors("PA2 Student Data/pa2-debug-f")