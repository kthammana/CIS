import numpy as np
import math
from Point3d import Point3d
from FileIO import read_calbody, read_calreadings
from Registration import registrationArunMethod

def bersteinPolynomial(N, k, v):
    return math.comb(N, k) * (1-v)**(N-k) * v**k

def calcDistortionCorrection(p, q):
    '''
    Parameters:
    p: ground-truth of 3D points, with no distortion
    q: corresponding 3D points read by navigation sensor, with distortion

    Returns:
    coef: coefficients for Bernstein polynomials to correct distortion
    q_min: lower bound for 3D points of bounding box for scaling
    q_max: upper bound for 3D points of bounding box for scaling
    '''

    # create a bounding box and scale measured values
    q_min = np.amin(q, axis=0)
    q_max = np.amax(q, axis=0)
    diff = q_max - q_min

    u = np.empty(q.shape)
    for i in range(q.shape[0]):
         u[i] = (q[i] - q_min)/diff
    # return coef, q_min, q_max


# def applyDistortionCorrection(X, coef):

#     return X_corr # return corrected C values

# calculating expected Cs
d, a, c = read_calbody("PA1 Student Data/pa1-debug-a-calbody.txt")
D, A, C = read_calreadings("PA1 Student Data/pa1-debug-a-calreadings.txt")
C_expected = np.zeros(C.shape)
for i in range(D.shape[0]):
    F_D = registrationArunMethod(d, D[i], "D")
    F_A = registrationArunMethod(a, A[i], "A")
    F_DA = F_D.inverse() * F_A
    for j in range(c.shape[0]):
        C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
        # gen_point = Point3d("DA", C0_exp[0][0], C0_exp[1][0], C0_exp[2][0])
        C_expected[i][j] = [C0_exp[0][0], C0_exp[1][0], C0_exp[2][0]]

calcDistortionCorrection(np.vstack(C_expected), np.vstack(C))
