import numpy as np
import math
from Point3d import Point3d
from FileIO import read_calbody, read_calreadings
from Registration import registrationArunMethod

def berstein(N, k, v):
    return math.comb(N, k) * (1-v)**(N-k) * v**k

def bersteinPolynomialF(N, q, q_min, q_max):
    diff = q_max - q_min
    u = np.empty(q.shape)
    for i in range(q.shape[0]):
         u[i] = (q[i] - q_min)/diff

    F = np.empty((u.shape[0], N**3))
    for x in range(u.shape[0]):
        index = 0
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    F[x][index] = berstein(N, i, u[x][0])*berstein(N, j, u[x][1])*berstein(N, k, u[x][2])
                    index += 1

    return F

def calcDistortionCorrection(p, q, N):
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

    # create F matrix of berstein polynomials based on measurements
    F = bersteinPolynomialF(N, q, q_min, q_max)
    print(F.shape)

    # singular value decomposition least squares
    U, S, Vt = np.linalg.svd(F, full_matrices=True)
    zeros = np.zeros((Vt.shape[0],U.shape[0]-Vt.shape[0]))
    Sinv = np.hstack((np.linalg.inv(np.diag(S)), zeros))
    y = Sinv @ U.transpose() @ p
    coef = Vt.transpose() @ y

    return coef, q_min, q_max

def correctDistortion(q, coef, q_max, q_min, N):
    F = bersteinPolynomialF(N, q, q_min, q_max)
    p = F @ coef

    return p # return values without distortion

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
        C_expected[i][j] = [C0_exp[0][0], C0_exp[1][0], C0_exp[2][0]]

coef, q_min, q_max = calcDistortionCorrection(np.vstack(C_expected), np.vstack(C), 3)

# print("C_expected: \n")
# print(C_expected)
# print("C with distortion: \n")
# print(correctDistortion(np.vstack(C), coef, q_min, q_max, 3))
