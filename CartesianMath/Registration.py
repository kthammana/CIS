import numpy as np
from FileIO import read_calbody, read_calreadings
from Frame import Frame
from Point3d import Point3d

def registrationArunMethod(a, b, frame_name):
    at = a.transpose()
    bt = b.transpose()
    a_c = np.mean(at, axis=1).reshape((-1,1))
    b_c = np.mean(bt, axis=1).reshape((-1,1))
    a = at - a_c
    b = bt - b_c
    H = np.zeros((3, 3))
    for i in range(a.shape[0]):
        H += [[a[0][i]*b[0][i], a[0][i]*b[1][i], a[0][i]*b[2][i]],
              [a[1][i]*b[0][i], a[1][i]*b[1][i], a[1][i]*b[2][i]],
              [a[2][i]*b[0][i], a[2][i]*b[1][i], a[2][i]*b[2][i]],]
    U, S, Vt = np.linalg.svd(H, full_matrices=True)
    R = np.matmul(np.transpose(Vt), np.transpose(U))
    T = b_c - np.matmul(R, a_c)
    F = Frame(frame_name, R, Point3d(frame_name, T[0][0], T[1][0], T[2][0]))
    return F

d, a, c = read_calbody("./PA1 Student Data/pa1-debug-c-calbody.txt")
D, A, C = read_calreadings("./PA1 Student Data/pa1-debug-c-calreadings.txt")
for i in range(D.shape[0]):
    F_D = registrationArunMethod(d, D[i], "D")
    # print(np.linalg.det(F_D.R))
    F_A = registrationArunMethod(a, A[i], "A")
    F_DA = F_D.inverse() * F_A
    # print(F_DA.R)
    for j in range(c.shape[0]):
        C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
        print(C0_exp)

# F_D = registrationArunMethod(d, D[0], "D")
# Why D[0] and not all of D (i.e., all frames)? 
# oh i was just testing it, it should be all of D
# print('D[0]:',D[0])
# print('d:',d)
# print('F_D:',F_D.R, F_D.p.coords)
# d_i = F_D^-1 * D_i
# F_A = registrationArunMethod(a, A[0], "A")
# print('F_A:',F_A.R, F_A.p.coords)
# F_DA = F_D.inverse() * F_A
# C0_exp = F_DA * C[0]
# print(C0_exp)
