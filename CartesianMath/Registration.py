import numpy as np
from FileIO import read_calbody, read_calreadings

def registrationArunMethod(a, b):
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
    return R, T

d, a, c = read_calbody("./PA1 Student Data/pa1-debug-a-calbody.txt")
D, A, C = read_calreadings("./PA1 Student Data/pa1-debug-a-calreadings.txt")
for i in range(D.shape[0]):
    registrationArunMethod(D[i], d)
print(D[0])
print(d)
