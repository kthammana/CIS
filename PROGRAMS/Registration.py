import numpy as np
from Frame import Frame
from Point3d import Point3d

def registrationArunMethod(a, b, frame_name):
    at = a.transpose()
    bt = b.transpose()

    # calculating average of all points
    a_c = np.mean(at, axis=1).reshape((-1,1))
    b_c = np.mean(bt, axis=1).reshape((-1,1))
    a = at - a_c
    b = bt - b_c

    # Arun's method
    H = np.zeros((3, 3))
    for i in range(a.shape[0]):
        H += [[a[0][i]*b[0][i], a[0][i]*b[1][i], a[0][i]*b[2][i]],
              [a[1][i]*b[0][i], a[1][i]*b[1][i], a[1][i]*b[2][i]],
              [a[2][i]*b[0][i], a[2][i]*b[1][i], a[2][i]*b[2][i]],]
        
    # singular value decomposition
    U, S, Vt = np.linalg.svd(H, full_matrices=True)
    R = np.matmul(np.transpose(Vt), np.transpose(U))
    T = b_c - np.matmul(R, a_c)
    F = Frame(frame_name, R, Point3d(frame_name, T[0][0], T[1][0], T[2][0]))
    return F
