import numpy as np
from Point3d import Point3d
from Frame import Frame
from FileIO import read_optpivot, read_output
from Registration import registrationArunMethod

def opticalCalibration(D, H):
    # Dt = D[0].transpose()
    # D0 = np.mean(Dt, axis=1).reshape((-1,1))
    # d = Dt - D0
    # F_D = registrationArunMethod(d.transpose(), D[0], "D")
    Ht = H[0].transpose()
    H0 = np.mean(Ht, axis=1).reshape((-1,1))
    h = Ht - H0
    F_H = registrationArunMethod(h.transpose(), H[0], "H")
    # print(np.linalg.det(F_G.R))
    negI = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
    for i in range(0, D.shape[0]):
        Dt = D[i].transpose()
        D0 = np.mean(Dt, axis=1).reshape((-1,1))
        d = Dt - D0
        F_D = registrationArunMethod(d.transpose(), D[i], "D")
        # Ht = H[i].transpose()
        # H0 = np.mean(Ht, axis=1).reshape((-1,1))
        # h = Ht - H0
        # F_H = registrationArunMethod(h.transpose(), H[i], "H")
        F_DH = F_D.inverse() * F_H
        # print('det:',np.linalg.det(F_DH.R))
        if i == 0:
            R_DHs = np.hstack((F_DH.R, negI))
            t_DHs = F_DH.p.coords.transpose()[..., np.newaxis]
        else:
            R_DH = np.hstack((F_DH.R, negI))
            R_DHs = np.vstack((R_DHs, R_DH))
            t_DHs = np.vstack((t_DHs, F_DH.p.coords.transpose()[..., np.newaxis]))
    t_DHs = t_DHs * -1
    U, S, Vt = np.linalg.svd(R_DHs, full_matrices=True)
    zeros = np.zeros((6,30))
    Sinv = np.hstack((np.linalg.inv(np.diag(S)), zeros))
    y = Sinv @ U.transpose() @ t_DHs
    x = Vt.transpose() @ y
    point = Point3d("EM", x[3][0],x[4][0],x[5][0])
    return point


D,H = read_optpivot("../PA1 Student Data/pa1-debug-a-optpivot.txt")
P_opt_exp = opticalCalibration(D,H)
print('Calculated output:', P_opt_exp.__str__())
__,__,P_opt = read_output("../PA1 Student Data/pa1-debug-a-output1.txt")
print('Expected output:',P_opt)
print('Error:',P_opt_exp.error(P_opt),'mm')