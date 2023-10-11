import numpy as np
from Point3d import Point3d
from Registration import registrationArunMethod

def opticalCalibration(d, D, H):
    # defining local coordinate system
    Ht = H[0].transpose()
    H0 = np.mean(Ht, axis=1).reshape((-1,1))
    h = Ht - H0
    F_H = registrationArunMethod(h.transpose(), H[0], "H")
    negI = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
    for i in range(0, D.shape[0]):

        # calculating F_D for each frame
        # to convert from optical tracker to EM coordinates
        F_D = registrationArunMethod(d, D[i], "D")
        F_H = registrationArunMethod(h.transpose(), H[i], "H")
        F_DH = F_D.inverse() * F_H
        if i == 0:
            R_DHs = np.hstack((F_DH.R, negI))
            t_DHs = F_DH.p.coords.transpose()[..., np.newaxis]
        else:
            R_DH = np.hstack((F_DH.R, negI))
            R_DHs = np.vstack((R_DHs, R_DH))
            t_DHs = np.vstack((t_DHs, F_DH.p.coords.transpose()[..., np.newaxis]))
    t_DHs = t_DHs * -1
    
    # singular value decomposition least squares
    U, S, Vt = np.linalg.svd(R_DHs, full_matrices=True)
    zeros = np.zeros((6,30))
    Sinv = np.hstack((np.linalg.inv(np.diag(S)), zeros))
    y = Sinv @ U.transpose() @ t_DHs
    x = Vt.transpose() @ y
    point = Point3d("EM", x[3][0],x[4][0],x[5][0])
    return point