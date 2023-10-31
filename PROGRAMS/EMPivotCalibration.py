import numpy as np
from Point3d import Point3d
from Registration import registrationArunMethod

def pivotCalibration(G):

    # defining a local coordinate system
    Gt = G[0].transpose()
    G0 = np.mean(Gt, axis=1).reshape((-1,1))
    g = Gt - G0

    # creating the [R -I][p_tip p_dimple]^T = [-t] matrix for each frame
    F_G = registrationArunMethod(g.transpose(), G[0], "G")
    negI = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
    R_Gs = np.hstack((F_G.R, negI))
    t_Gs = F_G.p.coords.transpose()[..., np.newaxis]
    for i in range(1, G.shape[0]):
        F_G = registrationArunMethod(g.transpose(), G[i], "G")
        R_G = np.hstack((F_G.R, negI))
        R_Gs = np.vstack((R_Gs, R_G))
        t_Gs = np.vstack((t_Gs, F_G.p.coords.transpose()[..., np.newaxis]))
    t_Gs = t_Gs * -1

    # singular value decomposition least squares
    U, S, Vt = np.linalg.svd(R_Gs, full_matrices=True)
    zeros = np.zeros((Vt.shape[0],U.shape[0]-Vt.shape[0]))
    Sinv = np.hstack((np.linalg.inv(np.diag(S)), zeros))
    y = Sinv @ U.transpose() @ t_Gs
    x = Vt.transpose() @ y

    # return p_dimple
    point = Point3d("EM", x[3][0],x[4][0],x[5][0])
    return point