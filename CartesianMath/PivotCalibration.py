import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from Point3d import Point3d
from Frame import Frame
from FileIO import read_empivot
from Registration import registrationArunMethod

def pivotCalibration(G):
    Gt = G[0].transpose()
    G0 = np.mean(Gt, axis=1).reshape((-1,1))
    g = Gt - G0
    F_G = registrationArunMethod(g.transpose(), G[0], "G")
    print(np.linalg.det(F_G.R))
    negI = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
    R_Gs = np.hstack((F_G.R, negI))
    t_Gs = F_G.p.coords.transpose()[..., np.newaxis]
    for i in range(1, G.shape[0]):
        Gt = G[i].transpose()
        # G0 = np.mean(Gt, axis=1).reshape((-1,1))
        g = Gt - G0
        F_G = registrationArunMethod(g.transpose(), G[i], "G")
        # print(F_G.R)
        # print(np.matmul(F_G.R, g[:,0]) + F_G.p.coords)
        # print(Gt[:,0])
        print(np.linalg.det(F_G.R))
        R_G = np.hstack((F_G.R, negI))
        R_Gs = np.vstack((R_Gs, R_G))
        t_Gs = np.vstack((t_Gs, F_G.p.coords.transpose()[..., np.newaxis]))
    t_Gs = t_Gs * -1
    # x = np.matmul((np.matmul(np.matmul(R_Gs.transpose(), R_Gs), R_Gs.transpose())), t_Gs)
    x = np.matmul(np.matmul(np.linalg.inv(np.matmul(R_Gs.transpose(), R_Gs)), R_Gs.transpose()), t_Gs)
    print(x)
    # print(np.matmul(F_G.R,x[0:3]) + F_G.p.coords.transpose()[..., np.newaxis])
    # print(x[3:6])


G = read_empivot("./PA1 Student Data/pa1-debug-a-empivot.txt")
pivotCalibration(G)