import numpy as np
from Point3d import Point3d
from Registration import registrationArunMethod

def GtoEM(G, p_tip):
    B = np.empty((G.shape[0], 3))

    # defining a local coordinate system
    Gt = G[0].transpose()
    G0 = np.mean(Gt, axis=1).reshape((-1,1))
    g = Gt - G0

    # creating the [R -I][p_tip p_dimple]^T = [-t] matrix for each frame
    for i in range(G.shape[0]):
        F_G = registrationArunMethod(g.transpose(), G[i], "G")
        B[i] = F_G.R.dot(p_tip.coords) + F_G.p.coords
        
    return B

# def EMToCT():
