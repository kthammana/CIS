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
        print("F_G rot:\n", F_G.R)
        print("F_G trans: ", F_G.p)
        print("p_tip: ", p_tip)
        B[i] = F_G.R @ p_tip.coords + F_G.p.coords
        print("Bi: ", B[i])
        
    return B, g

def EMtoCT(G, p_tip, g):
    B = np.empty((G.shape[0], 3))

    # creating the [R -I][p_tip p_dimple]^T = [-t] matrix for each frame
    for i in range(G.shape[0]):
        F_G = registrationArunMethod(g.transpose(), G[i], "G")
        print("F_G rot:\n", F_G.R)
        print("F_G trans: ", F_G.p)
        print("p_tip: ", p_tip)
        B[i] = F_G.R @ p_tip.coords + F_G.p.coords
        print("Bi: ", B[i])
        
    return B

# def EMToCT():
