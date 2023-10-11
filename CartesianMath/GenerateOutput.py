import numpy as np
from Point3d import Point3d
from Frame import Frame
from FileIO import read_calbody, read_calreadings, read_empivot, read_optpivot, read_output
from Registration import registrationArunMethod

# output file:
    # N_C , N_frames, NAME-OUTPUT1.TXT
    # P_EM
    # P_opt
    # C_1 through C_N_frames
d, a, c = read_calbody("../PA1 Student Data/pa1-unknown-h-calbody.txt")
D, A, C = read_calreadings("../PA1 Student Data/pa1-unknown-h-calreadings.txt")
G = read_empivot("../PA1 Student Data/pa1-unknown-h-empivot.txt")
D,H = read_optpivot("../PA1 Student Data/pa1-unknown-h-optpivot.txt")
f = open("h-outputs.txt", "x")

for i in range(D.shape[0]): # N_frames
    F_D = registrationArunMethod(d, D[i], "D")
    F_A = registrationArunMethod(a, A[i], "A")
    F_DA = F_D.inverse() * F_A
    for j in range(c.shape[0]): # N_C
        C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
        
f.close()
