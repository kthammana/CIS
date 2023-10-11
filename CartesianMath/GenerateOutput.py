import numpy as np
from Point3d import Point3d
from FileIO import read_calbody, read_calreadings, read_empivot, read_optpivot
from Registration import registrationArunMethod
from PivotCalibration import pivotCalibration
from OpticalCalibration import opticalCalibration

# output file:
    # N_C , N_frames, NAME-OUTPUT1.TXT
    # P_EM
    # P_opt
    # C_1 through C_N_frames
filename = "pa1-unknown-k-output1.txt" # output file name

# Calculate P_EM
d, a, c = read_calbody("./PA1 Student Data/pa1-unknown-k-calbody.txt")
G = read_empivot("./PA1 Student Data/pa1-unknown-k-empivot.txt")
P_EM = pivotCalibration(G)

# Calculate P_opt
D,H = read_optpivot("./PA1 Student Data/pa1-unknown-k-optpivot.txt")
P_opt = opticalCalibration(d,D,H)

# Calculate C_exp and write all variables to output file
D, A, C = read_calreadings("./PA1 Student Data/pa1-unknown-k-calreadings.txt")
f = open("./OUTPUT/"+filename, "w")
N_frames = D.shape[0]
N_C = c.shape[0]
f.write(str(N_C) + ", " + str(N_frames) + ", " + filename + "\n")
f.write(P_EM.__str__() + "\n")
f.write(P_opt.__str__() + "\n")
for i in range(N_frames): # N_frames
    F_D = registrationArunMethod(d, D[i], "D")
    F_A = registrationArunMethod(a, A[i], "A")
    F_DA = F_D.inverse() * F_A
    for j in range(N_C): # N_C
        C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
        point = Point3d("EM",C0_exp[0][0], C0_exp[1][0], C0_exp[2][0])
        f.write(point.__str__() + "\n")
f.close()
