import numpy as np
from FileIO import read_calbody, read_calreadings, read_empivot, read_optpivot, read_output
from Registration import registrationArunMethod
from PivotCalibration import pivotCalibration
from OpticalCalibration import opticalCalibration

d, a, c = read_calbody("../PA1 Student Data/pa1-debug-c-calbody.txt")
D, A, C = read_calreadings("../PA1 Student Data/pa1-debug-c-calreadings.txt")
for i in range(D.shape[0]):
    F_D = registrationArunMethod(d, D[i], "D")
    F_A = registrationArunMethod(a, A[i], "A")
    F_DA = F_D.inverse() * F_A
    for j in range(c.shape[0]):
        C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
    # do we want to output error? take average of C over dim 0, compare to C0_exp at each iteration, add difference and divide by 27?
    # that would be good

G = read_empivot("./PA1 Student Data/pa1-debug-g-empivot.txt")
P_em_exp = pivotCalibration(G)
print('Calculated output:', P_em_exp.__str__())
__,P_em,__ = read_output("./PA1 Student Data/pa1-debug-g-output1.txt")
print('Expected output:',P_em)
print('Error:',P_em_exp.error(P_em),'mm')

D,H = read_optpivot("./PA1 Student Data/pa1-debug-g-optpivot.txt")
d,_,_ = read_calbody("./PA1 Student Data/pa1-debug-g-calbody.txt")
P_opt_exp = opticalCalibration(d,D,H)
print('Calculated output:', P_opt_exp.__str__())
__,__,P_opt = read_output("./PA1 Student Data/pa1-debug-g-output1.txt")
print('Expected output:',P_opt)
print('Error:',P_opt_exp.error(P_opt),'mm')