import numpy as np
from FileIO import read_calbody, read_calreadings, read_empivot, read_optpivot, read_output1, read_ctfiducials, read_emfiducials, read_emnav, read_output2, read_mesh, read_output3, read_samplereadings, read_probbody
from Registration import registrationArunMethod
from EMPivotCalibration import pivotCalibration, GtoEM
from OpticalPivotCalibration import opticalCalibration
from Point3d import Point3d
from DistortionCorrection import bernstein, calcDistortionCorrection, correctDistortion

def testBernstein():
    # b_2,5(x) = 10x^2*(1 - x)^3 --> n = 5, k = 2, v = x
    x = 5
    expected = 10 * x**2 * (1-x)**3
    actual = bernstein(5,2,x)
    print('Bernstein error:',abs(expected-actual))

def testRegistration():
    print("Registration Test:")
    d, a, c = read_calbody("PA1 Student Data/pa1-debug-g-calbody.txt")
    # use given d (from PA1 data set g) and artificial D
    artiD = np.array([[0.3404, 0.5853, 0.2238],[-100.2885, -80.1432, 76.7532], 
                [111.0513, -62.0671, 79.7082],[10.4224, -142.7955, 156.2376],
                [11.1531, -109.2218, -101.3906],[-89.4758, -189.9503, -24.8612],
                [121.8640, -171.8742, -21.9062],[21.2351, -252.6026, 54.6232]])
    artiF_D = registrationArunMethod(d, artiD, "D")
    # artificial R and p calculated on MATLAB
    artiR = np.array([[0.0721, 0.7381, -0.6709],[-0.7320, -0.4177, -0.5382],[-0.6774, 0.5299, 0.5102]])
    print('Calculated R:', artiF_D.R)
    print('Expected R:',artiR)
    print('Calculated p:',artiF_D.p.__str__())
    print('Expected p: 0.3404, 0.5853, 0.2238')

def testDistortionReconstructsCexp(dataset):
    print("Distortion Test:")
    # calculating expected Cs
    d, a, c = read_calbody(dataset+"-calbody.txt")
    D, A, C = read_calreadings(dataset+"-calreadings.txt")
    C_expected = np.empty(C.shape)
    for i in range(D.shape[0]):
        F_D = registrationArunMethod(d, D[i], "D")
        F_A = registrationArunMethod(a, A[i], "A")
        F_DA = F_D.inverse() * F_A
        for j in range(c.shape[0]):
            C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
            C_expected[i][j] = [C0_exp[0][0], C0_exp[1][0], C0_exp[2][0]]

    C_expected_stacked = np.vstack(C_expected)
    C_stacked = np.vstack(C)
    C_errors = np.zeros(C_stacked.shape[0])
    coef, q_min, q_max = calcDistortionCorrection(C_expected_stacked, C_stacked, 3)
    C_undistorted = correctDistortion(C_stacked, coef, q_min, q_max, 3)
    for i in range(C_stacked.shape[0]): # N_C
        P_Cexp = Point3d("C", C_expected_stacked[i])
        C_errors[i] = P_Cexp.error(C_undistorted[i])
    print('Average calibration error per point:', np.mean(C_errors), 'mm')

def printPA2OutputErrors(dataset):
    print("PA2 Output Errors:")
    C_exp,P_em,_ = read_output1(dataset+"-output1.txt")

    # calculating expected Cs
    d, a, c = read_calbody(dataset+"-calbody.txt")
    D, A, C = read_calreadings(dataset+"-calreadings.txt")
    C_expected = np.zeros(C.shape)
    for i in range(D.shape[0]):
        F_D = registrationArunMethod(d, D[i], "D")
        F_A = registrationArunMethod(a, A[i], "A")
        F_DA = F_D.inverse() * F_A
        for j in range(c.shape[0]):
            C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
            C_expected[i][j] = [C0_exp[0][0], C0_exp[1][0], C0_exp[2][0]]

    coef, q_min, q_max = calcDistortionCorrection(np.vstack(C_expected), np.vstack(C), 5)
    C_errors = np.zeros(C.shape[0:2]) # stores error of each frame
    for i in range(C.shape[0]): # N_frames
        for j in range(C.shape[1]): # N_C
            P_Cexp = Point3d("C", C_exp[i][j])
            C_errors[i][j] = P_Cexp.error(C_expected[i][j])
    print('Average calibration error per point:', np.mean(C_errors), 'mm')

    # EM pivot calibration
    G = read_empivot(dataset+"-empivot.txt")
    G_corr = np.empty(G.shape)
    for i in range(G_corr.shape[0]):
        G_corr[i] = correctDistortion(G[i], coef, q_min, q_max, 5)
    P_em_exp, P_tip, g = pivotCalibration(G_corr)
    print('EM Pivot Error:',P_em_exp.error(P_em),'mm')

    # Fiducial point registration
    b = read_ctfiducials(dataset+"-ct-fiducials.txt")
    G_EM = read_emfiducials(dataset+"-em-fiducialss.txt")
    G_EM_corr = np.empty(G_EM.shape)
    for i in range(G_EM_corr.shape[0]):
        G_EM_corr[i] = correctDistortion(G_EM[i], coef, q_min, q_max, 5)
    B = GtoEM(G_EM_corr, P_tip, g)

    # Calculate F_reg
    F_reg = registrationArunMethod(B, b, "EM")

    # Compute tip location w.r.t. CT image
    G_nav = read_emnav(dataset+"-EM-nav.txt")
    G_nav_corr = np.empty(G_nav.shape)
    for i in range(G_nav_corr.shape[0]):
        G_nav_corr[i] = correctDistortion(G_nav[i], coef, q_min, q_max, 5)
    V = GtoEM(G_nav_corr, P_tip, g)
    v_exp = np.empty([G_nav.shape[0],3])
    for i in range(V.shape[0]):
        v_exp[i] = F_reg.R.dot(V[i]) + F_reg.p.coords

    # Calculate average v error
    v = read_output2(dataset+"-output2.txt")
    errors = np.empty(v.shape[0])
    for i in range(v.shape[0]):
        p_v = Point3d("CT",v[i])
        errors[i] = p_v.error(v_exp[i])
    print('Avg v error:', np.mean(errors))

def printPA1OutputErrors(dataset):
    print("PA1 Output Errors:")
    C_exp,P_em,P_opt = read_output1(dataset+"-output1.txt")

    # calculating expected Cs
    d, a, c = read_calbody(dataset+"-calbody.txt")
    D, A, C = read_calreadings(dataset+"-calreadings.txt")
    error = 0
    for i in range(D.shape[0]):
        F_D = registrationArunMethod(d, D[i], "D")
        F_A = registrationArunMethod(a, A[i], "A")
        F_DA = F_D.inverse() * F_A
        for j in range(c.shape[0]):
            C0_exp = np.matmul(F_DA.R, c[j].transpose()[..., np.newaxis]) + F_DA.p.coords.transpose()[..., np.newaxis]
            gen_point = Point3d("DA", C0_exp[0][0], C0_exp[1][0], C0_exp[2][0])
            error += gen_point.error(C[i][j])
    error = error/(C_exp.shape[0]*C_exp.shape[1])
    print('Average Expected C Error:',error)

    # EM pivot calibration
    G = read_empivot(dataset+"-empivot.txt")
    P_em_exp, P_tip, _ = pivotCalibration(G)
    print('EM Calculated output:', P_em_exp.__str__())
    print('EM Expected output:',P_em)
    print('EM Pivot Error:',P_em_exp.error(P_em),'mm')

    # Optical tracker pivot calibration
    D,H = read_optpivot(dataset+"-optpivot.txt")
    d,_,_ = read_calbody(dataset+"-calbody.txt")
    P_opt_exp = opticalCalibration(d,D,H)
    print('OPT Calculated output:', P_opt_exp.__str__())
    print('OPT Expected output:',P_opt)
    print('OPT Pivot Error:',P_opt_exp.error(P_opt),'mm')
    
def printPA3OutputErrors(dataset):
    print("PA3 Output Errors:")
    
    # test I/O functions
    Y_A, t_A = read_probbody("PA345 Student Data/Problem3-BodyA.txt")
    Y_B, t_B = read_probbody("PA345 Student Data/Problem3-BodyB.txt")
    V, i, n = read_mesh("PA345 Student Data/Problem3MeshFile.sur")
    M = read_samplereadings(dataset+"-SampleReadingsTest.txt")
    d, c, mag = read_output3(dataset+"-Output.txt")