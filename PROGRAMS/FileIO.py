import numpy as np

# calbody.txt describes the calibration object
def read_calbody(filename):
    '''
    Returns:
    d: N_D x 3 (N_D points)
    a: N_A x 3 (N_A points)
    c: N_C x 3 (N_C points)

    Definitions:
    N_D: number of optical markers on EM base
    N_A: number of optical markers on calibration objects
    N_C: number of EM markers on calibration object
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    d = np.empty((int(params[0]), 3))
    a = np.empty((int(params[1]), 3))
    c = np.empty((int(params[2]), 3))
    for i in range(int(params[0])):
        point = file.readline().replace(' ','').split(',')
        d[i] = point
    for i in range(int(params[1])):
        point = file.readline().replace(' ','').split(',')
        a[i] = point
    for i in range(int(params[2])):
        point = file.readline().replace(' ','').split(',')
        c[i] = point
    return d, a, c

# calreadings.txt provides the values read by the sensor
def read_calreadings(filename):
    '''
    Returns:
    D: N_frames x N_D x 3 (N_frames x N_D points)
    A: N_frames x N_A x 3 (N_frames x N_A points)
    C: N_frames x N_C x 3 (N_frames x N_C points)

    Definitions:
    N_D: number of optical markers on EM base
    N_A: number of optical markers on calibration objects
    N_C: number of EM markers on calibration object
    N_frames: number of data frames
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    D = np.empty((int(params[3]), int(params[0]), 3))
    A = np.empty((int(params[3]), int(params[1]), 3))
    C = np.empty((int(params[3]), int(params[2]), 3))
    for j in range(int(params[3])):
        for i in range(int(params[0])):
            point = file.readline().replace(' ','').split(',')
            D[j][i] = point
        for i in range(int(params[1])):
            point = file.readline().replace(' ','').split(',')
            A[j][i] = point
        for i in range(int(params[2])):
            point = file.readline().replace(' ','').split(',')
            C[j][i] = point
    return D, A, C

# empivot.txt provides values read by the sensor
def read_empivot(filename):
    '''
    Returns:
    G: N_frames x N_G x 3 (N_frames x N_G points)

    Definitions:
    N_G: number of EM markers on probe
    N_frames: number of data frames
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    G = np.empty((int(params[1]), int(params[0]), 3))
    for j in range(int(params[1])):
        for i in range(int(params[0])):
            point = file.readline().replace(' ','').split(',')
            G[j][i] = point
    return G

# optpivot.txt provides values read by the sensor
def read_optpivot(filename):
    '''
    Returns:
    D: N_frames x N_D x 3 (N_frames x N_D points)
    H: N_frames x N_H x 3 (N_frames x N_H points)

    Definitions:
    N_D: number of optical markers on EM base
    N_H: number of optical markers on probe
    N_frames: number of data frames
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    D = np.empty((int(params[2]), int(params[0]), 3))
    H = np.empty((int(params[2]), int(params[1]), 3))
    for j in range(int(params[2])):
        for i in range(int(params[0])):
            point = file.readline().replace(' ','').split(',')
            D[j][i] = point
        for i in range(int(params[1])):
            point = file.readline().replace(' ','').split(',')
            H[j][i] = point
    return D, H

# output1.txt provides output file for problem 1
def read_output1(filename):
    '''
    Returns:
    C_exp: N_frames x N_C x 3 (N_frames x N_C points)
    expected C coordinates
    P_em: Px, Py, Pz
    estimated post position with EM probe pivot calibration
    E_opt: Px, Py, Px
    estimated post position with optical probe pivot calibration

    Definitions:
    N_C: number of EM markers on cal object
    N_frames: number of data frames
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    C_exp = np.empty((int(params[1]), int(params[0]), 3))
    P_em = file.readline().replace("\n","").replace(' ','').split(',')
    P_opt = file.readline().replace("\n","").replace(' ','').split(',')
    for j in range(int(params[1])):
        for i in range(int(params[0])):
            point = file.readline().replace("\n","").replace(' ','').split(',')
            C_exp[j][i] = point
    return C_exp, P_em, P_opt

# ct-fiducials.txt describes CT fiducial coordinates b_j
def read_ctfiducials(filename):
    '''
    Returns:
    b: N_B x 3 (N_B points)
    coordinates of CT fiducials

    Definitions:
    N_B: number of CT fiducials
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    b = np.empty((int(params[0]), 3))
    for i in range(int(params[0])):
        point = file.readline().replace("\n","").replace(' ','').split(',')
        b[i] = point
    return b

# em-fiducials.txt describes frames of data in which the probe is in contact
# with the corresponding CT fiducials
def read_emfiducials(filename):
    '''
    Returns:
    G: N_frames x N_G x 3 (N_frames x N_G points)
    coordinates of CT fiducials

    Definitions:
    N_G: number of EM markers on probe
    N_frames: number of data frames
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    G = np.empty((int(params[1]), int(params[0]), 3))
    for j in range(int(params[1])):
        for i in range(int(params[0])):
            point = file.readline().replace("\n","").replace(' ','').split(',')
            G[j][i] = point
    return G

# em-nav.txt describes frames of data defining test points - able to find
# corresponding positions of the probe tip with respect to CT coordinates
def read_emnav(filename):
    '''
    Returns:
    G: N_frames x N_G x 3 (N_frames x N_G points)
    coordinates of CT fiducials

    Definitions:
    N_G: number of EM markers on probe
    N_frames: number of data frames
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    G = np.empty((int(params[1]), int(params[0]), 3))
    for j in range(int(params[1])):
        for i in range(int(params[0])):
            point = file.readline().replace("\n","").replace(' ','').split(',')
            G[j][i] = point
    return G

# output2.txt gives positions of prove tip in CT coordinates, corresponding to
# frames of data in em-nav.txt
def read_output2(filename):
    '''
    Returns:
    b: N_frames x 3 (N_frames points)
    coordinates of CT fiducials

    Definitions:
    N_frames: number of frames of data
    '''
    file = open(filename, 'r')
    params = file.readline().replace(' ','').split(',')
    v = np.empty((int(params[0]), 3))
    for i in range(int(params[0])):
        point = file.readline().replace("\n","").replace(' ','').split(',')
        v[i] = point
    return v