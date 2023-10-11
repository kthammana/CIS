import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from Point3d import Point3d
from Frame import Frame

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
    # d, a, c = [], [], []
    d = np.empty((int(params[0]), 3))
    a = np.empty((int(params[1]), 3))
    c = np.empty((int(params[2]), 3))
    for i in range(int(params[0])):
        # point = Point3d(file.readline().replace(' ','').split(','))
        # d.append(point)
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
    # print(params)
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
    # print(params)
    G = np.empty((int(params[1]), int(params[0]), 3))
    for j in range(int(params[1])):
        for i in range(int(params[0])):
            point = file.readline().replace(' ','').split(',')
            G[j][i] = point
    # print(G.shape)
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
    print(params)
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

# output.txt provides output file for problem 1
def read_output(filename):
    '''
    Returns:
    C_exp: N_frames x N_C x 3 (N_frames x N_D points)
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
    # print(params)
    C_exp = np.empty((int(params[1]), int(params[0]), 3))
    P_em = file.readline().replace(' ','').split(',')
    P_opt = file.readline().replace(' ','').split(',')
    for j in range(int(params[1])): # changed from params[2]
        for i in range(int(params[0])):
            point = file.readline().replace(' ','').split(',')
            C_exp[j][i] = point
    return C_exp, P_em, P_opt

# read_calbody("./PA1 Student Data/pa1-debug-a-calbody.txt")
# read_calreadings("./PA1 Student Data/pa1-debug-a-calreadings.txt")
# read_empivot("./PA1 Student Data/pa1-debug-a-empivot.txt")