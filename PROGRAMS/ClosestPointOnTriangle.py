import numpy as np

def projectOnSegment(c, x, y):
    # print("closest point on triangle")
    c = c.transpose()
    l = (((c - x).dot(y - x))/((y - x).dot(y - x)))[0]
    l_seg = max(0, min(l, 1))
    c = x + l_seg*(y - x)
    return c

def findClosestPointOnTriangle(a, v_coors):
    p, q, r = v_coors
    b = (a - p).transpose()[..., np.newaxis]
    A = np.array([q-p, r-p]).transpose()
    x = np.linalg.inv(A.transpose() @ A) @ (A.transpose() @ b)
    l = x[0][0]
    u = x[1][0]
    c = p + l*(q - p) + u*(r - p)
    c = c.transpose()[..., np.newaxis]
    if l >= 0 and u >= 0 and l + u <= 1:
        # print("closest point in triangle")
        return c[:,0]
    elif l < 0:
        return projectOnSegment(c, r, p)
    elif u < 0:
        return projectOnSegment(c, p, q)
    elif l + u > 0:
        return projectOnSegment(c, q, r)


a = np.array([1.75,0.5,1])
p = np.array([1,0,0])
q = np.array([2,0,0])
r = np.array([2,1,0])
v_coords = np.asarray([p, q, r])
c = findClosestPointOnTriangle(a, v_coords)
# print(c)