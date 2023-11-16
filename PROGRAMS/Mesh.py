import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from FileIO import read_mesh, read_output3, read_samplereadings, read_probbody
from Registration import registrationArunMethod
from ClosestPointOnTriangle import findClosestPointOnTriangle

class Mesh(object):
    
    def __init__(self, v_coors, v_inds, n_inds):
        """
            
        Parameters
        ----------
        v_coors : numpy array (N_vertices x 3)
            Coordinates of the triangles' vertices
        v_inds : numpy array (N_triangles x 3)
            Indices of the triangles' vertices
        n_inds : numpy array (N_triangles x 3)
            Indices of the (opposite) neighboring triangles' vertices
    
        Returns
        -------
        None.
    
        """
        self.v_coors = v_coors # V
        self.v_inds = v_inds # i
        self.n_inds = n_inds # n

    def getVerticesOfTriangle(self, t_idx):
        t_v = self.v_inds[t_idx]
        return [self.v_coors[int(t_v[0])].transpose(), self.v_coors[int(t_v[1])].transpose(), self.v_coors[int(t_v[2])].transpose()]
    
    def calcCentroid(self, vertices):
        x = 0
        y = 0
        z = 0
        for i in range(3):
            x += vertices[i][0]
            y += vertices[i][1]
            z += vertices[i][2]
        x /= 3
        y /= 3
        z /= 3
        return np.array([x,y,z])

# The following code is contributed by Prajwal Kandekar, as available on
# GeeksforGeeks and has been modified as necessary to accomodate the Mesh class

# Number of dimensions
k = 3

'''
Using a generic KD tree starter code, I modified it so that the node's point
is the centroid of its triangle and gave it a triangle variable to then hold
the vertices' coordinate information. To find the nearest node, the given point
is compared to each triangle's centroid to find the nearest one, then I call
findClosestPointOnTriangle(...) to calculate the true distance of the point
from the triangle, rather than distance from the centroid.
'''
# A structure to represent node of kd tree
class Node:
	def __init__(self, point, triangle, idx):
            self.point = point # centroid of triangle
            self.triangle = triangle
            self.left = None
            self.right = None
            self.idx = idx

# Inserts a new node and returns root of modified tree
# The parameter depth is used to decide axis of comparison
def insert(root, point, triangle, depth, idx):
    if not root: # Empty tree
        return Node(point, triangle, idx)

    cd = depth % k

    if point[cd] < root.point[cd]:
        root.left = insert(root.left, point, triangle, depth + 1, idx)
    else:
        root.right = insert(root.right, point, triangle, depth + 1, idx)

    return root

def distance_squared(point1, point2):
    return sum((p1 - p2) ** 2 for p1, p2 in zip(point1, point2))

def searchTree(root, point, depth, best_node, best_distance):
    if not root or root is None:
        return best_node, best_distance

    cd = depth % k
    next_best = best_node # None
    next_best_distance = best_distance

    if point[cd] < root.point[cd]:
        next_best, next_best_distance = searchTree(root.left, point, depth + 1, next_best, next_best_distance)
    else:
        next_best, next_best_distance = searchTree(root.right, point, depth + 1, next_best, next_best_distance)

    current_distance = distance_squared(root.point, point)

    if current_distance < next_best_distance:
        next_best = root
        next_best_distance = current_distance

    other_child = root.right if (point[cd] < root.point[cd]) else root.left
    if other_child:
        other_child_distance = (point[cd] - root.point[cd]) ** 2
        if other_child_distance < next_best_distance:
            next_best, next_best_distance = searchTree(other_child, point, depth + 1, next_best, next_best_distance)

    return next_best, next_best_distance

# Uses 
def search(root, point):
    nearest_node, __ = searchTree(root, point, 0, None, float('inf'))
    return nearest_node

def calcDistance(x, y):
    return np.sqrt((x[0] - y[0])**2 + (x[1] - y[1])**2 + (x[2] - y[2])**2)

if __name__ == '__main__':
    print("PA3 Output Errors:")
    
    dataset = "PA345 Student Data/PA3-F-Debug"
    
    # test I/O functions
    A, A_tip, N_A = read_probbody("PA345 Student Data/Problem3-BodyA.txt")
    B, B_tip, N_B = read_probbody("PA345 Student Data/Problem3-BodyB.txt")
    V, ind, n = read_mesh("PA345 Student Data/Problem3MeshFile.sur")
    mesh = Mesh(V, ind, n)
    a, b = read_samplereadings(dataset+"-SampleReadingsTest.txt", N_A, N_B)
    d_exp, c_exp, mag = read_output3(dataset+"-Output.txt")
    d_error = 0
    c_error = 0
    mag_error = 0

    d_k = np.empty((a.shape[0], 3))
    A_tip = A_tip.transpose()[..., np.newaxis]
    for i in range(a.shape[0]): # N_samples:
        F_Ai = registrationArunMethod(A, a[i], "A")
        F_Bi = registrationArunMethod(B, b[i], "B")
        F_BA = F_Bi.inverse() * F_Ai
        d_k[i] = (F_BA.R @ A_tip)[:,0] + F_BA.p.coords
        d_error += calcDistance(d_exp[i], d_k[i])

    # for PA3, F_reg = 1
    # TO-DO: construct octree/kdtree for searching for points
    # linear search to find the closest points to d_k
    root = None
    
    # Run time should now be O(n) instead of O(n^2)
    for i in range(ind.shape[0]):
        verts = mesh.getVerticesOfTriangle(i)
        root = insert(root, mesh.calcCentroid(verts), verts, 0, i)
    
    c_k = np.empty((a.shape[0], 3))
    for i in range(d_k.shape[0]):
        nearest_node = search(root, d_k[i])
        c_k[i] = findClosestPointOnTriangle(d_k[i], nearest_node.triangle)
        shortest_dist = calcDistance(d_k[i], c_k[i])
        # print(nearest_node.idx)
        c_error += calcDistance(c_exp[i], c_k[i])
        mag_error += (np.abs(mag[i]-shortest_dist))
    
    # Slight difference in these values
    print("d_k error: ", d_error/d_k.shape[0])
    print("c_k error: ", c_error/c_k.shape[0])
    print("mag error: ", mag_error/c_k.shape[0])
