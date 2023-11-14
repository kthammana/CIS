import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
    
# can add functions as necessary
# can change variable names to given V,i,n if that is easier\

    # WIP NEED TO MODIFY
# The following code is contributed by Prajwal Kandekar, as available on
# GeeksforGeeks and has been modified as necessary to accomodate the Mesh class

# Number of dimensions
k = 3

# A structure to represent node of kd tree
class Node:
	def __init__(self, point):
		self.point = point
		self.left = None
		self.right = None

# A method to create a node of K D tree
def newNode(point):
	return Node(point)

# Inserts a new node and returns root of modified tree
# The parameter depth is used to decide axis of comparison
def insertRec(root, point, depth):
	# Tree is empty?
	if not root:
		return newNode(point)

	# Calculate current dimension (cd) of comparison
	cd = depth % k

	# Compare the new point with root on current dimension 'cd'
	# and decide the left or right subtree
	if point[cd] < root.point[cd]:
		root.left = insertRec(root.left, point, depth + 1)
	else:
		root.right = insertRec(root.right, point, depth + 1)

	return root

# Function to insert a new point with given point in
# KD Tree and return new root. It mainly uses above recursive
# function "insertRec()"
def insert(root, point):
	return insertRec(root, point, 0)

# A utility method to determine if two Points are same
# in K Dimensional space
def arePointsSame(point1, point2):
	# Compare individual coordinate values
	for i in range(k):
		if point1[i] != point2[i]:
			return False

	return True

# Searches a Point represented by "point[]" in the K D tree.
# The parameter depth is used to determine current axis.
def searchRec(root, point, depth):
	# Base cases
	if not root:
		return False
	if arePointsSame(root.point, point):
		return True

	# Current dimension is computed using current depth and total
	# dimensions (k)
	cd = depth % k

	# Compare point with root with respect to cd (Current dimension)
	if point[cd] < root.point[cd]:
		return searchRec(root.left, point, depth + 1)

	return searchRec(root.right, point, depth + 1)

# Searches a Point in the K D tree. It mainly uses
# searchRec()
def search(root, point):
	# Pass current depth as 0
	return searchRec(root, point, 0)

def plot_kd_tree_3d(node, ax, xmin, xmax, ymin, ymax, zmin, zmax, depth=0):
    if node is not None:
        cd = depth % k

        if cd == 0:
            ax.plot([node.point[0], node.point[0]], [ymin, ymax], [zmin, zmax], color='black', linestyle='-', linewidth=1)
            plot_kd_tree_3d(node.left, ax, xmin, node.point[0], ymin, ymax, zmin, zmax, depth + 1)
            plot_kd_tree_3d(node.right, ax, node.point[0], xmax, ymin, ymax, zmin, zmax, depth + 1)
        elif cd == 1:
            ax.plot([xmin, xmax], [node.point[1], node.point[1]], [zmin, zmax], color='black', linestyle='-', linewidth=1)
            plot_kd_tree_3d(node.left, ax, xmin, xmax, ymin, node.point[1], zmin, zmax, depth + 1)
            plot_kd_tree_3d(node.right, ax, xmin, xmax, node.point[1], ymax, zmin, zmax, depth + 1)
        else:
            ax.plot([xmin, xmax], [ymin, ymax], [node.point[2], node.point[2]], color='black', linestyle='-', linewidth=1)
            plot_kd_tree_3d(node.left, ax, xmin, xmax, ymin, ymax, zmin, node.point[2], depth + 1)
            plot_kd_tree_3d(node.right, ax, xmin, xmax, ymin, ymax, node.point[2], zmax, depth + 1)


def plot_points_3d(ax, points):
    points = np.array(points)
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color='red', marker='o', label='Points')

def main():
    points = [[3, 6, 1], [17, 15, 0], [13, 15, 6], [6, 12, 4], [9, 1, 2], [2, 7, 3], [10, 19, 5]]
    root = None

    for point in points:
        root = insertRec(root, point, 0)

    xmin, xmax = min(p[0] for p in points) - 1, max(p[0] for p in points) + 1
    ymin, ymax = min(p[1] for p in points) - 1, max(p[1] for p in points) + 1
    zmin, zmax = min(p[2] for p in points) - 1, max(p[2] for p in points) + 1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    plot_kd_tree_3d(root, ax, xmin, xmax, ymin, ymax, zmin, zmax)
    plot_points_3d(ax, points)

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    ax.set_title('3D k-d Tree Visualization')
    ax.legend()
    plt.show()

def calcDistance(x, y):
    return np.sqrt((x[0] - y[0])**2 + (x[1] - y[1])**2 + (x[2] - y[2])**2)

if __name__ == '__main__':
    main()
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
    c_k = np.empty((a.shape[0], 3))
    for i in range(d_k.shape[0]):
        shortest_dist = np.infty
        for j in range(ind.shape[0]):
            c = findClosestPointOnTriangle(d_k[i], mesh.getVerticesOfTriangle(j))
            dist = calcDistance(d_k[i], c)
            if dist <= shortest_dist:
                closest_point = c
                shortest_dist = dist
        c_k[i] = closest_point
        c_error += calcDistance(c_exp[i], c_k[i])
        mag_error += (np.abs(mag[i]-shortest_dist))

    print("d_k error: ", d_error/d_k.shape[0])
    print("c_k error: ", c_error/c_k.shape[0])
    print("mag error: ", mag_error/c_k.shape[0])
