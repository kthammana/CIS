import numpy as np

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

# Driver program to test above functions
if __name__ == '__main__':
	root = None
	points = [[3, 6, 1], [17, 15, 0], [13, 15, 6], [6, 12, 4], [9, 1, 2], [2, 7, 3], [10, 19, 5]]

	n = len(points)

	for i in range(n):
		root = insert(root, points[i])

	point1 = [10, 19, 5]
	if search(root, point1):
		print("Found")
	else:
		print("Not Found")

	point2 = [12, 19, 5]
	if search(root, point2):
		print("Found")
	else:
		print("Not Found")

