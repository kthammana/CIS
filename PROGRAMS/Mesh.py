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
        return [self.v_coors[t_v[0]].transpose(), self.v_coors[t_v[1]].transpose(), self.v_coors[t_v[2]].transpose()]
    

# can add functions as necessary
# can change variable names to given V,i,n if that is easier