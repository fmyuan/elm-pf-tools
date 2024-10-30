#!/usr/bin/env python
import numpy as np

#---
# nearest_neibour for earth surface using kdtree
def nearest_using_kdtree(data, latname='Latitude', lonname='Longitude',kpt=2):
    import scipy.spatial as spatial
    #"Based on https://stackoverflow.com/q/43020919/190597"
    R = 6367000.0 # meters
    def dist_to_arclength(chord_length):
        """
        https://en.wikipedia.org/wiki/Great-circle_distance
        Convert Euclidean chord length to great circle arc length
        """
        central_angle = 2*np.arcsin(chord_length/(2.0*R)) 
        arclength = R*central_angle
        return arclength

    phi = np.deg2rad(data[latname])
    theta = np.deg2rad(data[lonname])
    data['x'] = R * np.cos(phi) * np.cos(theta)
    data['y'] = R * np.cos(phi) * np.sin(theta)
    data['z'] = R * np.sin(phi)
    points = list(zip(data['x'],data['y'],data['z']))
    tree = spatial.KDTree(points)
    distance, index = tree.query(points, k=kpt)
    
    #return nearest points other than itself
    return dist_to_arclength(distance[:, 1:]), index[:,1:]

#--------------------------------------------------------------------

