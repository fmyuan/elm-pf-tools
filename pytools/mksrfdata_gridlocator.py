#!/usr/bin/env python
import os, sys, csv, time, math
from optparse import OptionParser
import numpy as np
from scipy import interpolate
from cmath import nan, inf
from builtins import int
from scipy.ndimage import interpolation
from numpy import double
from _operator import index

#----- 
#----- nearest_neibour for earth surface using kdtree
#----- latitude/longitude
def nearest_pts_latlon_kdtree(allpoints, excl_points=0, latname='Latitude', lonname='Longitude'):
    import scipy.spatial as spatial
    from math import radians, cos, sin, asin, sqrt
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

    def kpt_default(data, latname='Latitude', lonname='Longitude',kpt=2):
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

    # nearest search, when data points are scattered (NOT in grided)
    # because mixed points, must make sure nearest points are NONE of first no. of self_points
    AllNearest=False
    Kpts = 2
    idx_void = None
    while not AllNearest:
        dist,idx=kpt_default(allpoints, latname='Latitude', lonname='Longitude',kpt=Kpts)

        if excl_points<=1:
            AllNearest=True
            return dist, idx
            
        else:
            #for lat/lon real index, must adjust 'idx' in mixed points 
            if idx_void is None:
                idx0  = np.copy(idx)   # saved for non-excluded points
                dist0 = np.copy(dist)
                idx_excl = idx[0:excl_points]
                dist_excl= dist[0:excl_points]
                # check if contains exclusive points
                idx_void= np.where(idx_excl<excl_points)
                if (len(idx_void[0])<=0): AllNearest = True

            else:
                # fill in the negative idx only
                idx_excl[idx_void] = (idx[...,0:excl_points])[idx_void]
                dist_excl[idx_void] = (dist[...,0:excl_points])[idx_void]
                
                # check if still exclusive points
                idx_void= np.where(idx_excl<excl_points)
                if (len(idx_void[0])>0):
                    Kpts = Kpts+1
                    # this will re-do the nearest search but with 1 more Kpts
                else:
                    AllNearest=True # exit 'while' block
                    # by now, idx_excl should be all non-exclusive and return data
                    idx0[0:excl_points] = idx_excl
                    dist0[0:excl_points] = dist_excl
                
            #
        #
    #
    return dist0, idx0
#
#----- projected geox/geoy ---------------------------------------------------

def nearest_pts_geoxy_kdtree(data, xname='geox', yname='geoy',kpt=2):
    import scipy.spatial as spatial
    
    data['x'] = data[xname]
    data['y'] = data[yname]
    points = list(zip(data['x'],data['y']))
    tree = spatial.KDTree(points)
    distance, index = tree.query(points, k=kpt)
    
    #return nearest points other than itself
    return distance[:, 1:], index[:,1:]

# ------
#
def pointLocator(lons_all, lats_all, lons_pt, lats_pt, is2d):
    
    xindxpts = []
    yindxpts = []
    
    # convert 'longitude' in format of 0-360, otherwise not works
    for lon in np.nditer(lons_all, op_flags=['readwrite']):
        if(lon<0.0): lon[...] = lon[...]+360.0
    
    if(lons_pt.strip()==''): # not-lon-specified point(s) 
        lons_pt=lons_all
    else:                
        if ('str' in str(type(lons_pt))):#values with operators
            nondecimal = re.compile(r'[^\d.,:;]+')
            v_pt = nondecimal.sub("",lons_pt)
            v_pt = np.float64(re.split("[,:;]",v_pt))
            
            # convert 'longitude' in format of 0-360 (sama as 'lons_all', otherwise not works
            for v in np.nditer(v_pt, op_flags=['readwrite']):
                if(v<0.0): v[...] = v[...]+360.0
            
            if("~" in lons_pt):#nearest point(s)
                pt_val=[]
                for v in v_pt:
                    pt_nearest=np.argmin(abs(lons_all-v))
                    pt_nearest=np.unravel_index(pt_nearest,lons_all.shape)
                    # if lons_all are multiple points, i.e. max. distance can be calculated,
                    # then 'pt_nearest' cannot be beyond that max. distance 
                    # (This is for 'lons_all' are from multiple files or datasets)
                    if(len(lons_all)>1):
                        max_dist = abs(np.diff(sorted(lons_all)))
                        max_dist = np.amax(max_dist)
                    else:
                        max_dist = abs(lons_all[pt_nearest]-v)
                    if(abs(lons_all[pt_nearest]-v)<=max_dist):
                        pt_val=np.append(pt_val,lons_all[pt_nearest])
            else:
                if("<=" in lons_pt): 
                    pt_val=lons_all[lons_all<=max(v_pt)]
                elif("<" in lons_pt):
                    pt_val=lons_all[lons_all<max(v_pt)]
            
                if(">=" in lons_pt): 
                    pt_val=lons_all[lons_all>=min(v_pt)]
                elif(">" in lons_pt):
                    pt_val=lons_all[lons_all>min(v_pt)]
 
            if pt_val is None:
                print('Error:  invalid ranged value expression - '+lons_pt)
                sys.exit()
            else:
                lons_pt = np.copy(pt_val)

        else: # exactly point(s)

            if(lons_all.dtype=='float64'):
                lons_pt = np.float64(lons_pt)
            else:
                lons_pt = np.float(lons_pt)
        
    #   
    if(lats_pt.strip()==''): # not-lat-specified point(s)
        lats_pt = lats_all
    else:
        if ('str' in str(type(lats_pt))):#values with operators
            nondecimal = re.compile(r'[^\d.,:;-]+')
            v_pt = nondecimal.sub("",lats_pt)
            v_pt = np.float64(v_pt.split("[,:;]"))
            if("~" in lats_pt):#nearest point(s)
                pt_val=[]
                for v in v_pt:
                    pt_nearest=np.argmin(abs(lats_all-v))
                    pt_nearest=np.unravel_index(pt_nearest,lats_all.shape)
                    # if lats_all are multiple points, i.e. max. distance can be calculated,
                    # then 'pt_nearest' cannot be beyond that max. distance 
                    # (This is for 'lats_all' are from multiple files or datasets)
                    if(len(lats_all)>1):
                        max_dist = abs(np.diff(sorted(lats_all)))
                        max_dist = np.amax(max_dist)
                    else:
                        max_dist = abs(lats_all[pt_nearest]-v)
                    if(abs(lats_all[pt_nearest]-v)<=max_dist):
                        pt_val=np.append(pt_val,lats_all[pt_nearest])
            else:
                if("<=" in lats_pt): 
                    pt_val=lats_all[lats_all<=max(v_pt)]
                elif("<" in lats_pt):
                    pt_val=lats_all[lats_all<max(v_pt)]
            
                if(">=" in lats_pt): 
                    pt_val=lats_all[lats_all>=min(v_pt)]
                elif(">" in lats_pt):
                    pt_val=lats_all[lats_all>min(v_pt)]
 
            if pt_val is None:
                print('Error:  invalid ranged value expression - '+lats_pt)
                sys.exit()
            else:
                lats_pt = np.copy(pt_val)

        else: # exactly point(s)
            if(lats_all.dtype=='float64'):
                lats_pt = np.float64(lats_pt)
            else:
                lats_pt = np.float(lats_pt)
        
                     
    if np.size(lons_pt)>0:  
        xind=np.where(np.isin(lons_all,lons_pt)) # paired indice
        xind=np.ravel_multi_index(xind,lons_all.shape) # flatten indice
    else:
        xind=[]
    if np.size(lats_pt)>0:
        yind=np.where(np.isin(lats_all,lats_pt))
        yind=np.ravel_multi_index(yind,lats_all.shape)
    else:
        yind=[]

    # must have both x and y index
    if(np.size(xind)>0 and np.size(yind)>0):
        
        # starting index
        if(is2d):
            xstart = xind
            ystart = yind
        else:
            # here we're picking points in 'xind' if it's in 'yind'. It shall be same if picking 'yind' if it's in 'xind'
            ij=np.where(np.isin(xind, yind))        
            ijindx=np.unravel_index(xind[ij], lons_all.shape)
            xstart = ijindx[0]
            ystart = ijindx[0]

   
        #unique and countinuing indices count
        x= np.array(xstart,dtype=int)
        if(np.size(x)>1):
            x=np.unique(x)
            xdiff=np.insert(np.diff(x),0,-1)
            xi=x[np.where(xdiff!=1)] # unique and non-continuing index
            xloc=np.where(np.isin(x,xi))[0] # index of xi in x, from which continuing indices can be counted
            xloc=np.append(xloc,np.size(x))
            xpts=np.diff(xloc)
        elif(np.size(x)==1):
            xi=np.copy(x)
            xpts=np.array([1],dtype=int)
        xindxpts = np.column_stack((xi,xpts))
     
        y= np.array(ystart,dtype=int)
        if(np.size(y)>1):
            y=np.unique(y)
            ydiff=np.insert(np.diff(y),0,-1)
            yi=y[np.where(ydiff!=1)] # unique and non-continuing index
            yloc=np.where(np.isin(y,yi))[0] # index of yi in y, from which continuing indices can be counted
            yloc=np.append(yloc,np.size(y))
            ypts=np.diff(yloc)
        elif(np.size(y)==1):
            yi=np.copy(y)
            ypts=np.array([1],dtype=int)
        yindxpts = np.column_stack((yi,ypts))
    
    #
    return xindxpts, yindxpts

#--------------------------------------------------------------------
"""


parser = OptionParser()

parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../../../../ccsm_inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--nco_path", dest="nco_path", default="", \
                     help = 'NCO bin PATH, default "" ')
(options, args) = parser.parse_args()


ccsm_input = os.path.abspath(options.ccsm_input)

#----
if datain2D:
    if interp_pft:
        finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='cubic')
    else:
        # nearest
        finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='linear')
        
        if UNSTRUCTURED_DOMAIN: 
            temp =[finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
            temp = np.squeeze(np.asarray(temp))
        else:
            temp = finterp_src1(lon_new, lat_new)

else:
    # nearest search, when data points are scattered (NOT in grided)
    len_new = len(lon_new)
    allpoints={}
    allpoints['Latitude']=np.append(lat_new, vlat)
    allpoints['Longitude']=np.append(lon_new, vlon)
    
    # because mixed points, must make sure nearest points are NONE of lat_new/lon_new
    AllNearest=False
    Kpts = 2
    idx_void = None
    while not AllNearest:
        dist,idx=nearest_using_kdtree(allpoints, latname='Latitude', lonname='Longitude',kpt=Kpts)
        if idx_void is None:
            idx_new = idx[0:len_new]-len_new # for vlat/vlon real index, must adjust 'idx' in mixed points 
        else:
            # fill in the negative idx only
            idx_new[idx_void] = (idx[...,0:len_new]-len_new)[idx_void]
            idx_void= np.where(idx_new<0)
            if (len(idx_void[0])>0):
                Kpts = Kpts+1
                # this will re-do the nearest search but with 1 more Kpts
            else:
                AllNearest=True # exit 'while' block
                # by now, idx_new should be all non-negative
                temp = vdata[i,idx_new.flatten()]
                temp[temp<0.0001]=0.0
                            
                vdata_new[i,] = np.copy(temp)
                del temp
        fdata_src1.close()
        #
        
        
        
#-----------------------------------------------
"""
