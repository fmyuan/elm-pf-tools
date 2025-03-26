#!/usr/bin/env python
import os, math
import numpy as np

#----- 
#----- nearest_neibour for earth surface using kdtree
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

def grids_nearest_or_within(src_grids={}, masked_pts={}, remask_approach = 'nearest', \
                     unlimit_xmin=False, unlimit_xmax=False, unlimit_ymin=False, unlimit_ymax=False,\
                     out2d=False, reorder_src=False, keep_duplicated=False):
    '''
    re-masking ELM domain grids, by providing centroids which also may be masked
        src_grids      - list of source grids, with required paired of 'xc' and 'yc', 
                         and/or its vertices 'xv', 'yv' for within-grid search
        masked_pts     - list of np arrays of 'xc', 'yc', 'mask' of user-provided grid-centroids or points
        remask_approach- 'nearest', or points in polygon if 'xv'/'yv' known
        unlimit_x(y)min(max) - x/y axis bounds are from source domain, when True; otherwise, from masked_pts
        out2d          - output data in 2D or flatten (so-called unstructured)
        reorder_src    - trunked source grids will re-order following that in masked_pts 
        keep_duplicated- don't sort or unique source grids. 
    '''

    import geopandas as geopd
    from shapely.geometry import Point
    from geopandas.tools import sjoin, sjoin_nearest
    from shapely.geometry import MultiLineString, LineString
    from shapely.ops import polygonize

    if len(masked_pts['xc'])<1 or len(masked_pts['yc'])<1: 
        print('Error: no masked paired xc/yc centroids provided! ')
        return    
    if len(src_grids['xc'])<1 or len(src_grids['yc'])<1: 
        print('Error: no source grids paired xc/yc centroids provided! ')
        return    

    # masked points, paired xc/yc (NOT blocked or box of [[xc],[yc]])
    mask_xc= masked_pts['xc']
    mask_yc= masked_pts['yc']
    LONGXY360=True
    if np.any(mask_xc<0.0): LONGXY360=False
    mask_id = np.indices(mask_xc.flatten().shape)[0]
        
    # 
    # Open the source domain or surface nc file, 
    # directly from an ELM domain nc file
    # TIP: in this case, both domain and land surface grid-systems are EXACTLY matching with each other
    src_yc = src_grids['yc']
    src_xc = src_grids['xc']
    if LONGXY360: 
        src_xc[src_xc<0.0]=360+src_xc[src_xc<0.0]
    else:
        src_xc[src_xc>180.0]=-360+src_xc[src_xc>180.0]
    if 'xv' in src_grids.keys() \
        and 'yv' in src_grids.keys():
        src_xv = src_grids['xv']
        if LONGXY360: 
            src_xv[src_xv<0.0]=360+src_xv[src_xv<0.0]
        else:
            src_xv[src_xv>180.0]=-360+src_xv[src_xv>180.0]
            # need to re-adj for cross -180/180 line's grid
            
        src_yv = src_grids['yv']    
    
    # 
    # 1.5 times of grid-size edges 
    # to allows one-grid edges of range of latitude/longitude boxes
    x_edge = 0.0; y_edge = 0.0
    if (1 in src_xc.shape):
        #probably unstructed grids
        if 'xv' in src_grids.keys() and \
            'yv' in src_grids.keys():
            xdiff = np.squeeze(abs(src_xv[...,0]-src_xv[...,1]))
            xdiff = xdiff[np.where(xdiff<=180.0)] # removal those cross line 0/360 or -180/180
            x_edge = np.nanmean(xdiff)*1.5
            y_edge = np.nanmean(np.squeeze(abs(src_yv[...,0]-src_yv[...,2])))*1.5        
    else:
        xdiff = np.diff(src_xc[0,:])
        xdiff = xdiff[np.where(abs(xdiff)<=180.0)] # removal those cross line 0/360 or -180/180
        x_edge = np.nanmean(xdiff)*1.5
        y_edge = np.nanmean(np.diff(src_yc[:,0]))*1.5
    # [mask_xc,mask_yc] box
    x_min = min(mask_xc)-x_edge
    if unlimit_xmin: x_min = min(src_xc[0,:])-x_edge # but may not want to do trunking
    x_max = max(mask_xc)+x_edge
    if unlimit_xmax: x_max = max(src_xc[0,:])+x_edge
    y_min = min(mask_yc)-y_edge
    if unlimit_ymin: y_min = min(src_yc[:,0])-y_edge
    y_max = max(mask_yc)+y_edge
    if unlimit_ymax: y_max = max(src_yc[:,0])+y_edge
    
    boxed_idx = \
        np.where((src_yc>=y_min) &  \
                 (src_yc<=y_max) &  \
                 (src_xc>=x_min) &  \
                 (src_xc<=x_max))
    xc = src_xc[boxed_idx]
    yc = src_yc[boxed_idx]
    poly_id = np.indices(src_xc.flatten().shape)[0]
    poly_id = poly_id.reshape(src_xc.shape)[boxed_idx]
    # this is the indices of trunked but in original paired [yc,xc] of source  

    if 'xv' in src_grids.keys() and \
        'yv' in src_grids.keys():    
        yv1 = src_yv[...,0]; xv1 = src_xv[...,0]
        yv2 = src_yv[...,1]; xv2 = src_xv[...,1]
        yv3 = src_yv[...,2]; xv3 = src_xv[...,2]
        yv4 = src_yv[...,3]; xv4 = src_xv[...,3]        
        yv1=yv1[boxed_idx]; xv1=xv1[boxed_idx]
        yv2=yv2[boxed_idx]; xv2=xv2[boxed_idx]
        yv3=yv3[boxed_idx]; xv3=xv3[boxed_idx]
        yv4=yv4[boxed_idx]; xv4=xv4[boxed_idx]
    
    # trunked source domain's xc,yc, i.e. centroids of xv/yv which formed polygons above, 
    # to be marked in X_axis, Y_axis, and re-indexing
    # (TIPs: this is import for visualizing data in map or transfrom btw 1D and 2D)
    X_axis = xc
    Y_axis = yc     
    [X_axis, i] = np.unique(X_axis, return_inverse=True)
    [Y_axis, j] = np.unique(Y_axis, return_inverse=True)
    YY, XX = np.meshgrid(Y_axis, X_axis, indexing='ij')
    xid = np.indices(XX.shape)[1]
    yid = np.indices(XX.shape)[0]
    xyid = np.indices(XX.flatten().shape)[0]
    xyid = np.reshape(xyid, xid.shape)

    if (1 not in src_xc.shape):
        # structured grids yc[nj,ni], xc[nj,ni]
        
        # new indices of [yc,xc] in grid-mesh YY/XX
        # note: here don't override xc,yc, so xc=X_axis[i], yc=Y_axis[j]
        polygonized = np.isin(XX, X_axis[i]) & np.isin(YY, Y_axis[j])
        ji = np.nonzero(polygonized) 
        sub_xid = xid[ji]     
        sub_yid = yid[ji]    
        sub_xyid= xyid[ji]  

    else:
        # unstructured  yc[1,ni], xc[1,ni]

        # new indices of paired [yc,xc] in grid-mesh YY/XX
        sub_xid = xid[j,i]     
        sub_yid = yid[j,i] 
        sub_xyid= xyid[j,i]        
    #

    # masked points
    points = geopd.GeoDataFrame({"pts_indices":mask_id, \
                                 "x":mask_xc,"y":mask_yc})
    points['geometry'] = points.apply(lambda p: Point(p.x, p.y), axis=1)

    #
    #remask_approach = 'nearest'
    if remask_approach=='nearest':
        print('nearest point searching')

        # re-meshed XX/YY centroids
        grid_centroids = geopd.GeoDataFrame({'pid':poly_id,"xyid":sub_xyid,"xid":sub_xid,"yid":sub_yid, \
                                             "xc":xc,"yc":yc})
        grid_centroids['geometry'] = grid_centroids.apply(lambda p: Point(p.xc, p.yc), axis=1)

        pointsInGrids = sjoin_nearest(points, \
                                      grid_centroids[['pid','xyid','xid','yid','geometry']], how='inner')
        
        remask_pid = pointsInGrids.pid  # the original (untrunked) grid/polygon flatten id
        remask_xyid = pointsInGrids.xyid  # the grid/polygon flatten id in new XX/YY mesh (re-ordered)
        remask_maskid = pointsInGrids.pts_indices # the orginal (masking) point id

    
    elif 'xv' in src_grids.keys() \
        and 'yv' in src_grids.keys():
        print('points in polygon searching')
        
        # box trunked source domain's xv,yv to form polygons
        vpts1=[(x,y) for x,y in zip(xv1,yv1)]
        vpts2=[(x,y) for x,y in zip(xv2,yv2)]
        vpts3=[(x,y) for x,y in zip(xv3,yv3)]
        vpts4=[(x,y) for x,y in zip(xv4,yv4)]
        lines = [(vpts1[i], vpts2[i], vpts3[i], vpts4[i], vpts1[i]) for i in range(len(vpts1))]
        lines_str = [LineString(lines[i]) for i in range(len(lines))]
        polygons = list(polygonize(MultiLineString(lines_str)))
        # if regular x/y mesh may be do like following: 
        #x = np.linspace(src_xv[0,0], xv[0,-1], len(src_xc[1]))         
        #y = np.linspace(src_yv[0,0], yv[0,-1], len(src_yc[0]))    
        #hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(x[:-1], x[1:]) for yi in y]
        #vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(y[:-1], y[1:]) for xi in x]
        #lines_str = vlines+hlines
        #polygons = list(polygonize(MultiLineString(lines_str)))
        
        grids = geopd.GeoDataFrame({'pid':poly_id,"xyid":sub_xyid,"xid":sub_xid,"yid":sub_yid,"geometry":polygons})
    
        pointsInPolys = sjoin(points, grids[['pid','xyid','xid','yid','geometry']], how='inner')
        remask_pid = pointsInPolys.pid.values  # the original (untrunked) grid/polygon flatten id
        remask_xyid = pointsInPolys.xyid.values  # the grid/polygon flatten id in new XX/YY mesh (re-ordered)        
        remask_maskid = pointsInPolys.pts_indices.values # the orginal (masking) point id
        
    #
        
    # need to obtain all original masked_pts id for each srcpts_unique_pid
    # so that, could do some maths within a srcpts_unique_pid
    srcpolys_unique_pid, pts_startidx, pts_count  = np.unique(remask_pid, \
                            return_index=True, return_counts=True)
    sorted_ptsingrids = np.column_stack((remask_pid, remask_maskid))
    sorted_ptsingrids = sorted_ptsingrids[sorted_ptsingrids[:,0].argsort()]      
    srcpolys_contained_maskid = {}
    for i in range(len(srcpolys_unique_pid)):
        # sorted unique grid id as key, for which masked_pid contained can be saved
        srcpolys_contained_maskid[srcpolys_unique_pid[i]] = \
            sorted_ptsingrids[pts_startidx[i]:pts_startidx[i]+pts_count[i],1]

    # re-masking grids, which actually conaining masked_xc/_yc/_id in new XX/YY mesh
    sub_remasked = np.isin(sub_xyid, remask_xyid)

    if keep_duplicated or reorder_src:
        
        sorted_ptsidx = remask_xyid # not really sorted here
        if keep_duplicated:
            xc = xc[sorted_ptsidx]
            yc = xc[sorted_ptsidx]
        elif reorder_src:
            # shuffle indices or data array by data order in mask
            sorted_ptsingrids = np.column_stack((remask_maskid, remask_pid, remask_xyid))
            sorted_ptsingrids = sorted_ptsingrids[sorted_ptsingrids[:,0].argsort()]
            sorted_ptsidx = sorted_ptsingrids[:,2] # i.e. sorted 'remask_xyid' by 'remask_maskid'

            # using mask_pts x/y, but still with src_xv/yv below
            xc = mask_xc[remask_maskid]
            yc = mask_yc[remask_maskid]
                    
        if 'xv' in src_grids.keys() and \
            'yv' in src_grids.keys():    
            yv1=yv1[sorted_ptsidx]; xv1=xv1[sorted_ptsidx]
            yv2=yv2[sorted_ptsidx]; xv2=xv2[sorted_ptsidx]
            yv3=yv3[sorted_ptsidx]; xv3=xv3[sorted_ptsidx]
            yv4=yv4[sorted_ptsidx]; xv4=xv4[sorted_ptsidx]
        
        # re-do indices of [yc,xc] in [YY,XX] meshgrid
        boxed_idx = (boxed_idx[0][sorted_ptsidx], boxed_idx[1][sorted_ptsidx])
        sub_xyid = sub_xyid[sorted_ptsidx]
        sub_xid= sub_xid[sorted_ptsidx]
        sub_yid= sub_yid[sorted_ptsidx]
        sub_remasked = sub_remasked[sorted_ptsidx]
        

    if out2d and (keep_duplicated or reorder_src):
        # to 2D but masked    
        grid_id_arr = np.reshape(sub_xyid,xyid.shape)
        grid_xid_arr = np.reshape(sub_xid,xyid.shape)
        grid_yid_arr = np.reshape(sub_yid,xyid.shape)
        grid_remasked_arr = np.reshape(sub_remasked,xyid.shape)
        
        xc_arr = np.reshape(xc,xyid.shape)
        yc_arr = np.reshape(yc,xyid.shape)
        xv_arr = np.stack((np.reshape(xv1,xyid.shape), \
                           np.reshape(xv2,xyid.shape), \
                           np.reshape(xv3,xyid.shape), \
                           np.reshape(xv4,xyid.shape)), axis=2)
        yv_arr = np.stack((np.reshape(yv1,xyid.shape), \
                           np.reshape(yv2,xyid.shape), \
                           np.reshape(yv3,xyid.shape), \
                           np.reshape(yv4,xyid.shape)), axis=2)
    else:
    
        grid_id_arr = sub_xyid
        grid_xid_arr = sub_xid
        grid_yid_arr = sub_yid
        grid_remasked_arr = sub_remasked
        xc_arr = xc
        yc_arr = yc
        xv_arr = np.stack((xv1,xv2,xv3,xv4), axis=1)
        yv_arr = np.stack((yv1,yv2,yv3,yv4), axis=1)
 
    # 
    # output arrays
    subdomain = {}
    subdomain['X_axis'] = X_axis
    subdomain['Y_axis'] = Y_axis
    subdomain['XX'] = XX
    subdomain['YY'] = YY
            
    subdomain['gridID']  = grid_id_arr
    subdomain['gridXID'] = grid_xid_arr
    subdomain['gridYID'] = grid_yid_arr
    subdomain['xc'] = xc_arr
    subdomain['yc'] = yc_arr
    subdomain['xv'] = xv_arr
    subdomain['yv'] = yv_arr
    subdomain['remasked'] = grid_remasked_arr
        
    return subdomain, boxed_idx, srcpolys_unique_pid, srcpolys_contained_maskid
#
