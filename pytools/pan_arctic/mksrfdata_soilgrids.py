#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset
from pyproj import Transformer
from pyproj import CRS

import geopandas as gpd
import rioxarray

try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False

#---
def download_geotiff_soilgrids(Range_XLONG=[], Range_YLATI=[], outputpath='./original', \
                                soilvars=['ocd','bdod','sand','silt','clay'], \
                                value='mean'):
    
    from owslib.wcs import WebCoverageService
    
    # SoilGrids map's CRS
    crs_wcs = "http://www.opengis.net/def/crs/EPSG/0/152160"
    crs_proj = CRS.from_proj4('+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')
    # standard horizons
    horizons = ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']   
    # mean
    #value='mean'
     
    # coverage of rectangle-boxed area
    #Range_XLONG = [-150.0, -149.0]
    #Range_YLATI =[68.5, 69.0]
    
    # EPSG: 4326
    # Proj4: +proj=longlat +datum=WGS84 +no_defs
    lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
    Tlonlat2xy = Transformer.from_proj(lonlatProj, crs_proj, always_xy=True)

    X1,Y1 = Tlonlat2xy.transform(Range_XLONG[0],Range_YLATI[0]) #left-bottom
    X2,Y2 = Tlonlat2xy.transform(Range_XLONG[0],Range_YLATI[1]) #right-bottom
    X3,Y3 = Tlonlat2xy.transform(Range_XLONG[1],Range_YLATI[1]) #right-top
    X4,Y4 = Tlonlat2xy.transform(Range_XLONG[1],Range_YLATI[0]) #left-top

    # when projection-transformed, not rectangle anymore, so need to re-do min/max (otherwise, coverage may be incompleted)
    subsets = [('X', min(X1,X2,X3,X4), max(X1,X2,X3,X4)), ('Y', min(Y1,Y2,Y3,Y4), max(Y1,Y2,Y3,Y4))]
    
    # organic carbon density: ocd, hg/m3 (aka 0.1kg/m3)
        
    for ivar in soilvars:
        soilgrids_wcs = WebCoverageService('http://maps.isric.org/mapserv?map=/map/'+ivar+'.map',
                         version='2.0.1')
        #infos for checking
        

        for iz in horizons:
            svar_horizon_id = ivar+'_'+iz+'_'+value
            svar_horizon = soilgrids_wcs.contents[svar_horizon_id]
            
            response = soilgrids_wcs.getCoverage(identifier=[svar_horizon_id], 
                                   crs=crs_wcs,
                                   subsets=subsets, 
                                   resx=250, resy=250, 
                                   format=svar_horizon.supportedFormats[0])  

            svar_horizon_value_tif = outputpath+'/'+svar_horizon_id+'.tif'
            with open(svar_horizon_value_tif, 'wb') as file:
                file.write(response.read())
    

#---
def surfnc_soils_from_geotiffsoilgrids(surfnc_soils_template='./', \
                                    inputpath='./original', \
                                    outputpath='./'):
    
    from shapely.geometry import MultiLineString
    from shapely.ops import polygonize
    #from osgeo import gdal

    horizons = ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']        
    soilvars=['ocd','bdod','sand','silt','clay']

    # the following assumes all data are with same extent and resolution    
    
    alldata = {}
    for ivar in range(len(soilvars)):
        var = soilvars[ivar]
        alldata[var] = {}
        for iz in range(len(horizons)):
            zstr = horizons[iz]
            file  = inputpath+'/'+var+'_'+zstr+'_mean.tif'
            
            fileout = outputpath+'/'+var+'_'+zstr+'_mean-latlon.tif'
            
            '''
            # the following didn't work, always got zeros for band data. If this works, will save a lot effort.
            gdal.Warp(fileout, file, \
                            format='GTiff', srcSRS='ESRI:54052', \
                            dstNodata=None, resampleAlg='average', \
                            xRes=0.0025, yRes=0.0025, \
                            outputBounds=[-150.0,68.5,-149.0,69.0], \
                            dstSRS='EPSG:4326')
            '''           
            rxdata = rioxarray.open_rasterio(file)
            # geometry only needs to do once
            if ivar==0 and iz==0:
            
                x = rxdata.x.data
                y = rxdata.y.data
                resx = np.diff(x)
                resy = np.diff(y)
                xv = np.linspace(x[0]-resx[0]/2.0, x[-1]+resx[-1]/2.0, len(x)+1)         
                yv = np.linspace(y[0]-resy[0]/2.0, y[-1]+resy[-1]/2.0, len(y)+1)    
                hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(xv[:-1], xv[1:]) for yi in yv]
                vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(yv[:-1], yv[1:]) for xi in xv]
                lines_str = vlines+hlines
                gridcells = list(polygonize(MultiLineString(lines_str)))
            
                layer0 = rxdata[0].data.flatten()
                gdf0 = gpd.GeoDataFrame(layer0,geometry=gridcells,crs=rxdata.rio.crs)
                gdf1 = gpd.GeoDataFrame.to_crs(gdf0, crs='EPSG:4326')
                
    # clip data into new lat/lon grid-system
    lon=gdf1.total_bounds()
    lat=gdf1.y
    res_lon = 0.01

    '''
    alldata['x'] = xx
    alldata['y'] = yy
    alldata['z'] = [0.0, 0.05, 0.15, 0.30, 0.60, 1.00, 2.00]
    XC, YC = np.meshgrid(xx, yy)
    lon2d,lat2d = Txy2lonlat.transform(XC,YC)
    alldata['lon2d']= lon2d
    alldata['lat2d']= lat2d
    alldata['crs_res'] = crs_res
    alldata['crs_wkt'] = crs.wkt
    '''
    return alldata
    

# 


#---
def mksrfdata_updatesoils(fsurfnc_in, user_srf_data={}, user_srfnc_file='', user_srf_vars='', OriginPFTclass=True):
    
    print('#--------------------------------------------------#')
    print("Replacing values in surface data by merging user-provided dataset")
    fsurfnc_out ='./'+fsurfnc_in.split('/')[-1]+'-merged'
    
    
    #--------------------------------------
    UNSTRUCTURED = False    
    if not user_srfnc_file !=' ':
        print('read data from: ', user_srfnc_file)
        f=Dataset(user_srfnc_file)
        if 'gridcell' in f.dimensions.items(): UNSTRUCTURED = True
    elif not len(user_srf_data)<=0:
        if len(user_srf_data['LATIXY'].shape)==1:
            UNSTRUCTURED = True
            
        if user_srf_vars=='':
            user_vname = user_srf_data.keys()
        else:
            user_vname = user_srf_vars.split(',')
        user_srf = user_srf_data

'''    
    if 'PCT_NAT_PFT' in user_srf.keys():
        if OriginPFTclass:
            # need to aggregate PCT_NAT_PFT from new one to default ELM categories
            user_natpft = user_srf['PCT_NAT_PFT']
            pct_natpft = np.zeros((len(natpft),)+user_natpft[0,...].shape)
            for i in range(len(user_pfts['pftnum'])):
                ip = user_pfts['pftnum'][i]
                pct_natpft[ip,...] = pct_natpft[ip,...]+user_natpft[i,...]
            #
            user_srf['PCT_NAT_PFT'] = pct_natpft
        
        # tiny fraction, which may be caught in ELM and crash model
        user_srf['PCT_NAT_PFT'][np.where(user_srf['PCT_NAT_PFT']<1.e-8)] = 0.0
        sumpft = np.sum(user_srf['PCT_NAT_PFT'], axis=0)
        for i in range(len(user_pfts['pftnum'])):
            ip = user_pfts['pftnum'][i]
            user_srf['PCT_NAT_PFT'][ip,...] = user_srf['PCT_NAT_PFT'][ip,...]/sumpft*100.0
        
        
    #--------------------------------------
    #                    
    # write into nc file
    with Dataset(fsurfnc_in,'r') as src, Dataset(fsurfnc_out, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # explicitly add a global attr for visualizing in gis tools
        dst.Conventions = "CF-1.0"
            
        # new surfdata dimensions
        for dname, dimension in src.dimensions.items():
            if dname == 'natpft':
                len_dimension = len(natpft)            # dim length from new data
            elif dname == 'gridcell':
                len_dimension = user_srf['LATIXY'].flatten().size
            elif dname in ['lsmlat','lat']:
                len_dimension = user_srf['LATIXY'].shape()[0]
            elif dname in ['lsmlon', 'lon']:
                len_dimension = user_srf['LONGXY'].shape()[1]                
            else:
                len_dimension = len(dimension)
            dst.createDimension(dname, len_dimension if not dimension.isunlimited() else None)
            #
            
        # create variables and write to dst
        for vname, variable in src.variables.items():

            if UNSTRUCTURED:
                vdim = variable.dimensions
                if 'gridcell' not in vdim and \
                    ('lsmlat' in vdim and 'lsmlon' in vdim):
                    vdim = vdim.replace('lsmlon', 'gridcell')
                    vdim = vdim.remove('lsmlat')
            else:
                vdim = variable.dimensions
    
            # create variables, but will update its values later 
            # NOTE: here the variable value size updated due to newly-created dimensions above
            dst.createVariable(vname, variable.datatype, vdim)
            # copy variable attributes all at once via dictionary after created
            dst[vname].setncatts(src[vname].__dict__)
                  
            # values
            print(vname, vdim)

            # dimension length may change, so need to 
            if vname=='natpft':
                varvals = dst[vname][...]
                varvals[...] = natpft
            
            elif vname in user_vname:
                varvals = dst[vname][...]
                varvals[...] = user_srf[vname]
                #
            else:
                varvals = src[vname][...]
                #
            
            #                                
            dst[vname][...] = varvals
                    
        # end of variable-loop        
                
        
        print('user surfdata merged and nc file created and written successfully!')
        
    #
#

'''
        
#--------------------------------------------------------------------
def test():

    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """  
    input_path  = './'
    output_path = './'

    # download geotiff soils from SoilGrids
    #download_geotiff_soilgrids(Range_XLONG=[-150.0, -149.0], Range_YLATI=[68.5, 69.0])
    
    # do necessay data processing for soils, and produce a ready-to-merge soil-only surfdata nc file
    surfnc_soils_from_geotiffsoilgrids(surfnc_soils_template='../../../surfdata_0.5x0.5_simyr1850_c240308_TOP.nc')
    
    #mksrfdata_updatesoils('./surfdata_0.0025deg.1D_simyr1850_c240308_TOP_TFSarcticpfts-default.nc', \
    #                    #user_srf_data=surf_fromcavm_jk, user_srf_vars='', OriginPFTclass=True)
    #                    user_srf_data=surf_fromcavm_jk, user_srf_vars='', OriginPFTclass=False)
    
      
if __name__ == '__main__':
    test()