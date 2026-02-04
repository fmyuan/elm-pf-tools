#!/usr/bin/env python
import os
import numpy as np
from pyproj import Transformer
from pyproj import CRS

import xarray as xr
import rasterio

#--- # downloading a boxed SoilGrid data (v2.0.1)
# 
# Poggio, L., de Sousa, L. M., Batjes, N. H., Heuvelink, G. B. M., Kempen, B., Ribeiro, E., and Rossiter, D.: 
# SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty, 
# SOIL, 7, 217–240, https://doi.org/10.5194/soil-7-217-2021, 2021.

# Following a notebook script found at: https://git.wur.nl/isric/soilgrids/soilgrids.notebooks.git
#
def download_geotiff_soilgrids(Range_XLONG=[], Range_YLATI=[], outputpath='./original', \
                                soilvars=['ocd','bdod','sand','silt','clay'], \
                                value='mean'):
    
    from owslib.wcs import WebCoverageService
    
    # SoilGrids map's CRS
    crs_wcs = "http://www.opengis.net/def/crs/EPSG/0/152160"
    # from above def, its equvalent is:
    crs_proj4 = CRS.from_proj4('+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')
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
    Tlonlat2xy = Transformer.from_proj(lonlatProj, crs_proj4, always_xy=True)

    X1,Y1 = Tlonlat2xy.transform(Range_XLONG[0],Range_YLATI[0]) #left-bottom
    X2,Y2 = Tlonlat2xy.transform(Range_XLONG[0],Range_YLATI[1]) #right-bottom
    X3,Y3 = Tlonlat2xy.transform(Range_XLONG[1],Range_YLATI[1]) #right-top
    X4,Y4 = Tlonlat2xy.transform(Range_XLONG[1],Range_YLATI[0]) #left-top

    # when projection-transformed, not rectangle anymore, so need to re-do min/max (otherwise, coverage may be incompleted)
    subsets = [('X', min(X1,X2,X3,X4), max(X1,X2,X3,X4)), ('Y', min(Y1,Y2,Y3,Y4), max(Y1,Y2,Y3,Y4))]
    
    # organic carbon density: ocd, hg/m3 (aka 0.1kg/m3)
    
    # obtain a full geotiff profile for tiff writing
    template_profile = rasterio.open(outputpath+'/soilgrids_template_withcrs.tif').profile.copy()
    
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
            with open('./tmp.tif', 'wb') as file:
                # better to save to xarray, but didn't figure out how to from response.read() -> xarray
                file.write(response.read())
            
            # ideally, the above file in tiff should have 'CRS' info but not
            # so we need to redo
            rdata=rasterio.open('./tmp.tif')
            newprofile = rdata.profile.copy()
            newprofile['crs'] = template_profile['crs']            
            svar_horizon_value_tif = outputpath+'/'+svar_horizon_id+'.tif'
            with rasterio.open(
                svar_horizon_value_tif,
                'w',
                **newprofile,
            ) as file:
                file.write(rdata.read(1), 1)
            
            if os.path.exists('./tmp.tif'): os.remove('./tmp.tif')
            
        # layer loop
    #variable loop
    

#--- #
# interplating layered soil data to match with ELM soil vertical structure
def mksrfdata_soilcolumn_interp(srf_soildata=np.empty((0)), srf_soilnode=np.empty((0)), \
                                nlevsoi=10, fill_method='zero'):
    
    from pytools.pan_arctic.arctic_domain import soilcolumn
    
    if srf_soildata.size!= srf_soilnode.size \
        or srf_soildata.size<=0: return np.empty((0))
    
    # default ELM soil column
    zisoi, dzsoi, zsoi = soilcolumn()
    # by default, ELM soil layer no  is 10, the rest 5 layers are as rock or alike 
    zisoi = zisoi[0:nlevsoi+1]  
    dzsoi = dzsoi[1:nlevsoi+1]
    zsoi = zsoi[1:nlevsoi+1]  # so-called node-depth, within a soil layer (but not middle)    
    znodes = xr.Dataset({'z': ("points", zsoi)})
    
    if srf_soilnode.size<=1:
        # but not sure if 2 data points works with 'data.interp()' function
        interp_vert = np.zeros(len(zsoi))
        if fill_method == 'extrapolate':
            interp_vert[:] = srf_soildata
        else:
            for iz in range(len(zsoi)):
                if zsoi[iz]<= srf_soilnode:
                    interp_vert[iz] = srf_soildata
                else:
                    continue
                           
        return interp_vert
    
    else:
        data = xr.DataArray(srf_soildata, 
                        dims=("z"), 
                        coords={"z": srf_soilnode})
        
        if fill_method == 'zero':
            interp_vert = data.interp(znodes, method='linear', \
                                  kwargs={'fill_value': 0})
        else:
            interp_vert = data.interp(znodes, method='linear', \
                                  kwargs={'fill_value': fill_method})        
        return interp_vert.as_numpy()
            
#--- # 
# 

#--------------------------------------------------------------------
def test(surf_from_soilgrid={}, surf_vars='PCT_SAND, PCT_CLAY, ORGANIC'):
    
    import pytools.pan_arctic.mksrfdata_updatevals as srfupdate
    
    
    input_path  = './'
    
    output_path = \
        os.path.join(('./', \
                        'surfdata_soilgrid_test.nc'))
    
    srfupdate.updatevals(output_path, \
                         user_srf_data=surf_from_soilgrid, \
                         user_srf_vars=surf_vars)
    
      
#if __name__ == '__main__':
#    test()




