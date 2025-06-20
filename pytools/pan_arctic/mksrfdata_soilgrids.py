#!/usr/bin/env python
import os
import numpy as np
from pyproj import Transformer
from pyproj import CRS

import xarray as xr

#--- # downloading a boxed SoilGrid data (v2.0.1)
# 
# Poggio, L., de Sousa, L. M., Batjes, N. H., Heuvelink, G. B. M., Kempen, B., Ribeiro, E., and Rossiter, D.: 
# SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty, 
# SOIL, 7, 217–240, https://doi.org/10.5194/soil-7-217-2021, 2021.

# Following a notebook script found at: https://git.wur.nl/isric/soilgrids/soilgrids.notebooks.git
#
soilgrids_vars     = ['ocd','bdod','sand','silt','clay']
soilgrids_horizons = ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']   
soilgrids_znodes   = np.array([0.0, 0.025, 0.10, 0.225, 0.45, 0.80, 1.50, 2.0])  #mid-horizon + 2-ends (from top-bottom range of data interpolation)


def download_geotiff_soilgrids(Range_XLONG=[], Range_YLATI=[], outputpath='./soilgrids',
                                soilvars=soilgrids_vars,
                                horizons=soilgrids_horizons, 
                                value='mean'):
    
    from owslib.wcs import WebCoverageService
    import rasterio
    import pathlib
    
    # SoilGrids map's CRS
    crs_wcs = "http://www.opengis.net/def/crs/EPSG/0/152160"
    # from above def, its equvalent is:
    crs_proj4 = CRS.from_proj4('+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')
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
    tif_withcrs = str(pathlib.Path(__file__).parent.resolve())+'/soilgrids_template_withcrs.tif'
    template_profile = rasterio.open(tif_withcrs).profile.copy()
    
    os.system('mkdir -p '+outputpath)
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
    
    #
    return outputpath
    

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

def srf_soils_from_soilgrid_tiff(soilvars=soilgrids_vars, horizons=soilgrids_horizons, 
                                soilgrids_dir='./', reproj_crs="EPSG:4326", 
                                elm_domainnc='', nlevsoi=10, outdata=True):

    import rioxarray
    import pytools.commons_utils.gridlocator as gridlocator
    
    #
    srf_soils = {}
    for ivar in soilvars:

        if ivar=='ocd': 
            elmvar = 'ORGANIC'
        elif ivar=='sand':
            elmvar = 'PCT_SAND'
        elif ivar=='clay':
            elmvar = 'PCT_CLAY'
        elif ivar=='silt':
            elmvar = 'PCT_SILT'
        else:
            print('NO corresponded ELM surface soil data: ', ivar)
            continue
        

        temp_data = []
        jx=[]
        for iz in range(len(horizons)):
            zstr=horizons[iz]
            
            # vertical nodes
            j = soilgrids_horizons.index(zstr)+1
            if len(jx)==0:
                jx=[j-1,j]          # top node inserted
            elif iz==len(horizons)-1:
                jx.extend([j,j+1])  # bottom node appened
            else:
                jx.extend([j])
            
            svar_horizon_id = ivar+'_'+zstr+'_mean'
            svar_horizon_value_tif = soilgrids_dir+'/'+svar_horizon_id+'.tif'
            with rioxarray.open_rasterio(svar_horizon_value_tif) as src:
                if reproj_crs=='':
                    svar = src[0]
                else:
                    svar = src[0].rio.reproject(reproj_crs)
 
                #matching with ELM domain
                if elm_domainnc!='':
                    #
                    df = gridlocator.grids_data_from_xrio(svar, elm_domainnc=elm_domainnc)
                    svar = df['g_mean']
                    if not 'LONGXY' in srf_soils.keys(): srf_soils['LONGXY'] = df.centroid.x.to_numpy()
                    if not 'LATIXY' in srf_soils.keys(): srf_soils['LATIXY'] = df.centroid.x.to_numpy()
                #
                else:
                    if not 'LONGXY' in srf_soils.keys(): srf_soils['LONGXY'] = svar.x.data.to_numpy()
                    if not 'LATIXY' in srf_soils.keys(): srf_soils['LATIXY'] = svar.y.data.to_numpy()
                temp_data.append(svar.to_numpy())
            #
        #
        temp_data = np.asarray(temp_data)
        temp_data = temp_data.reshape(temp_data.shape[0],-1)

        #mid-horizon + 2-ends (from top-bottom range of data interpolation)
        znodes = soilgrids_znodes[jx]
        ngrid = temp_data[0].size
            
        # ELM soil column data has 10 layers down to about 4.2 m
        #nlevsoi = 10 by default
        
        srf_soils[elmvar] = np.empty((nlevsoi, ngrid))
        for ig in range(ngrid):
            zdata = np.hstack([temp_data[0,ig],
                                    temp_data[:,ig],
                                    temp_data[-1,ig]])
            if ivar in ['bdod','sand','silt','clay']:
                srf_soils[elmvar][:,ig] = \
                    mksrfdata_soilcolumn_interp( \
                        srf_soildata=zdata, \
                        srf_soilnode=znodes, \
                        nlevsoi=nlevsoi, fill_method="extrapolate")
                if ivar in ['sand','silt','clay']: 
                    srf_soils[elmvar][:,ig] = srf_soils[elmvar][:,ig]*0.1 #but not sure if is in % or g/kg ??? it's very confusing of mapping unit or downloaded data unit
            elif ivar in ['ocd']:
                srf_soils[elmvar][:,ig] = \
                    mksrfdata_soilcolumn_interp( \
                        srf_soildata=zdata, \
                        srf_soilnode=znodes, \
                        nlevsoi=nlevsoi)
                srf_soils[elmvar][:,ig] = srf_soils[elmvar][:,ig]*0.1 #but not sure if 'ocd' in kgC or kgSOM ???
        
            
        #
    #
    if outdata:
        return srf_soils 


#--------------------------------------------------------------------
def test(surf_from_soilgrid={}, surf_vars='PCT_SAND,PCT_CLAY,ORGANIC'):
    
    import pytools.pan_arctic.mksrfdata_updatevals as srfupdate
    
    # boxed boundary
    with xr.open_dataset('./domain.lnd.r0125_IcoswISC30E3r5.250918_cavm1d.nc') as domain_xr:
        xlmt = domain_xr['xv'].to_numpy()
        ylmt = domain_xr['yv'].to_numpy()
        xrange = [np.nanmin(xlmt), np.nanmax(xlmt)]
        yrange = [np.nanmin(ylmt), np.nanmax(ylmt)]
        
    # SoilGrids download into ./soilgrids/*.tif
    soilgrids_dir=download_geotiff_soilgrids(Range_XLONG=xrange, Range_YLATI=yrange)  # if already downloaded, can skip and just use it
    #soilgrids_dir='./soilgrids' # if already downloaded
    
    # Generated soil data consistent with ELM surfdata structure/names from SoilGrids' tiff
    surf_from_soilgrid = srf_soils_from_soilgrid_tiff(
                                soilvars=soilgrids_vars, 
                                horizons=soilgrids_horizons, 
                                soilgrids_dir=soilgrids_dir, 
                                elm_domainnc='domain.lnd.r0125_IcoswISC30E3r5.250918_cavm1d.nc', 
                                nlevsoi=10, outdata=True)
    
    
    
    # put data into surfdata.nc
    output_path = \
        os.path.join('./', \
                        'surfdata_0.125x0.125_simyr1850_c250910_TOP_cavm1d.nc')
    
    srfupdate.updatevals(output_path, \
                         user_srf_data=surf_from_soilgrid, \
                         user_srf_vars=surf_vars)
    
      
if __name__ == '__main__':
    test()




