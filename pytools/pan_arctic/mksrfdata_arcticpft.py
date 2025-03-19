#!/usr/bin/env python
import os
import numpy as np
from netCDF4 import Dataset
from cmath import nan, inf
from pyproj import Transformer
from pyproj import CRS

from pytools.commons_utils.gridlocator import grids_nearest_or_within

try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False

#---
def arctic_veged_fromraster_jkumaretal(inputpath='./', \
                                    rasterfiles=['./graminoid_toolik.tif'], 
                                    normallized=False):
    #
    from pytools.commons_utils.Modules_netcdf import geotiff2nc
    
    # the following data are not normalized over land area, and only really-vegetated
    rasterfiles = ['lichen_toolik.tif', 
                   'bryophyte_toolik.tif', 
                   'forb_toolik.tif', 
                   'graminoid_toolik.tif', 
                   'evergreen_shrub_toolik.tif', 
                   'deciduous_shrub_toolik.tif']

    # the following are, for data above (exactly matching up with), standared PFT names and order in 14 arcticpft physiology parameters
    # a full list is: (0)not_vegetated, (1)arctic_lichen, (2)arctic_bryophyte,
    #                 (3)arctic_needleleaf_tree, (4)arctic_broadleaf_tree,
    #                 (5)arctic_evergreen_shrub_dwarf, (6)arctic_evergreen_shrub_tall,
    #                 (7)arctic_deciduous_shrub_dwarf, (8)arctic_deciduous_srhub_low, (9)arctic_deciduous_shrub_tall, (10)arctic_deciduous_shrub_alder,
    #                 (11)arctic_forb, (12)arctic_dry_graminoid,(13)arctic_wet_graminoid
    bandinfos={'bands':["arctic_lichen",
                    "arctic_bryophyte",
                    "arctic_forb",
                    "arctic_graminoid",
                    "arctic_evergreen_shrub",
                    "arctic_deciduous_shrub"],
                'pftnum': [1,2,11,12,5,8]
               };
    
    alldata = {}
    for i in range(len(bandinfos['pftnum'])):
        file = inputpath+rasterfiles[i]
        xx, yy, crs_res, crs_wkt, data = geotiff2nc(file, outdata=True)
        alldata[bandinfos['pftnum'][i]] = data[0] # only 1 banded data
    
    if not normallized:
        for i in range(len(bandinfos['pftnum'])):
            if i==0:
                sumall = alldata[bandinfos['pftnum'][i]]
            else:
                sumall = sumall+alldata[bandinfos['pftnum'][i]]
        zerosum_indx = np.where(((sumall==0.0) & (~sumall.mask)))
        nonzero_indx = np.where(((sumall>0.0) & (~sumall.mask)))
        for i in bandinfos['pftnum']: 
            alldata[i][nonzero_indx] = alldata[i][nonzero_indx] \
                                        /sumall[nonzero_indx]
            alldata[i][zerosum_indx] = 0.0  # those probably is fraction of water bodies
    
    alldata['bandinfos'] = bandinfos
    alldata['xx'] = xx
    alldata['yy'] = yy
    alldata['crs_res'] = crs_res
    alldata['crs_wkt'] = crs_wkt
    
    return alldata
    

# 


def arctic_landunit_natpft_fromcavm_jkumaretal(cavm_domainnc='./domain.lnd.pan-arctic_CAVM.1km.1d.c241018.nc', \
                                            vegedonly_pftdata={}, \
                                            outnc_surf='', \
                                            outdata=False):
    #
    # CAVM provides freshwater, barens fractions
    cavm_alldata = Dataset(cavm_domainnc,'r')
    landcov = cavm_alldata['grid_code'][...]
    xc = cavm_alldata['xc'][...]
    yc = cavm_alldata['yc'][...]
    xv = cavm_alldata['xv'][...]
    yv = cavm_alldata['yv'][...]
    area_km = cavm_alldata['area_LAEA'][...]  # appears I got something wrong in orginal data
    area_km[...] = 1.0
    
    cavm_alldata.close()

    # fresh-water (open water), as PCT_LAKE in surfdata
    pct_water = np.zeros_like(landcov)
    pct_water[np.where(landcov==91)] = 1.0

    # glaciers, as PCT_GLACIER in surfdata
    pct_glacier = np.zeros_like(landcov)
    pct_glacier[np.where(landcov==93)] = 1.0
    
    # partial barren, assumed as category 0 in natural PFTs, i.e. not_vegetated
    # (NOTE: this will have to merged into veged only, and re-normalled to 100%)    
    pct_barren = np.zeros_like(landcov)
    #codes are: 1-B1, 2-B2a, 3-B3, 4-B4, 5-B2b, 
    # assuming 25% veged for B1, and 50% for other barren, if non-zeros existed in veged data
    pct_barren[np.where(landcov==1)] = 0.25
    pct_barren[np.where(landcov==2)] = pct_barren[np.where(landcov==2)]+0.50
    pct_barren[np.where(landcov==3)] = pct_barren[np.where(landcov==3)]+0.50
    pct_barren[np.where(landcov==4)] = pct_barren[np.where(landcov==4)]+0.50
    pct_barren[np.where(landcov==5)] = pct_barren[np.where(landcov==5)]+0.50

    pct_vegedonly = np.ones_like(landcov)
    pct_vegedonly = pct_vegedonly - pct_water - pct_glacier - pct_barren

    # pct_barren + pct_vegedonly, assumed to be PCT_NATVEG in surfdata
    pct_natveg = pct_barren + pct_vegedonly
    pct_barren = np.zeros_like(pct_natveg)
    pct_barren[np.where(pct_natveg>0.0)] = pct_barren[np.where(pct_natveg>0.0)] \
                                            /pct_natveg[np.where(pct_natveg>0.0)]
    pct_vegedonly[np.where(pct_natveg>0.0)] = pct_vegedonly[np.where(pct_natveg>0.0)] \
                                            /pct_natveg[np.where(pct_natveg>0.0)]

    # grids containing veged_natpft data (ormalize PFT fractions) from JKumar & TQ Zhang et al.
    # note that the resolution is about 20m, but in regular lat/lon mesh
    masked_pts = {}
    if 'xc' in vegedonly_pftdata.keys() and 'yc' in vegedonly_pftdata.keys():
        # fully paired centroid xc/yc, unstructured or flatten 
        masked_pts['xc'] = vegedonly_pftdata['xc']
        masked_pts['yc'] = vegedonly_pftdata['yc']
    elif 'xx' in vegedonly_pftdata.keys() and 'yy' in vegedonly_pftdata.keys():
        # xx/yy meshed, need to flatten for easier point search
        x_axis=vegedonly_pftdata['xx']
        y_axis=vegedonly_pftdata['yy']
        yy,xx=np.meshgrid(y_axis, x_axis, indexing='ij')
        masked_data = vegedonly_pftdata[vegedonly_pftdata['bandinfos']['pftnum'][0]]
        masked_pts['xc']=xx[~masked_data.mask]
        masked_pts['yc']=yy[~masked_data.mask]
    #        
    cavm_grids={}
    cavm_grids['xc'] = xc
    cavm_grids['yc'] = yc
    cavm_grids['xv'] = xv
    cavm_grids['yv'] = yv
    cavm_grids_sub, grids_uid, grids_maskptsid = grids_nearest_or_within( \
                                            src_grids=cavm_grids, masked_pts=masked_pts, \
                                            remask_approach = 'within', \
                                            keep_duplicated=False)
    # only need trunked data
    xc = xc.flatten()[grids_uid]
    yc = yc.flatten()[grids_uid]
    xv = xv.flatten()[grids_uid]
    yv = yv.flatten()[grids_uid]
    area_km = area_km.flatten()[grids_uid]
    pct_water = pct_water.flatten()[grids_uid]
    pct_barren = pct_barren.flatten()[grids_uid]
    pct_natveg = pct_natveg.flatten()[grids_uid]
    
    pct_vegedonly = pct_vegedonly.flatten()[grids_uid]
    indx_vegedonly = np.where(pct_vegedonly>0.0)    
    # the following are, for data above (exactly matching up with), standared PFT names and order in 14 arcticpft physiology parameters
    # a full list is: (0)not_vegetated, (1)arctic_lichen, (2)arctic_bryophyte,
    #                 (3)arctic_needleleaf_tree, (4)arctic_broadleaf_tree,
    #                 (5)arctic_evergreen_shrub_dwarf, (6)arctic_evergreen_shrub_tall,
    #                 (7)arctic_deciduous_shrub_dwarf, (8)arctic_deciduous_srhub_low, (9)arctic_deciduous_shrub_tall, (10)arctic_deciduous_shrub_alder,
    #                 (11)arctic_forb, (12)arctic_dry_graminoid,(13)arctic_wet_graminoid

    pct_nat_pft = np.zeros(np.append(14, pct_natveg.shape), dtype=float)
    pct_nat_pft[0,...] = pct_barren
    for iv in range(1,14):
        pct_iv = np.zeros_like(pct_vegedonly)
        
        # get real veged PFT fraction from  JKumar & TQ Zhang et al.
        # vegedonly_pftdata, from JKumar and T-Q Zhang etal
        bandinfos = vegedonly_pftdata['bandinfos']
        for ip in bandinfos['pftnum']:
            for i in range(len(grids_uid)):
                igrid = grids_uid[i]
                newpftdata = vegedonly_pftdata[vegedonly_pftdata['bandinfos']['pftnum'][ip]]
                newpftdata = newpftdata[~newpftdata.mask]
                pct_iv[i] = np.nanmean(newpftdata[grids_maskptsid[igrid]])
            #
            # after all grids_uid done            
            pct_iv[indx_vegedonly] = pct_iv[indx_vegedonly]/pct_vegedonly[indx_vegedonly]
        
        #
        pct_nat_pft[iv,...] = pct_iv
    

    if outnc_surf!='':
        print('write to land unit and natpft to surfdata: ', outnc_surf)

    if outdata:
        surfdata = {}
        surfdata['LONGXY'] = xc
        surfdata['LATIXY'] = yc
        surfdata['AREA']   = area_km
        surfdata['PCT_GLACIER'] = pct_glacier
        surfdata['PCT_LAKE']    = pct_water
        surfdata['PCT_NATVEG']  = pct_natveg
        surfdata['PCT_NAT_PFT'] = pct_nat_pft
        
        return surfdata
    

# 
def mksrfdata_arctic(fsurfnc_all, fsrfnc_user='', vsrfnc_user='', redo=False, OriginType=False):
    
    print('#--------------------------------------------------#')
    print("Creating surface data by merging user-provided dataset")
    fsurf_arcticpft ='./surfdata_arcticpft.nc'
    
    
    #--------------------------------------
    #
    # datasets in original arctic pft tiff or nc file (with 12 arctic PFTs + 2 additional tree PFTs)
    bandinfos={'bands':["arctic_lichen",
                    "arctic_bryophyte",
                    "arctic_forb",
                    "arctic_graminoid",
                    "arctic_wet_graminoid",
                    "arctic_evergreen_shrub",
                    "arctic_evergreen_tall_shrub",
                    "arctic_deciduous_dwarf_shrub",
                    "arctic_deciduous_low_shrub",
                    "arctic_low_to_tall_willowbirch_shrub",
                    "arctic_low_to_tall_alder_shrub",
                    "arctic_needleleaf_tree",
                    "arctic_broadleaf_tree",
                    "non_vegetated",
                    "water"],
                'pftnum': [1,2,11,12,13,5,6,7,8,9,10,3,4,0,-1]
               };
    
    if OriginType:
        # lichen as not_vegetated (0), moss/forb/graminoids as c3 arctic grass (12),
        # evergreen shrub(9), deci. boreal_shrub(11),
        # evergreen boreal tree(2), deci boreal tree (3)
        bandinfos['pftnum'] = [0,12,12,12,12,9,9,11,11,11,11,2,3,0,-1]
        natpft = np.asarray(range(17))
    else:
        natpft = np.asarray(range(max(bandinfos['pftnum'])+1)) # this is the real arcticpft order number 
    
    #--------------------------------------
    # 
    dnames_elm=['gridcell']
    dnames_lndtopo=['x','y'] #should be projected [x,y]
    if (redo or os.path.isfile(fsurf_arcticpft)==False):
        src1=Dataset(fsrfnc_user,'r')

        
        pct_pft_orig = {}
        for v in src1.variables['pftname'][0:]:
            try:
                ib = bandinfos['bands'].index(v)
                iv = np.where(src1.variables['pftname'][0:]==v)
                pct_pft_orig[ib] = src1.variables['pftfrac'][iv][0]
            except ValueError:
                iv = -9999
            
            if v=='water':
                pct_water = src1.variables['pftfrac'][iv][0]

        x = np.asarray(src1.variables[dnames_lndtopo[0]])
        y = np.asarray(src1.variables[dnames_lndtopo[1]])
        src1.close()
       
        # assign data to elm-ordered array 
        # AND have to make sure all PCTs summed to 100%
        pct_nat_pft = np.zeros((len(natpft),len(y),len(x)),dtype=float)
        for ip in natpft: 
            iv = np.where(bandinfos['pftnum']==ip)[0]
            for i in iv: # in case having multiple classes
                if i>=0 and i in pct_pft_orig.keys(): 
                    pct_nat_pft[ip,...] = pct_nat_pft[ip,...]+pct_pft_orig[i]
        sum_pct = np.sum(pct_nat_pft,0)
        nonlnd_idx = np.where(sum_pct<=0.0)  # useful to mask later, either 'water' or 'non-land'
        lnd_idx = np.where(sum_pct>0.0)
        for ip in natpft: 
            pct_nat_pft[ip][lnd_idx] = pct_nat_pft[ip][lnd_idx]/sum_pct[lnd_idx]*100.0
        
        # water fraction redoing (TODO)
        
        
        
        # x/y projection to lat/lon 
        # NAD83/Alaska Albers, for AK Seward Peninsula, ifsar 5-m DEM data
        #Proj4 = +proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
        geoxy_proj_str = "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
        geoxyProj = CRS.from_proj4(geoxy_proj_str)
        # EPSG: 4326
        # Proj4: +proj=longlat +datum=WGS84 +no_defs
        lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
        Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)

        xx_new, yy_new = np.meshgrid(x, y)
        lon_new,lat_new = Txy2lonlat.transform(xx_new,yy_new)
        ij=np.where(lon_new<0.0)
        if(len(ij[0])>0): lon_new[ij]=lon_new[ij]+360.0 # for convenience, longitude from 0~360
        
        # write into nc file
        with Dataset(fsrfnc_user,'r') as src1, Dataset(fsurfnc_all,'r') as src2, Dataset(fsurf_arcticpft, "w") as dst:
            
            # new surfdata dimensions
            for dname2, dimension2 in src2.dimensions.items():
                if dname2 in dnames_elm:
                    i=np.where(np.asarray(dnames_elm)==dname2)[0][0]
                    dimension1 = src1.dimensions[dnames_lndtopo[i]]
                    len_dimension2 = len(dimension1)            # dim length from new data
                    dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
                if dname2 == 'natpft':
                    len_dimension2 = len(natpft)            # dim length from new data
                    dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
            #
            
            # pft dim
            pdim = ('natpft')
            
            # 2-D structured grids (lat/lon)
            if 'lsmlat' in dnames_elm or 'lsmlon' in dnames_elm:
                vname = 'lsmlat'
                vdim = ('lsmlat')
                vtype = src1.variables['lat'].datatype
                laty=dst.createVariable(vname, vtype, vdim)
                laty.units = 'degrees_north'
                laty.standard_name = 'latitude'
                dst[vname][...] = np.copy(src1.variables['lat'][...])
                
                vname = 'lsmlon'
                vdim = ('lsmlon')
                vtype = src1.variables['lon'].datatype
                lonx=dst.createVariable(vname, vtype, vdim)
                lonx.units = 'degrees_east'
                lonx.standard_name = 'longitude'
                dst[vname][...] = np.copy(src1.variables['lon'][...])

                gdim = ('lsmlat','lsmlon')
                pgdim = ('natpft','lsmlat','lsmlon')
            
            else:
                #1d unstructured domain
                gdim = ('gridcell')
                len_dimension2 = lat_new.size
                dst.createDimension('gridcell', len_dimension2)

                pgdim = ('natpft','gridcell')

                
                xdim = ('geox')
                dst.createDimension('geox', x.size)
                vname = 'geox'
                vtype = x.dtype
                x1d=dst.createVariable(vname, vtype, xdim, fill_value=nan)
                x1d.units = 'm'
                x1d.standard_name = 'geo-projected coordinate x'
                dst[vname][...] = x
                
                ydim = ('geoy')
                dst.createDimension('geoy', y.size)
                vname = 'geoy'
                vtype = y.dtype
                y1d=dst.createVariable(vname, vtype, ydim, fill_value=nan)
                y1d.units = 'm'
                y1d.standard_name = 'geo-projected coordinate y'
                dst[vname][...] = y
                
                vname = 'gridcell_jy'
                vtype = np.int32
                yy=dst.createVariable(vname, vtype, gdim, fill_value=-9999)
                yy.units = '-'
                yy.standard_name = 'geo-projected coordinate y indices, 0-based, for gridcell'
                dst[vname][...] = np.unravel_index(range(yy_new.size), yy_new.shape)[0]
                
                vname = 'gridcell_ix'
                vtype = np.int32
                xx=dst.createVariable(vname, vtype, gdim, fill_value=-9999)
                xx.units = '-'
                xx.standard_name = 'geo-projected coordinate x indices, 0-based, for gridcell'
                dst[vname][...] = np.unravel_index(range(xx_new.size), xx_new.shape)[1]
            
            vname = 'LATIXY'
            vdim = gdim
            vtype = src2.variables['LATIXY'].datatype
            lat2d=dst.createVariable(vname, vtype, vdim, fill_value=nan)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = lat_new.flatten()
                dst[vname][...] = np.copy(vals)
            else:
                vals = np.moveaxis(vals,0,1)
                vals[:,...] = np.copy(src1.variables['lat'][...])
                dst[vname][...] = np.moveaxis(vals,1,0)
            del vals
            
            vname = 'LONGXY'
            vdim = gdim
            vtype = src2.variables['LONGXY'].datatype
            lon2d=dst.createVariable(vname, vtype, vdim, fill_value=nan)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = lon_new.flatten()
            else:
                vals[:,...] = lon_new
            dst[vname][...] = np.copy(vals)
            del vals
                        
            vname = 'natpft'
            vdim = pdim
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, vdim, fill_value=-99)
            dst[vname].setncatts(src2[vname].__dict__)
            vals = dst.variables[vname][...]
            vals[:] = natpft
            dst[vname][...] = np.copy(vals)
            del vals

            #double PCT_NAT_PFT(natpft, lsmlat, lsmlon) ;
            vname = 'PCT_NAT_PFT'
            vdim = pgdim
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, vdim, fill_value=nan)
            dst[vname].setncatts(src2[vname].__dict__)
            vals = dst.variables[vname][...]
            if 'gridcell' in vdim:
                for ip in natpft:
                    vals[ip,...] = pct_nat_pft[ip].flatten()
            else:
                vals[:,...] = pct_nat_pft
            dst[vname][...] = np.copy(vals)
            del vals

            #need to re-calculate 'water' fraction from above
            vals_sum = np.sum(pct_nat_pft, axis=0)
            
            # if 'sum' is less than 100, the residue is actually 'water' in original data
            # then need to adjust PCT_NATVEG in a gridcell
            sumpft = np.ones_like(vals_sum)*100.0
            idx=np.where(vals_sum<100.0) # excluding 100 and above
            if len(idx[0])>0:
                sumpft[idx]=vals_sum[idx]
                vname = 'PCT_NATVEG'
                vdim = gdim
                vtype = src2.variables[vname].datatype
                dst.createVariable(vname, vtype, vdim)
                dst[vname].setncatts(src2[vname].__dict__)
                vals = dst.variables[vname][...]
                if 'gridcell' in vdim:
                    vals[:,...] = sumpft.flatten()
                else:
                    vals[:,...] = sumpft
                dst[vname][...] = np.copy(vals)
                del vals
                
                # 'water' in original data shall be called 'lake' in ELM land units
                vname = 'PCT_LAKE'
                vtype = src2.variables[vname].datatype
                dst.createVariable(vname, vtype, vdim)
                dst[vname].setncatts(src2[vname].__dict__)
                vals = dst.variables[vname][...]
                if 'gridcell' in vdim:
                    vals[:,...] = 100.0-sumpft.flatten()
                else:
                    vals[:,...] = 100.0-sumpft
                dst[vname][...] = np.copy(vals)
                del vals
        # 

#--------------------------------------------------------------------
def main():

    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """  
    input_path  = './'
    output_path = './'

    # read-in geotiff data from Zhang and Jitu etal
    vegedonly_pftdata_toolik=arctic_veged_fromraster_jkumaretal(\
        inputpath='/Users/f9y/Desktop/NGEE-P4_sites/ToolikFieldStation_TFS/toolik_clip_2024_11_14/')

    arctic_landunit_natpft_fromcavm_jkumaretal( \
        cavm_domainnc='./domain.lnd.pan-arctic_CAVM.1km.1d.c241018.nc', \
        vegedonly_pftdata=vegedonly_pftdata_toolik, \
        outnc_surf='surfdata_test.nc')
      
if __name__ == '__main__':
    main()




