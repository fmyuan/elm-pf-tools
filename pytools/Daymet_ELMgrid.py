#!/usr/bin/env python

import os, sys, time, math
import glob
import re
import numpy as np
from datetime import datetime, date
from matplotlib.dates import date2num, num2date

from optparse import OptionParser

from numpy import long, int16
from netCDF4 import Dataset
from copy import deepcopy

import netCDF4
from Modules_CLMoutput_nc4 import CLM_NcRead_1simulation
from Modules_CLMoutput_nc4 import CLMvar_1Dtseries

from pyproj import Proj
from pyproj import Transformer
from pyproj import CRS

#------------------------------------------------------------------------------------------------------------------
def Daymet_ELM_gridmatching(Grid1_Xdim, Grid1_Ydim, Grid2_x, Grid2_y, \
                         Grid1ifxy=False, Grid2ifxy=True, Grid1_cells=()):
    
    # Match Grid2 within each Grid1, so return Grid2-xy index for each Grid1 cell
    # by default, (1) Grid1 (parent) in lon/lat (aka ifxy=False), Grid2 in geox/y (aka ifxy=True)
    #             (2) Grid1 XY is grid-mesh-nodes, Grid2xy is grid-centroids. Then is good for searching Grid2 in Grid1
    #             (3) all cells in Grid1 are assigned Grid2's cell-index - maybe only those indiced in 'Grid1_cells'
    
    # it's supposed that: Grid1 X/Y are in 1-D regularly-intervaled (may not evenly) nodes along axis
    # while, Grid2 might be either like Grid2 or in 2-D mesh.
    if (len(Grid2_x.shape)<2): 
        # Grid2 must be converted to 2D paired x/y mesh, if not
        Grid2_xx, Grid2_yy = np.meshgrid(Grid2_x, Grid2_y) # mid-points of grid
    elif (len(Grid2_x.shape)==2):
        # Grid2 grid-centroids are in paired x/y for each grid
        Grid2_xx = Grid2_x
        Grid2_yy = Grid2_y
        
    if (len(Grid1_Xdim.shape)==1): #  Grid1 mesh in TWO 1-D dimensional nodes
        Grid1_x = Grid1_Xdim
        Grid1_y = Grid1_Ydim
        Grid1_xx, Grid1_yy = np.meshgrid(Grid1_Xdim, Grid1_Ydim) # nodes of grid-mesh
    else:
        #Grid1 mesh in 2-D for X/Y axis 
        print ('TODO - matching range-Grid1 in 2D mesh')
        sys.exit()
    
    # For projection conversion
    #     short lambert_conformal_conic ;
    #    lambert_conformal_conic:grid_mapping_name = "lambert_conformal_conic" ;
    #    lambert_conformal_conic:longitude_of_central_meridian = -100. ;
    #    lambert_conformal_conic:latitude_of_projection_origin = 42.5 ;
    #    lambert_conformal_conic:false_easting = 0. ;
    #    lambert_conformal_conic:false_northing = 0. ;
    #    lambert_conformal_conic:standard_parallel = 25., 60. ;
    #    lambert_conformal_conic:semi_major_axis = 6378137. ;
    #    lambert_conformal_conic:inverse_flattening = 298.257223563 ;
    #Proj4: +proj=lcc +lon_0=-100 +lat_1=25 +lat_2=60 +k=1 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m  +no_defs
    geoxy_proj_str = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
    geoxyProj = CRS.from_proj4(geoxy_proj_str)

    # EPSG: 4326
    # Proj4: +proj=longlat +datum=WGS84 +no_defs
    lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
    
    # only if 2 grids are in different projections, do tansformation
    if (Grid2ifxy and not Grid1ifxy):
        Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)
        Grid2_gxx,Grid2_gyy = Txy2lonlat.transform(Grid2_xx,Grid2_yy)
        
        ij=np.where(Grid2_gxx<0.0)
        if(len(ij[0])>0): Grid2_gxx[ij]=Grid2_gxx[ij]+360.0 # for convenience, longitude from 0~360
        ij=np.where(Grid1_x<0.0)
        if(len(ij[0])>0): Grid1_x[ij]=Grid1_x[ij]+360.0 # for convenience, longitude from 0~360

    elif (not Grid2ifxy and Grid1ifxy):
        Tlonlat2xy = Transformer.from_proj(lonlatProj, geoxyProj, always_xy=True)
        Grid2_gxx,Grid2_gyy = Tlonlat2xy.transform(Grid2_xx,Grid2_yy)
    
    else:
        Grid2_gxx = Grid2_xx
        Grid2_gyy = Grid2_yy

    # DAYMET grids' index (Grid2) included in each ELM land-grid (Grid1)
    Grid2in1_indx = {}
    if (len(Grid1_cells)<=0): 
        Grid1_ij = np.where(~np.isnan(Grid1_xx[:-1,:-1])) # cell-index rather than mesh-line index
    else:
        Grid1_ij = Grid1_cells
        
    for indx in range(len(Grid1_ij[0])): # Grid1 grid-cell no.
        j = Grid1_ij[0][indx]  # ELM output data is in (t,elmy,elmx) dimensional-order
        i = Grid1_ij[1][indx]
        
        iwst = np.min(Grid1_x[i:i+2])
        iest = np.max(Grid1_x[i:i+2])
        jsth = np.min(Grid1_y[j:j+2])
        jnth = np.max(Grid1_y[j:j+2])
        ij = np.where( ((Grid2_gxx<=iest) & (Grid2_gxx>iwst)) & \
                       ((Grid2_gyy<=jnth) & (Grid2_gyy>jsth)) )
        Grid2in1_indx[str(indx)] = deepcopy(ij)
            
        if(len(ij[0])<1):
             # none of DAYMET cell centroid inside a ELM grid, find the close one instead
            closej  = np.where((Grid2_gyy<=jnth) & (Grid2_gyy>jsth)) # do lat/y first, due to evenly-intervaled along lat/y
            if closej[0].size<=0:
                closei  = np.where((Grid2_gxx<=iest) & (Grid2_gxx>iwst)) # do lon/x first
                if(closei[0].size>0):
                    closeiy = np.argmin(abs(Grid2_gyy[closei]-(jnth+jsth)/2.0))
                    closeij = (np.asarray(closei[0][closeiy]),np.asarray(closei[1][closeiy]))
                else:
                    closeij = deepcopy(closei)
            else:
                closejx  = np.argmin(abs(Grid2_gxx[closej]-(iwst+iest)/2.0))
                closeij = (np.asarray(closej[0][closejx]),np.asarray(closej[1][closejx]))
            Grid2in1_indx[str(indx)] = deepcopy(closeij)
    
    # done with all grids
    return Grid2in1_indx, Grid2_gxx, Grid2_gyy
 
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------

# Write to geo-referenced CF compliant nc file, if filename given

def Write1GeoNc(vars, vardatas, ptxy=[], ncfname='', newnc=True):
    # INPUTS: vars      - variable names, separated by ','. if 'all' means every var-key in 'vardatas'
    #         vardatas  - python np list, with dicts/data
    #    (optional) ptxy    - paired lon/lat (x/y), in [x/lon,y/lat] (only 1 point)
    #    (optional) ncfname - if not empty, write processed data into (geo)NC file
    #    (optional) newnc   - if not True, write into existed 'ncfname' nc file
    # OUTPUTS (optional): if NOT write to NC file and 'vars' only has 1 variable
    #                     Output year, doy, lon, lat, data for 'vars' in np.array  
    
    if ncfname !='':
        if not ncfname.endswith('.nc'): ncfname = ncfname+'.nc'
        if newnc:
            if os.path.isfile(ncfname): os.system('rm -rf '+ncfname)
            ncfile = netCDF4.Dataset(ncfname, mode='w',format='NETCDF4') 
            print('Create and Write NC file: '+ncfname)
        else:
            ncfile = netCDF4.Dataset(ncfname, mode='a',format='NETCDF4') 
            print('Write NC file: '+ncfname)


    # mid of day
    mid_day = vardatas['date']
    daynums = np.asarray([date2num(x) for x in mid_day])
    year    = np.asarray([x.year for x in mid_day])
    doy0    = np.asarray([date2num(date(x,1,1)) for x in year])
    doy     = daynums-doy0+1
    YEARLY  = False
    if(np.min(doy) == np.max(doy)): YEARLY = True
    try:
        nt = len(date2num(mid_day))
    except:
        nt = 1

    # Construct the grid in lat/lon.
    if 'lat' in vardatas.keys():
        xlon = vardatas['lon']
        xlat = vardatas['lat']
        PROJECTED = False
    elif 'geox' in vardatas.keys():
        xlon = vardatas['geox']
        xlat = vardatas['geoy']
        PROJECTED = True
        
    # extracting pts, if specified
    if(len(ptxy)>1):
        d = abs(xlon-ptxy[0])
        ix = np.where(d==np.amin(d))
        lon = xlon[ix]
        d = abs(xlat-ptxy[1])
        iy = np.where(d==np.amin(d))
        lat = xlat[iy]
    else:
        lon = xlon
        lat = xlat

    #write to nc file
    DONE_header = False
    if not newnc: DONE_header = True
    DONE_time = False
    
    if vars[0]=='all': 
        vars=vardatas.keys()
        
    for varname in vars:
        #varname = 'Day_CMG_Snow_Cover'

        # header only needs to be done once
        if not DONE_header:
            if ncfname !='':
                # dimensions for nc file
                time_dim= ncfile.createDimension('time', None)

                if PROJECTED:
                    lon_dim = ncfile.createDimension('geox',  lon.size)
                    lat_dim = ncfile.createDimension('geoy',  lat.size)

                    vlat = ncfile.createVariable('geoy', np.float32, ('geoy',))
                    vlat.units = 'meters'
                    vlat.long_name = 'Northing (Lambert Conformal Conic projection)'
                
                    vlon = ncfile.createVariable('geox', np.float32, ('geox',))
                    vlon.units = 'meters'
                    vlon.long_name = 'Easting (Lambert Conformal Conic projection)'

                    vproj = ncfile.createVariable('lambert_conformal_conic', np.int32)
                    vproj.grid_mapping_name = "lambert_conformal_conic"
                    vproj.longitude_of_central_meridian = -100.
                    vproj.latitude_of_projection_origin = 42.5
                    vproj.false_easting = 0.
                    vproj.false_northing = 0.
                    vproj.standard_parallel = 25., 60.
                    vproj.semi_major_axis = 6378137.
                    vproj.inverse_flattening = 298.257223563 

                else:
                    lon_dim = ncfile.createDimension('lon',  lon.size)
                    lat_dim = ncfile.createDimension('lat',  lat.size)
                    vlat = ncfile.createVariable('lat', np.float32, ('lat',))
                    vlat.units = 'degrees_north'
                    vlat.long_name = 'latitude'
                
                    vlon = ncfile.createVariable('lon', np.float32, ('lon',))
                    vlon.units = 'degrees_east'
                    vlon.long_name = 'longitude'

                vlat[:] = lat
                vlon[:] = lon
        
                # time, create only
                vdaysnum = ncfile.createVariable('daysnum', np.float32, ('time',))
                vdaysnum.units = 'days'
                vdaysnum.long_name = 'days since 1980-01-01 UTC + 1'

                vdate = ncfile.createVariable('date', np.unicode_, ('time',))
                vdate.units = ''
                vdate.long_name = 'date in standard python-datetime calendar'

                vdoy = ncfile.createVariable('doy', np.int16, ('time',))
                vdoy.units = 'doy'
                vdoy.long_name = 'day of year'

                vyear = ncfile.createVariable('year', np.int16, ('time',))
                vyear.units = 'year'
                vyear.long_name = 'year in format yyyy'

                # when using VISIT, it shows 'time' as cycle of multiple time-series
                vtime = ncfile.createVariable('time', np.int16, ('time',))
                if YEARLY:
                    vtime.units = 'year'
                    vtime.long_name = 'year in format yyyy'
                    
                else:
                    # DOY is much useful info in this case
                    vtime.units = 'doy'
                    vtime.long_name = 'day of year'

            
                #global attributes
                ncfile.data_source = ('Daymet Software Version 3.0,' +
                                      'Please see http://daymet.ornl.gov/ for current Daymet data citation information')
                if ('ELM' in ncfname):
                    ncfile.data_source2 = ('E3SM v1.1 Land Model offline simulations for No60 and above, ' +
                                           'For NGEE-Arctic Project sponsored by DOE Office of Science')

                if ('ELM' in ncfname):
                    ncfile.history = '2020-06-29: ELM forcing data assinged for each Daymet grid-cell.'
                else:
                    ncfile.history = '2020-06-29: Daymet data aggregating and re-projecting into ELM grid.'
                    
                ncfile.contact = 'F.-M. Yuan, CCSI/ESD-ORNL. yuanf@ornl.gov'
            # done if ncfname !='':
            
            DONE_header = True
        
        # write time
            
        if not DONE_time:
            
            
            if ncfname !='':
                if newnc:
                    vdaysnum[0:nt] = daynums
                    vdate[0:nt] = np.asarray([np.unicode_(x) for x in mid_day])
                    vdoy[0:nt]  = doy
                    vyear[0:nt] = year
                    if YEARLY:
                        vtime[0:nt] = year
                    else:
                        vtime[0:nt] = doy
                        
                else:
                    vdaysnum = ncfile.variables['daysnum']
                    prv_nt = len(vdaysnum)
                    vdaysnum[prv_nt:prv_nt+nt] = daynums
                
                    vdate = ncfile.variables['date']
                    vdate[prv_nt:prv_nt+nt] = np.asarray([np.unicode_(x) for x in mid_day])
                
                    vdoy = ncfile.variables['doy']
                    vdoy[prv_nt:prv_nt+nt] = doy

                    vyear = ncfile.variables['year']
                    vyear[prv_nt:prv_nt+nt] = year

                    vtime = ncfile.variables['time']
                    if YEARLY:
                        vtime[prv_nt:prv_nt+nt] = year
                    else:
                        vtime[prv_nt:prv_nt+nt] = doy

            # done write to nc (if ncfname !='':)
            
            DONE_time = True
        
        # 
        # appears vardatas are in S-N/E-W ordered, so must be flip over
        data = vardatas[varname]
        data = np.float32(data) # data type is 'uint8', convert to short (othwise cannot be read by Visit)
        
        if ncfname !='':
            if newnc:
                if PROJECTED:
                    vtemp = ncfile.createVariable(varname, np.float32, ('time','geoy','geox')) # note: unlimited dimension is leftmost
                else:
                    vtemp = ncfile.createVariable(varname, np.float32, ('time','lat','lon')) # note: unlimited dimension is leftmost
                
                vtemp.units = ''
                vtemp.standard_name = varname.strip() # this is a CF standard name
                
                if ('prcp' in varname or 'SNOW' in varname or 'RAIN' in varname):
                    vtemp.long_name = "daily total precipitation/SNOW/RAIN"
                    vtemp.units = "mm/day"
                    vtemp.missing_value = -9999.
                    vtemp.coordinates = "geoy geox"
                    vtemp.grid_mapping = "lambert_conformal_conic"
                    vtemp.cell_methods = "area: mean time: sum"
                    
                vtemp[0:nt,:,:] = data.astype(np.float32)
            else:
                vtemp = ncfile.variables[varname]
                vtemp[prv_nt:prv_nt+nt,:,:] = data.astype(np.float32)
    # done with 'for varname in vars'
    
    # end of for varname in vars:
    if ncfname !='': ncfile.close()
    
    
#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
parser.add_option("--daymetheader", dest="daymetheader", default="", \
                  help = "DAYMET Netcdf4 file header with path ")
parser.add_option("--elmheader", dest="elmheader", default="", \
                  help = "ELM output Netcdf file header with path but no .nc ")
parser.add_option("--daymet_varname", dest="daymet_varname", default="prcp", \
                  help = "daymet Netcdf file's variable name  to process")
parser.add_option("--elm_varname", dest="elm_varname", default="", \
                  help = "ELM output Netcdf file's variable name to process")
parser.add_option("--ptslon", dest="ptx", default=None, \
                  help = "point(s)/range longitude, in point or range(:) format, to extract data \
                          default 'None' ")
parser.add_option("--ptslat", dest="pty", default=None, \
                  help = "point(s)/range latitude, in point or range(:) format, to extract data \
                          default 'None' ")
parser.add_option("--clmout_timestep", dest="clmout_ts", default="daily", \
                  help="clm output variable timestep (default = 'daily', other option monthly)")
parser.add_option("--startyr", dest="startyr", default="1", \
                  help="clm output starting year to plot (default = 1, i.e. first year of simulation" \
                   " and can be user-defined)")
parser.add_option("--endyr", dest="endyr", default="", \
                  help="clm output ending year to plot (default = none, i.e. end of simulation)")
parser.add_option("--lookup_snowfreeseason", dest="lookup_snowfreeseason", default=False, \
                  help = " lookup snow ending/starting doy in a year and write a Netcdf4 file, default: FALSE ", action="store_true")

(options, args) = parser.parse_args()

#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)
cwdir = options.workdir

if (options.elmheader == '' and options.daymetheader == ''):
    print('MUST input file name header, including fullpath, by " --elmheader=??? or --daymetheader"')
    sys.exit()

# 
if options.ptx is not None:
    ptx = re.split(',|:|;| ',options.ptx)
    ptx = np.asarray(ptx,dtype=np.float)
else:
    ptx = []
if options.pty is not None:
    pty = re.split(',|:|;| ',options.pty)
    pty = np.asarray(pty,dtype=np.float)
else:
    pty = []
##

startdays = (int(options.startyr)-1)*365.0
if(options.clmout_ts=='daily'): 
    startdays=startdays+1.0
elif(options.clmout_ts=='monthly'): 
    startdays=startdays-1.0
    
enddays = -9999
if(options.endyr !=""): enddays = int(options.endyr)*365.0


if not cwdir.endswith('/'): cwdir = cwdir+'/'

daymet_varname = options.daymet_varname #'Days_snowfree'#'snowcov'
elm_varname = options.elm_varname #'Days_snowfree'#'FSNO'

# sea cell constant
FillValue_SEA = -999.0
FillValue_LND = 0.0

if('prcp' in daymet_varname or 'SNOW' in elm_varname): 
    # for better numerical range, set sea-cell filled value to -55
    FillValue_SEA = -22.0
    FillValue_LND = -11.0

#------------------------------------------------------------------------------

# ELM out data reading
if (options.elmheader != ""):
    elmpathfileheader = options.elmheader
    
    ftype = 'nc'

    alldirfiles = sorted(glob.glob("%s*.%s" % (elmpathfileheader, ftype)))
    if(len(alldirfiles)<=0):
        sys.exit("No file exists - %s*.%s IN %s" %(elmpathfileheader, ftype, cwdir))
    else:
        print('Total Files of ELM outputs: '+str(len(alldirfiles)))

    ncfileheader = elmpathfileheader.split('/')[-1]
    elm_odir   = elmpathfileheader.replace(ncfileheader,'')
    if(elm_odir.strip()==''):elm_odir='./'
    elmfincl = 'h0'
    if ('.h1.' in alldirfiles[0]): 
        elmfincl='h1'
    elif ('.h2.' in alldirfiles[0]): 
        elmfincl='h2'
    nonheader = '.clm2.'+elmfincl
    if ncfileheader.endswith('.'): nonheader=nonheader+'.'
    if (nonheader in ncfileheader):
        ncfileheader = ncfileheader.replace(nonheader,'') # need to remove '.clm2.h?.' for use in next line
    
    # read-in datasets from 1 simulation
    nx, ny, nlgrnd, nldcmp, ncolumn, npft, varsdata, varsdims, ttunits = \
        CLM_NcRead_1simulation(elm_odir, \
                           ncfileheader, \
                           elmfincl, \
                           False, \
                           [elm_varname], \
                           startdays, enddays, \
                           False)

    if2dgrid = True
    if('topo' in varsdata.keys()):
        if(len(varsdata['topo'].shape)==1): if2dgrid = False
    # names for variable and its time-axis
    vars_list = list(varsdata.keys())    # var_list is in format of 'h*_varname', in which 'varname' is the real variable names
    for hv in vars_list:
        if re.search(elm_varname, hv): 
            var_h = hv
            hinc  = hv.replace(elm_varname,'')
            var_t = '%stime' %hinc
            break
    vdata = varsdata[var_h]
    vdims = varsdims[var_h]
    tt = varsdata[var_t]   # time unit: days (default)
    
    # processing original data
    t, gdata, sdata, zdim_indx, pdim_indx = \
        CLMvar_1Dtseries(tt, vdims, vdata, nx, ny, -999, -999, if2dgrid, \
                    False, False)
    t = np.asarray(t)

    # ELM output usually is grid-wised
    if(if2dgrid):
        elm_vdata = np.reshape(gdata,(len(t),ny,nx))
    else:
        elm_vdata = gdata
    if(elm_varname == 'SNOW' or elm_varname == 'RAIN'):         
        ij=np.where(~np.isnan(elm_vdata))
        elm_vdata[ij] = elm_vdata[ij]*86400.0 # -> mm/day
    
    # lon/lat of ELM output
    elmx = varsdata['lon']
    elmy = varsdata['lat']
    elmx_res = 0.50
    elmy_res = 0.50
    if ('landmask' in varsdata.keys()):
        landmask = varsdata['landmask'] # ELM 2-D data is in [y,x] order
    else:
        landmask = np.full((ny,nx),0.0)
        landmask[np.where((elm_vdata[0]>=0.0) & (~np.isnan(elm_vdata[0])))] = 1.0

    # extract point/range data, if any
    if (len(ptx)==1): # point
        if ptx<0.0: ptx=ptx+360.0
        pti = np.squeeze(np.argmin(abs(elmx-ptx)))
    elif(len(ptx)==2): # range
        ptx[np.where(ptx<0.0)] = ptx[np.where(ptx<0.0)]+360.0
        pti = np.squeeze(np.where((elmx>=min(ptx)) & (elmx<=max(ptx))))
    else:
        pti = np.asarray([])
    if pti.size>0:
        if pti.size==1:
            elmx = elmx[pti,None]
            elm_vdata = elm_vdata[:,:,pti,None]
            landmask = landmask[:,pti,None]
        else:
            elmx = elmx[pti]
            elm_vdata = elm_vdata[:,:,pti]
            landmask = landmask[:,pti]
    
    if (len(pty)==1): # point
        ptj = np.squeeze(np.argmin(abs(elmy-pty)))
    elif(len(pty)==2): # range
        ptj = np.squeeze(np.where((elmy>=min(pty)) & (elmy<=max(pty))))
    else:
        ptj = np.asarray([])
    if ptj.size>0:
        if ptj.size==1:
            elmy = elmy[None,ptj]
            elm_vdata = elm_vdata[:,None,ptj,:]
            landmask = landmask[None,ptj,:]
        else:
            elmy = elmy[ptj]
            elm_vdata = elm_vdata[:,ptj,:]
            landmask = landmask[ptj,:]

    
    # by default, time is in unit of days since '1850-01-01 00:00:00', without leap-year (SO CANNOT use date/time functions from python)
    if ('days' in ttunits):
        elm_t0 = 1850*365.0
        elm_daynums = t+elm_t0
    elif('year' in ttunits):
        elm_daynums = t*365.0 # first day of year t (i.e. year)
    else:
        elm_daynums = t
        
    
#-------------------------------------------------------------------------
# Daymet variable reading from daily NC dataset
if (options.daymetheader != ""):
    daymetpathfileheader = options.daymetheader

    alldirfiles = sorted(glob.glob("%s*.%s" % (daymetpathfileheader, ftype)))
    if(len(alldirfiles)<=0):
        sys.exit("No file exists - %s*.%s IN %s" %(daymetpathfileheader, ftype, cwdir))
    else:
        print('Total Files of DAYMET NC datasets: '+str(len(alldirfiles)))

    ncfileheader = daymetpathfileheader.split('/')[-1]
    daymet_odir   = daymetpathfileheader.replace(ncfileheader,'')
    if(daymet_odir.strip()==''):daymet_odir='./'

    for ncfile in alldirfiles:
        print ('Processing - ', ncfile)
        f = Dataset(ncfile,'r')
        
        if ncfile == alldirfiles[0]:
            daymetx = np.asarray(f.variables['x']) # 1-D geox
            daymety = np.asarray(f.variables['y']) # 1-D geoy

            daymet_lon = np.asarray(f.variables['lon']) #2-D longitude
            ij=np.where(daymet_lon<0.0)
            if(len(ij[0]>0)): daymet_lon[ij]=daymet_lon[ij]+360.0
            daymet_lat = np.asarray(f.variables['lat']) #2-D latitude
            
            val_filling = f.variables[daymet_varname]._FillValue
            val_missing = f.variables[daymet_varname].missing_value
            
        # read-in datasets from one nc file 
        if ('time' in f.variables.keys()):
            tt = np.asarray(f.variables['time']) #units - days since date/time 1980-01-01 00:00:00
            if ('days since' in f.variables['time'].units):
                # note: Daymet uses leap_year system, but remove doy of 366 data to keep same length of days in a yearr
                tt_t0 = date2num(date(1980,1,1))
                tt = tt + tt_t0
        vdata = np.asarray(f.variables[daymet_varname])
        ij=np.where(vdata==val_filling)
        if (vdata[ij].size>0): vdata[ij]=np.nan
        ij=np.where(vdata==val_missing)
        if (vdata[ij].size>0): vdata[ij]=np.nan

        # ---------------------------------------------------------------------
        # need to merge ELM  datasets for comparison
        if (options.elmheader !=""):
            
            # timing variables for counting/recording when ELM-DAYMET time-matching up with each other
            # multiple files and time-steps
            daynums_all   = np.empty((0),np.float32)
            daymet_it_all = np.empty((0),np.int)
            elm_it_all    = np.empty((0),np.int)

            # ELM centroid lon/lat <==> DAYMET  centroid geox/geoy
            if(ncfile==alldirfiles[0]): #only need to do once

                # the DAYMET grid centroids in lon/lat mesh
                daymet_gx = daymet_lon # mid-points of grid
                daymet_gy = daymet_lat
                
                # ELM grid-mesh nodes: elmx/y are grid-centroids
                if elmx.size>1:
                    halfx = np.mean(np.diff(elmx))/2.0
                    elmnodex = elmx + 0.50*np.hstack((np.diff(elmx),2.0*halfx))
                    elmnodex = np.hstack((elmx[0]-halfx, elmnodex))
                elif elmx.size == 1:
                    halfx = elmx_res/2.0
                    elmnodex = np.asarray([elmx[0]-halfx, elmx[0]+halfx])
                if elmy.size>1:
                    halfy = np.mean(np.diff(elmy))/2.0
                    elmnodey = elmy + 0.50*np.hstack((np.diff(elmy),2.0*halfy))
                    elmnodey = np.hstack((elmy[0]-halfy, elmnodey))
                elif elmy.size == 1:
                    halfy = elmy_res/2.0
                    elmnodey = np.asarray([elmy[0]-halfy, elmy[0]+halfy])
                
                # truncating ELM space (may save some computing time)
                daymet_lonmax=np.max(np.max(daymet_lon,0))
                daymet_lonmin=np.min(np.min(daymet_lon,0))
                i = np.asarray(np.where((elmnodex>=daymet_lonmin) & (elmnodex<=daymet_lonmax)))[0]
                elmnodex  = elmnodex[i] #1-D
                elmx      = elmx[i[:-1]] # mid-point less 1 element
                landmask  = landmask[:,i[:-1]] # 2-D
                elm_vdata = elm_vdata[:,:,i[:-1]] #2-D

                daymet_latmax=np.max(np.max(daymet_lat,0))+1.0
                daymet_latmin=np.min(np.min(daymet_lat,0))
                j = np.asarray(np.where((elmnodey>=daymet_latmin) & (elmnodey<=daymet_latmax)))[0]
                elmnodey  = np.asarray(elmnodey[j]) #1-D
                elmy      = elmy[j[:-1]] # mid-point less 1 element
                landmask  = landmask[j[:-1],:] #2-D
                elm_vdata = elm_vdata[:,j[:-1],:] #2-D

                elm_seaij = np.where(landmask==0)
                elm_lndij = np.where(landmask==1)

                # matching daymet grid-centroids within ELM grid-mesh
                daymetinelm_indx, daymet_lon, daymet_lat = \
                    Daymet_ELM_gridmatching(elmnodex, elmnodey, \
                                                daymet_gx, daymet_gy, \
                                                Grid1ifxy=True, Grid2ifxy=True, \
                                                Grid1_cells=elm_lndij)

                    
                # need ALL time-values from ELM output, but only calculated once
                yr_elm   = np.int32(elm_daynums/365.0)
                doy_elm  = elm_daynums-np.int32(elm_daynums/365.0)*365.0+1.0
            # done with if (ncfile is the first of alldirfiles)
            
            for daymet_it in range(tt[0:].size):
                # match YEAR/DOY btw 'vdata' and 'elm_vdata' (NOT date/time due to no_leap in ELM time)
                # note: Daymet uses leap_year system, but remove doy of 366 data to keep same length of days in a yearr
                print ('Time for '+ncfile+' : '+str(daymet_it))
                day_tt = num2date(tt[daymet_it]).date()
                yr_tt  = day_tt.year
                doy_tt = np.floor(tt[daymet_it] - date2num(date(yr_tt,1,1)) + 1)
                elm_it = np.squeeze(np.where((yr_elm==yr_tt) & (doy_elm==doy_tt)))
                if(elm_it.size>0):
                    elm_it_all    = np.hstack((elm_it_all,elm_it))         # timer count 
                    daymet_it_all = np.hstack((daymet_it_all,daymet_it))          # timer count 
                    daynums_all   = np.hstack((daynums_all,tt[daymet_it])) # timer
                
                # data in DAYMET grid-cells
                vdata_it = np.float32(vdata[daymet_it,])
                vdata_fromelm=np.full(np.shape(vdata_it),np.float32(FillValue_SEA)) # to hold from one-time ELM data (pre-filled as sea) in DAYMET cells
                vdata_fromelm[np.where(~np.isnan(vdata_it))] = np.float32(FillValue_LND) # pre-filling DAYMET land cell with LND fillingvalue (aka beyond data coverage)

                # data in ELM grid-cells
                vdata_1delmcell = np.full((len(elm_lndij[0])),np.float32(FillValue_LND)) # to hold 1-D aggregated DAYMET 'vdata' by ELM lnd-grid
                vdata_std_1delmcell = np.full((len(elm_lndij[0])),np.float32(FillValue_LND)) # to hold 1-D aggregated DAYMET 'vdata' by ELM lnd-grid
                vdata_n_1delmcell = np.full((len(elm_lndij[0])),np.float32(FillValue_LND)) # to hold 1-D aggregated DAYMET 'vdata' by ELM lnd-grid
                for idx in range(len(elm_lndij[0])):
                    ij = daymetinelm_indx[str(idx)] #  paired-tuple index of daymet cells in ONE elm-grid
                    # assign DAYMET averaged to elm land-cell
                    if(ij[0].size>0):
                        i = daymetinelm_indx[str(idx)][0]
                        j = daymetinelm_indx[str(idx)][1]
                        if (vdata_it[ij].size>0):
                            vdata_1delmcell[idx] = np.nanmean(vdata_it[ij]) # DAYMET ==> ELM
                            vdata_std_1delmcell[idx] = np.nanstd(vdata_it[ij]) # DAYMET ==> ELM
                            vdata_n_1delmcell[idx] = np.count_nonzero(~np.isnan(vdata_it[ij])) # DAYMET ==> ELM
                            #vdata_1delmcell[idx] = np.mean(daymet_lon[ij]) # DAYMET ==> ELM, for testing if daymetx/y correctly mapped into ELM grid
                            #vdata_1delmcell[idx] = np.mean(daymety[j]) # DAYMET ==> ELM. for test if daymetx/y correctly mapped into ELM grid

                    # assign ELM output to DAYMET' land grids, if time matches
                    if(elm_it.size>0):
                        temp = elm_vdata[elm_it]
                        #if(elm_it.size)>1: temp = np.mean(temp,axis=0) # in-case 'daymet_it' contains multiple 'elm_it (e.g. sub-daily')
                        j = elm_lndij[0][idx]  # ELM output data is in (t,elmy,elmx) dimensional-order
                        i = elm_lndij[1][idx]
                        daymetlnd = np.where(~np.isnan(vdata_it[ij]))[0] # extract sub-set of DAYMET land-cells and to be filled in them below
                        if daymetlnd.size>0:
                            if ij[0].size<=1:
                                vdata_fromelm[(ij[0],ij[1])] = temp[(j,i)] # ELM ==> DAYMET
                            else:
                                vdata_fromelm[(ij[0][daymetlnd],ij[1][daymetlnd])] = temp[(j,i)] # ELM ==> DAYMET
                            #vdata_fromelm[ij] = elmy[j] # ELM ==> DAYMET, for testing if latitude from ELM mapping to DAYMET grids
                            #vdata_fromelm[ij] = elmx[i] # ELM ==> DAYMET, for testing if longitude from ELM mapping to DAYMET grids
                # done with all elm landcells 
                
                # when all elm land-grids are done 
                # reshape vdata_1delmcell according to ELM grids
                temp = np.full((elmy.size,elmx.size),np.float32(FillValue_SEA))
                temp[elm_lndij] = vdata_1delmcell
                vdata_2delmcell = deepcopy(temp)
                
                temp = np.full((elmy.size,elmx.size),np.float32(FillValue_SEA))
                temp[elm_lndij] = vdata_std_1delmcell
                vdata_std_2delmcell = deepcopy(temp)
                
                temp = np.full((elmy.size,elmx.size),np.float32(FillValue_SEA))
                temp[elm_lndij] = vdata_n_1delmcell
                vdata_n_2delmcell = deepcopy(temp)
                
                #----------------------------------------------------------------------------
                
                # after done with 1 time-step, need to save matched datasets in either ELM or DAYMET grids
                
                # into ELM grid system
                if (elm_it.size>0):
                    # DAYMET data aggregated into ELM cells, if sucessfully calculated
                    varname = elm_varname+'_daymet'
                    ncdata = {}
                    ncdata['date'] = [num2date(tt[daymet_it]).date()]
                    ncdata['lon']  = deepcopy(elmx)
                    ncdata['lat']  = deepcopy(elmy)

                    lndij = np.where(vdata_n_2delmcell>=0) #  ELM-land cell
                    ij = np.where(vdata_n_2delmcell>0) # cell numbers from DAYMET in a ELM-grid

                    ncdata[elm_varname] = np.full(vdata_2delmcell.shape,np.float32(FillValue_SEA))
                    ncdata[elm_varname][lndij] = FillValue_LND
                    ncdata[elm_varname][ij] = deepcopy(elm_vdata[elm_it][ij])
                    
                    ncdata[varname] = np.full(vdata_2delmcell.shape,np.float32(FillValue_SEA))
                    ncdata[varname][lndij] = FillValue_LND
                    ncdata[varname][ij] = deepcopy(vdata_2delmcell[ij])

                    ncdata[varname+'_std'] = np.full(vdata_2delmcell.shape,np.float32(FillValue_SEA))
                    ncdata[varname+'_std'][lndij] = FillValue_LND
                    ncdata[varname+'_std'][ij] = deepcopy(vdata_std_2delmcell[ij])
                    
                    ncdata[varname+'_n'] = np.full(vdata_2delmcell.shape,np.float32(FillValue_SEA))
                    ncdata[varname+'_n'][lndij] = FillValue_LND
                    ncdata[varname+'_n'][ij] = deepcopy(vdata_n_2delmcell[ij])

                    # difference btw ELM and DAYMET
                    temp = np.full(vdata_2delmcell.shape,np.float32(FillValue_SEA)) # shape as DAYMET data, but will be filled with data from ELM next-line
                    temp = ncdata[elm_varname] - ncdata[varname] # ELM-DAYMET
                    # since 'diff' could be negative, better to mark value-filled-land in either dataset as ZERO
                    ij = np.where( (ncdata[elm_varname]==FillValue_LND)  | \
                                   (np.isnan(ncdata[elm_varname])) | \
                                   (np.isinf(ncdata[elm_varname])) )
                    temp[ij] = 0.0 
                    ij = np.where( (ncdata[varname]==FillValue_LND)  | \
                                   (np.isnan(ncdata[varname])) | \
                                   (np.isinf(ncdata[varname])) )
                    temp[ij] = 0.0
                    # for best visual effect, mask sea-cell as nan
                    ij = np.where( (ncdata[elm_varname]==FillValue_SEA)  | \
                                   (ncdata[varname]==FillValue_SEA) )
                    temp[ij] = np.nan
                    ncdata[elm_varname+'_diff'] = temp

                    ncfname = ncfile.split('/')[-1] # remove directory name if any
                    ncfname = './ELM_forcing_from_'+ncfname
                    newfile = True
                    if(elm_it!=daymet_it_all[0]): newfile = False
                    Write1GeoNc([elm_varname, elm_varname+'_diff', varname, varname+'_std', varname+'_n'], \
                            ncdata, ptxy=[], ncfname=ncfname, newnc=newfile)
                
                # into DAYMET grid system
                if (len(elm_it_all)>0): 
                    # ELM simulation assigned into DAYMET grids
                    varname = daymet_varname+'_elm'
                    ncdata = {}
                    ncdata['date'] = [num2date(tt[daymet_it]).date()]
                    ncdata['geox']  = deepcopy(daymetx)
                    ncdata['geoy']  = deepcopy(daymety)
                    vdata_ij = deepcopy(vdata[daymet_it,])
                    daymetlndmask = ~np.isnan(vdata_ij) # mask land-cells
                    vdata_ij[np.where(np.isnan(vdata_ij))] = FillValue_SEA # fill sea-cell value
                
                    # mask daymet data beyond ELM boundary
                    ijout = np.where( ((daymet_lon<np.min(elmnodex)) | (daymet_lon>np.max(elmnodex)) |
                                       (daymet_lat<np.min(elmnodey)) | (daymet_lat>np.max(elmnodey))) &
                                      (daymetlndmask) ) # land-cells only
                    vdata_ij[ijout] = np.float32(FillValue_LND)
                
                    ncdata[daymet_varname] = deepcopy(vdata_ij) 
                    ncdata[varname] = deepcopy(vdata_fromelm)

                    # difference btw ELM and DAYMET
                    temp = np.full(ncdata[daymet_varname].shape,np.float32(FillValue_SEA)) # shape as Daymet data, but filled with data from ELM next-line
                    temp = ncdata[varname] - ncdata[daymet_varname] # ELM - DAYMET
                    # since 'diff' could be negative, better to mark value-filled-land in either dataset as ZERO
                    ij = np.where( (ncdata[daymet_varname]==FillValue_LND) | \
                                   (np.isnan(ncdata[daymet_varname])) | \
                                   (np.isinf(ncdata[daymet_varname])) )
                    temp[ij] = 0.0 
                    ij = np.where( (ncdata[varname]==FillValue_LND) | \
                                   (np.isnan(ncdata[varname])) | \
                                   (np.isinf(ncdata[varname])) )
                    temp[ij] = 0.0
                    # for best visual effect, mask sea-cell as nan
                    ij = np.where( (ncdata[daymet_varname]==FillValue_SEA)  | \
                                   (ncdata[varname]==FillValue_SEA) )
                    temp[ij] = np.nan
                    ncdata[daymet_varname+'_diff'] = temp

                    # write NC file(s)
                    newfile = True
                    if(elm_it!=daymet_it_all[0]): newfile = False
                    ncfname = ncfile.split('/')[-1] # remove directory name if any
                    ncfname = './ELM_forcing_for_'+ncfname
                    Write1GeoNc([daymet_varname, daymet_varname+'_diff', varname], \
                                ncdata, ptxy=[], ncfname=ncfname, newnc=newfile)
                
            # done with 'for it in range(tt[0:].size):'
        # done with 'if (options.elmheader !=""):'
    # ---------------------------------------------------------------------
    #
    # done with all daymet nc files in directory


#------------------------------------------------------------------------------------------------