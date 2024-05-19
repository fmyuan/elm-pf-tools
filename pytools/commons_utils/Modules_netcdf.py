#Python utilities for operating a netcdf file using available tools, python netcdf4, nco
import os, sys
import numpy as np
from netCDF4 import Dataset
from numpy import newaxis

#----------------------------------------------------------------------------             
def getvar(ncfile, varname):
    odata  = {}
    odata_dims = {}
    odata_attr = {}

    try:
        f = Dataset(ncfile,'r')
        if varname=="" or varname[0]=="": 
            print('FILE: '+ncfile+' Variable(s) is(are): ')
            print(f.variables.keys()) 
            return odata, odata_dims, odata_attr

    except ValueError:
        exit()
                
    # If key is in varname list and in nc file, uniquely add to a dictionary    

    if varname=="ALL" or varname[0]=="ALL": 
        for key in f.variables.keys():
            odata[key]      = np.asarray(f.variables[key])                
            odata_dims[key] = f.variables[key].dimensions
            odata_attr[key] = f.variables[key].ncattrs()

    else:
        for key in varname:           
            if(key in f.variables.keys()):
                odata[key]      = np.asarray(f.variables[key])                
                odata_dims[key] = f.variables[key].dimensions
                odata_attr[key] = f.variables[key].ncattrs()
    
    # in case no-data for 'varname'
    if(len(odata)<1):
        print('FILE: '+ncfile+' Variable(s) is(are): ')
        print(f.variables.keys())

    #
    f.close()
    return odata, odata_dims, odata_attr
 
#----------------------------------------------------------------------------             
def putvar(ncfile, varname, varvals, varatts=''):
    # varatts: 'varname::att=att_val; varname::att=att_val; ...'
    try:
        f = Dataset(ncfile,'r+')
        #
        if varname=="": 
            print('FILE: '+ncfile+' Variable(s) is(are): ')
            for key in f.variables.keys(): 
                print (key)
            return -1
        #
        else:
            if varatts!='':
                varatts=varatts.split(";")
                for varatt in varatts:
                    v=varatt.split('::')[0].strip()
                    att=varatt.split('::')[1].strip().split('=')[0]
                    att_val=varatt.split('::')[1].strip().split('=')[1]
                    if (att=='add_offset' or att=='scale_factor'):
                        # this must be done before data written, 
                        # otherwise the written would use the original add_offset and scale_factor
                        # TIP: when data written, the input valule is in unpacked and the program will do packing
                        f.variables[v].setncattr(att,float(att_val))
                    else:
                        f.variables[v].setncattr(att,att_val)

            
            if isinstance(varvals, dict): #multiple dataset for corresponding same numbers of varname
                for v in varname:
                    if v in f.variables.keys():
                        if f.variables[v].dtype!=float:
                            val=np.ma.masked_where(np.isnan(varvals[v]), varvals[v])
                        else:
                            val=np.copy(varvals[v])
                        f.variables[v][...]=val
            else:                        #exactly 1 dataset in np.array for 1 varname
                v=varname[0]
                if v in f.variables.keys():
                    if f.variables[v].dtype!=float or \
                      (f.variables[v].dtype!=np.float64 or f.variables[v].dtype!=np.float32):
                        # nan is not dtype of integers
                        val=np.ma.masked_where(np.isnan(varvals), varvals)
                    else:
                        val=np.copy(varvals)
                    f.variables[v][...]=val
            #
            f.sync()
            f.close()
            
            return 0

    except:
        return -2

#----------------------------------------------------------------------------             
# duplicate and/or expand all variables along named dimension by multiple times
def dupexpand(ncfilein, ncfileout,dim_name='',dim_len=-999, dim_multipler=1):
    # dim_len -999, no-change; dim_len 0, unlimited; dime_len positive, re-size
    # but dim_multipler must be 1 (obviously)
    if type(dim_name)==str and (dim_len>=0 and dim_multipler>1):
        print('Error: cannot have new-dim_len and multiplied len -', dim_len, dim_multipler)
        os.sys.exit(-1)
    
    with Dataset(ncfilein) as src, Dataset(ncfileout, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
   
        if type(dim_name)==str:
            # the following will allow multiple name/leng change 
            dim_name=[dim_name]
            dim_multipler=[dim_multipler]
            dim_len=[dim_len]
            

        # copy, resize or multiply dimensions
        for name, dimension in src.dimensions.items():
            len_dimension = len(dimension)
            if name in dim_name:
                indx=dim_name.index(name)
                if (dim_len[indx]>0):
                    len_dimension = dim_len[indx]
                elif (dim_multipler[indx]>1):
                    len_dimension = len(dimension)*dim_multipler[indx]
                elif (dim_len[indx]==0):
                    dimension=None
            
            dst.createDimension(name, len_dimension if not dimension.isunlimited() else None)
        #   
    
        # copy all file data except for matched-up variables, which instead multiply copied
        for name, variable in src.variables.items():

            # create variables, but will update its values later 
            # NOTE: here the variable value size updated due to newly-created dimensions above
            dst.createVariable(name, variable.datatype, variable.dimensions)
            
            # copy variable attributes all at once via dictionary after created
            dst[name].setncatts(src[name].__dict__)

            #
            varvals = np.copy(src[name][...])
            for dim in dim_name:
                if dim in variable.dimensions:
                    if dim_len[dim_name.index(dim)]>0:
                        dim_indx = variable.dimensions.index(dim)
                        resizer = dim_len[dim_name.index(dim)]
                        if(resizer>np.size(varvals, axis=dim_indx)):
                            tmp=np.zeros(dst.variables[name].shape)
                            # the following is slow, if length too big
                            tmp = varvals
                            resizer = resizer - np.size(varvals, axis=dim_indx)
                            while resizer>1:
                                varval = np.take(varvals, [0], axis=dim_indx)
                                tmp=np.concatenate((tmp,varval), axis=dim_indx)
                                resizer-=1
                        else:
                            tmp = np.take(varvals, range(resizer), axis=dim_indx)
                    
                    elif dim_multipler[dim_name.index(dim)]>1:
                        dim_indx = variable.dimensions.index(dim)
                        multipler = dim_multipler[dim_name.index(dim)]
                        tmp=varvals
                        while multipler>1:
                            tmp=np.concatenate((tmp,varvals), axis=dim_indx)
                            multipler-=1
                    #
                    varvals=tmp
            #           
                
            dst[name][...] = np.copy(varvals)
        
        #
            
    #            
    

#----------------------------------------------------------------------------             
# merge all variables along 1 named dimension
# by append one file into another, with both exactly same defined dimensions and variables
def mergefilesby1dim(ncfilein1, ncfilein2, ncfileout, dim_name):
    with Dataset(ncfilein1) as src1, Dataset(ncfilein2) as src2, Dataset(ncfileout, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src1.__dict__)
   
        if type(dim_name)==str: 
            dim_name=[dim_name]
        elif (len(dim_name)>1): # make sure only 1 dimension allowed
            print("Error: only 1 dimension allowed")
            sys.exit(-1)

        # changing dimension length for dst
        for name, dimension in src1.dimensions.items():
            if name in dim_name:
                dimension2 = src2.dimensions[name]
                d12 = []
                if name in src1.variables.keys() and name in src2.variables.keys(): # overlay or merge
                    d1 = np.asarray(src1[name])
                    d2 = np.asarray(src2[name])
                    d12, idx1, idx2  = np.intersect1d(d1, d2, return_indices=True)
                    if len(d12)>0:
                        for d in d12: 
                            d1 = d1[d1!=d]
                            d2 = d2[d2!=d]
                    len_dimension = len(d1)+len(d12)+len(d2)
                else:
                    len_dimension = len(dimension) + len(dimension2)
            else:
                len_dimension = len(dimension)
            
            dst.createDimension(name, len_dimension if not dimension.isunlimited() else None)
        #   
    
        # copy all file data except for matched-up variables, which instead multiply copied
        for name, variable in src1.variables.items():

            # create variables, but will update its values later 
            # NOTE: here the variable value size updated due to newly-created dimensions above
            dst.createVariable(name, variable.datatype, variable.dimensions)
            
            # copy variable attributes all at once via dictionary after created
            dst[name].setncatts(src1[name].__dict__)

            #
            varvals1 = np.copy(src1[name][...])
            tmp = varvals1
            dim=dim_name[0]  #only 1 dim supported
            if dim in variable.dimensions:
                dim_indx = variable.dimensions.index(dim)
                varvals2=np.copy(src2[name][...])
                if len(d12)>0:
                    # overlay cells, if any
                    varvals12 = (np.take_along_axis(varvals1,idx1,axis=dim_indx) + \
                                 np.take_along_axis(varvals2,idx2,axis=dim_indx))/2.0
                    
                    varvals1 = np.put_along_axis(varvals1, idx1, varvals12, axis=dim_indx)
                    varvals2 = np.delete(varvals2, idx2, axis=dim_indx)
                tmp=np.concatenate((varvals1,varvals2), axis=dim_indx)
                
            dst[name][...] = np.copy(tmp)
        
        #
        
    #            
    #print('done!')

#----------------------------------------------------------------------------             

#----------------------------------------------------------------------------             
# merge/overlay 2 ncfiles with same dimension names and resolution but might not be same length or range
# 
def overlayfiles(ncfilein1, ncfilein2, ncfileout):
    with Dataset(ncfilein1) as src1, Dataset(ncfilein2) as src2, Dataset(ncfileout, "w") as dst:
   
        # dimensions in common to be merge/overlay
        idx1={}
        idx2={}
        udims = {}
        for name, dimension in src1.dimensions.items():
            if (name in src2.dimensions.keys()):
                dimension2 = src2.dimensions[name]
                if name in src1.variables.keys() and name in src2.variables.keys(): # overlay or merge
                    d1 = np.asarray(src1[name])
                    d2 = np.asarray(src2[name])
                    d12 = np.union1d(d1, d2)
                    idx1[name] = np.asarray(np.where(np.isin(d12, d1))[0]) # indexing of d1 in d12, but converted to integer array
                    idx2[name] = np.asarray(np.where(np.isin(d12, d2))[0])
                    len_dimension = len(d12)
                else:
                    len_dimension = len(dimension) + len(dimension2)   # concatenating
                    idx1[name] = np.asarray(range(len(dimension)))
                    idx2[name] = len(dimension)+np.asarray(range(len(dimension2)))
                    
            else:
                len_dimension = len(dimension)
                idx1[name] = np.asarray(range(len(dimension)))
                idx2[name] = []

            if dimension.isunlimited(): # save current unlimited dim length for using below
                udims[name] = len_dimension
            
            dst.createDimension(name, len_dimension if not dimension.isunlimited() else None)
        
        #   
        # copy all file data except for matched-up variables, which do merge/overlay
        for name, variable in src1.variables.items():

            # create variables, but will update its values later 
            # NOTE: here the variable value size updated due to newly-created dimensions above
            dst.createVariable(name, variable.datatype, variable.dimensions)
            
            # copy variable attributes all at once via dictionary after created
            dst[name].setncatts(src1[name].__dict__)

            #
            varvals1 = np.copy(src1[name][...])
            tmp = varvals1
            
            
            if (name in src2.variables.keys()) and (not dst[name].dtype==np.str):
                varvals2=np.copy(src2[name][...])

                # nan-filled 'dst[name]' array for putting data from src1/src2
                dst_shape = np.asarray(dst[name].shape)
                for udim in udims.keys():
                    if udim in variable.dimensions:
                        udim_idx = variable.dimensions.index(udim)
                        # unlimited length in 'dst_shape' from dst is default as 0
                        # must assign current length, otherwise 'tmp1'/'tmp2' below is incorrect in shape
                        dst_shape[udim_idx] = udims[udim]
                temp1 = np.full(dst_shape, np.nan, dtype=np.double) # nan-filled array for src1
                temp2 = np.full(dst_shape, np.nan, dtype=np.double) # nan-filled array for src2
                
                vdim = np.asarray(variable.dimensions)
                temp_indx1 = np.asarray(np.where(np.isnan(temp1)))
                temp_indx2 = np.asarray(np.where(np.isnan(temp2))) # whole-set multi-tuples indices, with 0 dim as indice for 'vdim'
                for i in range(len(vdim)):
                    if vdim[i] in idx1.keys():  # dimension to be merge/overlay
                        idx = np.asarray(np.where(np.isin(temp_indx1[i],idx1[vdim[i]])))
                        temp_indx1 = temp_indx1[:,np.squeeze(idx)]
                        idx = np.asarray(np.where(np.isin(temp_indx2[i],idx2[vdim[i]])))
                        temp_indx2 = temp_indx2[:,np.squeeze(idx)]
                #
                vdim_indx1=np.ravel_multi_index(temp_indx1,temp1.shape)
                vdim_indx2=np.ravel_multi_index(temp_indx2,temp2.shape)
                np.put(temp1, vdim_indx1, varvals1)
                np.put(temp2, vdim_indx2, varvals2)
                temp1 = np.expand_dims(temp1, axis=0)
                temp2 = np.expand_dims(temp2, axis=0)
                tmp  = np.nanmean(np.concatenate((temp1, temp2), axis=0), axis=0)
                #
            #
            dst[name][...] = np.copy(tmp)
        #
        
    #
    #print('done!')

#----------------------------------------------------------------------------             

#----------------------------------------------------------------------------             
#Using NCO to extract all variables along named dimensions
def nco_extract(ncfilein, ncfileout,dim_names,dim_starts, dim_lens, ncksdir=""):
    
    if(not os.path.isfile(ncfilein)): 
        print('Error: invalid input NC file: '+ncfilein)
        os.sys.exit()

    if(os.path.isfile(ncfileout)): 
        os.system('rm -rf '+ncfileout)
    
    #single string/integers --> list, so that multiple dims can be done
    if (not isinstance(dim_names, (list,tuple))):
        dim_names = [dim_names]
    if (not isinstance(dim_starts, (list,tuple))):
        dim_starts = [dim_starts]
    if (not isinstance(dim_lens, (list,tuple))):
        dim_lens = [dim_lens]

    
    dim_string = ''
    for idim in range(0,len(dim_names)):
        dim_string = dim_string + ' -d ' \
                                +dim_names[idim].strip()+',' \
                                +str(dim_starts[idim])+',' \
                                +str(dim_starts[idim]+dim_lens[idim]-1)

    #
    if((ncksdir!="") & (not ncksdir.endswith('/'))):
        ncksdir = ncksdir.strip()+'/'
    os.system(ncksdir+'ncks --no_abc -O '+ dim_string+ \
          ' '+ncfilein+' '+ncfileout)

    #
    
#----------------------------------------------------------------------------             
# convert nc to csv
def nc2csv(file, varnames):
    #file = 'test01.nc'
    #varnames = ['LATIXY','LONGXY']
    odata, odata_dims, odata_attr = getvar(file,varnames)
    for key in odata.keys():
        v=np.asarray(odata[key]).ravel()
        if(key is odata.keys()[0]):
            t_v = np.hstack((key, v))
            same_dims = odata_dims[key]
        else:
            if(odata_dims[key]==same_dims):
                t_v = np.vstack( (t_v,np.hstack((key, v)) ) )
            else:
                print(key+' has different dimensions, so will NOT included in .csv file')
                            
    t_v = np.swapaxes(t_v, 1, 0)
    np.savetxt(file+'.csv', t_v, delimiter=',', fmt='%s')
    #print('done!')

#----------------------------------------------------------------------------             
# convert geotiff to ncfile
def geotiff2nc(file, bandinfos):
    # file: file name
    # bandinfos in 2-D strings: tiff file' bands and its NC name, units, long_name, standard_name
    import datetime as dt
    import rasterio
    import re

    f = rasterio.open(file, mode='r')
    alldata = f.read()
    nband,ny,nx = np.shape(alldata)
        

    # create NetCDF output file
    ncof = Dataset(file+'.nc','w',clobber=True)
    
    # geox/y coordinates in 1-D (axis), centroid
    geox = np.arange(nx)*f.transform[0]+f.transform[2]+f.transform[0]/2.0
    geoy = np.arange(ny)*f.transform[4]+f.transform[5]+f.transform[4]/2.0
        
    
   
    # create dimensions, variables and attributes:
    ncof.createDimension('geox',nx)
    ncof.createDimension('geoy',ny)

    geoxo = ncof.createVariable('geox',np.double,('geox'))
    geoxo.units = 'degrees_east (m) or lat'
    geoxo.standard_name = 'geo-reference x coordinates'

    geoyo = ncof.createVariable('geoy',np.double,('geoy'))
    geoyo.units = 'degrees_north (m) or lat'
    geoyo.standard_name = 'geo-reference y coorindates'

    georef = ncof.createVariable('crs','int')
    georef.CRS = f.crs.wkt

    bounds = ncof.createVariable('bounds','int')
    bounds.box=f.bounds

    res = ncof.createVariable('resolution','int')
    res.box=f.res

    # create variables
    if('bands' in bandinfos.keys()):
        nvars = len(bandinfos['bands'])
    else:
        exit('Must have bands name in bandinfos!')
        
    varo = {}
    for iv in range(0,nvars):
        v = bandinfos['bands'][iv]
        varo[v] = ncof.createVariable(v, np.double,  ('geoy', 'geox'), fill_value=-1.0e20)
        if 'units' in bandinfos.keys():
            varo[v].units = bandinfos['units'][iv]
        if 'long_name' in bandinfos.keys():
            varo[v].long_name = bandinfos['long_name'][iv]
        if 'standard_name' in bandinfos.keys():
            varo[v].standard_name = bandinfos['standard_name'][iv]

        # write variable
        varo[v][:,:] = alldata[iv]
        
    #write x,y
    geoxo[:]=geox
    geoyo[:]=geoy

    ncof.close()
#
#----------------------------------------------------------------------------             
# average specific variable(s) along named dimension and write back for all
def varmeanby1dim(ncfilein, ncfileout,dim_name,var_name='ALL',var_excl=''):
    dim_name = dim_name.strip()
    if dim_name=='':
        print('Error: please provide a valid dim_name')
        sys.exit()
    
    with Dataset(ncfilein) as src, Dataset(ncfileout, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
   
        if type(var_name)==str:
            # in case variable(s) not in a string array
            if var_name=='ALL':
                var_name = list(src.variables.keys())
            else:
                var_name=var_name.strip().split(',')
        else:
            print('Error in ',var_name)
            sys.exit()
        if var_excl!='' and type(var_excl)==str:
            var_remv = var_excl.strip().split(',')
            for v in var_remv:
                if v in var_name: 
                    var_name.remove(v)
                else:
                    print('Warning: ',v,' NOT in ',var_name)
        

        # copy dimensions
        for name, dimension in src.dimensions.items():
            len_dimension = len(dimension)
            
            dst.createDimension(name, len_dimension if not dimension.isunlimited() else None)
        #   
    
        # copy all file data except for matched-up variables, which instead averaged along specified dimension
        for name, variable in src.variables.items():

            # create variables, but will update its values later 
            dst.createVariable(name, variable.datatype, variable.dimensions)
            
            # copy variable attributes all at once via dictionary after created
            dst[name].setncatts(src[name].__dict__)

            #
            varvals = np.copy(src[name][...])
            if (dim_name in variable.dimensions) and (name in var_name):
                dim_indx = variable.dimensions.index(dim_name)
                tmp = np.mean(varvals, axis=dim_indx)
                if dim_indx==0:
                    varvals[:,...]=tmp
                elif dim_indx==1:
                    varvals[:,:,...]=tmp[:,newaxis]
                elif dim_indx==2:
                    varvals[:,:,:,...]=tmp[:,:,newaxis]
                elif dim_indx==3:
                    varvals[:,:,:,:,...]=tmp[:,:,:,newaxis]
                else:
                    print('Error:  dimension over 4 not supported')
                    exit(-1)
            #
                
            dst[name][...] = np.copy(varvals)
        
        #
            
    #            
    


#####################################################################
# module testing

#varname=['PCT_NATVEG','natpft','PCT_NAT_PFT']
#varname=['PCT_NAT_PFT','PCT_NATVEG','PCT_LAKE','PCT_WETLAND','PCT_GLACIER','PCT_URBAN','PCT_CROP']
#odata, odata_dims, odata_attr = getvar('surfdata.nc',varname)
#print('WAITING for a while to check ...')

#os.system('cp test01.nc test02.nc')
#vardata=odata
#putvar('test02.nc',varname,vardata)
#os.system('rm -rf test02.nc')

# duplicately expanding nc files along named dimension(s)
#dupexpand('surfdata.nc', 'surfdata_5x1pt.nc', dim_name=['lsmlat','lsmlon'], dim_len=[1,5])
#dupexpand('surfdata.pftdyn.nc', 'surfdata.pftdyn_5x1pt.nc', dim_name=['lsmlat','lsmlon'], dim_len=[1,5])
#dupexpand('surfdata.nc', 'surfdata_3x1pt.nc', dim_name='gridcell', dim_len=3)
#dupexpand('domain.nc', 'domain_3x1pt.nc', dim_name='ni', dim_multipler=3)
#dupexpand('domain.lnd.1x1pt_Oakharbor-GRID_navy.nc', 'domain.lnd.5x1pt_Oakharbor-GRID_navy.nc', dim_name=['nj','ni'], dim_len=[1,3])
#dupexpand('landuse.timeseries_1x1pt_kougarok-NGEE_simyr1850-2015_c181015m64.nc', 'test6x1.nc', ['lsmlat','lsmlon'], [1,6])

# merge 2 files along 1-only named-dimension
#mergefilesby1dim('f1.nc', 'f2.nc', 'fout.nc', 'n')
#files=['qbot1.nc','qbot2.nc','qbot3.nc','qbot4.nc','qbot5.nc']
#f0='tmp1.nc'
#fout='GSWP3_daymet4_QBOT_1980-2014_z01.nc'
#dimname='n'
#os.system('cp -rf '+ files[0]+' '+f0)
#for f in files[1:]:
#    mergefilesby1dim(f0, f, fout, dimname)
#    os.system('cp -rf '+ fout+' '+f0)
#os.system('rm -rf '+ f0)

#
#nco_extract('test02.nc', 'test03.nc', ['lsmlat','lsmlon'], [1,1], [1,1],'/usr/local/nco/bin')

#print('good!')

# convert nc to csv
#file = 'surfdata_51x63pt_kougarok-NGEE_simyr1850.nc'
#varname = ['LONGXY','LATIXY']
#nc2csv(file,varname)

# convert geotiff to nc
#file = 'land_cover.tif'
#bandinfos={'bands': ['tree','shrub','herbaceous','barren'], 
#           'units': ['%','%','%','%'], 
#           'long_name': ['land cover percentage of trees',
#                         'land cover percentage of shrubs',
#                         'land cover percentage of herbaceous',
#                         'land cover percentage of barren'
#                         ],
#           'standard_name': ['percentage of trees', 
#                             'percentage of shrubs',
#                             'percentage of herbaceous',
#                             'percentage of barren']
#           }
#file='ext100x100m_topo_seward.tif'
#bandinfos={'bands':['aspect','esl','slope']}
#file='ext100x100m_arcticpft_seward.tif'
#bandinfos={'bands':["arctic_lichen",
#                    "arctic_bryophyte",
#                    "arctic_forb",
#                    "arctic_graminoid",
#                    "arctic_evergreen_shrub",
#                    "arctic_deciduous_dwarf_shrub",
#                    "arctic_deciduous_low_shrub",
#                    "arctic_low_to_tall_willowbirch_shrub",
#                    "arctic_low_to_tall_alder_shrub",
#                    "arctic_needleleaf_tree",
#                    "arctic_broadleaf_tree",
#                    "non_vegetated",
#                    "water"]}
#geotiff2nc(file, bandinfos)

#varmeanby1dim - examples for surface data generation to check single variable
#varmeanby1dim('surfdata_pft.nc', 'surfdata_VIC.nc','gridcell', \
#              var_name='ALL', \
#              var_excl='LONGXY,LATIXY,AREA,Ds,Dsmax,Ws,binfl')
#              #var_excl='LONGXY,LATIXY,AREA,PCT_CROP,PCT_GLACIER,PCT_LAKE,PCT_NATVEG,PCT_NAT_PFT,PCT_URBAN,PCT_WETLAND,PFTDATA_MASK')
#              #var_excl='LONGXY,LATIXY,AREA,F0,FMAX,P3,ZWT0')
#              #var_excl='LONGXY,LATIXY,AREA,F0,FMAX,P3,ZWT0')
#              #var_excl='LONGXY,LATIXY,AREA,TOPO,SLOPE,STD_ELEV')
#              #var_excl='LONGXY,LATIXY,AREA,TOPO,SLOPE,STD_ELEV')
#              #var_excl='LONGXY,LATIXY,AREA,PCT_CLAY,PCT_SAND,ORGANIC,SOIL_COLOR,SOIL_ORDER')
            #
            #var_name='ALL', # all variables averaged, but excluding in 'var_excl'
            #var_excl='PCT_CROP,PCT_GLACIER,PCT_LAKE,PCT_NATVEG,PCT_NAT_PFT,PCT_URBAN,PCT_WETLAND,PFTDATA_MASK') # only those averaged

# overlayfiles - examples for joint 2 daymet nc files
#ncfilein1 = 'ELM_sim_for_NSIDC_yearly_dayssnowfree_1997-2019_1WestAlaskaArctic.nc'
#ncfilein2 = 'ELM_sim_for_NSIDC_yearly_dayssnowfree_1997-2019_2NorthAlaskaArctic.nc'
#ncfileout = 'merged.nc'
#overlayfiles(ncfilein1, ncfilein2, ncfileout)

# read 1 single variable from a nc file
#ncfile = 'surfdata.nc'
#varname = 'PCT_NATVEG'
#odata, odata_dims, odata_attr = getvar(ncfile, varname)

