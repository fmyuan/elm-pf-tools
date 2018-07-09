#Python utilities for operating a netcdf file using available tools, python netcdf4, nco
import os
import numpy as np
from netCDF4 import Dataset
from netCDF4 import Variable
from array import array
ncopath = '/usr/local/nco/bin'

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
def putvar(ncfile, varname, varvals):
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
            if isinstance(varvals, dict): #multiple dataset for corresponding varname           
                for v in varname:
                    val=np.asarray(varvals[v])
                    if v in f.variables.keys():
                        f.variables[v][...]=np.copy(val)
            else:                        #exactly 1 dataset in np.array for 1 varname
                val=np.asarray(varvals)
                v=varname[0]
                if v in f.variables.keys():
                    f.variables[v][...]=np.copy(val)
                
                
            f.sync()
            f.close()
            
            f = Dataset(ncfile,'r')
            return 0

    except:
        return -2

#----------------------------------------------------------------------------             
# duplicately expand all variables along named dimension by multiple times
def dupexpand(ncfilein, ncfileout,dim_name,dim_multipler):
    with Dataset(ncfilein) as src, Dataset(ncfileout, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
   
        if type(dim_name)==str: 
            dim_name=[dim_name]
            dim_multipler=[dim_multipler]

        # copy or multiply dimensions
        for name, dimension in src.dimensions.items():
            if name in dim_name:
                indx=dim_name.index(name)
                len_dimension = len(dimension)*dim_multipler[indx]
            else:
                len_dimension = len(dimension)
            
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
                if dim in variable.dimensions and dim_multipler[dim_name.index(dim)]>1:
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
    print('done!')

#----------------------------------------------------------------------------             
#Using NCO to extract all variables along named dimensions
def nco_extract(ncfilein, ncfileout,dim_names,dim_starts, dim_lens, ncksdir=ncopath):   
    
    if(not os.path.isfile(ncfilein)): 
        print('Error: invalid input NC file: '+ncfilein)
        sys.exit()

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
    if(not ncksdir.endswith('/')):
        ncksdir = ncksdir.strip()+'/'
    os.system(ncksdir+'ncks -a -O '+ dim_string+ \
          ' '+ncfilein+' '+ncfileout)

    #
    
#----------------------------------------------------------------------------             

#####################################################################
# module testing
#varname=['PCT_NATVEG','natpft','PCT_NAT_PFT']
#odata, odata_dims, odata_attr = getvar('test01.nc',varname)

#os.system('cp test01.nc test02.nc')
#vardata=odata
#putvar('test02.nc',varname,vardata)
#os.system('rm -rf test02.nc')

#
#dupexpand('test01.nc', 'test02.nc', ['lsmlat','lsmlon'], [2,3])

#
#nco_extract('test02.nc', 'test03.nc', ['lsmlat','lsmlon'], [1,1], [1,1],'/usr/local/nco/bin')

#print('good!')

  
