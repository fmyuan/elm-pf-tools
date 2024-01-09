
# data_partition module for batch processing
# based on array_split and function definition

import os,sys 
import netCDF4 as nc
import numpy as np
from time import process_time
from datetime import datetime

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m%d%Y')

def forcing_save_1dNA(input_path, file, var_name, period, tlength, output_path):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    source_file = input_path + '/'+ file
    src = nc.Dataset(source_file, 'r', format='NETCDF4')

    total_rows = src.dimensions['x'].size
    total_cols = src.dimensions['y'].size
    total_times= src.dimensions['time'].size

    #print('total timesteps is :' + str(total_timesteps))
    if tlength == -1:
        tlength = total_times

    data = src[var_name][0:tlength, :, :] # read (time, y, x) format
    time = src['time'][0:tlength]
    latxy= src['lat'][:,:]
    lonxy= src['lon'][:,:]

    #create land mask
    mask = data[0]    # data is in (time, Y, X) format
    mask = np.where(~np.isnan(mask), 1, np.nan)

    #create gridIDs
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
    grid_ids = grid_ids.reshape(total_cols,total_rows)

    # create an flattened list of land gridID and reduce the size of gridIDs array
    grid_ids = np.multiply(mask,grid_ids)
    grid_ids = grid_ids[~np.isnan(grid_ids)]

    # extract the data over land gridcells
    landcells = len(grid_ids)
    data = np.multiply(mask, data)
    data = data[~np.isnan(data)]
    data = np.reshape(data,(tlength,landcells))

    latxy = np.multiply(mask,latxy)
    latxy = latxy[~np.isnan(latxy)]
    lonxy = np.multiply(mask,lonxy)
    lonxy = lonxy[~np.isnan(lonxy)]

    # convert local grid_id_lists into an array
    grid_id_arr = np.array(grid_ids)

    #data_arr = np.array(FSDS_list[i])
    data_arr = np.array(data)
    lonxy_arr= np.array(lonxy)
    latxy_arr= np.array(latxy)
    
    
    dst_name = output_path + '/clmforc.Daymet4.1km.1d.' + var_name + '.' + period +'.nc'

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    dst = nc.Dataset(dst_name, 'w', format='NETCDF4')
    dst.title = var_name + '('+period+') created from '+ input_path +' on ' +formatted_date

    # create the gridIDs, lon, and lat variable
    x = dst.createDimension('ni', grid_id_arr.size)
    x = dst.createDimension('nj', 1)
    x = dst.createDimension('time', tlength)

    w_nc_var = dst.createVariable('gridID', np.int32, ('nj','ni'))
    #  set variable attributes
    w_nc_var.long_name = "gridId in the NA domain" ;
    w_nc_var.decription = "Covers all land and ocean gridcells, with #0 at the upper left corner of the domain" ;
    dst.variables['gridID'][...] = grid_id_arr.reshape(grid_id_arr.size,1)

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():
        if (name == var_name):
            w_nc_var = dst.createVariable(var_name, np.float32, ('time', 'nj', 'ni'))
            dst.variables[var_name][...] =data_arr.reshape(tlength,grid_id_arr.size,1)
            for attr_name in variable.ncattrs():
                dst[name].setncattr(attr_name, variable.getncattr(attr_name))
        
        if (name == 'time'):
            dvname = 'time'
            w_nc_var = dst.createVariable(dvname, np.float32, ('time'))
            dst.variables[dvname][...] = time
            for attr_name in variable.ncattrs():
                dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

        if (name == 'lat'):
            dvname = 'LATIXY'
            w_nc_var = dst.createVariable(dvname, np.float64, ('nj','ni'))
            dst.variables[dvname][...] = latxy_arr
            for attr_name in variable.ncattrs():
                dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

        if (name == 'lon'):
            dvname = 'LONGXY'
            w_nc_var = dst.createVariable(dvname, np.float64, ('nj','ni'))
            dst.variables[dvname][...] = lonxy_arr
            for attr_name in variable.ncattrs():
                dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

    src.close()  # close the source file 
    dst.close()  # close the new file        
    
def get_files(input_path):
    print(input_path)
    files = os.listdir(input_path) 

    files.sort() 

    file_no =0

    files_nc = [f for f in files if (f[-2:] == 'nc')] 
    print("total " + str(len(files_nc)) + " files need to be processed")
    return files_nc

def main():
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    time = int(args[2])
    
    files_nc = get_files(input_path)

    for f in files_nc: 
        if (not f.startswith('clmforc')): continue
        var_name = f[20:-11]
        period = f[-10:-3]
        print('processing '+ var_name + '(' + period + ') in the file ' + f )
        start = process_time() 
        forcing_save_1dNA(input_path, f, var_name, period, time, output_path)
        end = process_time()
        print("Generating 1D forcing data takes {}".format(end-start))

if __name__ == '__main__':
    main()
    
