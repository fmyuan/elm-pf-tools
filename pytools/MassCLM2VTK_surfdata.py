#!/usr/bin/env python

import sys
from optparse import OptionParser
from netCDF4 import Dataset
import numpy as np
import h5py as h5
import glob
import math

pfts=["not_vegetated",
      "arctic_lichen",
      "arctic_bryophyte",
      "needleleaf_evergreen_temperate_tree",
      "needleleaf_evergreen_boreal_tree",
      "needleleaf_deciduous_boreal_tree",
      "broadleaf_evergreen_tropical_tree",
      "broadleaf_evergreen_temperate_tree",
      "broadleaf_deciduous_tropical_tree",
      "broadleaf_deciduous_temperate_tree",
      "broadleaf_deciduous_boreal_tree",
      "broadleaf_evergreen_shrub",
      "broadleaf_deciduous_temperate_shrub",
      "broadleaf_deciduous_boreal_shrub",
      "evergreen_arctic_shrub",
      "deciduous_arctic_shrub",
      "c3_arctic_sedge",
      "c3_arctic_forb",
      "c3_arctic_grass",
      "c3_non-arctic_grass",
      "c4_grass","c3_crop",
      "c3_irrigated",
      "corn",
      "irrigated_corn",
      "spring_temperate_cereal",
      "irrigated_spring_temperate_cereal",
      "winter_temperate_cereal",
      "irrigated_winter_temperate_cereal",
      "soybean",
      "irrigated_soybean"]


def GetPflotranMeshData(filename, zscale, nztrunc):
    """This routine assumes that the mesh was created in layers, that is
    a 2D surface mesh extruded in the vertical dimension.
    """
    # read info
    mesh  = h5.File(filename,"r")
    xyz   = np.asarray(mesh["Domain/Vertices"]) # ordered vertices coordinates in (x,y,z), indexed from 1
    cells = np.asarray(mesh["Domain/Cells"],dtype=int);# ordered cells vertices: no. of vertices, vertices node labeling no. 
    cells[:,1:] -= 1
        
    topo  = int((cells.shape[1]-1)/2) # 
    
    
    # usually soil z ranging within meters, but XY could be in much larger scales
    # in this case, user can scale Z to have a good visualization
    if zscale != 1:
        xyz[:,2] *= zscale
    
    # determine layering
    c1   = int(cells.shape[0]/2)
    toppest = cells[c1,1:(topo+1)]
    c21 = cells[:,-topo:]
    c22 = c21-toppest
    c23 = np.prod(c22==0,axis=1)
    c24 = np.where(c23==1)
    c2 = c24[0][0]
    
    cpl  = c1-c2 # cells per layer

    # extract the top topology and remap vertex ids
    scells = np.copy(cells[:cpl,:(topo+1)]); 
    scells[:,0] = topo
    inds   = np.unique(scells[:,1:])
    cmap   = inds.searchsorted(range(cells.shape[0]))
    scells[:,1:] = cmap[scells[:,1:]]
    sxyz   = np.copy(xyz[inds,:]); sxyz[:,-1] = 0

    # may need to truncate z layers of cells, but keep the xyz not modified
    nz = cells.shape[0]/cpl
    if(nz!=int(nz)):
        print("soil layers NOT an integer!")
        exit        
    elif nz > nztrunc:
        z_trunc = np.cumsum(np.ones((nz-nztrunc)*cpl,dtype=int))-1
        cells = np.delete(cells, z_trunc, axis=0)

    return xyz,cells,sxyz,scells

def WriteCLMDataToVTK(filehead,xyz,cells,clmdata, updown, nztrunc):
    if(len(clmdata["ORGANIC"].shape)==3):  # data in [nz, nj, ni]
        nx= clmdata["ORGANIC"].shape[2]
        ny= clmdata["ORGANIC"].shape[1]
    else:
        nx= 1
        ny= clmdata["ORGANIC"].shape[1]  # data in [nz, gridcell]        
    nxy = nx*ny
    nz  = clmdata["ORGANIC"].shape[0]
    nxyz= nxy*nz    
    npft   = clmdata["PCT_NAT_PFT"].shape[0]   #data in [npft, nj,ni], or [npft, gridcell]
    
    ncells = cells.shape[0]

    # build another dictionary with references to this time steps variables
    var = {}
    for key in clmdata.keys():
        data = clmdata[key]
        xdim = len(data.shape)
            
        shp  = data.shape
        if len(shp) == 0: continue   # skip because this is only a 0D array
                   
        if np.prod(shp) == ncells:
            # shape of data matches the number of cells in the VTK mesh
            
            if updown and (np.prod(shp) %nxyz ==0): #if VTK 3-D PF mesh is upside-down, then need to do so for CLM data from defined bottom layer
                if len(shp) == 2: 
                    zdata= data[nz-1::-1,:]
                    data = zdata #np.swapaxes(zdata, axis1=1, axis2=0)   # CLM: z first, then xy(gridcell); PF: xy first, then z 
                elif len(shp) == 3: 
                    zdata= data[nz-1::-1,:,:]
                    data = zdata #np.swapaxes(zdata, axis1=2, axis2=1, axis3=0)          

            var[key] = data.reshape((cells.shape[0]))
            
        elif np.prod(shp) % ncells == 0:
            # shape of data matches a multiple of the number of cells, but not nz times
            if np.prod(shp) %npft == 0 and key == "PCT_NAT_PFT":
                sub = data.reshape((-1, nxy))
                for j in range(sub.shape[0]):
                    var["%s_%s" % (key,pfts[j])] = sub[j,:]
            else:
                if(np.prod(shp)%nz==0):
                    continue
                else:
                    var[key] = data

    #write var to file
    if len(var.keys())==0: return # if no variables, no need to plot
    WriteVTK("%s-clmsurfdata.vtk" % filehead,xyz,cells,False, cell_variables=var)

def WriteVTK(outfile,xyz,cells,shift,cell_variables={}):
    #print "Writing %s..." % outfile
    f = open(outfile,'w')
    f.write('# vtk DataFile Version 2.0\n')
    f.write('converted\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS %d float\n' % (xyz.shape[0]))
    sft = np.zeros(3)
    if shift:
        sft[0] = xyz[:,0].min()
        sft[1] = xyz[:,1].min()   
           
    (xyz-sft).tofile(f,sep=' ',format='%s')
        
    f.write('\nCELLS %d %d\n' % (cells.shape[0],cells.size))   
    cells.tofile(f,sep=' ',format='%s')
    
    f.write('\nCELL_TYPES %d\n' % (cells.shape[0]))
    ((cells[:,0] == 3)*5+
     (cells[:,0] == 4)*9+
     (cells[:,0] == 6)*13+
     (cells[:,0] == 8)*12).tofile(f,sep=' ',format='%s')
    if cell_variables.keys(): f.write("\nCELL_DATA %d" % (cells.shape[0]))
    for key in cell_variables.keys():
        f.write("\nSCALARS %s float 1\nLOOKUP_TABLE default\n" % key)
        if (cell_variables[key].shape == 2):
            var_rot=np.rot90(cell_variables[key],1)
            var_rot.tofile(f,sep=' ',format="%s")           
        else:
            cell_variables[key].tofile(f,sep=' ',format="%s")




#--------------------------------------------------------------------------------------
parser = OptionParser()

parser.add_option("--clmsurf_dir", dest="clmsurf_dir", default="./", \
                  help="clm surface data directory (default = ./, i.e., under current directory)")
parser.add_option("--clmsurf_file", dest="clmsurf_file", default="", \
                  help = "clm surface file name in *.nc")
parser.add_option("--pflotran_meshfile", dest="pf_meshfile", default="", \
                  help="pflotran mesh file including its path (default = ./*.h5)")
parser.add_option("--nzmax", dest="nztrunc", default="", \
                  help="max. z depth numbers (default = 10)")
parser.add_option("--zscale", dest="zscale", default="1", \
                  help="z-axis scaling (default = 1)")
(options, args) = parser.parse_args()



if (options.nztrunc==""):
    nzmax = 15
else:
    nzmax = int(options.nztrunc)

# parse the PFLOTRAN mesh file and extract surface mesh
if (options.pf_meshfile==""):
    print('MUST provide PFLOTRAN mesh file name including its path! Sorry but Exit -- ')
    sys.exit()
else:
    xyz,cells,sxyz,scells = GetPflotranMeshData(options.pf_meshfile, int(options.zscale), int(nzmax))

# setup which CLM files to parse
if (options.clmsurf_dir=="./"):
    print('Warning: CLM surface data directory IS the current ! ')

if (options.clmsurf_file==""):
    print('MUST provide CLM surface data file name, usually as "surfdata_*x*pt_US-Brw_simyr1850.nc"! ')
    sys.exit()
else:
    clmsurf = options.clmsurf_dir+'/'+options.clmsurf_file
    

#--------------------------------------------------------------------------------------

# which variables would you like to keep?
keep    = ['LATIXY', 'LONGXY', 'TOPO', 'AREA', 'ORGANIC', 'PCT_NAT_PFT']

variables = {}

try:
    f = Dataset(clmsurf,'r')
except:
    print('CLM surface data file NOT exists! ')
    sys.exit()

# If key is on the keep list, uniquely add to a dictionary
for key in f.variables.keys():
    if key not in keep: continue
    if key not in variables.keys(): 
        variables[key] = np.asarray(f.variables[key])

WriteCLMDataToVTK("surface",sxyz,scells, variables, True, nzmax)
WriteCLMDataToVTK("volume" , xyz, cells, variables, True, nzmax)


