import sys
from optparse import OptionParser
from netCDF4 import Dataset
import numpy as np
import h5py as h5
import glob
import math
from numpy import integer

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

#-------------------------------------------------------------------------------
def GetPflotranMeshData(filename, zscale):
    """This routine assumes that the mesh was created in layers, that is
    a 2D surface mesh extruded in the vertical dimension.
    """
    # read info
    mesh  = h5.File(filename,"r")
    xyz   = np.asarray(mesh["Domain/Vertices"])
    cells = np.asarray(mesh["Domain/Cells"],dtype=int);
    cells[:,1:] -= 1
    topo  = int((cells.shape[1]-1)/2)
    
    if zscale != 1:
        xyz[:,2] *= zscale
    
    # determine layering
    c1   = int(cells.shape[0]/2)
    test = cells[c1,1:(topo+1)]
    #c2   = np.where(np.prod((cells[:,-topo:]-test)==0,axis=1)==1)[0][0]
    c21 = cells[:,-topo:]
    c22 = c21-test
    c23 = np.prod(c22==0,axis=1)
    c24 = np.where(c23==1)
    c2 = c24[0][0]
    
    cpl  = c1-c2 # cells per layer

    # extract the top topology and remap vertex ids
    tcells = np.copy(cells[:cpl,:(topo+1)]); tcells[:,0] = topo
    inds   = np.unique(tcells[:,1:])
    cmap   = inds.searchsorted(range(cells.shape[0]))
    tcells[:,1:] = cmap[tcells[:,1:]]
    txyz   = np.copy(xyz[inds,:]); txyz[:,-1] = 0

    return xyz,cells,txyz,tcells # temporary hack

#-------------------------------------------------------------------------------
def Grid1dTo2d(grids1d_ij, grids1d_data,grids2d_shape):
    # Mapping 1D format of grided data to 2d-Grided data   

    ij2d=np.ravel_multi_index(np.array(grids1d_ij)-1, grids2d_shape)
    
    if(len(grids1d_data.shape)>0):
        shp2d = np.append(grids2d_shape,np.array(grids1d_data.shape[1:],int))
        shp1d = np.append(grids2d_shape[0]*grids2d_shape[1],np.array(grids1d_data.shape[1:],int))
        grids2d_data = np.full(shp2d, fill_value=np.nan, dtype=grids1d_data.dtype)
        grids2d_1d = np.reshape(grids2d_data,shp1d)
    
        grids2d_1d[ij2d,...] = grids1d_data
        grids2d_data = np.reshape(grids2d_1d, shp2d)
        
    else:
        grids2d_data = []
    #
    return grids2d_data

#-------------------------------------------------------------------------------
def Grid1dAggregating(subgrids1d_ij, subgrids1d_wt, subgrids1d_data):
    #aggregating data from pft/col/lun, in which XY is in 1D, to 1d-Grided data   
    
    #subgrids1d_ijxy: 2 sets of grid-index(i/j) for 2D grid (x/y)
    ix=subgrids1d_ij[0]
    jy=subgrids1d_ij[1]

    if(len(subgrids1d_data.shape)>1):
        wtsubdata = np.multiply(np.transpose(subgrids1d_data),subgrids1d_wt)
    else:
        wtsubdata = np.multiply(subgrids1d_data,subgrids1d_wt)
        
    ijxy=np.ravel_multi_index(np.array(subgrids1d_ij)-1,(max(ix),max(jy)))
    j=np.where(np.bincount(ijxy,subgrids1d_wt)>0.) #bincount() outputs counts for '0' as well, which should be excluded
    if (len(wtsubdata.shape)>1):
        grid1d_data = []
        for i in range(wtsubdata.shape[0]):
            if(i>0):
                grid1d_data = np.vstack((grid1d_data, np.bincount(ijxy, wtsubdata[i,:])[j]))
            else:
                grid1d_data = np.bincount(ijxy, wtsubdata[i,:])[j]
        grid1d_data = np.transpose(grid1d_data)# need to transpose back to the original
                
    elif(len(wtsubdata.shape)==1):
        grid1d_data = np.bincount(ijxy, wtsubdata)[j]
        
    else:
        grid1d_data = []
                
    #        
    return grid1d_data
#-------------------------------------------------------------------------------   
def WriteVTK(outfile,xyz,cells,shift,cell_variables={}):
    print ("Writing %s..." % outfile)
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
            
#-------------------------------------------------------------------------------
def WriteCLMDataToVTK(filehead,xyz,cells,clmdata, updown, nztrunc):
    if(len(clmdata["litr1c_vr"].shape)==3):
        nxy    = clmdata["litr1c_vr"].shape[0]*clmdata["litr1c_vr"].shape[1]
        nz     = min(nztrunc, clmdata["litr1c_vr"].shape[2])
    elif(len(clmdata["litr1c_vr"].shape)==2): 
        nxy    = clmdata["litr1c_vr"].shape[0]
        nz     = min(nztrunc, clmdata["litr1c_vr"].shape[1])
    
    ncells = cells.shape[0]

    curdate  = np.asscalar(clmdata["timemgr_rst_curr_ymd"])   # "time" in clm output *.r.*.nc IS in yyyymmdd
    if(type(curdate) is int):
        ymd= curdate
    else:
        ymd=curdate[0]
             
    tyr = int(math.floor(ymd/10000))
    tmd = int(ymd-tyr*10000)  

    # build another dictionary with references to this time steps variables
    var = {}
    for key in clmdata.keys():
        data = clmdata[key]
        shp  = clmdata[key].shape
        if len(shp) <= 1: continue   # skip because this is only a 1D or 0D array
        if len(shp)==3:
            data = data.reshape([shp[0]*shp[1],shp[2]])
            shp = data.shape

        if (key=='T_SOISNO' or key=='H2OSOI_LIQ' or key=='H2OSOI_ICE'):  # layers include 5 snow-layers
            data = data[...,5:]
            shp = data.shape
        
        if updown: #if VTK 3-D PF mesh is upside-down, then need to do so for CLM data from defined bottom layer
            if len(shp) == 2: 
                zdata=data[:,nz-1::-1]
                data = np.swapaxes(zdata, axis1=1, axis2=0)   # CLM: z first, then xy(gridcell); PF: xy first, then z 
            
        if np.prod(data.shape) == ncells:
            # shape of data matches the number of cells in the VTK mesh
            var[key] = data.reshape((cells.shape[0]))
        elif np.prod(data.shape) % ncells == 0:
            # shape of data matches a multiple of the number of cells, but not nz times
            if np.prod(data.shape) / ncells == nz: continue
            sub = data.reshape((nxy,-1))
            for j in range(sub.shape[1]):
                if np.allclose(sub[:,j],1e36): continue
                var["%s_%s" % (key,str(j))] = sub[:,j]

    if len(var.keys())==0: return # if no variables, no need to plot
    WriteVTK("%s-yyyymmdd-%s%s.vtk" % (filehead,str.zfill(str(tyr),4),str.zfill(str(tmd),4)),
                 xyz,cells,False, cell_variables=var)
  
#--------------------------------------------------------------------------------------
parser = OptionParser()

parser.add_option("--clmout_dir", dest="clm_odir", default="./", \
                  help="clm output directory (default = ./, i.e., under current directory)")
parser.add_option("--clmfile_head", dest="clm_filehead", default="", \
                  help = "clm output file name header, usually the portion before *.clm2.r.*.nc")
parser.add_option("--pflotran_meshfile", dest="pf_meshfile", default="", \
                  help="pflotran mesh file including its path (default = ./*.h5)")
parser.add_option("--nzmax", dest="nztrunc", default="", \
                  help="max. z depth numbers (default = 10)")
parser.add_option("--zscale", dest="zscale", default="1", \
                  help="z-axis scaling (default = 1)")
parser.add_option("--adspinup", dest="adspinup", action="store_true", default=False, \
                  help="whether results of an ad_spinup run (default = False)")
(options, args) = parser.parse_args()

# parse the PFLOTRAN mesh file and extract surface mesh
if (options.pf_meshfile==""):
    print('MUST provide PFLOTRAN mesh file name including its path! Sorry but Exit -- ')
    sys.exit()
else:
    xyz,cells,sxyz,scells = GetPflotranMeshData(options.pf_meshfile, int(options.zscale))

# setup which CLM files to parse
if (options.clm_odir=="./"):
    print('Warning: CLM output data directory IS the current ! ')

if (options.clm_filehead==""):
    print('MUST provide CLM output data file Header, which usually is the portion usually the portion before "*.clm2.r.*.nc"! ')
    sys.exit()
else:
    clmhead = options.clm_odir+'/'+options.clm_filehead+'.clm2'
    
if (options.nztrunc==""):
    nzmax = 15
else:
    nzmax = int(options.nztrunc)
    
#--------------------------------------------------------------------------------------

include = ['r']
allfile = glob.glob("%s*.nc" % clmhead)
chunks  = []
for filename in allfile:
    filename = filename.split(".")
    if chunks.count(filename[-2]) == 0: chunks.append(filename[-2])

# which variables would you like to keep?
keep    = ['litr1c_vr', 'litr2c_vr', 'litr3c_vr',
           'soil1c_vr', 'soil2c_vr', 'soil3c_vr', 'soil4c_vr',
           'smin_no3_vr', 'smin_nh4_vr']
keep.append('totsomc_vr')
keep.append('totlitr_vr')
keep.append('T_SOISNO')
keep.append('H2OSOI_LIQ')
keep.append('H2OSOI_ICE')

ind_totsomc = 9  # the index in 'keep' array
ind_totlitr = 10
totsomc = ['soil1c_vr', 'soil2c_vr', 'soil3c_vr', 'soil4c_vr']
ad_factor =[1,1,10,100]
totlitr = ['litr1c_vr', 'litr2c_vr', 'litr3c_vr']

#
keep.append("timemgr_rst_curr_ymd")    # need this for clm timing

# process CLM files
for chunk in chunks:

    # for each chunk of time slices, append variables to a single dictionary
    variables = {}
    variables_dims = {}
    for inc in include:
        
        # try reading all files
        filename = "%s.%s.%s.nc" % (clmhead,inc,chunk)
        try:
            f = Dataset(filename,'r')
        except:
            continue

        #grids mapping
        if 'grid1d_lon' in f.variables.keys() and 'grid1d_lat' in f.variables.keys():
            grid1d_xy=[np.asarray(f.variables['grid1d_lon']),np.asarray(f.variables['grid1d_lat'])]
        else:
            grid1d_xy=[]
        if 'grid1d_ixy' in f.variables.keys() and 'grid1d_jxy' in f.variables.keys():
            grid1d_ij=[np.asarray(f.variables['grid1d_ixy']),np.asarray(f.variables['grid1d_jxy'])]
        else:
            grid1d_ij=[]
        #landunits mapping to grids
        if 'land1d_wtxy' in f.variables.keys():
            land1d_wt=np.asarray(f.variables['land1d_wtxy'])
        else:
            land1d_wt=[]
        if 'land1d_ixy' in f.variables.keys() and 'land1d_jxy' in f.variables.keys():
            land1d_ij=[np.asarray(f.variables['land1d_ixy']),np.asarray(f.variables['land1d_jxy'])]
        else:
            land1d_ij=[]
        #columns mapping to grids
        if 'cols1d_wtxy' in f.variables.keys():
            cols1d_wt=np.asarray(f.variables['cols1d_wtxy'])
        else:
            cols1d_wt=[]
        if 'cols1d_ixy' in f.variables.keys() and 'cols1d_jxy' in f.variables.keys():
            cols1d_ij=[np.asarray(f.variables['cols1d_ixy']),np.asarray(f.variables['cols1d_jxy'])]
        else:
            cols1d_ij=[]
        #pfts mapping to grids
        if 'pfts1d_wtxy' in f.variables.keys():
            pfts1d_wt=np.asarray(f.variables['pfts1d_wtxy'])
        else:
            pfts1d_wt=[]
        if 'pfts1d_ixy' in f.variables.keys() and 'pfts1d_jxy' in f.variables.keys():
            pfts1d_ij=[np.asarray(f.variables['pfts1d_ixy']),np.asarray(f.variables['pfts1d_jxy'])]
        else:
            pfts1d_ij=[]
            
            
        
        # If key is on the keep list, uniquely add to a dictionary
        for key in f.variables.keys():
            if key not in keep: continue
            if key not in variables.keys(): 
                variables[key] = np.asarray(f.variables[key])
                variables_dims[key] = f.variables[key].dimensions
                
                #aggregating data into gridcell-level
                gdata=[]
                if('landunit' in variables_dims[key] and len(land1d_ij)>0):
                    gdata = Grid1dAggregating(land1d_ij, land1d_wt, variables[key])
                    variables_dims[key]=np.core.defchararray.replace(variables_dims[key],'landunit','gridcell')
                if('column' in variables_dims[key] and len(cols1d_ij)>0):
                    gdata = Grid1dAggregating(cols1d_ij, cols1d_wt, variables[key])
                    variables_dims[key]=np.core.defchararray.replace(variables_dims[key],'column','gridcell')
                if('pft' in variables_dims[key] and len(pfts1d_ij)>0):
                    gdata = Grid1dAggregating(pfts1d_ij, pfts1d_wt, variables[key])
                    variables_dims[key]=np.core.defchararray.replace(variables_dims[key],'pft','gridcell')
                if(len(gdata)>0): variables[key] = gdata
                
                
                #1d --> 2d
                if('gridcell' in variables_dims[key]):
                    variables[key] = Grid1dTo2d(grid1d_ij, variables[key], [720, 360])
                    variables_dims[key] = np.core.defchararray.replace(variables_dims[key],'gridcell',['lon','lat'])
                
                
                #ad-spinup run
                if (options.adspinup and key in totsomc):
                    key_index = totsomc.index(key)
                    variables[key] = variables[key]*ad_factor[key_index]
                
                # summing up of total liter C
                if key in totlitr:
                    if key==totlitr[0]:
                        variables[keep[ind_totlitr]] = variables[key]
                    else:
                        variables[keep[ind_totlitr]] = variables[keep[ind_totlitr]] + variables[key]

                # summing up of total som C
                if key in totsomc:
                    if key==totsomc[0]:
                        variables[keep[ind_totsomc]] = variables[key]
                    else:
                        variables[keep[ind_totsomc]] = variables[keep[ind_totsomc]] + variables[key]

    WriteCLMDataToVTK("surf_"+options.clm_filehead,sxyz,scells, variables, True, nzmax)
    WriteCLMDataToVTK("vol_"+options.clm_filehead , xyz, cells, variables, True, nzmax)

