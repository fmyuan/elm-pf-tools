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

def WriteVTK(outfile,xyz,cells,shift,cell_variables={}):
    print "Writing %s..." % outfile
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

def WriteCLMDataToVTK(filehead,xyz,cells,clmdata, updown, nztrunc):
    nxy    = clmdata["litr1c_vr"].shape[0]
    nz     = min(nztrunc, clmdata["litr1c_vr"].shape[1])
    ncells = cells.shape[0]

    curdate  = np.asscalar(clmdata["timemgr_rst_curr_ymd"])   # "time" in clm output *.r.*.nc IS in yyyymmdd
    ymd = curdate[0]
    tyr = int(math.floor(ymd/10000))
    tmd = int(ymd-tyr*10000)  

    # build another dictionary with references to this time steps variables
    var = {}
    for key in clmdata.keys():
        data = clmdata[key]
        shp  = clmdata[key].shape
        if len(shp) == 1: continue   # skip because this is only a 1D array

        if (key=='T_SOISNO' or key=='H2OSOI_LIQ' or key=='H2OSOI_ICE'):  # layers include 5 snow-layers
            data = clmdata[key][:,5:]
            shp = data.shape
        
        if updown: #if VTK 3-D PF mesh is upside-down, then need to do so for CLM data from defined bottom layer
            if len(shp) == 2: 
                zdata=data[:,nz-1::-1]
                data = np.swapaxes(zdata, axis1=1, axis2=0)   # CLM: z first, then xy(gridcell); PF: xy first, then z 
            if len(shp) == 3: 
                zdata=data[:,:,nz-1:-1]
                data = np.swapaxes(zdata, axis1=2, axis2=1, axis3=0)          
            
        if np.prod(data.shape) == ncells:
            # shape of data matches the number of cells in the VTK mesh
            var[key] = data.reshape((cells.shape[0]))
        elif np.prod(data.shape) % ncells == 0:
            # shape of data matches a multiple of the number of cells, but not nz times
            if np.prod(data.shape) / ncells == nz: continue
            pfts1d_wtgcell = clmdata['pfts1d_wtxy'].reshape((nxy,-1))
            sub = data.reshape((nxy,-1))
            for j in range(sub.shape[1]):
                if np.allclose(sub[:,j],1e36): continue
                var["%s_%s" % (key,pfts[j])] = sub[:,j]
                var["%s_total" % key] = np.apply_along_axis(np.sum,1,pfts1d_wtgcell*sub)

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
    nzmax = 10
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

# I need to add a few to get some info for plotting
keep.append("timemgr_rst_curr_ymd")    # need this for clm timing
keep.append("pfts1d_wtxy")

# process CLM files
for chunk in chunks:

    # for each chunk of time slices, append variables to a single dictionary
    variables = {}
    for inc in include:
        
        # try reading all files
        filename = "%s.%s.%s.nc" % (clmhead,inc,chunk)
        try:
            f = Dataset(filename,'r')
        except:
            continue

        # If key is on the keep list, uniquely add to a dictionary
        for key in f.variables.keys():
            if key not in keep: continue
            if key not in variables.keys(): 
                variables[key] = np.asarray(f.variables[key])
                
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

