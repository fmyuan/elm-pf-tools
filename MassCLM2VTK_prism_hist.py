import sys
from optparse import OptionParser
from netCDF4 import Dataset
import numpy as np
import h5py as h5
import glob
import math
from numpy.f2py.crackfortran import endifpattern

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

#----------------------------------------------------------------------------------------------------------
# local functions

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
    print("Writing %s..." %outfile)
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

def WriteCLMDataToVTK(filehead,xyz,cells,clmdata, clmdata_dims, updown, nztrunc):
    steps  = np.asarray(clmdata["nstep"],dtype=int)
    nstep  = steps.shape[0]
    nx     = clmdata["topo"].shape[0]
    if (len(clmdata["topo"].shape)==1):
        ny     = 1
    elif(len(clmdata["topo"].shape)==2):
        ny     = clmdata["topo"].shape[1]
    
    
    nz1    = clmdata["levgrnd"].shape[0]
    nz2    = clmdata["levdcmp"].shape[0]
    nz     = min(nztrunc, min(nz1,nz2))
    ncells = cells.shape[0]

    times  = np.asarray(clmdata["time"],dtype=int)   # "time" in clm output *.nc IS in days from starting time of the year

    # loop over all time steps
    for i in range(nstep):
        # build another dictionary with references to this time steps variables
        var = {}
        tyr = math.floor(times[i]/365)+starttime          # converting time unit from days to year-doy
        tdoy= int(times[i]-(tyr-starttime)*365)           # converting time unit from days to year-doy
        for key in clmdata.keys():
            data = clmdata[key]
            shp  = clmdata[key].shape           
            dim_list = clmdata_dims[key]
            
            if ("time" not in dim_list): continue      # skip if NOT time-series data array

            # do having vertical dimension?                
            if ("levgrnd" in dim_list or "levdcmp" in dim_list):
                if ("levgrnd" in dim_list):
                    ind_zdim = dim_list.index("levgrnd")
                elif ("levdcmp" in dim_list):
                    ind_zdim = dim_list.index("levdcmp")
                
                # truncating data along Z-dimension if needed
                if (nz<shp[ind_zdim]):
                    if(ind_zdim==0): 
                        data = clmdata[key][abs(nz-shp[ind_zdim]):,]
                    elif(ind_zdim==1):
                        data = clmdata[key][:,abs(nz-shp[ind_zdim]):,]
                    elif(ind_zdim==2):
                        data = clmdata[key][:,:,abs(nz-shp[ind_zdim]):,]
                    elif(ind_zdim==3):
                        data = clmdata[key][:,:,:,abs(nz-shp[ind_zdim]):,]
                    else:
                        print('Warning: Z-dimension beyond 4-D: please check ! ')
                        sys.exit()
                        
                    shp = data.shape
       
                     
                #if VTK 3-D PF mesh is upside-down, then need to do so for CLM data from defined bottom layer
                if updown: 
                    if(ind_zdim==0): 
                        data = data[nz-1::-1,]
                    elif(ind_zdim==1): 
                        data = data[:,nz-1::-1,]
                    elif(ind_zdim==2):
                        data = data[:,:,nz-1::-1,]
                    elif(ind_zdim==3):
                        data = data[:,:,:,nz-1::-1,]
                    else:
                        print('Warning: Z-dimension beyond 4-D: please check ! ')
                        sys.exit()
            
            if np.prod(data[i,...].shape) == ncells:
                # shape of data matches the number of cells in the VTK mesh
                var[key] = data[i,...].reshape((cells.shape[0]))
            elif np.prod(data[i,...].shape) % ncells == 0:
                # shape of data matches a multiple of the number of cells, but not nz times
                if np.prod(data[i,...].shape) / ncells == nz: continue
                if(clmdata['pfts1d_wtgcell'].exists()): continue
                pfts1d_wtgcell = clmdata['pfts1d_wtgcell'].reshape((nxy,-1))
                sub = data[i,...].reshape((nxy,-1))
                for j in range(sub.shape[1]):
                    if np.allclose(sub[:,j],1e36): continue
                    var["%s_%s" % (key,pfts[j])] = sub[:,j]
                var["%s_total" % key] = np.apply_along_axis(np.sum,1,pfts1d_wtgcell*sub)
        if len(var.keys())==0: return # if no variables, no need to plot
        WriteVTK("%s-yyyydoy-%s%s.vtk" % (filehead,str.zfill(str(int(tyr)),4),str.zfill(str(int(tdoy)),3)),
                 xyz,cells,False,cell_variables=var)


#--------------------------------------------------------------------------------------
parser = OptionParser()

parser.add_option("--clmout_dir", dest="clm_odir", default="./", \
                  help="clm output directory (default = ./, i.e., under current directory)")
parser.add_option("--clmfile_head", dest="clm_filehead", default="", \
                  help = "clm output file name header, usually the portion before *.clm2.h[0-5].*.nc")
parser.add_option("--pflotran_meshfile", dest="pf_meshfile", default="", \
                  help="pflotran mesh file including its path (default = ./*.h5)")
parser.add_option("--nzmax", dest="nztrunc", default="", \
                  help="max. z depth numbers (default = 10)")
parser.add_option("--zscale", dest="zscale", default="1", \
                  help="z-axis scaling (default = 1)")
parser.add_option("--adspinup", dest="adspinup", action="store_true", default=False, \
                  help="whether results of an ad_spinup run (default = False)")
parser.add_option("--starting_time",dest="starttime", default="1", \
                  help="clm run starting time (default=1 for spinup; for transient it should be 1850")
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
    print('MUST provide CLM output data file Header, which usually is the portion usually the portion before "*.clm2.[h0-h2].*.nc"! ')
    sys.exit()
else:
    clmhead = options.clm_odir+'/'+options.clm_filehead+'.clm2'

starttime = int(options.starttime)

if (options.nztrunc==""):
    nzmax = 10
else:
    nzmax = int(options.nztrunc)

#--------------------------------------------------------------------------------------

include = ['h0','h1','h2','h3','h4','h5']
allfile = glob.glob("%s*.nc" % clmhead)
chunks  = []
for filename in allfile:
    filename = filename.split(".")
    if chunks.count(filename[-2]) == 0: chunks.append(filename[-2])

# which variables would you like to keep?
keep    = ['TOTSOMC','TOTLITC','TOTVEGC','TLAI','H2OSOI','SOILLIQ','SOILICE','TSOI']   # 8

keep.append('litr1c_vr')
keep.append('litr2c_vr')
keep.append('litr3c_vr')
keep.append('soil1c_vr')
keep.append('soil2c_vr')
keep.append('soil3c_vr')
keep.append('soil4c_vr')
keep.append('litr1n_vr')
keep.append('litr2n_vr')
keep.append('litr3n_vr')
keep.append('soil1n_vr')
keep.append('soil2n_vr')
keep.append('soil3n_vr')
keep.append('soil4n_vr')
keep.append('smin_no3_vr')
keep.append('smin_nh4_vr')   # 24
keep.append('totsomc_vr')    # 25
keep.append('totlitc_vr')
keep.append('totsomn_vr')
keep.append('totlitn_vr')

ind_totsomc = 25  # the index in 'keep' array
ind_totlitc = 26
totsomc = ['soil1c_vr', 'soil2c_vr', 'soil3c_vr', 'soil4c_vr']
ad_factor =[1,1,10,100]
totlitc = ['litr1c_vr', 'litr2c_vr', 'litr3c_vr']

# I need to add a few to get some info for plotting
keep.append("nstep")   # need this for time step numbers
keep.append("time")    # need this for clm timing
keep.append("topo")    # need this for number of surface cells
keep.append("levgrnd") # need this for the number of layers 
keep.append("levdcmp") # need this for the number of layers 
keep.append("pfts1d_wtgcell")

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

        # If key is on the keep list, uniquely add to a dictionary
        for key in f.variables.keys():
            if key not in keep: continue
            if key not in variables.keys(): 
                variables[key] = np.asarray(f.variables[key])
                variables_dims[key] = f.variables[key].dimensions

            #ad-spinup run
            if (options.adspinup and key in totsomc):
                key_index = totsomc.index(key)
                variables[key] = variables[key]*ad_factor[key_index]
            
            # summing up of total liter C
            if key in totlitc:
                if key==totlitc[0]:
                    variables[keep[ind_totlitc]] = variables[key]
                    variables_dims[keep[ind_totlitc]] = variables_dims[key]                    
                else:
                    variables[keep[ind_totlitc]] = variables[keep[ind_totlitc]] + variables[key]
        
            # summing up of total som C
            if key in totsomc:
                if key==totsomc[0]:
                    variables[keep[ind_totsomc]] = variables[key]
                    variables_dims[keep[ind_totsomc]] = variables_dims[key]                    
                else:
                    variables[keep[ind_totsomc]] = variables[keep[ind_totsomc]] + variables[key]


    WriteCLMDataToVTK("VAR_surf2d_"+options.clm_filehead,sxyz,scells, variables, variables_dims, True, nzmax)
    WriteCLMDataToVTK("VAR_soil3d_"+options.clm_filehead , xyz, cells, variables, variables_dims, True, nzmax)

