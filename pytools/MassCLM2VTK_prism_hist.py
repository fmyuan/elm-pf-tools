#!/usr/bin/env python

import sys
import glob
import math
import numpy as np
import h5py as h5
from Modules_CLM_nc4 import CLM_NcRead_1simulation
from optparse import OptionParser

#
inityr = 1

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
    nl   = int(cells.shape[0]/cpl)

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

#
def WriteCLMDataToVTK(filehead, xyz, cells, timekey, clmdata, clmdata_dims, updown, nztrunc):
    
    vhdr = timekey.split("_")[0]
    
    steps  = np.asarray(clmdata[vhdr+"_nstep"],dtype=int)
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

    times  = np.asarray(clmdata[vhdr+"_time"],dtype=int)   # "time" in clm output *.nc IS in days from starting time of the year

    # loop over all time steps
    for i in range(nstep):
        # build another dictionary with references to this time steps variables
        var = {}
        tyr = math.floor(times[i]/365.0)                    # converting time unit from days to year
        tdoy= int(times[i]-tyr*365.0)           # converting time unit from days to doy
        for key in clmdata.keys():
            
            if vhdr not in key: continue   # skip if not same time-series data
                        
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
                        data = data[:nz,]
                    elif(ind_zdim==1):
                        data = data[:,:nz,]
                    elif(ind_zdim==2):
                        data = data[:,:,:nz,]
                    elif(ind_zdim==3):
                        data = data[:,:,:,:nz,]
                    else:
                        print('Warning: Z-dimension beyond 4-D: please check ! ')
                        sys.exit()
                    
                    shp = data.shape           

                #if VTK 3-D PF mesh is upside-down, then need to do so for CLM data from defined bottom layer
                # (this MUST be done after truncation of bottom layers above)
                if updown: 
                    if(ind_zdim==0): 
                        data = data[shp[ind_zdim]::-1,]
                    elif(ind_zdim==1): 
                        data = data[:,shp[ind_zdim]::-1,]
                    elif(ind_zdim==2):
                        data = data[:,:,shp[ind_zdim]::-1,]
                    elif(ind_zdim==3):
                        data = data[:,:,:,shp[ind_zdim]::-1,]
                    else:
                        print('Warning: Z-dimension beyond 4-D: please check ! ')
                        sys.exit()
                                 
            if np.prod(data[i,...].shape) == ncells:
                # shape of data matches the number of cells in the VTK mesh
                var[key] = data[i,...].reshape((cells.shape[0]))
        
        WriteVTK("%s-yyyydoy-%s%s.vtk" % (filehead,str.zfill(str(int(tyr+inityr)),4),str.zfill(str(int(tdoy)),3)),
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
                  help="max. z depth numbers (default = 15)")
parser.add_option("--zscale", dest="zscale", default="1", \
                  help="z-axis scaling (default = 1)")
parser.add_option("--adspinup", dest="adspinup", action="store_true", default=False, \
                  help="whether results of an ad_spinup run (default = False)")
parser.add_option("--upsidedown", dest="updown", action="store_true", default=False, \
                  help="whether flip-over data vertically (upside-down) (default = False)")
parser.add_option("--startyr", dest="startyr", default="1", \
                  help="clm run starting year (default = 1, this is for spinup; for transient it should be 1850; " \
                   " and can be user-defined)")
parser.add_option("--endyr", dest="endyr", default="", \
                  help="clm run ending year (default = none, i.e. end of simulation)")
parser.add_option("--varname", dest="vars", default="ALL", \
                  help = "variable name(s) (default: ALL) to be reading, separated by comma ")
(options, args) = parser.parse_args()

#
if (options.nztrunc==""):
    nzmax = 15
else:
    nzmax = int(options.nztrunc)

#
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

if (options.vars == ''):
    varnames = 'ALL'
else:
    varnames = options.vars.split(':')  

startyyyy = 1
if(int(options.startyr) > 1): startyyyy = int(options.startyr)
endyyyy   = 9999
if(options.endyr !=""):   endyyyy = int(options.endyr)

#--------------------------------------------------------------------------------------
fincludes = ['h0','h1','h2','h3','h4','h5']

# read-in datasets from 1 simulation year by year
tmax = 0
for iyr in range(startyyyy,endyyyy):
    startdays = (int(iyr)-1)*365.0+1.0
    enddays   = (int(iyr)+0)*365.0

    nx, ny, nlgrnd, nldcmp, npft, varsdata, varsdims = \
        CLM_NcRead_1simulation(options.clm_odir, \
                           options.clm_filehead, \
                           False, \
                           varnames, \
                           startdays, enddays, \
                           options.adspinup)
    tt=[]
    for finc in fincludes:
        timekey = finc+"_time"          
        if timekey in varsdata: 
            tt = np.append(tt,varsdata[timekey])       
        
            WriteCLMDataToVTK("VAR_surf2d_"+options.clm_filehead+"_"+finc,sxyz,scells, \
                                  timekey, varsdata, varsdims, options.updown, nzmax)
            WriteCLMDataToVTK("VAR_soil3d_"+options.clm_filehead+"_"+finc, xyz, cells, \
                                  timekey, varsdata, varsdims, options.updown, nzmax)

    
    if(len(tt)==0 or tmax>=max(tt)): 
        exit("DONE!") # end of all times in the simulation
    else: 
        tmax = max(tt) #updating max. time

