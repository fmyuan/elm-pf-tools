#!/usr/bin/env python

from conda_build._link import SITE_PACKAGES
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob
import sys
from optparse import OptionParser

# ---------------------------------------------------------------
# commonly used for all

# CLM pft names
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
      "c4_grass",
      "c3_crop",
      "c3_irrigated",
      "corn",
      "irrigated_corn",
      "spring_temperate_cereal",
      "irrigated_spring_temperate_cereal",
      "winter_temperate_cereal",
      "irrigated_winter_temperate_cereal",
      "soybean",
      "irrigated_soybean"]

# some needed variables NOT output from CLM outputs, but can sum-up from others
totsomc_vr = ['SOIL1C_vr', 'SOIL2C_vr', 'SOIL3C_vr', 'SOIL4C_vr']
totsomn_vr = ['SOIL1N_vr', 'SOIL2N_vr', 'SOIL3N_vr', 'SOIL4N_vr']
totlitc_vr = ['LITR1C_vr', 'LITR2C_vr', 'LITR3C_vr']
totlitn_vr = ['LITR1N_vr', 'LITR2N_vr', 'LITR3N_vr']

# factors to adjust actual values of SOMC, if ad-spinup output
ad_factor =[1,1,10,100]

# -------------------------------------------------------------------
# read 1 CLM nc file into chunk
def CLM_NcRead(ncfile, varnames_print, keep_vars, chunk_keys):
    nx = 0
    ny = 0
    nldcmp = 0
    nlgrnd = 0
    npft   = 0
    chunk  = {}

    try:
        f = Dataset(ncfile,'r')
        if varnames_print: 
            print('FILE: '+ncfile+' ------- ')

    except:
        return nx,ny,nldcmp,nlgrnd,npft,chunk
        
        
    # If key is on the keep list and in nc file, uniquely add to a dictionary
    
    for key in f.variables.keys():
            
        # print out all variables contained in nc files
        if options.varnames_print:
            print (key)
            
        if key not in keep_vars: continue  #cycle the loop
            
        if (len(chunk_keys)<=0) or (key not in chunk_keys): # only needs to read data once, if more than one available
            chunk[key] = np.asarray(f.variables[key])                

            if key == 'topo':
                [nx, ny] = chunk[key].shape
            if key == 'levgrnd':
                nlgrnd = len(chunk[key])
            if key == 'levdcmp':
                nldcmp = len(chunk[key])                
            if key == 'pft':
                npft = len(chunk[key])

        #ad-spinup run, recal. each component of totsomc_vr, if needed
        if options.adspinup:
            if key in totsomc_vr:
                sub_indx = totsomc_vr.index(key)
                chunk[key] = chunk[key]*ad_factor[sub_indx]

            if key in totsomn_vr:
                sub_indx = totsomn_vr.index(key)
                chunk[key] = chunk[key]*ad_factor[sub_indx]
            
    
    # summing up of total liter/som C/N, if needed and available
    if 'TOTLITC_vr' in keep_vars:
        indx = keep_vars.index('TOTLITC_vr')
        subs = totlitc_vr
        for isub in subs:
            if isub in chunk.keys():
                if isub==subs[0]:
                    chunk[keep_vars[indx]] = chunk[isub]
                else:
                    chunk[keep_vars[indx]] = chunk[keep[indx]] + chunk[isub]
            else:
                continue

    if 'TOTLITN_vr' in keep_vars:
        indx = keep_vars.index('TOTLITN_vr')
        subs = totlitn_vr
        for isub in subs:
            if isub in chunk.keys():
                if isub==subs[0]:
                    chunk[keep_vars[indx]] = chunk[isub]
                else:
                    chunk[keep_vars[indx]] = chunk[keep[indx]] + chunk[isub]
            else:
                continue
        
    if 'TOTSOMC_vr' in keep_vars:
        indx = keep_vars.index('TOTSOMC_vr')
        subs = totsomc_vr
        for isub in subs:
            if isub in chunk.keys():
                if isub==subs[0]:
                    chunk[keep_vars[indx]] = chunk[isub]
                else:
                    chunk[keep_vars[indx]] = chunk[keep[indx]] + chunk[isub]
            else:
                continue

    if 'TOTSOMN_vr' in keep_vars:
        indx = keep_vars.index('TOTSOMN_vr')
        subs = totsomn_vr
        for isub in subs:
            if isub in chunk.keys():
                if isub==subs[0]:
                    chunk[keep_vars[indx]] = chunk[isub]
                else:
                    chunk[keep_vars[indx]] = chunk[keep[indx]] + chunk[isub]
            else:
                continue

    # out datasets
    return nx, ny, nldcmp, nlgrnd, npft, chunk


# Plot soil data with layers
def SoilLayeredVarPlotting(plt, nrow, ncol, isubplot, t, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    if(varname == 'SOILPSI'):
        sdata = -sdata
        ax.set_yscale("log", nonposy='clip')
    
    plt.plot(t, sdata[:,0], 'b-.', t, sdata[:,1], 'g-.', \
                  t, sdata[:,2], 'm-.', t, sdata[:,3], 'c-.', \
                  t, sdata[:,4], 'r-.', \
                  t, sdata[:,5], 'b-', t, sdata[:,6], 'g-', \
                  t, sdata[:,7], 'm-', t, sdata[:,8], 'c-', \
                  t, sdata[:,9], 'r-', \
                  t, sdata[:,10], 'b--', \
                  t, sdata[:,11], 'g--',  t, sdata[:,12], 'm--', \
                  t, sdata[:,13], 'c--',  t, sdata[:,14], 'r--')
    plt.legend(('Layer 15','Layer 14','Layer 13','Layer 12','Layer 11', \
            'Layer 10','Layer 9','Layer 8','Layer 7','Layer 6', \
            'Layer 5','Layer 4','Layer 3','Layer 2','Layer 1'), \
            loc=0, fontsize=12)
    plt.xlabel('Years')
    plt.ylabel(varname)

    lx = 0.40
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

# Plot Grided data
def GridedVarPlotting(plt, nrow, ncol, isubplot, t, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    
    plt.plot(t, sdata)
    
    #plt.legend(loc=0, fontsize=12)
    plt.xlabel('Years')
    plt.ylabel(varname)

    lx = 0.40
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--clmout_dir", dest="clm_odir", default="./", \
                  help="clm output directory (default = ./, i.e., under current directory)")
parser.add_option("--clmfile_head", dest="ncfileheader", default="", \
                  help = "clm output file name header, usually the portion before *.clm2.h[0-5].*.nc")
parser.add_option("--varname_help", dest="varnames_print", default=False, \
                  help = "print out VARIABLES available ")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4) to be reading/plotting, separated by comma ")
parser.add_option("--starting_time", dest="starttime", default="1", \
                  help="clm run starting time (default = 1, this is for spinup; for transient it should be 1850)")
parser.add_option("--adspinup", dest="adspinup", action="store_true", default=False, \
                  help="whether results of an ad_spinup run (default = False)")
parser.add_option("--Xindex", dest="xindex", default=0, \
                  help = " X direction grid index to be reading/plotting, default 0 ")
parser.add_option("--Yindex", dest="yindex", default=0, \
                  help = " Y direction grid index to be reading/plotting, default 0 ")

(options, args) = parser.parse_args()

#
if (options.clm_odir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.clm_odir)

if (options.ncfileheader == ''):
    print('MUST input clm file header by " --ncfileheader=??? "' //
     ' , which usually is the portion before "*.clm2.[h0-h5].*.nc"! ')

    sys.exit()
else:
    files = options.ncfileheader.split(':')  
    nfiles = len(files)
    if(nfiles>2):
        print('Currently ONLY support at most 2 sets of nc data series to be differenced')
        sys.exit()
    else:  
        clmhead = options.clm_odir+'/'+files[0]+'.clm2'
        print('nc file set 1 : '+ clmhead + '*.nc')
        if(nfiles==2):
            clmhead2 = options.clm_odir+'/'+files[1]+'.clm2'
            print('nc file set 2 : '+ clmhead2 + '*.nc')

if (options.vars == ''):
    print('No variable name by " --varname=??? "; So print out ALL variable names')
    varnames =[]
    options.varnames_print = True
else:
    varnames = options.vars.split(':')  
    nvars = len(varnames)
    if(nvars>4):
        print('Currently ONLY support at most 4 variables to be plotted')
        sys.exit()


starttime = int(options.starttime)
#

ix=int(options.xindex);
iy=int(options.yindex);

#--------------------------------------------------------------------------------------

# variables name dictionary
keep = []

# a few variables common to all for plotting
keep.append("nstep")   # need this for time step numbers
keep.append("time")    # need this for clm timing
keep.append("topo")    # need this for shape of surface cells
keep.append("levgrnd") # need this for the number of layers  for PHY
keep.append("levdcmp") # need this for the number of layers for BGC
keep.append("pft")
keep.append("pfts1d_wtgcell")

keep.extend(varnames)

# for some variables, they are requiring a group of variables in CLM output to sum up
if 'TOTSOMC_vr' in varnames:
    keep.extend(totsomc_vr)

if 'TOTSOMN_vr' in varnames:
    indx_totsomn_vr = keep.index('TOTSOMN_vr')
    keep.extend(totsomn_vr)

if 'TOTLITC_vr' in varnames:
    indx_totlitc_vr = keep.index('TOTLITC_vr')
    keep.extend(totlitc_vr)

if 'TOTLITN_vr' in varnames:
    indx_totlitn_vr = keep.index('TOTLITN_vr')
    keep.extend(totlitn_vr)

# build all nc files into one or two arrays
fincludes = ['h0','h1','h2','h3','h4','h5']
allfile = glob.glob("%s.h*.*.nc" % clmhead)
fchunks  = []
for filename in allfile:
    filename = filename.split(".")
    if fchunks.count(filename[-2]) == 0: fchunks.append(filename[-2])
#
if nfiles == 2:
    allfile = glob.glob("%s.h*.*.nc" % clmhead2)
    fchunks2  = []
    for filename in allfile:
        filename = filename.split(".")
        if fchunks.count(filename[-2]) == 0: fchunks.append(filename[-2])
    
# final datasets initialization
tt = []
vardata0 = []
vardata1 = []
vardata2 = []
vardata3 = []

nx     = 0
ny     = 0
nldcmp = 0
nlgrnd = 0
npft   = 0          

# process CLM files
for chunk in fchunks:  # actually one time-series in a file

    # for each chunk of time slices, append variables to a single dictionary
    chunkdata = {}
    for finc in fincludes:        
        # try reading all files one by one
        filename = "%s.%s.%s.nc" % (clmhead,finc,chunk)
        
        nxi,nyi,nldcmpi,nlgrndi,npfti,chunkdatai=CLM_NcRead(filename, options.varnames_print, keep, chunkdata.keys())
        if nxi>0 and nx<=0: nx = nxi
        if nyi>0 and ny<=0: ny = nyi
        if nldcmpi>0 and nldcmp<=0: nldcmp = nldcmpi
        if nlgrndi>0 and nlgrnd<=0: nlgrnd = nlgrndi
        if npfti>0 and npft<=0: npft = npfti
        
        if len(chunkdatai)>0: 
            if (len(chunkdata)<=0):
                chunkdata = chunkdatai
            else:
                for ikey in chunkdatai:
                    if ikey in chunkdata:
                        chunkdata[ikey].extend(chunkdatai[ikey])
                    else:
                        chunkdata[ikey] = chunkdatai[ikey]

    #append chunks
    if options.varnames_print: 
        sys.exit("printing out variable names is DONE")
    else:
        #time
        tt.extend(chunkdata['time'])
        
        indx=varnames[0]
        vardata0.extend(chunkdata[indx])
        
        if(nvars>=2):
            indx=varnames[1]
            vardata1.extend(chunkdata[indx])
            
        if(nvars>=3):
            indx=varnames[2]
            vardata2.extend(chunkdata[indx])

        if(nvars==4):
            indx=varnames[3]
            vardata3.extend(chunkdata[indx])
            
    

# data-sets sorted by time-series
t  = sorted(tt)
it = sorted(range(len(tt)), key=lambda k: tt[k])
nt = len(tt)
nl = max(nldcmp, nlgrnd)

# plotting

nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
if(nvars>=2):
    nrow = 2
if(nvars>=3):
    ncol = 2

sdata = np.zeros((nt,nl))   #temporary data holder
gdata = np.zeros((nt,nx*ny))

# plot 1
nztrunc = 1
if(len(vardata0[0].shape)==3): nztrunc=vardata0[0].shape[0]
for i in range(len(tt)):
    if(nztrunc==nldcmp or nztrunc==nlgrnd):
        sdata[i,:nztrunc] = vardata0[it[i]][:, iy, ix]
    elif(nztrunc==1):
        gdata[i,:]=vardata0[it[i]].reshape(nx*ny)

if(nztrunc==nldcmp or nztrunc==nlgrnd):
    SoilLayeredVarPlotting(plt, nrow, ncol, 1, t, sdata, varnames[0], '(a) Grid ( '+str(ix)+', '+str(iy)+')')
elif(nztrunc==1):
    GridedVarPlotting(plt, nrow, ncol, 1, t, gdata, varnames[0], '(a) All Grids')
    

# plot 2
if (nvars >= 2):
    nztrunc = 1
    if(len(vardata1[0].shape)==3): nztrunc=vardata1[0].shape[0]
    for i in range(len(tt)):
        if(nztrunc==nldcmp or nztrunc==nlgrnd):
            sdata[i,:] = vardata1[it[i]][:, iy, ix]    
        elif(nztrunc==1):
            gdata[i,:]=vardata1[it[i]].reshape(nx*ny)
    
    if(nztrunc==nldcmp or nztrunc==nlgrnd):
        SoilLayeredVarPlotting(plt, nrow, ncol, 2, t, sdata, varnames[1], '(b) Grid ( '+str(ix)+', '+str(iy)+')')
    elif(nztrunc==1):
        GridedVarPlotting(plt, nrow, ncol, 1, t, gdata, varnames[1], '(b) All Grids')


# plot 3
if (nvars >= 3):
    nztrunc = 1
    if(len(vardata2[0].shape)==3): nztrunc=vardata2[0].shape[0]
    for i in range(len(tt)):
        if(nztrunc==nldcmp or nztrunc==nlgrnd):
            sdata[i,:] = vardata2[it[i]][:, iy, ix]    
        elif(nztrunc==1):
            gdata[i,:]=vardata2[it[i]].reshape(nx*ny)
    
    if(nztrunc==nldcmp or nztrunc==nlgrnd):
        SoilLayeredVarPlotting(plt, nrow, ncol, 2, t, sdata, varnames[2], '(c) Grid ( '+str(ix)+', '+str(iy)+')')
    elif(nztrunc==1):
        GridedVarPlotting(plt, nrow, ncol, 1, t, gdata, varnames[2], '(c) All Grids')

# plot 4
if (nvars == 4):
    nztrunc = 1
    if(len(vardata3[0].shape)==3): nztrunc=vardata3[0].shape[0]
    for i in range(len(tt)):
        if(nztrunc==nldcmp or nztrunc==nlgrnd):
            sdata[i,:] = vardata3[it[i]][:, iy, ix]    
        elif(nztrunc==1):
            gdata[i,:]=vardata3[it[i]].reshape(nx*ny)
    
    if(nztrunc==nldcmp or nztrunc==nlgrnd):
        SoilLayeredVarPlotting(plt, nrow, ncol, 2, t, sdata, varnames[3], '(d) Grid ( '+str(ix)+', '+str(iy)+')')
    elif(nztrunc==1):
        GridedVarPlotting(plt, nrow, ncol, 1, t, gdata, varnames[3], '(d) All Grids')


# printing plot in PDF
ofname = 'Figure_CLM.pdf'
fig = plt.gcf()
fig.set_size_inches(12, 15)
plt.savefig(ofname)
plt.show()

#with PdfPages(ofname) as pdf:
#    fig = plt.figure()
#    plt.show()
#    pdf.savefig(fig)

plt.close('all')


