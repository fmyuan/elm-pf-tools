#!/usr/bin/env python

import os, sys, time
import glob
import re

from optparse import OptionParser

#------------------------------------------------------------------------------------------------------------------
def filematching(path, case_sensitive=True, \
                 header='', ending='', pattern='*'):
    
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
    
    #
    return fullpath, rel_pathfile
 
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------

    
    
#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "work directory (default = ./, i.e., under current dir)")
parser.add_option("--srcheader", dest="srcheader", default="", \
                  help = "ELM src file header with path but no .F90 ")

(options, args) = parser.parse_args()

#
if (options.workdir == './'):
    print('file directory is the current')
else:
    print('file directory: '+ options.workdir)
cwdir = options.workdir

if (options.srcheader == ''):
    print('MUST input file name header, including fullpath, by " --srcheader=???"')
    sys.exit()

#

if not cwdir.endswith('/'): cwdir = cwdir+'/'

#------------------------------------------------------------------------------

# ELM src code files
if (options.srcheader != ""):
    elmpathfileheader = options.srcheader
    
    ftype = '.F90'

    alldirfiles = sorted(glob.glob("%s*.%s" % (elmpathfileheader, ftype)))
    if(len(alldirfiles)<=0):
        sys.exit("No file exists - %s*.%s IN %s" %(elmpathfileheader, ftype, cwdir))
    else:
        print('Total Files of ELM src: '+str(len(alldirfiles)))

    srcfileheader = elmpathfileheader.split('/')[-1]
    elm_odir   = elmpathfileheader.replace(srcfileheader,'')
    if(elm_odir.strip()==''):elm_odir='./'
    nonheader = '.f90'
    if srcfileheader.endswith('.'): nonheader=nonheader+'.'
    if (nonheader in srcfileheader):
        srcfileheader = srcfileheader.replace(nonheader,'')
    


#------------------------------------------------------------------------------------------------