#!/usr/bin/env python

#----------------------
# Extract point(s) data from ELM outputs
#
import os, sys
import glob
import re
import math
from optparse import OptionParser
import numpy as np

import netcdf_modules as nfmod
from builtins import int, float

#-------------------Local functions --------------------------------------------
def Daymet_pixel_extract(mapfile, lonlat_pts=[], geoxy_pts=[]):

    # read-in mapping file
    #mapfile = options.gridmap.strip()
    with open(mapfile, 'r') as f:
        # remove header txts
        next(f)
        try:
            data = [x.strip().split() for x in f]
            f.close()
        except Exception as e:
            print(e)
            print('Error in reading - '+mapfile)
            sys.exit(-1)
    data = np.asarray(data,float)
    lon=data[:,0]
    lat=data[:,1]
    geox=data[:,2]
    geoy=data[:,3]
    xidx=np.asanyarray(data[:,4],int)-1  # in mappings.txt, those index are 1-based
    yidx=np.asanyarray(data[:,5],int)-1
    gidx=np.asanyarray(data[:,6],int)-1
    
    # extrac pixel either by paired lat/lon or geox/geoy
    f = open('Extracted_'+mapfile, 'w')
    fheader='   lon          lat            geox            geoy        i     j     g '
    
    pts_idx=[]
    if len(lonlat_pts)>0:
        #
        
        f.write(fheader+'\n')
        for ipt in range(len(lonlat_pts)):
            if lonlat_pts[ipt][0]<0:
                xx = np.abs(lonlat_pts[ipt][0]+360.0 - lon)
            else:
                xx = np.abs(lonlat_pts[ipt][0] - lon)
            yy = np.abs(lonlat_pts[ipt][1] - lat)
            dd = xx**2+yy**2
            ig = np.argmin(dd)
            
            #print(ipt, lonlat_pts[ipt], ig, dd[ig], xidx[ig], yidx[ig], gidx[ig])
            # 
            if dd[ig]>=0:
                pts_idx.append([xidx[ig],yidx[ig],gidx[ig]])
                
                #'(f12.5,1x,f12.6,1x, 2(f15.1, 1x),3(I5,1x))'
                f.write('%12.5f ' % lon[ig] )
                f.write('%12.6f ' % lat[ig] )
                f.write('%15.1f ' % geox[ig] )
                f.write('%15.1f ' % geoy[ig] )
                f.write('%5d ' % (xidx[ig]+1) )  #x/yidx were 0-based, but need to 1-based for the mapping file
                f.write('%5d ' % (yidx[ig]+1) )
                f.write('%5d ' % (gidx[ig]+1) )
                f.write('\n')
        f.close()
    #
    elif len(geoxy_pts)>0:
        #TODO
        f.close()
        
    else:
        print('No point read in!')
        f.close()

    
    return pts_idx


#------------------------------ --------------------------------------------
parser = OptionParser()

parser.add_option("--outdir", dest="outdir", default="", \
                  help = 'output directory')
parser.add_option("--ofheader", dest="header", default="none", \
                  help = 'header of output files either in fullpath or filename in default path')
parser.add_option("--newdata_affix", dest="newdata_affix", default="", \
                  help = 'new dataset output with affix')
parser.add_option("--lons", dest="ptlon", default="", \
                      help = "point or ranged(<,<=,>=,>,~) longitude(s), separated with',:;'")
parser.add_option("--lats", dest="ptlat", default="", \
                      help = "point or ranged (<,<=,>=,>,~) latitude(s), separated with',:;'")
parser.add_option("--daymet_elm_mapfile", dest="mapfile", default="daymet_elm_mappings.txt", \
                  help = 'daymet_elm_mapfile name')
parser.add_option("--ncobinpath", dest="ncobinpath", default="", \
                      help = "NCO bin path if not in $PATH")


(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------
currdir = os.getcwd()

#--------------- nco bin path -----------
if(not options.ncobinpath.strip()==''):
    ncopath = options.ncobinpath
else:
    ncopath = '/usr/local/nco/bin'

###########################################################################
if(options.outdir.startswith('/')):
    outdir = os.path.abspath(options.outdir)
    if(not outdir.endswith('/')): outdir = outdir+'/'
else:
    outdir = './'


#------------extract all data as needed -----------------------------------
if (os.path.isdir(outdir.strip())):
    lonlat_pts=[]
    optionx = options.ptlon.split(':')
    optiony = options.ptlat.split(':')
    for ip in range(len(optionx)):
        lonlat_pts.append([float(optionx[ip]),float(optiony[ip])])

    pts_index = Daymet_pixel_extract(options.mapfile, lonlat_pts)
    if len(pts_index)<=0:
        print('Error: NO point found in daymet_elm_map -', options.mapfile)
        system.exit(-1)

    #checking where is the original data
    dirfiles = os.listdir(outdir)
    for dirfile in dirfiles:
        filehead = options.header
        filehead_new = filehead+'_'+options.newdata_affix
        
        # in outdir directory
        if(os.path.isfile(outdir+'/'+dirfile)): # file names
            if(filehead in dirfile):
            
                print('\n file: '+dirfile)
    
                outfile_old = outdir+'/'+dirfile
                #
                outfile_temp = outdir+'/temp.nc'
                os.system('cp -f ' + outfile_old + ' '+ outfile_temp)
                
                for ipt in range(len(pts_index)):
                    wdigit=int(math.log10(len(pts_index)))
                    dirfile_new = dirfile.replace(filehead,filehead_new+str(ipt).zfill(wdigit))
                    outfile_new = outdir+'/'+dirfile_new

                    print('INFO: extracting clm output - \n', outfile_old, '\n -->', outfile_new)
                    
                    dim_string = ' -d geox,'+str(pts_index[ipt][0])+','+str(pts_index[ipt][0]) \
                                +' -d geoy,'+str(pts_index[ipt][1])+','+str(pts_index[ipt][1])
                    os.system(ncopath+'./ncks --no_abc -O'+ dim_string+ \
                          ' '+outfile_temp.strip()+' -o '+outfile_new.strip())

                #
                os.system('rm -rf '+outfile_temp.strip())
                
                #merge file
                ierr = os.system('find '+outdir+'/ -name "'+filehead_new+'*'+ \
                                 '" | xargs ls | sort | '+ncopath+'./ncecat -O -h -o ' \
                                 +outdir+'/temp0.nc')
                if(ierr!=0): 
                    raise RuntimeError('Error: ncecat ')
                else:
                    os.system('find '+outdir+'/ -name "'+filehead_new+'*'+ \
                                 '" -exec rm {} \; ')
                
                #point dim 'record' from 'ncecat' swap to 'geox' position
                #print(ncopath+'./ncpdq -h -O -a geox,record '+outdir+'/temp0.nc -o '+outdir+'/temp1.nc')
                ierr = os.system(ncopath+'./ncpdq -h -O -a geox,record '+outdir+'/temp0.nc -o '+outdir+'/temp1.nc')
                if(ierr!=0): raise RuntimeError('Error: ncpdq ')
                
                # remove 'geox/geoy' dims
                ierr = os.system(ncopath+'ncwa -h -O -a geox '+outdir+'/temp1.nc -o '+outdir+'/temp2.nc')
                if(ierr!=0): raise RuntimeError('Error: ncwa ')
                ierr = os.system(ncopath+'ncwa -h -O -a geoy '+outdir+'/temp2.nc -o '+outdir+'/temp1.nc')
                if(ierr!=0): raise RuntimeError('Error: ncwa ')
                # rename 'record' dim to 'gridcell'
                ierr = os.system(ncopath+'ncrename -h -O -d record,gridcell '+outdir+'/temp1.nc '+outdir+'/temp2.nc')
                if(ierr!=0): raise RuntimeError('Error: ncrename ')
                #rename temp filename to new file name
                ierr = os.system(ncopath+'./ncks --no_abc -O -x -v geox,geoy '+outdir+'/temp2.nc -o '+outdir+'/multi_'+dirfile.replace(filehead,filehead_new))
                if(ierr!=0): 
                    raise RuntimeError('Error: ncks ')
                else:
                    os.system('rm -f '+outdir+'/temp*.nc')


            # (END) if-file contains 'fileheader'
        # (END) if-it's a file
    #(END) for-loop of all files in outdir directory

    
#------------END of pointELM_daymet_out -----------------------------------------------------------------------
    

