#!/usr/bin/env python

import os, sys, time, math
import re
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser


#-------------------PLotting submodule-----------------------------------------------
# plotting 1 graph with at most 4 sub-plots 

def SinglePlot(tt, time_unit, varnames, varunits, vardatas, figno=None):
    nvars = len(varnames)
    
    nrow = 1   # sub-plot vertically arranged number (row no.)
    ncol = 1   # sub-plot horizontally arranged number (column no.)
    if(nvars>=2):
        nrow = 2
    if(nvars>=3):
        ncol = 2


    fig = plt.figure(figsize=(7.5,9.5))
     
    # plot 1, or subplot 1 if more than 1 subplot
    if (nvars >= 2):
        sdata = vardatas[0][:]
    else:
        sdata = vardatas
       
    ax0=plt.subplot(nrow, ncol, 1)
    if (not figno is None): plt.suptitle('FIGURE '+figno)
    plt.plot(tt, sdata)
    plt.xlabel('Time ('+time_unit+')')
    plt.ylabel(varnames[0]+varunits[0])

    lx = 0.05
    ly = 0.95
    if(nvars>=2):
        plt.text(lx, ly, '(a) ', transform=ax0.transAxes)
    else:
        plt.text(lx, ly, '', transform=ax0.transAxes)
        
    # subplot 2
    if (nvars >= 2):
        sdata = vardatas[1][:]
        
        ax1=plt.subplot(nrow, ncol, 2)
        plt.plot(tt, sdata)
        plt.xlabel('Time ('+time_unit+')')
        plt.ylabel(varnames[1]+varunits[1])
        lx = 0.05
        ly = 0.95
        plt.text(lx, ly, '(b) ', transform=ax1.transAxes)

    # subplot 3
    if (nvars >= 3):
        sdata = vardatas[2][:]

        ax2=plt.subplot(nrow, ncol, 3)
        plt.plot(tt, sdata)
        plt.xlabel('Time ('+time_unit+')')
        plt.ylabel(varnames[2]+varunits[2])
        lx = 0.05
        ly = 0.95
        plt.text(lx, ly, '(c) ', transform=ax2.transAxes)

    # sub-plot 4
    if (nvars >= 4):
        sdata = vardatas[3][:]

        ax3=plt.subplot(nrow, ncol, 4)
        plt.plot(tt, sdata)
        plt.xlabel('Time ('+time_unit+')')
        plt.ylabel(varnames[3]+varunits[3])
        lx = 0.05
        ly = 0.95
        plt.text(lx, ly, '(d) ', transform=ax3.transAxes)

    #
    ofname = 'Figure_pflotran.pdf'
    if (not figno is None): ofname = 'Figure_pflotran-'+figno+'.pdf'
    plt.savefig(ofname)
    plt.show()

    plt.close('all')



#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
parser.add_option("--pfh5file", dest="pfh5file", default="", \
                  help = "pflotran output h5 file name without .h5 ")
parser.add_option("--time_unit", dest="tunit", default="", \
                  help = "time unit, default UNKNOWN ")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4 in one PLOT, or, with ALL for each var in one PLOT) to be reading/plotting, separated by comma ")
parser.add_option("--Xindex", dest="xindex", default=0, \
                  help = " X direction mesh index to be reading/plotting, default 0 ")
parser.add_option("--Yindex", dest="yindex", default=0, \
                  help = " Y direction mesh cell index to be reading/plotting, default 0 ")
parser.add_option("--Zindex", dest="zindex", default=0, \
                  help = " Z direction mesh cell index to be reading/plotting, default 0 ")

(options, args) = parser.parse_args()

#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)

if (options.pfh5file == ''):
    print('MUST input file name by " --pfh5file=??? "')
    sys.exit()
else:
    files = options.pfh5file.split(':')  
    nfiles = len(files)
    if(nfiles>2):
        print('Currently ONLY support at most 2 h5 datasets to be plotted')
        sys.exit()
    else:  
        filename1=options.workdir+'/'+files[0]+'.h5'
        print('h5file 1 : '+filename1)
        if(nfiles==2):
            filename2=options.workdir+'/'+files[1]+'.h5'
            print('h5file 2: '+filename2)

if (options.vars == ''):
    print('MUST input at least one variable name by " --varname=??? "')
    sys.exit()
elif (options.vars == 'ALL'):
    print('WILL plot each Variable One by One - could be hugh!')
    varnames = []    
else:
    varnames = options.vars.split(':')  
    nvars = len(varnames)
    if(nvars>4):
        print('Currently ONLY support at most 4 variables to be plotted')
        sys.exit()

if (options.tunit == ''):
    time_unit =''
    print('You may input time unit by " --time_unit=??? "')
else:
    time_unit = options.tunit
    print('Time Unit: '+time_unit)

#
pts = [None]*3
str_pts = [options.xindex, options.yindex, options.zindex]
for i in np.arange(len(str_pts)):
    str_pt = str_pts[i]
    if(str_pt==0):
        pts[i] = 0
    else:
        if(str_pt.strip()==':'):
            pts[i]=1
        else:                
            if ('str' in str(type(str_pt))):#values with operators
                nondecimal = re.compile(r'[^\d.,:]+')
                v_pt = nondecimal.sub("",str_pt)
                if(':' in str_pt):
                    if sys.version_info[0] < 3:
                        v_pt = v_pt.split(':')
                    else:
                        v_pt = np.split(v_pt,':')
                    pts[i] = np.arange(int(v_pt[0]),int(v_pt[1])+1, dtype=int)
            
                elif(',' in str_pt):
                    pts[i] = np.fromstring(v_pt,sep=',')
                else:
                    pts[i] = np.int(v_pt)                                                  
                   
            else: # exactly point(s)
                pts[i] = np.int(str_pt)

xpts = pts[0]
ypts = pts[1]
zpts = pts[2]
#
f0 = h5.File(filename1, 'r')
if(nfiles==2):
    f1 = h5.File(filename2, 'r')

tt = []
varunits = []
vardatas = {}

for i in f0.keys():
    title = i.split()
    # PFLOTRAN h5 data grouped by 'Time' steps, so needs to concat. by time for each variable 
    if title[0] == 'Time:':
        tt.append(float(title[1]))	
        group0 = f0[i]
        if time_unit=='': time_unit = title[2]
        for j in group0.keys():
            h5vars = j.split()

            if(nfiles==1):
                vdata = group0[j][xpts,ypts,zpts]
            elif(nfiles==2):
                group1 = f1[i]
                vdata = group0[j][xpts,ypts,zpts]-group1[j][xpts,ypts,zpts]
                        
            # vars individually
            varname = h5vars[0]
                
            if(options.vars == 'ALL'): # all variables, with 1 variable in 1 plot in 1 single figure
                if varname not in varnames: 
                    varnames.append(varname) # build a variable list
                    if len(h5vars)>1: 
                        varunits.append(h5vars[1])
                    else:
                        varunits.append('')

                    vardatas[varname] = vdata
                else:
                    vardatas[varname]=np.vstack([vardatas[varname],vdata])
                
            else: # user-picked variable(s) at most 4 for (sub-)plotting in 1 single figure               
                if varname in varnames: 
                    if varname not in vardatas.keys():
                        vardatas[varname] = vdata
                        if len(h5vars)>1: 
                            varunits.append(h5vars[1])
                        else:
                            varunits.append('')
                    else:
                        vardatas[varname]=np.vstack([vardatas[varname],vdata])
                        
                    


# CANNOT assume data-sets are time-sorted, so must sort them by time at first
t  = sorted(tt)
it = sorted(range(len(tt)), key=lambda k: tt[k])

nvars = len(varnames)
nt = len(tt)
nl = np.size(xpts)*np.size(ypts)*np.size(zpts)
if(options.vars == 'ALL'): # all variables, with 1 variable in 1 plot
    sdata = np.zeros((nt,nl))
    for ivar in range(0,nvars):
        # plotting 'varname' one by one
        vardata = vardatas[varnames[ivar]]
        varname = varnames[ivar]
        varunit = varunits[ivar]
        for i in range(len(tt)):
            sdata[i,:] = np.squeeze(vardata[it[i]])
        SinglePlot(t, time_unit, [varname], [varunit], sdata, figno=str(ivar)+' - '+varname)
        
                
else: # user-picked variable(s) at most 4 for (sub-)plotting in 1 single plot               
    sdata = np.zeros((nvars,nt,nl))
    ivar = -1
    for varname in varnames:
        # plotting one or more sub-plots in 1 single figure
        ivar = ivar + 1
        vardata = vardatas[varname]
        for i in range(len(tt)):
            sdata[ivar,i,:] = np.squeeze(vardata[it[i]])

    # (sub-)plotting in 1 single plot
    SinglePlot(t, time_unit, varnames, varunits, sdata)

f0.close()
if(nfiles==2):
    f1.close()

