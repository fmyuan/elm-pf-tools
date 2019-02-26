#!/usr/bin/env python

import os, sys, time, math
import re
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
parser.add_option("--pfh5file", dest="pfh5file", default="", \
                  help = "pflotran output h5 file name without .h5 ")
parser.add_option("--time_unit", dest="tunit", default="", \
                  help = "time unit, default UNKNOWN ")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4) to be reading/plotting, separated by comma ")
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
vardata0 = []
vardata1 = []
vardata2 = []
vardata3 = []

for i in f0.keys():
    title = i.split()
    if title[0] == 'Time:':
        tt.append(float(title[1]))	
        group0 = f0[i]
            
        for j in group0.keys():
            h5vars = j.split()

            if(nfiles==1):
                vdata = group0[j][xpts,ypts,zpts]
            elif(nfiles==2):
                group1 = f1[i]
                vdata = group0[j][xpts,ypts,zpts]-group1[j][xpts,ypts,zpts]
            
            varname0 = varnames[0]
            if (h5vars[0] == varname0 and varname0 != ''):
                vardata0.append(vdata)

            if(nvars >= 2):
                varname1 = varnames[1]
                if (h5vars[0] == varname1 and varname1 != ''):
                    vardata1.append(vdata)

            if(nvars >= 3):
                varname2 = varnames[2]
                if (h5vars[0] == varname2 and varname2 != ''):
                    vardata2.append(vdata)

            if(nvars == 4):
                varname3 = varnames[3]
                if (h5vars[0] == varname3 and varname3 != ''):
                    vardata3.append(vdata)

# data-sets
t  = sorted(tt)
it = sorted(range(len(tt)), key=lambda k: tt[k])
nt = len(tt)
nl = np.size(xpts)*np.size(ypts)*np.size(zpts)

# plotting

nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
if(nvars>=2):
    nrow = 2
if(nvars>=3):
    ncol = 2

# plot 1
sdata = np.zeros((nt,nl))
for i in range(len(tt)):
    sdata[i,:] = np.squeeze(vardata0[it[i]])

ax0=plt.subplot(nrow, ncol, 1)
if(varname0 == 'Liquid_Pressure'):
    sdata = -sdata
    ax0.set_yscale("log", nonposy='clip')
    
plt.plot(t, sdata)
plt.xlabel('Time ('+time_unit+')')
plt.ylabel(varname0)

lx = 0.05
ly = 0.9
plt.text(lx, ly, '(a) ', transform=ax0.transAxes)

# plot 2
if (nvars >= 2):
    sdata = np.zeros((nt,nl))
    for i in range(len(tt)):
        sdata[i,:] = np.squeeze(vardata1[it[i]])

    ax1=plt.subplot(nrow, ncol, 2)
    if(varname1 == 'Liquid_Pressure'):
        sdata = -sdata
        #ax1.set_yscale("log", nonposy='clip')
    
    plt.plot(t, sdata)

    plt.xlabel('Time ('+time_unit+')')
    plt.ylabel(varname1)
    lx = 0.05
    ly = 0.9
    plt.text(lx, ly, '(b) ', transform=ax1.transAxes)

# plot 3
if (nvars >= 3):
    sdata = np.zeros((nt,nl))
    for i in range(len(tt)):
        sdata[i,:] = np.squeeze(vardata2[it[i]])

    ax2=plt.subplot(nrow, ncol, 3)
    if(varname2 == 'Liquid_Pressure'):
        sdata = -sdata
        #ax2.set_yscale("log", nonposy='clip')
    
    plt.plot(t, sdata)
    plt.xlabel('Time ('+time_unit+')')
    plt.ylabel(varname2)
    lx = 0.05
    ly = 0.9
    plt.text(lx, ly, '(c) ', transform=ax2.transAxes)

# plot 4
if (nvars >= 4):
    sdata = np.zeros((nt,nl))
    for i in range(len(tt)):
        sdata[i,:] = np.squeeze(vardata3[it[i]])

    ax3=plt.subplot(nrow, ncol, 4)
    if(varname3 == 'Liquid_Pressure'):
        sdata = -sdata
        #ax3.set_yscale("log", nonposy='clip')
    
    plt.plot(t, sdata)
    plt.xlabel('Time ('+time_unit+')')
    plt.ylabel(varname3)
    lx = 0.05
    ly = 0.9
    plt.text(lx, ly, '(d) ', transform=ax3.transAxes)

#
ofname = 'Figure_pflotran.pdf'
fig = plt.gcf()
fig.set_size_inches(12, 15)
plt.savefig(ofname)
plt.show()

#with PdfPages(ofname) as pdf:
#    fig = plt.figure()
#    plt.show()
#    pdf.savefig(fig)

plt.close('all')

f0.close()
if(nfiles==2):
    f1.close()

