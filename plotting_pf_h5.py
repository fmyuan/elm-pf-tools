#!/usr/bin/env python

from conda_build._link import SITE_PACKAGES
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
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4) to be reading/plotting, separated by comma ")
parser.add_option("--Xindex", dest="xindex", default=0, \
                  help = " X direction grid index to be reading/plotting, default 0 ")
parser.add_option("--Yindex", dest="yindex", default=0, \
                  help = " Y direction grid index to be reading/plotting, default 0 ")

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

#

x=int(options.xindex);
y=int(options.yindex);

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
                vdata = group0[j][x,y,:]
            elif(nfiles==2):
                group1 = f1[i]
                vdata = group0[j][x,y,:]-group1[j][x,y,:]
            
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
nl = len(vardata0[0])

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
    sdata[i,:] = vardata0[it[i]]

ax0=plt.subplot(nrow, ncol, 1)
if(varname0 == 'Liquid_Pressure'):
    sdata = -sdata
    ax0.set_yscale("log", nonposy='clip')
    
plt.plot(t, sdata)
plt.xlabel('Years')
plt.ylabel(varname0)

lx = 0.05
ly = 0.9
plt.text(lx, ly, '(a) ', transform=ax0.transAxes)

# plot 2
if (nvars >= 2):
    sdata = np.zeros((nt,nl))
    for i in range(len(tt)):
        sdata[i,:] = vardata1[it[i]]

    ax1=plt.subplot(nrow, ncol, 2)
    if(varname1 == 'Liquid_Pressure'):
        sdata = -sdata
        ax1.set_yscale("log", nonposy='clip')
    
    plt.plot(t, sdata)

    plt.xlabel('Time')
    plt.ylabel(varname1)
    lx = 0.05
    ly = 0.9
    plt.text(lx, ly, '(b) ', transform=ax1.transAxes)

# plot 3
if (nvars >= 3):
    sdata = np.zeros((nt,nl))
    for i in range(len(tt)):
        sdata[i,:] = vardata2[it[i]]

    ax2=plt.subplot(nrow, ncol, 3)
    if(varname2 == 'Liquid_Pressure'):
        sdata = -sdata
        ax2.set_yscale("log", nonposy='clip')
    
    plt.plot(t, sdata)
    plt.xlabel('Time')
    plt.ylabel(varname2)
    lx = 0.05
    ly = 0.9
    plt.text(lx, ly, '(c) ', transform=ax2.transAxes)

# plot 4
if (nvars >= 4):
    sdata = np.zeros((nt,nl))
    for i in range(len(tt)):
        sdata[i,:] = vardata3[it[i]]

    ax3=plt.subplot(nrow, ncol, 4)
    if(varname3 == 'Liquid_Pressure'):
        sdata = -sdata
        ax3.set_yscale("log", nonposy='clip')
    
    plt.plot(t, sdata)
    plt.xlabel('Time')
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

