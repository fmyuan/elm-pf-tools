#!/usr/bin/env python

#----------------------
# Converting from ELM outputs to CMIP format
#

import os, sys, time, math
from optparse import OptionParser
import re
from netCDF4 import Dataset

# the nc file variable(s) always to be included together with specific variable for concatentation, if existed
# NOTE: 'time','time_bounds' will automatically be included for each file.
incl_vars = ['lat','lon','mcdate','dtime']

#-------------------Parse options-----------------------------------------------
parser = OptionParser()

parser.add_option("--outdir", dest="outdir", default="", \
                  help = 'output directory')
parser.add_option("--ofheader", dest="header", default="none", \
                  help = 'header of output files either in fullpath or filename in default path')
parser.add_option("--ELM1var_Concat", dest="ELM1var", default=False, \
                  help = 'To indicate ELM outputs concatented into a single nc file', action="store_true")
parser.add_option("--cmorncheader", dest="cmorlist", default="", \
                      help = "CMOR variable NC header, containing basic information on conversion")
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
if(not ncopath.endswith('/')): ncopath = ncopath.strip()+'/'

###########################################################################
SINGLE_VAR_NCFILE = False
if(options.ELM1var): SINGLE_VAR_NCFILE = True

if(options.outdir.startswith('/')):
    outdir = os.path.abspath(options.outdir)
    if(not outdir.endswith('/')): outdir = outdir+'/'
else:
    outdir = './'

#------------get all nc output files if it's what needed -----------------------------------
filehead = options.header
if (os.path.isdir(outdir.strip())):
    #checking where is the original data
    dirfiles = os.listdir(outdir)
    ncfiles = []
    for dirfile in dirfiles:                
        if(dirfile.startswith(filehead) and dirfile.endswith('.nc')): 
            ncfiles.append(dirfile)

    #
    ncfiles.sort(key = lambda x: x.split('.')[-2])
                
#------------get all CMOR variable NC header files-----------------------------------
ftxt = []
if(options.cmorlist!=''):  # varlist txt file is checked first, so over-ride 'options.varlist'
    if (os.path.isdir(options.cmorlist.strip())):
        fdir = options.cmorlist.strip()
        if(not fdir.endswith('/')):
            fdir = fdir + '/'
        # a bunch of ascii files in the directory
        flist=os.listdir(fdir)
        for ifile in flist:
            if(ifile.endswith('.txt')): ftxt.append(fdir+ifile)
    
    elif (os.path.isfile(options.cmorlist.strip())):
        # a single ascii file
        ftxt.append(options.cmorlist.strip())
        
    else:
        print ('\n '+options.cmorlist+'is not a txt file or directory to contain CMOR variable NC header')
        sys.exit()
        
          
else:
    print ('\n MUST provide a txt file or directory to contain CMOR variable NC header')
    sys.exit()


vars_name={}
vars_att={}
for i in range(0,len(ftxt)):
    vars_att[i]=[]
    vars_name[i]=[]
    with open(ftxt[i]) as f:
        lns = f.read().splitlines()
        lns = list(filter(None,lns))
        for l in lns:
            l = l.strip()
            if (':' in l):
                l = re.split('[:;]+',l,1)
                vars_att[i].append(l[1].strip()) # [0] is the variable name
            else:
                l = re.split('[,;\s\(\)]+',l)
                l = list(filter(None,l))
                vars_name[i].append(l[1].strip()) # the [0] is the type
                
                print('\nCMOR-var: '+str(vars_name[i]))
              
    f.close()
              
    
#------------pick and merge vars in the NC file list -----------------------------------
CMOR_varname={} # dict for holding informations on CMOR variables actually available
CMOR_varatt={}
CMOR_ncfile={}
n_var = -1

#
if (SINGLE_VAR_NCFILE):
    # each variable has already been concatenated into a single nc file
    #-------------------------------------------------------------
    for ivar in range(0,len(vars_name)):
        var_new = vars_name[ivar]
        var_new = var_new[0] #list to one single string
                         
        for iatt in vars_att[ivar]:
            if ('original_name' in iatt): 
            # ELM's output variable name(s), 
            #  e.g. tsl:original_name = "TSOI" ;
            #       mrso:original_name = "SOILICE, SOILLIQ" ;
                varname_orig = re.split('[=,;\s\"]+',iatt.strip())
                varname_orig = varname_orig[1:]
                varname_orig = list(filter(None,varname_orig))
                                                 
        #-------------------------------------------------------------
        # from 1 or more ELM variable --> 1 CMOR variable
        for var_old in varname_orig:
            for nf in ncfiles:                
                outfile_old = outdir+'/'+nf
                f = Dataset(outfile_old,'r')
                allvars = f.variables.keys()
                f.close()
            
                if (var_old in allvars):
                    if (var_old==varname_orig[0]):
                        n_var = n_var + 1
                        
                        outfile_new = './CMOR_'+var_new+'_'+nf
                        if(var_old in outfile_new):
                            outfile_new=outfile_new.replace(var_old, '')
                            
                        print('\n Copying file: '+nf + '---> '+ outfile_new)

                        CMOR_varname[n_var]=var_new
                        CMOR_varatt[n_var]=vars_att[ivar]
                        CMOR_ncfile[n_var]=outfile_new # save the new files for using next-step. NOTE: must be in exactly-order as 'ivar' in vars_name and vars_att
            
                        # (1) if the 1st ELM variable --> 1 CMOR variable, copying into the 'outfile_new'
                        vstring = ' -v time,'+var_old
                        for v in incl_vars:
                            if (v not in vstring and v in allvars): 
                                vstring += (','+str(v).strip())
            
                        os.system(ncopath+'ncks -O '+vstring + ' '+ outfile_old + ' -o '+ outfile_new) 
                                      
                    else:
                        # (2) if more than 1 ELM variable --> 1 CMOR variable, appending into the 'outfile_new'
                        print('\n Appending file: '+nf + '---> '+ outfile_new)
                        os.system(ncopath+'ncks --append -v '+ var_old+' ' + outfile_old + ' -o '+ outfile_new)                    
                    
                    # end if-block of 'var_old == varname_orig[0]'
                    break # the for-loop of ncfile look-up (because it's not in other file)
                
                #end if-block of 'var_old in allvars in ncfile'    
            # end of for-loop of ncfile look-up
        
        #end of for-loop of 'var_old in varname_orig'
    
    #end of for-loop of 'ivar in range(0, len(varname))'

else:
    # one ELM nc file contains all variables, which may be in one nc file or multiple file for each time-period
    
    for nf in ncfiles:            
    
        outfile_old = outdir+'/'+nf

        f = Dataset(outfile_old,'r')
        allvars = f.variables.keys()
        f.close()
                   
        #-------------------------------------------------------------
        for ivar in range(0,len(vars_name)):
            var_new = vars_name[ivar]
            var_new = var_new[0] #list to one single string
         
            outfile_new = './CMOR_'+var_new+'_'+nf
                
            for iatt in vars_att[ivar]:
                if ('original_name' in iatt): 
                    # ELM's output variable name(s), 
                    #  e.g. tsl:original_name = "TSOI" ;
                    #       mrso:original_name = "SOILICE, SOILLIQ" ;
                    varname_orig = re.split('[=,;\s\"]+',iatt.strip())
                    varname_orig = varname_orig[1:]
                    varname_orig = list(filter(None,varname_orig))
                                                 
            #-------------------------------------------------------------
            # from 1 or more ELM variable --> 1 CMOR variable
            YESinNCFILE = False
            if(len(varname_orig)>=1):              
                #checking if any is in the nc file
                for v in varname_orig:
                    if (v in allvars):
                        YESinNCFILE = True
                        break
            else:
                continue
           
            # if variable not in current ELM nc file, jump to the next file
            if(not YESinNCFILE): 
                continue
            else:
                n_var = n_var + 1
                print('\n Now Processing file: '+nf + '---> '+ outfile_new)

                CMOR_varname[n_var]=var_new
                CMOR_varatt[n_var]=vars_att[ivar]
                CMOR_ncfile[n_var]=outfile_new # save the new files for using next-step. NOTE: must be in exactly-order as 'ivar' in vars_name and vars_att
       
            #-------------------------------------------------------------
            for var_old in varname_orig:
                if (var_old == varname_orig[0]):
                    #(1) copy the original ELM variable for the first one, with necessary timing/lat/lon variables 
                    vstring = ' -v time,'+var_old
                    for v in incl_vars:
                        if (v not in vstring and v in allvars): 
                            vstring += (','+str(v).strip())
            
                    os.system(ncopath+'ncks -O '+vstring + ' '+ outfile_old + ' -o '+ outfile_new) 
                                      
                elif(var_old in allvars):
                    # (2) if more than 1 ELM variable --> 1 CMOR variable, appending into the 'outfile_new'
                    os.system(ncopath+'ncks --append -v '+ var_old+' ' + outfile_old + ' -o '+ outfile_new)                                                     
                # end of if-block                
            #end of for-loop (var_old in varname_orig)
        
        # (END) for-loop for all vars (new) in CMOR nc header files                   
    #(END) for-loop of all original NC files in outdir directory
    
#(END) if-block (SINLGE_VAR_NCFILE)
           
# a little bit clean-up to avoid confusion
del allvars
del vars_att
del vars_name
del outfile_new
del outfile_old
del var_new
del var_old
del varname_orig
                                                  
#------------any arithmatic operations in the CMOR NC files -----------------------------------

for ivar in range(0,len(CMOR_varname)):
    var_new = CMOR_varname[ivar]
        
    # multiplier for converting from ELM variable to CMOR
    multiplier = 1.0
        
    # summed over soil layers, if required
    summed_soillayers = False
        
    for iatt in CMOR_varatt[ivar]:
        if ('original_name' in iatt): 
            # ELM's output variable name(s), 
            #  e.g. tsl:original_name = "TSOI" ;
            #       mrso:original_name = "SOILICE, SOILLIQ" ;
            varname_orig = re.split('[=,;\s\"]+',iatt.strip())
            varname_orig = varname_orig[1:]
            varname_orig = list(filter(None,varname_orig))
                
        elif ('comment' in iatt):
            # in CMOR nc headers, arithematic operations are described in variable attribute: comment
            # e.g. evspsbl:comment = "QSOIL+QVEGE+QVEGT" ; or hfls:comment = "EFLX_LH_TOT unchanged" ;
            if('convert g to kg' in iatt): 
                multiplier = 0.001
                
            if('reversed sign' in iatt): 
                multiplier = multiplier*(-1.0)

            if('convert to percentage' in iatt): 
                multiplier = 100.0
            
            if('summed over soil layer depths' in iatt):
                #mrfso:comment = "SOILICE summed over soil layer depths" ;
                summed_soillayers = True

    # multiplier for a specific variable: burntArea
    if(var_new=='burntArea'):
    #burntArea:comment = "FAREA_BURNED times days per month and 86400*100 to convert to percentage" ;
        multiplier=365.0/12.0*86400.0*100.0 # not exactly (TODO)
                      
    # multiplier for a specific variable: snw
    if(var_new=='snw'):
        # convert from H2OSNO in kg/m2 to meters water
        multiplier=0.001
                  
    #-------------------------------------------------------------
    var_old = varname_orig[0]
    
    #-------------------------------------------------------------
    outfile_new = CMOR_ncfile[ivar]      
    f = Dataset(outfile_new,'r')
    allvars = f.variables.keys()
    f.close()
        
    print('\n Now Arithematic Operation of File: '+ outfile_new)

    # add-up operation.
    # NOTE: if other than 'sum' of involved variables, please modify this (TODO)    
    for iv in range(1,len(varname_orig)):
        vstring = var_old + '='+var_old+'+'+varname_orig[iv]
        os.system(ncopath+"ncap2 --append -s '"+ vstring +";' "+ outfile_new + " -o "+ outfile_new)
        # remove no-more-needed original variable
        os.system(ncopath+'ncks -O -x -v '+ varname_orig[iv] +' '+outfile_new + ' -o '+ outfile_new)             
        # (END) for-loop for vars from 1 ~ end in varname_orig

    # rename variable, after clean-up original variable(s)
    os.system(ncopath+'ncrename -O -v '+var_old+','+var_new+' '+ outfile_new + ' -o '+ outfile_new)    
    
    # because the final nc file is to be used by ILAMB, have to switch timing variable, if any
    if('dtime' in allvars):
        # in this case, 'dtime' is the original 'time' (days since the simulation year first day)
        # but when doing concatentating, it may be renamed as 'dtime' while 'mcdate' renamed as 'time' for using in VISIT
        if('time' in allvars): # do this first
            os.system(ncopath+'ncrename -O -v time,mcdate ' + outfile_new + ' -o '+ outfile_new)    
        os.system(ncopath+'ncrename -O -v dtime,time '+ outfile_new + ' -o '+ outfile_new)    
            
              
    # summed over soil layers, dim of 'levgrnd' usually the second slice of 'time,levgrnd,lat,lon'
    if(summed_soillayers):
        os.system(ncopath+"ncwa -O -a levgrnd " + outfile_new + " -o "+ outfile_new) 
        vstring = var_new + '='+var_new+'*15.0'
        os.system(ncopath+"ncap2 --append -s '"+ vstring +";' "+ outfile_new + " -o "+ outfile_new)
        
    # multiplier, if any other than 1
    if(multiplier != 1.0):
        vstring = var_new + '='+var_new+'*'+str(multiplier)
        os.system(ncopath+"ncap2 --append -s '"+ vstring +";' "+ outfile_new + " -o "+ outfile_new)

    # rename/add attributes for varialbe
    for iatt in CMOR_varatt[ivar]:
        att=re.split('[=;]+',iatt)
        att=list(filter(None,att))
        
        att_name=att[0].strip()
        att_val = att[1].strip()
        if('\"' in att_val):
            att_type='c'
        elif('.' in att_val):
            if(att_val.endswith('f')):
                att_type='f'
                att_val = float(att_val.replace('f',''))
            else:
                att_type='d'
                att_val = float(att_val)
            
        else:
            att_type='i'
            att_val = int(att_val)
    
        #delete old attribute if any
        att_string = " -a "+att_name+","+var_new+",d,, "+outfile_new
        os.system(ncopath+"ncatted "+att_string) 
                
        #always create a new one
        ed_mode = 'c'
        att_string =  att_name+','+var_new+','+ed_mode+','+att_type+','+str(att_val) #att_name,var_name,mode,att_type,att_value        
        att_string = " -a "+att_string+" "+outfile_new
        os.system(ncopath+"ncatted -O "+att_string) 
    #(END) for-loop of renaming attributes
    
                
#(END) for-loop of 'var_new in vars_name of CMOR' 
                                        
#------------END of CLM_outvars_ToCMOR -----------------------------------------------------------------------
