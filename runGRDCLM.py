#!/usr/bin/env python

import os, sys
import glob
from optparse import OptionParser

#DMR 4/16/13
#
# FMY 12/6/2015
# modified to work for global/regional simulation with PFLOTRAN coupled, used by NGEE-Arctic project

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--caseroot", dest="caseroot", default='./', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--jobname", dest="jobname", default="", \
                  help="userdefined job name, default is the casename")
parser.add_option("--runroot", dest="runroot", default="../runs", \
                  help="Directory where the run would be created")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  help = "input data directory for CESM (required)")
parser.add_option("--cesmdir", dest="cesmdir", default='..', \
                  help = "cesm directory (default = .., i.e., upper directory of this script)")
parser.add_option("--compset", dest="compset", default='I1850CLM45CN', \
                  help = "component set to use (required)")
parser.add_option("--res", dest="res", default = '', \
                  help = "grid resolution ---- \n"  
                         "default = '' \n"
                         "options = checking by ./create_newcase -list grid \n "
                         "NOTE: make sure the option you chose available or already prepared")
parser.add_option("--machine", dest="machine", default = '', \
                  help = "machine to use ---- \n"  
                         "default = '' \n"
                         "options = checking by ./create_newcase -list machines \n "
                         "NOTE: make sure the option you chose well-defined in config_machines.xml")
parser.add_option("--mach_specific", dest="mach_specific", default = '', \
                  help = "machine to use ---- \n"  
                         "default = '' \n"
                         "user-defined 'env_mach_specific' file to replace the pre-defined ")
parser.add_option("--osname", dest="osname", default = '', \
                  help = "name of machine OS ---- \n"
                         "   default = 'gnu', the default compiler for the chosen machine) \n "
                         "   options = intel,ibm, pgi,pathscale,gnu,cray,lahey,userdefined \n "
                         "NOTE: make sure the option you chose well-defined in config_compilers.xml")
parser.add_option("--compiler", dest="compiler", default = 'gnu', \
                  help = "compiler to use on machine ---- \n"
                         "   default = '', the default compiler for the chosen machine) \n "
                         "   options = intel,ibm, pgi,pathscale,gnu,cray,lahey,userdefined \n "
                         "NOTE: make sure the option you chose well-defined in config_compilers.xml")
parser.add_option("--mpilib", dest="mpilib", default = 'mpi-serial', \
                  help = "mpi library to use (default = mpi-serial)"
                         "options=openmpi,mpich,mpt,ibm,mpi-serial, BUT upon your system")
parser.add_option("--no_fire", dest="nofire", action="store_true", \
                  default=False, help="Turn off fire algorightms")
parser.add_option("--centbgc", dest="centbgc", default=False, \
                  help = 'To turn on CN with multiple soil layers, CENTURY C module', action="store_true")
parser.add_option("--vertsoilc", dest="vsoilc", default=False, \
                  help = 'To turn on CN with multiple soil layers (not necessary to bundle with centbgc)', action="store_true")
parser.add_option("--nitrif_denitrif", dest="kovenN", default=False, \
                  help = 'To turn on CN with Koven nitrif-denitrif (not necessary to bundle with centbgc)', action="store_true")
parser.add_option("--CH4", dest="CH4", default=False, \
                  help = 'To turn on CN with CLM4me (not necessary to bundle with centbgc)', action="store_true")
parser.add_option("--arcticpft", dest="arcticpft", default=False, \
                  help = 'To turn on Expanded Arctic PFTs flag (-DEXTENDED_PFT) in CLM4.5. Must provide --parm_file', action="store_true")
parser.add_option("--parm_file", dest="parm_file", default="", \
                  help = 'CLM user-defined physiological parameter file')
parser.add_option("--co2_file", dest="co2_file", default="fco2_datm_1765-2007_c100614.nc", \
                  help = 'CLM transient CO2 file for diagnostic option')
parser.add_option("--ad_spinup", action="store_true", \
                  dest="ad_spinup", default=False, \
                  help = 'Run accelerated decomposition spinup (note: exit-ad will run in the end as well')
parser.add_option("--nyears_ad_spinup", dest="ny_ad", default=600, \
                  help = 'number of years to run ad_spinup')
parser.add_option("--coldstart", dest="coldstart", default=False, \
                  help = "set cold start (mutually exclusive w/finidat)", \
                  action="store_true")
parser.add_option("--finidat_case", dest="finidat_case", default='', \
                  help = "case containing initial data file to use" \
                  +" (should be in your runroot directory)")
parser.add_option("--finidat", dest="finidat", default='', \
                  help = "initial data file to use" \
                  +" (should be in the runroot directory)")
parser.add_option("--finidat_year", dest="finidat_year", default=1, \
                  help = "model year of initial data file (default is" \
                  +" last available)")
parser.add_option("--run_units", dest="run_units", default='nyears', \
                  help = "run length units (ndays, nmonths, nyears, date)")
parser.add_option("--run_n", dest="run_n", default=600, \
                  help = "run length (in run units)")
parser.add_option("--restart_units", dest="restart_units", default='nyears', \
                  help = "restart length units for output restart file" \
                  + "(default: same as run_units)")
parser.add_option("--restart_n", dest="restart_n", default='', \
                  help = "restart length for outputfor output restart file" \
                  + "(default: same as run_n)")
parser.add_option("--rmold", dest="rmold", default=False, action="store_true", \
                  help = 'Remove old case directory with same name' \
                  +" before create/update new one")
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--hist_userdefined", dest="hist_file", default='', \
                  help = 'user=defined hist file')
parser.add_option("--clean_build", dest="clean_build", default=False, \
                  help = 'Perform clean build before building', \
                  action="store_true")
parser.add_option("--debug_build", dest="debug_build", default=False, \
                  help = 'Perform debug build', \
                  action="store_true")
parser.add_option("--clean_config", dest="clean_config", default=False, \
                  help = 'Perform clean setup before setting-up', \
                  action="store_true")
parser.add_option("--no_config", dest="no_config", default=False, \
                  help = 'do NOT configure case', action="store_true")
parser.add_option("--no_build", dest="no_build", default=False, \
                  help = 'do NOT build CESM', action="store_true")
parser.add_option("--no_submit", dest="no_submit", default=False, \
                  help = 'do NOT submit CESM to queue', action="store_true")
parser.add_option("--cleanlogs",dest="cleanlogs", help=\
                   "Removes temporary and log files that are created",\
                   default=False,action="store_true")
parser.add_option("--queue", dest="queue", default='', \
                  help = 'PBS submission queue')
parser.add_option("--np", dest="np", default=1, \
                  help = 'number of processors requested')
parser.add_option("--ppn", dest="ppn", default=1, \
                  help = 'processors per node, usually input with --np option above, with PBS -l nodes=np/ppn+1:ppn')
parser.add_option("--mppwidth", dest="mppwidth", default=1, \
                  help = 'processors range, which is --np option but using default ppn, e.g. on Hopper')
parser.add_option("--mppnodes", dest="mppnodes", default=1, \
                  help = 'number of nodes requested, with PBS -l nodes=mppnodes, e.g. on Titan')
parser.add_option("--walltime", dest="walltime", default="", \
                  help = 'user-requested walltime in format hh:mm:ss')
parser.add_option("--tstep", dest="tstep", default=0.5, \
                  help = 'CLM timestep (hours)')
parser.add_option("--surfdatafile", dest="surfdatafile", default="", \
                  help = 'CLM user-defined surface data file in inputdata/lnd/clm2/surfdata_map/ \n'
                  "the default file name IS: surfdata_?x?pt_sitename_simyr????.nc")
parser.add_option("--pftdynfile", dest="pftdynfile", default="", \
                  help = 'CLM user-defined pftdyn data file in inputdata/lnd/clm2/surfdata_map/ \n'
                  "the default file name IS: surfdata.pftdyn_?x?pt_sitename.nc")
#CLM-PFLOTRAN coupling options
parser.add_option("--clm_pflotran", action="store_true", dest="pflotran", default = False, \
                  help = "clm coupled with PFLOTRAAN (default = False, i.e. not coupled). Must provide --pflotran_dir")
parser.add_option("--clm_pf_colmode", action="store_true", dest="clmpfcolmode", default = False, \
                  help = "clm coupled with PFLOTRAAN by column-mode (default = False, i.e. fully in 3-D; otherwise, vertical-only )." )
parser.add_option("--pflotran_srcdir", dest="pflotran_srcdir", default = "", \
                  help = "PFLOTRAN source root directory, under which /src/clm-pflotran/libpflotran.a exists! ")
parser.add_option("--pflotran_indir", dest="pflotran_indir", default = "", \
                  help = "PFLOTRAN input root directory, under which input files for pflotran will be copied from (default: $cesminput/pflotran/ ")
parser.add_option("--petsc_dir", dest="petsc_dir", default = "", \
                  help = "PETSC_DIR petsc library root directory, if not set in env_specific_mach")
parser.add_option("--petsc_arch", dest="petsc_arch", default = "", \
                  help = "PETSC_ARCH petsc ARCH name, if not set in env_specific_mach")
#

(options, args) = parser.parse_args()

#------------------- arguments ------------------------------------------------------------

# current directory
runCLMdir = os.getcwd()

# cesm model directory
if (options.cesmdir==''):
    print('UNKNOWN cesm root directory: ')
    sys.exit()
else:
    csmdir = os.path.abspath(options.cesmdir)
    scriptsdir = csmdir+'/scripts'

# grid resolution
if (options.res==''):
    print('grid resolution option is required !')
    sys.exit()  

# machine options
if (options.machine==''):
    print('machine option is required !')
    sys.exit()
else:
    machineoptions = ' -mach '+options.machine
    if (options.machine == 'userdefined'):
        if (options.osname == '' or options.compiler == '' or options.mpilib == ''):
            print('"osname", "compiler", and "mpilib" options are required for " -mach userdefined"!')
            sys.exit()
        
        if (options.mach_specific == ''):
            print('please provide a specific mach name for setting up "env_mach_specific" if " -mach userdefined"! \n' + \
                  'OR, the script will look for "env_mach_spefic.'+options.osname+'_'+options.compiler+'"')        
            options.mach_specific = options.osname+'_'+options.compiler

if (options.compiler != ''):
    machineoptions += ' -compiler '+options.compiler

if (options.mpilib != ''):
    machineoptions += ' -mpilib '+options.mpilib

# case directory
if (options.caseroot == '' or \
(os.path.exists(options.caseroot) == False)):
    caseroot = csmdir+'/cases'
else:
    caseroot = os.path.abspath(options.caseroot)
print('CASE root directory: '+options.caseroot)

# case run root directory
if (options.runroot == '' or \
(os.path.exists(options.runroot) == False)):
    runroot = csmdir+"/runs"
else:
    runroot = os.path.abspath(options.runroot)
print('CASE RUN root directory: '+runroot)

#check for valid input data directory
if (options.ccsm_input == '' or \
(os.path.exists(options.ccsm_input) == False)):
    print('Error:  invalid input data directory')
    sys.exit()
else:
    ccsm_input = os.path.abspath(options.ccsm_input)

#check for valid compset
compset = options.compset
if (compset.startswith('I1850CLM45') == False \
    and compset.startswith('I1850CRUCLM45') == False \
    and compset.startswith('I20TRCLM45') == False \
    and compset.startswith('I20TRCRUCLM45') == False):
    print('Error:  please enter one of following valid options for compset:')
    print('        I1850(CRU)CLM45CN, I1850(CRU)CLM45BGC, I20TR(CRU)CLM45CN, I20TR(CRU)CLM45BGC')
    sys.exit()

# check if pflotran directory as an input
if (options.pflotran):
    if(options.pflotran_srcdir==""):
        print(" PFLOTRAN directories NOT defined, please provide one using '--pflotran_srcdir=' ! \n")
        sys.exit()
    elif(os.path.exists(options.pflotran_srcdir) == False):
        print(" PFLOTRAN directories NOT exist, please the directory provided using '--pflotran_srcdir=' ! \n")
        sys.exit()

#check consistency of options   
if (compset.startswith('I20TRCLM45') == True or compset.startswith('I20TRCRUCLM45') == True):
    #ignore spinup option if transient compset
    if (options.ad_spinup):
        print('Spinup options not available for transient compset.')
        sys.exit()
elif(options.ad_spinup):
    options.coldstart = True

if (options.arcticpft and options.parm_file == ''):  # must provide user-defined 'pft-physiology.???.nc'
    print('MUST provide user-defined parameter file! Exit \n')
    sys.exit()
              
#finidat file and finidat year
if (options.coldstart and (options.finidat != '' or options.finidat_case != '')):
        print('Error: Cannot have an finidat/finidat_case AND coldstart simultaneously! Exit \n')
        sys.exit()
    
if (options.finidat == '' and options.finidat_case == ''):   # not user-defined
    if (options.coldstart==False and compset.startswith('I1850')==True):
        if (options.mycaseid != ''):
            options.finidat_case = options.mycaseid+'_'+ \
                                  compset+'_ad_spinup'
        else:
            options.finidat_case = compset+'_ad_spinup'
    
        if (options.finidat_year == -1):
            options.finidat_year = int(options.ny_ad)+1
    
    if (compset.startswith('I20TRCLM45') == True):
               
        #finidat and finidat_year is required for transient compset
        if (os.path.exists(runroot+'/'+options.finidat_case) == False \
            or options.finidat_year == -1):
            print('Error:  must provide initial data file for I20TR(CRU)CLM45??? compset, OR, '+ \
                  runroot+'/'+options.finidat_case+' existed as refcase')
            sys.exit()
elif (options.finidat != ''):  # user-defined finidat file
    if (options.finidat.startswith('/')):  # full path and file names
        finidat = options.finidat
    else:  # otherwise, finidat is assummed under the $ccsm_input/lnd/clm2/inidata/
        finidat = ccsm_input+'/lnd/clm2/inidata/'+options.finidat
    
    if (options.finidat_year == -1):
        print('Error: must define the finidat_year if finidat defined! Exit \n')
        sys.exit()

finidat_year = int(options.finidat_year)
finidat_yst = str(finidat_year)
if (finidat_year >= 100 and finidat_year < 1000):
    finidat_yst = '0'+str(finidat_year)
if (finidat_year >= 10 and finidat_year < 100):
    finidat_yst = '00'+str(finidat_year)
if (finidat_year < 10):
    finidat_yst = '000'+str(finidat_year)

if (options.finidat_case != '' and options.finidat == ''):
    finidat = runroot+'/'+options.finidat_case+'/run/'+ \
        options.finidat_case+'.clm2.r.'+finidat_yst+'-01-01-00000.nc'

#get simyr for fixed surface data
mysimyr=1850
if (options.compset.startswith('I20TR') == True):
    mysimyr=2000

#----- Construct default casename
casename    = options.res+"_"+compset
if (options.mycaseid != ""):
    casename = options.mycaseid+'_'+casename
if (options.ad_spinup):
    casename = casename+'_ad_spinup'

#----- Construct case directory (clm4.5 case has two parts: directory+casename)
if (caseroot != "./"):
    casedir=caseroot+"/"+casename
else:
    casedir=casename
print ("CASE directory is: "+casedir+"\n")

if (os.path.exists(casedir)):
    print('Warning:  Case directory exists and --rmold not specified')
    var = raw_input('proceed (p), remove old (r), or exit (x)? ')
    if var[0] == 'r':
        os.system('rm -rf '+casedir)
    if var[0] == 'x':
        sys.exit()    

#construct case build and run directory
blddir=runroot+"/"+casename+'/bld'
print ("CASE bulid and exeroot is: "+blddir+"\n")
rundir=runroot+"/"+casename+"/run"
print ("CASE rundir is: "+rundir+"\n")

#pft parameter file
# default
pftfile = ccsm_input+'/lnd/clm2/paramdata/clm_params.'+casename+'.nc'
os.system('cp -f '+ccsm_input+'/lnd/clm2/paramdata/clm_params.c130821.nc ' \
              + pftfile)

# new or user-defined pft-phys file if desired
if (options.parm_file != ''):
    pftfile = ccsm_input+'/lnd/clm2/paramdata/' + \
                  options.parm_file.replace('.nc','') + '.' + casename + '.nc'
    os.system('cp -f '+ccsm_input+'/lnd/clm2/paramdata/'+ options.parm_file + \
               ' '+pftfile)

#set number of run years, if not user-defined
if (options.ny_ad != options.run_n and options.ad_spinup):
    options.run_n = options.ny_ad
if (compset.startswith('I20TRCLM45') == True and options.run_n == 600):    # 600 is the default (i.e., NOT user-defined)
    if (options.run_units == 'nyears'): 
        options.run_n = endyear - 1850 +1
    elif (options.run_units == 'date'):
        options.run_n = endyear + 1  # date in yyyymmdd, but here only needs year value

    
# set 'debug' mode
debugoption = ''
if (options.debug_build):
    debugoption = ' -confopts _D'
    print ("case build will be configured for debugging \n")

#----------------------------------------------------------------------------------------------------
if (os.path.isdir(scriptsdir)):

    os.chdir(scriptsdir)

# ------------------ IF no refcase, create, setup and build -------------
# (1) create a new case
    #create new case
    comps = options.compset
    if(compset == 'I20TRCLM45CN'):
        comps = 'I20TRCLM45BGC'   # this is a hack, since 'I20TRCLM45CN' is not default anymore
    
    print ('./create_newcase -case '+casedir+' '+machineoptions + \
                 ' -compset '+ comps +' -res '+options.res + \
                 debugoption)
    os.system('./create_newcase -case '+casedir+' '+machineoptions + \
                 ' -compset '+ comps +' -res '+options.res + \
                 debugoption + \
                  ' > create_newcase.log')
    if (os.path.isdir(casedir)):
        print(casename+' created.  See create_newcase.log for details')
        os.system('mv create_newcase.log ' + casedir +"/"+casename+"_case.log")
    else:
        print('failed to create case.  See create_newcase.log for details')

    # go to newly created case directory
    os.chdir(casedir)

# (2) env_build.xml modification ---------------------------
    
    # user-defined machine
    if (options.machine == "userdefined"):
        os.system('./xmlchange -file env_build.xml -id ' \
                  +'OS -val '+options.osname)
        
        os.system('./xmlchange -file env_build.xml -id ' \
                  +'COMPILER -val '+options.compiler)
         
        os.system('./xmlchange -file env_build.xml -id ' \
                  +'MPILIB -val '+options.mpilib)
        
        os.system('./xmlchange -file env_build.xml -id ' \
                  +'GMAKE -val make')

    if (options.runroot != '' or options.machine == "userdefined"):
        os.system('./xmlchange -file env_build.xml -id ' \
                  +'EXEROOT -val '+blddir) 
    
    # the following is new
    os.system('./xmlchange -file env_build.xml -id ' \
                  +'CESMSCRATCHROOT -val '+runroot) 

    # turn off rof module
    os.system('./xmlchange -file env_build.xml -id ' \
                  +'RTM_MODE -val NULL') 

    # clm4_5 config options (note: this configuration is re-designed, with some options moved to CLMNMLBLD)
    # i) basical options
    clmconfig_opts = "-phys clm4_5 -bgc cn"

    # ii) turn off fire module 
    if (options.nofire):
        clmconfig_opts += " -nofire"

    os.system('./xmlchange -file env_build.xml -id ' \
                  +'CLM_CONFIG_OPTS -val "'+clmconfig_opts+'"')
    print ("CLM configuration options: " + clmconfig_opts +"\n")

# (3) env_run.xml modification ------------------------------------
    # input/run/output directory
    if (options.runroot != '' or options.machine == "userdefined"):
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'RUNDIR -val '+rundir) 
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'DOUT_S -val TRUE') 
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'DOUT_S_ROOT -val '+runroot+'/archives/'+casename) 
    
    if (options.ccsm_input != '' or options.machine == "userdefined"):
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'DIN_LOC_ROOT -val '+ccsm_input) 
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'DIN_LOC_ROOT_CLMFORC -val '+ccsm_input+'/atm/datm7') 
    
    # datm options
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_MODE -val CLM1PT') 
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_START -val '+str(startyear))
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_END -val '+str(endyear))
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_ALIGN -val '+str(alignyear))

    # run timestep
    if (options.tstep != 0.5):
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'ATM_NCPL -val '+str(int(24/float(options.tstep))))
    
    # run-type adjusting -- needs checking ('rof' not working??) 
    if (options.ad_spinup==False and options.coldstart==False):
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_REFDATE -val '+finidat_yst+'-01-01')
   
     # starting date/time
    if (compset.startswith('I20TRCLM45') == True):
        #by default, 'I20TR' starting from 1850, but if not, then
        if(int(options.finidat_year) > 1850 and \
          (not ('I1850' in options.finidat or 'I1850' in options.finidat_case))):
            os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_STARTDATE -val '+finidat_yst+'-01-01')
        else:
            os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_STARTDATE -val 1850-01-01')                              
    else:
            os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_STARTDATE -val '+finidat_yst+'-01-01')

    #adds capability to run with transient CO2
    if (compset.startswith('I20TRCLM45') == True):
        os.system('./xmlchange -file env_run.xml -id ' \
                          +'CCSM_BGC -val CO2A')
        os.system('./xmlchange -file env_run.xml -id ' \
                          +'CLM_CO2_TYPE -val diagnostic')

    # user-defined running stop options
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_OPTION -val '+options.run_units)
    if (options.run_units == 'date'):
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_DATE -val '+str(options.run_n)+'0101')
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_N -val -9999')        
    else:
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_DATE -val -9999')
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_N -val '+str(options.run_n))        

    # user-defined restart options
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_OPTION -val '+options.restart_units)
    if (options.restart_n != ''):
        if (options.restart_units == 'date'):
            os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_DATE -val '+str(options.run_n)+'0101')
            os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_N -val -9999')
        else:
            os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_N -val '+str(options.restart_n))            
            os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_DATE -val -9999')
   
    # CLM build namelist options 
    # (NOTE: originally those options in configuration NOW moved here)
    stdout  = os.popen("./xmlquery -valonly -silent CLM_BLDNML_OPTS")
    clmbldnmls = stdout.read().rstrip( )   

    # i) CLM-CN or CLM-CENTURY
    if (options.centbgc):
        clmbldnmls += " -bgc bgc"   
        # this will turn-on the whole package of clm4me, nitrif-denitrif, century-decomp, vertical-soil
        print(" \n ======= CLM-BGC ==========================================" )
        print(" CLM bgc will be turned on with CENTURY BGC with clm4me, nitrif-denitrif ! " )
    else:
        clmbldnmls += " -bgc cn"
        print(" \n ======= CLM-CN ==========================================" )
        print(" CLM bgc will be turned on with classic CLM-CN ! " )
        # will hack 'vertical soil scheme' and 'clm4me' for CN later on    
               
    # ii) ad-spinup with exit-spinup included
    if (options.ad_spinup):
        clmbldnmls += " -bgc_spinup on"
        
    # iii) write into the xml valule
    os.system('./xmlchange -file env_run.xml -id ' \
                      +'CLM_BLDNML_OPTS -val "'+clmbldnmls+'"')
    print ("CLM namelist build options: " + clmbldnmls +"\n")

    #save intermit restart files (NOW, at most 3 sets are saved)    
    os.system('./xmlchange -file env_run.xml -id DOUT_S_SAVE_INTERIM_RESTART_FILES ' \
                  +'-val TRUE' )
    
        
# (4) env_mach_pes.xml modification ------------------------------------
    # normal 1 pt mode, not threaded
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 0')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val 0')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val 0')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val 0')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val 0')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val 0')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val 0')
    os.system('./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id TOTALPES -val 1')
    
    #if number of land instances > 1
    if (int(options.ninst) > 1):
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'NINST_LND -val '+options.ninst)
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'NTASKS_LND -val '+options.ninst)

    #if running with > 1 processor
    if(int(options.np)<=1):
        if(int(options.mppwidth)>1):     # this input is for HOPPER PBS
            options.np = options.mppwidth
            
        if( int(options.mppnodes)>1 or 
            (int(options.mppnodes)==1 and int(options.ppn)>1) ):     # this input is for Titan PBS
            options.np = str(int(options.mppnodes)*int(options.ppn))

    else:
        options.mppwidth = options.np
        options.mppnodes = str((int(options.np)-1)/int(options.ppn)+1)
        
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'NTASKS_ATM -val '+options.np)
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'NTASKS_LND -val '+options.np)
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'NTASKS_CPL -val '+options.np)
        if (options.machine == 'mira'):
            os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'MAX_TASKS_PER_NODE -val '+str(int(options.ppn)))
        else:
            os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'MAX_TASKS_PER_NODE -val '+str(min(int(options.ppn),int(options.np))))
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'TOTALPES -val '+options.np)    
                        
# (5) cesm setup -------------------------------------------
    #clean configure if requested prior to configure
    if (options.clean_config):
        os.system('./cesm_setup -clean')
        os.system('rm -f Macro')
        os.system('rm -f user-nl-*')

# (5a) user-defined machine settings ------
    if (options.no_config == False):
        if (options.mach_specific != ''):
            os.system('cp -f '+runCLMdir+'/userdefined_machines/env_mach_specific.'+options.mach_specific + \
                      ' env_mach_specific')
            os.system('cp -f '+runCLMdir+'/userdefined_machines/mkbatch.'+options.mach_specific + \
                      ' ./Tools/mkbatch.userdefined')
            if (os.path.isfile(runCLMdir+'/userdefined_machines/Macros.'+options.mach_specific)):
                os.system('cp -f '+runCLMdir+'/userdefined_machines/Macros.'+options.mach_specific + \
                      ' Macros')
            
# (5b) settings for clm coupled with pflotran, if requested ------
        if (options.pflotran):
            print(" \n ============== CLM-PFLOTRAN ===================================" )
            print(" NOTE: PFLOTRAN coupled CLM will be configured ! \n" )
            print(" make sure of libpflotran.a compiled, and exists in:\n")
            if (options.pflotran_srcdir !=""):
                print(" PFLOTRAN directory '" \
                      +options.pflotran_srcdir+"/src/clm-pflotran/'\n" )                             
                with open("env_mach_specific", "a") as myfile:
                    myfile.write("\n")
                    myfile.write("#pflotran settings\n")
                    myfile.write("setenv PFLOTRAN TRUE\n")
                    myfile.write("setenv PFLOTRAN_COUPLED_MODEL "+options.pflotran_srcdir+"\n")
                    if(options.petsc_dir!=""):
                        myfile.write("setenv PETSC_DIR "+options.petsc_dir+"\n")
                    if(options.petsc_arch!=""):
                        myfile.write("setenv PETSC_ARCH "+options.petsc_arch+"\n")
                    
                    if(os.path.exists(options.petsc_dir+"/lib")):
                        myfile.write("setenv PETSC_LIB ${PETSC_DIR}/lib \n")
                    else:
                        myfile.write("setenv PETSC_LIB ${PETSC_DIR}/${PETSC_ARCH}/lib \n")

                    myfile.write("\n")
                myfile.close()
            else:
                print("PFLOTRAN directory NOT defined! Exit!")
                sys.exit()

            if(options.clmpfcolmode):
                print(" CLM-PFLOTRAN coupled in COLUMN_MODE ! \n" )
                with open("env_mach_specific", "a") as myfile:
                    myfile.write("\n")
                    myfile.write("#CLM-PFLOTRAN coupled in COLUMN_MODE\n") 
                    myfile.write("setenv COLUMN_MODE TRUE\n")
                    myfile.write("\n") 
                    myfile.close()
                
# (5c) flags for option to turn on expanded arctic PFTs, if requested ------
        if (options.arcticpft):
            print(" \n ======= EXTENDED-PFT ==========================================" )
            print(" \n Expanded PFTs for Arctic Tundra in CLM4.5 will be turned on ! " )
            print(" NOTE: make sure of using right CLM4.5 code ! \n" )
            with open("env_mach_specific", "a") as myfile:
                myfile.write("\n")
                myfile.write("#expanded arctic PFT switch for CLM4.5\n") 
                myfile.write("setenv EXTENDED_PFT TRUE\n")
                myfile.write("\n") 
                myfile.close()

# (5e) setup ---------
        os.system('./cesm_setup > configure.log')
        
    else:
        print("Warning:  No case configure performed.  GRDCLM will not " \
                  +"make any requested modifications to env_*.xml files.  Exiting.")
        sys.exit()    

# (6) clm user-defined namelist modification ('user_nl_clm') -----
    output = open("user_nl_clm",'w')
    output.write('&clm_inparm\n')
        
    #(6a) user-defined initial data file ---
    if (options.coldstart == False):
        if (finidat != ''):
            output.write(" finidat = '"+finidat+"'\n")
    
    #(6b) surfdata - user-defined ----
    if(options.surfdatafile == ""):
        output.write(" fsurdat = '"+ccsm_input+"/lnd/clm2/surfdata_map/" + \
                 options.surfdatafile+"'\n")
    
    #(6c) pft dynamics file for transient run (TODO: name changed) ----
    if (compset.startswith('I20TRCLM45') == True):
        if(options.pftdynfile == ""):
            output.write(" flanduse_timeseries = '"+ccsm_input+"/lnd/clm2/surfdata_map/" + \
                          "surfdata.pftdyn"+str(numxpts)+"x"+str(numypts)+"pt_"+ \
                          options.site+"_hist_simyr1850-2005_clm4_5_pftgrd_c140204.nc'\n")
        else:
            output.write(" flanduse_timeseries = '"+ccsm_input+"/lnd/clm2/surfdata_map/" + \
                          options.pftdynfile+"'\n")
    
    #(6d) user-defined pft physiological file ----
    if (pftfile != ''):
        output.write(" paramfile=   '" + pftfile + "'\n")
        
    #(6e) clm output hist user-defined file ----
    if (options.hist_file != ''):
        histfile = runCLMdir+"/userdefined_outputs/"+options.hist_file
        hvars_file = open(histfile)
        output.write("\n")
        for s2 in hvars_file:
            myline = s2
            output.write(myline)
        output.write("\n")
        hvars_file.close()

    # (6f) flags for option to turn on 'VERTISOIL' and '-LCH4' for CLM45CN without CENTURY BGC, if requested ------
    # this is because configuration of 'VERTISOILC' and 'CLM4ME' is bundled with centurary BGC
    if (options.vsoilc and options.centbgc==False):
        print(" \n ======= VERTSOILC ==========================================" )
        print(" CLM-CN+ classic bgc will be turned on with vertical soil profiling! " )
        output.write(" use_vertsoilc = .true. \n")

    if (options.kovenN and options.centbgc==False):
        print(" \n ======= NITRIF-DENITRIF ====================================" )
        print(" CLM-CN+ classic bgc will be turned on with Koven' nitrif_denitrif option! " )
        output.write(" use_nitrif_denitrif = .true. \n")
    elif (options.centbgc==False):
        print(" CLM-CN+ classic bgc will be with simple inorganic N cycle modules ! " )            

    if (options.CH4 and options.centbgc==False):
        print(" \n ======= CH4 ==========================================" )
        print(" CLM-CN+ classic bgc will be turned on with clm4me option! " )
        output.write(" use_lch4 = .true. \n")

    #(6g) namelist options for PFLOTRAN coupling ----
    if (options.pflotran):
        output.write(" use_pflotran = .true.\n")
        output.write("/\n")

        output.write("&clm_pflotran_inparm\n")
        output.write(" pflotran_prefix = '"+ casename + "'\n")
        output.write("/\n")
    
    output.close()

#(7) copy user-defined sourcemods codes  ----
    if (options.srcmods_loc != ''):
        if (options.srcmods_loc.startswith('/')):
            options.srcmods_loc = os.path.abspath(options.srcmods_loc)
        else:   
            options.srcmods_loc = os.path.abspath(runCLMdir+'/'+options.srcmods_loc)
            
        if (os.path.exists(options.srcmods_loc) == False):
            print('Invalid srcmods directory.  Exiting')
            sys.exit()
        else:
            print('User-defined source codes will be copied from: '+options.srcmods_loc)
        os.system('cp -rf '+options.srcmods_loc+'/* ./SourceMods')   
   
# (8) datm namelist modifications (cycle met data streams - a bug since clm4.5.10) ---
# the issue is that climate data will not repeating (taxmode is 'extend'). 
    myoutput = open('user_nl_datm','w')
    myoutput.write("&shr_strdata_nml\n")
    if (compset.startswith('I20TRCLM45') == True):             
        myoutput.write(" taxmode = 'cycle', 'extend','extend'\n")
    else:
        myoutput.write(" taxmode = 'cycle', 'extend'\n")
    myoutput.write("/\n")
    myoutput.close()      

# (9) ------- build clm45 within cesm ---------------------------------------------- 
    #clean build if requested prior to build
    if (options.clean_build):
        os.system('./'+casename+'.clean_build')
        os.system('rm -rf '+rundir+'/*')
        os.system('rm -rf '+blddir+'/*')
            
    #compile cesm
    if (options.no_build == False):        
        os.system('./'+casename+'.build')        
        
        # if pflotran coupled, need to copy input files
        if (options.pflotran):
            pfindir = ""
            if (options.pflotran_indir.startswith('/')):  # full path
                pfindir = options.pflotran_indir
            else:
                if(options.pflotran_indir=="" or options.pflotran_indir=="./"):
                    pfindir = options.ccsm_input+'/pflotran/default'           
                else:
                    pfindir = options.ccsm_input+'/pflotran/'+options.pflotran_indir           
            
            if(os.path.isfile(pfindir+'/pflotran_clm.in')):
                os.system('cp '+pfindir+'/pflotran_clm.in '+rundir+'/'+casename+'.in')
            else:
                print('Error: must provide a "pflotran_clm.in" for CLM-PFLOTRAN in '+pfindir+'! Exit \n')
                sys.exit()
                
            if(glob.glob(pfindir+'/*.meshmap')):
                os.system('cp '+pfindir+'/*.meshmap '+rundir+'/')
            else:
                print('Error: must provide a set of "*.meshmap" for CLM-PFLOTRAN in '+pfindir+'! Exit \n')
                sys.exit()

            if(os.path.isfile(pfindir+'/hanford-clm.dat')):
                os.system('cp '+pfindir+'/hanford-clm.dat '+rundir+'/')         
            if(os.path.isfile(pfindir+'/CLM-CN_database.dat')):
                os.system('cp '+pfindir+'/CLM-CN_database.dat '+rundir+'/')
            else:
                if(os.path.isfile(pfindir+'/hanford-clm.dat') == False):
                    print('Waring: NO PFLOTRAN-bgc database file "handford-clm.dat" or "CLM-CN_database.dat" exists! \n')
            
            if(glob.glob(pfindir+'/*.h5')):      
                os.system('cp '+pfindir+'/*.h5 '+rundir+'/')
            else:
                print('Waring: NO PFLOTRAN *.h5 input file exists! -- be sure it is the case! \n')
            
           # decomposing pflotran mesh horizontally exactly as CLM  (Something NOT right here - temporarily OFF)
           # myfile = open(rundir+'/'+casename+'.in','r')
           # tmpfile = open(rundir+'/pfclm.in','w')
           # proc = False
           # for s in myfile:
           #     s2 = s.strip()
           #     if (s2.startswith('PROC')):
           #        tmpfile.write('PROC 1 '+options.np+' 1\n')
           #         proc = True
           #     else:
           #         tmpfile.write(s)
           # 
           # if(proc==False):
           #     tmpfile.write('PROC 1 '+options.np+' 1\n')      # not defined in pflotran-clm.in
           #     
           # myfile.close()
           # tmpfile.close()
           # os.system('mv '+rundir+'/pfclm.in '+ rundir+'/'+casename+'.in')
                  
# ---------------------------- Reference case set ------------------------------------------
# the following has not yet checked for CLM4.5
else:  
    print ('Not supported yet !')
    sys.exit()

# --------------------------- end of refcase ------------------------------------------------
    
    
# ----- copy rpointers and restart files to current run directory prior to run model ---
if (options.finidat_case != ''):
    os.system('cp -f '+runroot+'/'+options.finidat_case+'/run/' + \
              options.finidat_case+'.*'+finidat_yst+'* ' + rundir)
    os.system('cp -f '+runroot+'/'+options.finidat_case+'/run/'+ \
              'rpointer.* '+rundir)

# -------- make necessary modificaitons to run script for user-defined submit options ------------------------------
os.chdir(casedir)
myinput  = open("./"+casename+".run")
myoutput = open("./"+casename+".temprun",'w')
for s in myinput:  
    if (s[0:8]  == '#PBS -N '):
        if(options.jobname != ''):
            myoutput.write("#PBS -N "+options.jobname+"\n")
        else:
            myoutput.write(s)

    elif s[0:8] == '#PBS -q ':
        if(options.queue != ''):
            myoutput.write("#PBS -q "+options.queue+"\n")
        else:
            myoutput.write(s)
            
    elif s[0:14] =="#PBS -l nodes=":   
        if( (options.np !=1 or options.ppn != 1) 
            and (s.find(":ppn="))!=-1 ):   # general linux (?) there is a 'ppn' in PBS -l nodes= option
            myoutput.write("#PBS -l nodes="+str( (int(options.np)-1)/int(options.ppn)+1) + \
                                 ":ppn="+str(min(int(options.np),int(options.ppn)))+"\n")
        elif(options.mppnodes != 1):  # titan CNL - np and ppn not an option, rather using 'mppnodes'
            myoutput.write("#PBS -l nodes="+str(int(options.mppnodes))+"\n")
        else:
            myoutput.write(s)  

    elif s[0:17] == "#PBS -l mppwidth=": # for CNL with default 'ppn' on hopper
        myoutput.write("#PBS -l mppwidth="+str(int(options.mppwidth))+"\n")  
    
    elif s[0:17] == '#PBS -l walltime=':
        if(options.walltime != ''):
            myoutput.write("#PBS -l walltime="+options.walltime+"\n")
        else:
            myoutput.write(s)
    
    #elif s[0:16] == '#mpirun -np ':    # for linux
    #    myoutput.write("mpirun -np "+str(options.np)+" --hostfile $PBS_NODEFILE ./ccsm.exe >&! ccsm.log.$LID\n")
    
    #elif s[0:11] == '#aprun -n ':    # for CNL (TODO)
    #    myoutput.write(s)  
    
    else:
        myoutput.write(s)
    
myoutput.close()
myinput.close()
os.system("mv "+casename+".temprun "+casename+".run")  

# ----- submit job if requested ---
if (options.no_submit == False):    
    
    os.chdir(casedir)
    stdout  = os.popen("which qsub")
    stdout_val = stdout.read().rstrip( )   
    if (stdout_val == ""):
        os.system("csh ./"+casename+".run")
    else:   
        os.system("qsub ./"+casename+".run")

