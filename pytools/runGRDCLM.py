#!/usr/bin/env python

import os, sys
import glob
from optparse import OptionParser

#DMR 4/16/13
#python ./runPTCLM.py does the following:
#  1. Call routines to create point data (makepointdata.py, makemetdata.py) - removed
#  2. Set point and case-specific namelist options
#  2. configure case
#  3. build (compile) CESM with clean_build first if requested
#  4. apply patch for transient CO2 if transient run
#  6. apply user-specified PBS and submit information - removed
#  7. submit job to PBS queue if requested.
#
# FMY 12/6/2015
# modified to work for global/regional simulation with PFLOTRAN coupled, used by NGEE-Arctic project

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--caseroot", dest="caseroot", default='./', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--runroot", dest="runroot", default="../runs", \
                  help="Directory where the run would be created")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  help = "input data directory for CESM (required), but embedded in config_machines.xml for specified machine")
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
parser.add_option("--compiler", dest="compiler", default = '', \
                  help = "compiler to use on machine ---- \n"
                         "   default = '', the default compiler for the chosen machine) \n "
                         "   options = intel,ibm, pgi,pathscale,gnu,cray,lahey, ... \n "
                         "NOTE: make sure the option you chose well-defined in config_compilers.xml")
parser.add_option("--mpilib", dest="mpilib", default = 'mpi-serial', \
                  help = "mpi library to use (default = mpi-serial)"
                         "options=openmpi,mpich,mpt,ibm,mpi-serial, BUT upon your system")
#compsets and module options
parser.add_option("--no_fire", dest="nofire", action="store_true", \
                  default=False, help="Turn off fire algorightms")
parser.add_option("--centbgc", dest="centbgc", default=False, \
                  help = 'To turn on CN with multiple soil layers, CENTURY C module', action="store_true")
parser.add_option("--nitrif_denitrif", dest="kovenN", default=False, \
                  help = 'To turn on CN with Koven nitrif-denitrif (not necessary to bundle with centbgc)', action="store_true")
parser.add_option("--CH4", dest="CH4", default=False, \
                  help = 'To turn on CN with CLM4me (not necessary to bundle with centbgc)', action="store_true")
#parser.add_option("--extended_pft", dest="extended_pft", default=False, \
#                  help = 'To turn on Expanded (Arctic) PFTs flag (-DEXTENDED_PFT) in CLM. Must provide --parm_file', action="store_true")
parser.add_option("--parm_file", dest="parm_file", default="clm_params_c160822.nc", \
                  help = 'CLM user-defined physiological parameter file, with default: clm_params.c160822.nc')
parser.add_option("--co2_file", dest="co2_file", default="fco2_datm_1765-2007_c100614.nc", \
                  help = 'CLM transient CO2 file for diagnostic option')
parser.add_option("--ad_spinup", action="store_true", \
                  dest="ad_spinup", default=False, \
                  help = 'Run accelerated decomposition spinup (note: exit-ad will run in the end as well)')
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
parser.add_option("--restart_units", dest="restart_units", default='', \
                  help = "restart length units for output restart file" \
                  + "(default: same as run_units)")
parser.add_option("--restart_n", dest="restart_n", default='', \
                  help = "restart length for outputfor output restart file" \
                  + "(default: same as run_n)")
#user-defined codes/output lists
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--hist_userdefined", dest="hist_file", default='', \
                  help = 'user=defined hist file')
#build options
parser.add_option("--rmold", dest="rmold", default=False, action="store_true", \
                  help = 'Remove old case directory with same name' \
                  +" before create/update new one")
parser.add_option("--clean_build", dest="clean_build", default=False, \
                  help = 'Perform clean build before building', \
                  action="store_true")
parser.add_option("--debug_build", dest="debug_build", default=False, \
                  help = 'Perform debug build', \
                  action="store_true")
parser.add_option("--clean_config", dest="clean_config", default=False, \
                  help = 'Perform clean setup before setting-up', \
                  action="store_true")
parser.add_option("--cleanlogs",dest="cleanlogs", help=\
                   "Removes temporary and log files that are created",\
                   default=False,action="store_true")
#submit/run options
parser.add_option("--no_submit", dest="no_submit", default=False, \
                  help = 'do NOT submit CESM to queue', action="store_true")
parser.add_option("--jobname", dest="jobname", default="", \
                  help="userdefined job name, default is the casename")
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
#user-defined datm, if '--res=CLM_USR_DAT' or empty
parser.add_option("--align_year", dest="align_year", default="", \
                  help = 'Alignment year for datm data, if --res=CLM_USR_DAT')
parser.add_option("--start_year", dest="start_year", default="", \
                  help = 'Starting year for datm data, if --res=CLM_USR_DAT')
parser.add_option("--end_year", dest="end_year", default="", \
                  help = 'Ending year for datm data, if --res=CLM_USR_DAT')
parser.add_option("--usrdat-name", dest="usr_dat", default="", \
                  help = 'CLM_USR_DAT name for datm data and surf_data, if --res=CLM_USR_DAT')



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
    scriptsdir = csmdir+'/cime/scripts'

# grid resolution
if (options.res==''):
    print('grid resolution option is required !')
    sys.exit()  

if (options.machine==''):
    print('machine option is required !')
    sys.exit()
else:
    machineoptions += ' -mach '+options.machine
    if (options.machine == 'userdefined'):
        print('WARNING: must manually edit env_case.xml, env_build.xml, env_run.xml, env_mach_pes.xml ...! \n')        

if (options.compiler != ''):
    machineoptions += ' -compiler '+options.compiler
else:
    print('compiler is required !')
    sys.exit()

if (options.mpilib != ''):
    machineoptions += ' -mpilib '+options.mpilib
else:
    machineoptions += ' -mpilib mpi-serial'


#check for valid compset
compset = options.compset
if (compset.startswith('I1850CLM45') == False \
    and compset.startswith('I1850CRUCLM45') == False \
    and compset.startswith('I20TRCLM45') == False \
    and compset.startswith('I20TRCRUCLM45') == False):
    print('Error:  please enter one of following valid options for compset:')
    print('        I1850(CRU)CLM45CN, I1850(CRU)CLM45BGC, I20TR(CRU)CLM45CN, I20TR(CRU)CLM45BGC')
    sys.exit()

if (compset.startswith('I20TR') == True):
    #ignore spinup option if transient compset
    if (options.ad_spinup):
        print('Spinup options not available for transient compset.')
        sys.exit()
elif(options.ad_spinup):
    options.coldstart = True

# cases root directory
if (options.caseroot == ''):
    caseroot = csmdir+'/cases'
else:
    caseroot = os.path.abspath(options.caseroot)
print('CASE root directory: '+options.caseroot)
if(os.path.exists(options.caseroot) == False):
    os.system('mkdir -p '+caseroot)

#----- Construct default casename
casename    = options.site+"_"+compset
if (options.mycaseid != ""):
    casename = options.mycaseid+'_'+casename
if (options.ad_spinup):
    casename = casename+'_ad_spinup' 
    
#case directory
if (caseroot != "./"):
    casedir = caseroot+"/"+casename
else:
    casedir = casename
print ("CASE directory is: "+casedir+"\n")

#Check for existing case directory
if (os.path.exists(casedir)):
    print('Warning:  Case directory exists and --rmold not specified')
    var = raw_input('proceed (p), remove old (r), or exit (x)? ')
    if var[0] == 'r':
        os.system('rm -rf '+casedir)
    if var[0] == 'x':
        sys.exit()    


# cases run root directory
if (options.runroot == ''):
    runroot = csmdir+"/runs"
else:
    runroot = os.path.abspath(options.runroot)
print('CASE RUN root directory: '+runroot)
if(os.path.exists(options.runroot) == False):
    os.system('mkdir -p '+runroot)

blddir=runroot+"/"+casename+'/bld'
print ("CASE bulid and exeroot is: "+blddir+"\n")
rundir=runroot+"/"+casename+"/run"
print ("CASE rundir is: "+rundir+"\n")
 
           
#finidat file and finidat year
if (options.coldstart and (options.finidat != '' or options.finidat_case != '')):
        print('Error: Cannot have an finidat/finidat_case AND coldstart simultaneously! Exit \n')
        sys.exit()
    
if (options.finidat == '' and options.finidat_case == ''):   # not user-defined
    if (options.coldstart==False and compset.startswith('I1850')==True):
        options.finidat_case = casename+'_ad_spinup'
       
    if (compset.startswith('I20TR') == True):
        options.finidat_case = casename.replace('I20TR','I1850')
                
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

#default simyr
mysimyr=1850
if (options.compset.startswith('I20TR') == True):
    mysimyr=2000


#pft parameter file
# new or user-defined pft-phys file if desired
if (options.extended_pft and options.parm_file == ''):
    print('MUST provide user-defined parameter file! Exit \n')
    sys.exit()

if (options.parm_file != ''):
    pftfile = ccsm_input+'/lnd/clm2/paramdata/' + \
                  options.parm_file.replace('.nc','') + '.' + casename + '.nc'
    os.system('cp -f '+ccsm_input+'/lnd/clm2/paramdata/'+ options.parm_file + \
               ' '+pftfile)

#set number of run years, if user-defined
run_n = options.run_n
run_units = options.run_units

if (compset.startswith('I20TR') == True and options.run_n == 600):    # 600 is the default (i.e., NOT user-defined)
    if (run_units == 'nyears'): 
        run_n = endyear - 1850 +1
    elif (run_units == 'date'):
        run_n = endyear + 1  # date in yyyymmdd, but here only needs year value

#-------------------------------------------------------------

os.chdir(scriptsdir)

# ------------------ create, setup and build -----------------------------------------
#--- (1) create a new case
#create new case
comps = options.compset
    
print ('./create_newcase -case '+casedir +' '+machineoptions + \
    ' -compset '+ comps +' -res CLM_USRDAT ')
os.system('./create_newcase -case '+casedir+' '+machineoptions + \
                 ' -compset '+ comps +' -res CLM_USRDAT ' + \
                  ' > create_newcase.log')
if (os.path.isdir(casedir)):
    print(casename+' created.  See create_newcase.log for details')
    os.system('mv create_newcase.log ' + casedir +"/"+casename+"_case.log")
else:
    print('failed to create case.  See create_newcase.log for details')

# go to newly created case directory
os.chdir(casedir)

# (2) env_build.xml modification ---------------------------
if (options.runroot != ''):
    # the following is new
    os.system('./xmlchange -file env_build.xml -id ' \
                  +'CESMSCRATCHROOT -val '+runroot) 
    #build directory
    os.system('./xmlchange -file env_build.xml -id ' \
                  +'EXEROOT -val '+blddir) 
    
# turn off rtm module
os.system('./xmlchange -file env_build.xml -id ' \
                  +'RTM_MODE -val NULL') 

# turn off rtm flood module
os.system('./xmlchange -file env_build.xml -id ' \
                  +'RTM_FLOOD_MODE -val NULL') 

# clm4_5 config options (note: this configuration has re-designed, with most options moved to CLMNMLBLDOPT in env_run.xml)
# base physic options
clmconfig_opts = "-phys clm4_5"

os.system('./xmlchange -file env_build.xml -id ' \
    +'CLM_CONFIG_OPTS -val "'+clmconfig_opts+'"')
print ("CLM configuration options: " + clmconfig_opts +"\n")

# (3) env_run.xml modification ------------------------------------
# input/run/output directory
if (options.runroot != ''):
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'RUNDIR -val '+rundir) 
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DOUT_S -val TRUE') 
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DOUT_S_ROOT -val '+runroot+'/archives/'+casename) 
    
if (options.ccsm_input != ''):
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DIN_LOC_ROOT -val '+ccsm_input) 
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DIN_LOC_ROOT_CLMFORC -val '+ccsm_input+'/atm/datm7') 
    
if():
    # datm options
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_MODE -val CLM1PT') 
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_START -val '+str(startyear))
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_END -val '+str(endyear))

# run timestep
if (options.tstep != 0.5):
    os.system('./xmlchange -file env_run.xml -id ' \
                      +'ATM_NCPL -val '+str(int(24/float(options.tstep))))
    
# run-type adjusting -- needs checking ('rof' not working??) 
if (options.ad_spinup==False and options.coldstart==False):
    os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_REFDATE -val '+finidat_yst+'-01-01')
   
# run starting date/time
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
if (compset.startswith('I20TR') == True):
    os.system('./xmlchange -file env_run.xml -id ' \
                          +'CCSM_BGC -val CO2A')
    os.system('./xmlchange -file env_run.xml -id ' \
                          +'CLM_CO2_TYPE -val diagnostic')

# user-defined running stop options
os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_OPTION -val '+run_units)
if (options.run_units == 'date'):
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_DATE -val '+str(run_n)+'0101')
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_N -val -9999')        
else:
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_DATE -val -9999')
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_N -val '+str(run_n))        

# user-defined restart options
if (options.restart_units != ''):
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_OPTION -val '+options.restart_units)
if (options.restart_n != ''):
    if (options.restart_units == 'date'):
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_DATE -val '+str(options.restart_n)+'0101')
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_N -val -9999')
    else:
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_N -val '+str(options.restart_n))            
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'REST_DATE -val -9999')

#User-defined resolution data name   
os.system('./xmlchange -file env_run.xml -id CLM_USRDAT_NAME ' \
                  +'-val '+clm_usrdat_name)


# CLM build namelist options 
# (NOTE: clm-cn is by default
stdout  = os.popen("./xmlquery -valonly CLM_BLDNML_OPTS")
clmbldnmls = stdout.read().rstrip( )   

# i) CLM-CN or CLM-CENTURY
# NOTE: vertical-soil now is for both option
if (options.centbgc):
    clmbldnmls = " -bgc bgc"   
    # this will turn-on the whole package of clm4me, nitrif-denitrif, century-decomp
    print(" \n ======= CLM-BGC ==========================================" )
    print(" CLM bgc will be turned on with CENTURY BGC with clm4me, nitrif-denitrif ! " )

else:
    clmbldnmls = " -bgc cn"
    print(" \n ======= CLM-CN ==========================================" )
    print(" CLM bgc will be turned on with classic CLM-CN ! " )
    # option to turn on '-LCH4' for CLM45CN without CENTURY BGC, if requested ------
    # this is because configuration of 'CLM4ME' is bundled with centurary BGC
    if (options.kovenN):
        print(" \n ======= NITRIF-DENITRIF ====================================" )
        print(" CLM-CN+ classic bgc will be turned on with Koven' nitrif_denitrif option! " )
        clmbldnmls += " -nitrif_denitrif"
    #elif (options.centbgc==False):
    #print(" CLM-CN+ classic bgc will be with simple inorganic N cycle modules ! " )            

    if (options.CH4 and options.centbgc==False):
        print(" \n ======= CH4 =========================================" )
        print(" CLM-CN+ classic bgc will be turned on with clm4me option! " )
        clmbldnmls += " -methane"

               
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
    
if(int(options.np)<=1):
    if(int(options.mppwidth)>1):     # this input is for HOPPER PBS
        options.np = options.mppwidth
            
    if( int(options.mppnodes)>1 or 
        (int(options.mppnodes)==1 and int(options.ppn)>1) ):     # this input is for Titan PBS
            options.np = str(int(options.mppnodes)*int(options.ppn))

#if running with > 1 processor
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
    os.system('./case.setup -clean')
    os.system('rm -f user-nl-*')


# (5a)check if pflotran directory as an input
if (options.pflotran):
    if(options.pflotran_srcdir==""):
        print(" PFLOTRAN directories NOT defined, please provide one using '--pflotran_srcdir=' ! \n")
        sys.exit()
    elif(os.path.exists(options.pflotran_srcdir) == False):
        print(" PFLOTRAN directories NOT exist, please the directory provided using '--pflotran_srcdir=' ! \n")
        sys.exit()
            
# (5b) settings for clm coupled with pflotran, if requested ------
    else:
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
if (options.extended_pft):
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
os.system('./case.setup > configure.log')
        
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
                 "surfdata_"+str(numxpts)+"x"+str(numypts)+"pt_"+options.site+ \
                 "_simyr"+str(mysimyr)+".nc'\n")   
else:
    output.write(" fsurdat = '"+ccsm_input+"/lnd/clm2/surfdata_map/" + \
                 options.surfdatafile+"'\n")
    
#(6c) pft dynamics file for transient run ----
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
    histfile = runCLMdir+"/outputs-userdefined/"+options.hist_file
    hvars_file = open(histfile)
    output.write("\n")
    for s2 in hvars_file:
        myline = s2
        output.write(myline)
    output.write("\n")
    hvars_file.close()

# (6f) hacking 'nofire' namelist
if (options.nofire):
    print(" \n ======= FIRE ==========================================" )
    print(" Turn OFF fire option! " )
    output.write(" use_nofire = .true. \n")

#(6g) force namelist options for 'maxpatch_pft' changed if extended arctic pft ----
if (options.extended_pft):
    output.write(" maxpatch_pft = 23\n")

#(6h) namelist options for PFLOTRAN coupling ----
if (options.pflotran):
    output.write(" use_clm_interface = .true.\n")
    output.write(" use_pflotran = .true.\n")
    output.write("/\n")

    output.write("&clm_pflotran_inparm\n")
    output.write(" pflotran_prefix = '"+ casename + "'\n")
    output.write("/\n")
    
output.close()

# (7) copy user-defined sourcemods codes  ----
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
   
# (8) transient CO2 patch for transient run ----
if (compset.startswith('I20TR') == True):      
#   (8a) historical co2 stream data: globally 1 value ----
    os.system('cp '+csmdir+'/models/lnd/clm/doc/UsersGuide/co2_streams.txt ./')
    myinput  = open('co2_streams.txt')
    myoutput = open('co2_streams.txt.tmp','w')
    for s in myinput:
        s2 = s.strip()
        if (s2 =='<filePath>'):
            myoutput.write(s)
            myoutput.write('            '+ccsm_input+'/atm/datm7/CO2\n')
            next(myinput);
        elif (s2 =='<fileNames>'):
            myoutput.write(s)
            myoutput.write('            '+options.co2_file+'\n')
            next(myinput);
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()       
    os.system('mv co2_streams.txt.tmp co2_streams.txt')
        
# (8b) modifying default 'datm_atm.in' to include historical co2 stream data ---
    myinput  = open('./Buildconf/datmconf/datm_atm_in')
    myoutput = open('user_nl_datm','w')
    for s in myinput:
        s2 = s.strip()
        if (s2.startswith('dtlimit')):
            myoutput.write(' '+s2+',1.5\n')
        elif (s2.startswith('fillalgo')):
            myoutput.write(" "+s2+",'nn'\n")
        elif (s2.startswith('fillmask')):
            myoutput.write(" "+s2+",'nomask'\n")
        elif (s2.startswith('mapalgo')):
            myoutput.write(" "+s2+",'nn'\n")
        elif (s2.startswith('mapmask')):
            myoutput.write(" "+s2+",'nomask'\n")
        elif (s2.startswith('streams')):
            myoutput.write(" "+s2+",'datm.global1val.streams.co2.txt 1766 1766 2010'\n")
        elif (s2.startswith('taxmode')):
            myoutput.write(" taxmode = 'cycle', 'extend', 'extend'\n")
        elif (s2.startswith('tintalgo')):
            myoutput.write(" "+s2+",'linear'\n")
        else:
            myoutput.write(s)
 
    myinput.close()
    myoutput.close()       
        
# datm namelist modifications (cycle met data streams - a bug since clm4.5.10)
# the issue is that climate data will not repeating (taxmode is 'extend'). 
myoutput = open('user_nl_datm','w')
myoutput.write("&shr_strdata_nml\n")
if (compset.startswith('I20TR') == True):             
    myoutput.write(" taxmode = 'cycle', 'extend','extend'\n")
else:
    myoutput.write(" taxmode = 'cycle', 'extend'\n")
myoutput.write("/\n")
myoutput.close()      

# (9) ------- build clm45 within cesm ---------------------------------------------- 
#clean build if requested prior to build
if (options.clean_build):
    os.system('./case.build --clean-all')
    os.system('rm -rf '+rundir+'/*')
    os.system('rm -rf '+blddir+'/*')
            
#compile cesm     
os.system('./case.build')        
        
# note: *.build will sweep everything under ./Buildconf, but we need 'co2_streams.txt' in transient run ---
if (compset.startswith('I20TR') == True):
    os.system('cp -f co2_streams.txt ./Buildconf/datmconf/datm.global1val.streams.co2.txt')
    os.system('cp -f co2_streams.txt '+rundir+'/datm.global1val.streams.co2.txt')       
        
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

        if(os.path.isfile(pfindir+'/CLM-CN_database.dat')):
            os.system('cp '+pfindir+'/CLM-CN_database.dat '+rundir+'/')
        elif(os.path.isfile(pfindir+'/hanford-clm.dat')):
            os.system('cp '+pfindir+'/hanford-clm.dat '+rundir+'/')         
        else:
            if(os.path.isfile(pfindir+'/hanford-clm.dat') == False):
                print('Waring: NO PFLOTRAN-bgc database file "handford-clm.dat" or "CLM-CN_database.dat" exists! \n')
            
        if(glob.glob(pfindir+'/*.h5')):      
            os.system('cp '+pfindir+'/*.h5 '+rundir+'/')
        else:
            print('Warning: NO PFLOTRAN *.h5 input file exists! -- be sure it is the case! \n')
            
                 
    
    
# ----- copy rpointers and restart files to current run directory prior to run model ---
if (options.finidat_case != ''):
    os.system('cp -f '+runroot+'/'+options.finidat_case+'/run/' + \
              options.finidat_case+'.*'+finidat_yst+'* ' + rundir)
    os.system('cp -f '+runroot+'/'+options.finidat_case+'/run/'+ \
              'rpointer.* '+rundir)

# ----- submit job if requested ---
if (options.no_submit == False):    
    
    os.chdir(casedir)
    os.system("./case.submit")

