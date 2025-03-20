#!/usr/bin/env python

#######################################################################################################
#
# test modules
#
##################################################################################
"""
odata = \
    subsetDaymetReadNCfile('daymet_v3_prcp_1980_na.nc', lon_range=[-170.0,-141.0], lat_range=[60.0, 90.0], SUBSETNC=True)
"""

#----------------------------
"""

site,odata_header,odata = \
    singleDaymetReadCsvfile('daymet_kougarok-NGEE00.csv')

"""

#----------------------------
"""

# read-in metdata from CPL_BYPASS_FULL
cplbypass_dir='./cpl_bypass_full'
cplbypass_fileheader=''
cplbypass_mettype='GSWP3'
lon = 195.25
lat = 65.25  
met_cplbypass = clm_metdata_cplbypass_read(cplbypass_dir,cplbypass_fileheader, cplbypass_mettype, 'TBOT', lon, lat)

"""

#----------------------------
"""

site,odata_header,odata = \
    singleNCDCReadCsvfile('2016100_initproc.csv')

"""

#----------------------------
"""

odata_header,odata = \
    singleReadCsvfile('5_minute_data_v21.csv')

#odata_header,odata = \
#    singleReadCsvfile('Met_Final_forELM.csv')

"""

#----------------------------
"""

clm_metdata_CRUJRA('./rawdata/crujra.v2.3.5d')

"""

#----------------------------
"""

#clm_metdata_cplbypass_extraction('./', 'GSWP3_daymet4', 203.1241, 70.5725,ncopath='/usr/local/nco/bin/') #BEO
#clm_metdata_cplbypass_extraction('./', 'GSWP3_daymet4', -157.4089, 70.4696,ncopath='/usr/local/nco/bin/')  #ATQ
##clm_metdata_cplbypass_extraction('./', 'CRUJRAV2.3.c2023.0.5x0.5', -97.0287, 27.9798, ncopath='/software/user_tools/current/cades-ccsi/nco/nco-5.1/bin/')  #test

#multiple_cplbypass_extraction('info_14sites.txt')
#multiple_cplbypass_extraction('info_PIE3sites.txt')

"""

#----------------------------
"""

#NGEE-p4 sites: BEO(TILE14412), C71(TILE13869), CHARS, K64(TILE13868), KFC(TILE13868), QHI(TILE14241)), T27(TILE13868), T47(TILE13868), TFS(TILE14236), TVC(TILE14244)

# clm_metdata_cplbypass_extraction('./TILE14258/cpl_bypass_full/', \
#                                  'GSWP3_daymet4', \
#                                  -105.0415, 69.1198, \
#                                  ncopath= \
#'/sw/baseline/spack-envs/base/opt/linux-rhel8-x86_64/gcc-8.5.0/nco-5.1.5-scmfzcuxfmldbn3lwwcjuesd22eumhwk/bin/')

"""

'''
from pytools.metdata_processing.met_utils import multiple_cplbypass_extraction
multiple_cplbypass_extraction('info_TFS_meq2_sites.txt')
'''
