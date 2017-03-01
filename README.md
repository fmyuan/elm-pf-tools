# **README**

**A package of Tools, including scripts, settings, test-suites, user-guides, etc. for ALM-PFLOTRAN models.**

## **I. In Briefing**

### I-1. Directory Tree
```
|
| --machine-userdefined   (a few sets of user-defined machine files for CLM4.5.35)
  --userdefined_output    (a  few examples of user-defined outputs from CLM: it's copied to 'lnd_in' under CLM run directory)
  --link-tree             (a script to soft-link CLM inputdata to user's inputdata directory. It was in scripts/ before 4.5.81, but removed after and appears useful)

  --runPTCLM.py           (a python script to run PoinT-mode CLM)
  --runGRDCLM.py          (a python script to run Global-mode CLM)
  
  --alm_runPTCLMpy_mymac.sh   (a bash example script to link/copy data and call 'runPTCLM.py' on MAC OSX to configure, build, and run off-line CLM)
  --alm_runPTCLMpy_cades.sh   (a bash example script to link/copy data and call 'runPTCLM.py' on CADES-OR-CONDO to configure, build, and run off-line CLM)
  


```

### I-2.  


## **II. Wiki**
[Wiki Page](https://code.ornl.gov/alm-pflotran/clm-pf-tools/wikis/home#guides-for-coupling-alm-and-pflotran)

### II-1. HOW-TOs



### II-2. for ALM users

[MAC OSX users](https://code.ornl.gov/alm-pflotran/clm-pf-tools/wikis/Alm%20on%20mac:%20environments,%20building,%20and%20running)


[CCSI-CADES users](https://code.ornl.gov/alm-pflotran/clm-pf-tools/wikis/Alm%20on%20cades,%20or%20condo:%20environments,%20building,%20and%20running)


### II-3. for PFLOTRAN users



## III. ALM-PFLOTRAN Developments


UPDATED: June-06-2016.
@AUTHORS