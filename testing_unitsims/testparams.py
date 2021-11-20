import os

PWD=os.getcwd()

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# CONFIGURATION

# if True, will scan log files for last known completed state and run from there
continueRun = False

# stages:
#   1 : extract redshift slices from data
#   2 : void extraction using zobov
#   3 : removal of small voids and voids near the edge 
startCatalogStage = 1
endCatalogStage   = 3

#changed the directories-chanun
# directory for the input simulation files
catalogDir = PWD+"/testoutputs/"

# void catalog output directory
voidOutputDir = PWD+"/testoutputs/voids/"

# output directory for log files
logDir = PWD+"/testoutputs/logs/"

# output directory for figures
figDir = PWD+"/testoutputs/figs/"

# where to place the pipeline scripts
scriptDir = PWD+"/testoutputs/scripts"

# don't change
dataType = "simulation"


#unsure but the UNIT sims uses gadget so i changed it - chanun
# available formats for simulation: gadget, sdf, multidark
dataFormat = "gadget"

#i think we can keep this-chanun
# units of position in Mpc/h
dataUnit = 1

#probably false? -chanun
# place particles on the lightcone (z-axis in sims)?
useLightCone = False

#probably false?-chanun
# add peculiar velocities?
doPecVel = False

#left this as is -chanun
# optimization: maximum number of parallel threads to use
numZobovThreads = 2
# optimization: number of subdivisions of the box
numZobovDivisions = 2

#let's not merge -chanun
# maximum density for merging voids
#   0 (equivalent to infinitely large value) -> Merge everything (no threshold)
#   1e-9 (or smaller != 0) -> Do not merge anything
mergingThreshold = 1e-9

# prefix to give all outputs
prefix = "test_"

#we can keep it at 1 for xyz -chanun
# how many independent slices along the z-axis?
numSlices = 1

# how many subdivisions along the x- and y- axis?
#   ( = 2 will make 4 subvolumes for each slice, = 3 will make 9, etc.)
numSubvolumes = 1


###############################################################################
# Particles
#didnt chnage anything here-chanun

# common filename of particle files
particleFileBase = "example_simulation_NNNN.dat"


# this flag will be replacedui by values in fileNums list below
particleFileDummy = 'NNNN'

# list of file numbers for the particle files
fileNums = ["z0.0"]

# redshift of each file in the above fileNums list
redshifts = ["0.0"]

# list of desired subsamples - these are in unts of h Mpc^-3!
subSamples = [1.0]

# if True, do the subsampling in preparation (available for sdf and multidark)
doSubSamplingInPrep = False


# if 'absolute', subSamples are given in particles per cubic Mpc/h
# if 'relative', subSamples are given as a fraction of input particles
subSampleMode = "relative"

# shift the z-coord of sims with redshift
shiftSimZ = False

###############################################################################
# Halos
#didnt change anything though maybe we should since we have halos -chanun

# common filename of halo files, leave blank to ignore halos
haloFileBase = ""
#haloFileBase = "mf_4s_1G_1k_bgc2_NNNNN.sdf"

# this flag will be replaced by values in fileNums list above
haloFileDummy = ''
#haloFileDummy = 'NNNNN'

# minimum halo mass cuts to apply for the halo catalog
#   use "none" to get all halos
minHaloMasses = []
#minHaloMasses = ["none", 1.2e13]

# locations of data in the halo catalog
haloFileMCol  = 6    # mass
haloFileXCol  = 0    # x
haloFileYCol  = 1    # y
haloFileZCol  = 2    # z
haloFileVXCol = 3    # v_x
haloFileVYCol = 4    # v_y
haloFileVZCol = 5    # v_z
haloFileColSep = ',' # separator
haloFileNumComLines = 0 # number of comments before data


###############################################################################
# simulation information
#i think these params are right based off the sim config on desi -chanun

numPart = 4096*4096*4096
lbox = 1000 # Mpc/h
omegaM = 0.3089
hubble = 0.6774 # h_0/100


###############################################################################
# HOD
#did not change anything -chanun

# each of the HOD sets will be applied to each halo catalog defined above
hodParmList = [
  #{'name'       : "LowRes", #BOSS: Manera et al. 2012, eq. 26
  # 'Mmin'       : 0.0,
  # 'M1'         : 1.e14,
  # 'sigma_logM' : 0.596,
  # 'alpha'      : 1.0127,
  # 'Mcut'       : 1.19399e13,
  # 'galDens'    : 0.0002,
  #},
]

# END CONFIGURATION
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
