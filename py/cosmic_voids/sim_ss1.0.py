#!/usr/bin/env/python
import os
from vide.backend.classes import *

continueRun = False # set to True to enable restarting aborted jobs
startCatalogStage = 1
endCatalogStage   = 3
               
regenerateFlag = False

mergingThreshold = 1e-9

dataSampleList = []
           
setName = "sim_ss1.0"

workDir = "/global/u2/c/chanun/cosmic_voids/py/vide-test/examples/example_simulation//sim_ss1.0/"
inputDataDir = "/global/u2/c/chanun/cosmic_voids/py/vide-test/examples/"
figDir = "/global/u2/c/chanun/cosmic_voids/py/vide-test/figs/example_simulation//sim_ss1.0/"
logDir = "/global/u2/c/chanun/cosmic_voids/py/vide-test/logs/example_simulation//sim_ss1.0/"

numZobovDivisions = 2
numZobovThreads = 2
               
newSample = Sample(dataFile = "example_simulation_z0.0.dat",
                   dataFormat = "multidark",
                   dataUnit = 1,
                   fullName = "sim_ss1.0_z0.00_d00",
                   nickName = " SS 1.0, z = 0.00",
                   dataType = "simulation",
                   zBoundary = (0.00, 0.36),
                   zRange    = (0.00, 0.36),
                   zBoundaryMpc = (0.00, 999.98),
                   shiftSimZ = False,
                   omegaM    = 0.2847979853038958,
                   minVoidRadius = 1,
                   profileBinSize = "auto",
                   includeInHubble = True,
                   partOfCombo = False,
                   boxLen = 999.983,
                   usePecVel = False,
                   numSubvolumes = 1,
                   mySubvolume = "00",
                   useLightCone = False,
                   subsample = "1.0")
dataSampleList.append(newSample)
  