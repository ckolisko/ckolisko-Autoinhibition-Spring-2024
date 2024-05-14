from pathlib import Path
import sys
sys.path.insert(1, str(Path('../')))

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import AssayAnalysis as AA
import ArrayPlots as AP
import GeneralUse as GU
import GeneletSystemModel as GSM

fileName = "C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230808_G1S1_ProductionTest_Optimize_Files/20230808_G1S1_ProductionTest.txt"
wellPlate1 = AA.WellPlate(fileName)

wellCoords = ["K8", "K9", "K10", "K11", "L8", "L9", "L10", "L11"]
wellNames = []
heatCorrectionTimes = []

displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2000/60, upperBound = 4000/60) # Find the "low" fluorescence
# TODO, see is lower up is based right after spike.
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 34330/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)

# Alway put in the actual time bounds from the data (just before big spike)
timeBds = [4000/60, 33880/60] # [52100/60, None]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest
wellsCursedVersion = [
    [
        [[WDI["K8"]]],
        [[WDI["K9"]]],
        [[WDI["K10"]]],
        [[WDI["K11"]]]
    ],
    [
        [[WDI["L8"]]],
        [[WDI["L9"]]],
        [[WDI["L10"]]],
        [[WDI["L11"]]]
    ]
    
]

# wellsActualStuff = [
#     [
#         [[WDI["L8"]]],
#         [[WDI["L9"]]],
#         [[WDI["L10"]]],
#         [[WDI["L11"]]]
#     ],
#     [
#         [],
#         [],
#         [],
#         []
#     ]
# ]

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    localXTitle="Time (min)",
    localYTitle="Fractional Reporter Displaced",
    colorOffset = 2.5,
    plotFileName="C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230808_G1S1_ProductionTest_Optimize_Files/20230808_G1S1_ProductionTest4",
    spacer=2,
    specialOpts = [
        ["Legend", [[["-G1S1, -A1"], ["-A1"], ["-RNAP"], ["+RNase H"]], [["2 nM"], ["4 nM"], ["8 nM (DEFAULT)"], ["16 nM"]]]]
    ],
    fixedY=[-0.05, 1.05]
)

# AP.generatePlots(timePtsOfInterest, wellsActualStuff,
#     localXTitle="Time (hrs)",
#     localYTitle="Fractional Reporter Displaced",
#     colorOffset = 2.5,
#     plotFileName="P:/PhD Work/Assets/Fluorescence Plots/G1 Activity/DELETE_TEST_12",
#     spacer=2
# )