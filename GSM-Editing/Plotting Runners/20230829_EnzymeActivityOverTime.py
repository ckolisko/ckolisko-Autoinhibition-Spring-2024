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

fileName = "P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230829_EnzymeActivityOverTime.txt"
wellPlate1 = AA.WellPlate(fileName)

wellCoords = ["J6", "J7", "J8", "K6", "K7", "K8", "L6", "L7", "L8", "M6", "M7", "M8", "N6", "N7", "N8", "O6", "O7", "O8"]
wellNames = []
heatCorrectionTimes = [1700, 3500, 5300, 7400, 9000, 10700, 12500, 14400]

displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 0/60, upperBound = 1700/60) # Find the "low" fluorescence
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 40500/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

wellPlate1.adjustHeatCorrection(0.08, 12)

# if len(heatCorrectionTimes)>0:
#     AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [0/60, 40500/60]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest
wellsCursedVersion = [
    # [
    #     [[WDI["J6"], WDI["J7"], WDI["J8"], WDI["K6"], WDI["K7"], WDI["K8"], WDI["L6"], WDI["L7"], WDI["L8"]]],
    # ],
    [
        [[WDI["M6"],
        WDI["M7"],
        WDI["M8"],
        WDI["N6"],
        WDI["N7"],
        WDI["N8"],
        WDI["O6"],
        WDI["O7"],
        WDI["O8"]]]
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
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230829_EnzymeActivityOverTime_1UuL_320nM",
    spacer=2,
    specialOpts = [
        ["Legend", [[["0 hr", "0.5 hr", "1 hr", "1.5 hr", "2 hr", "2.5 hr", "3 hr", "3.5 hr", "4 hr"]]]]
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