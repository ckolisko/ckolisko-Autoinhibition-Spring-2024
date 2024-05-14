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

fileName = "P:/PhD Work/Assets/Fluorescence Plots/G1 Activity/06262023_G1Y1_3PrimeQTH_Test_MgCl.txt"
wellPlate1 = AA.WellPlate(fileName)

wellCoords = ["B8", "B9", "B10", "B11", "C8", "C9", "C10", "C11"]
wellNames = []
heatCorrectionTimes = []

displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 1500/60, upperBound = 2700/60) # Find the "low" fluorescence
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 66000/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [40, 300] # [52100/60, None]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest
wellsCursedVersion = [
    [
        [[WDI["B8"]]],
        [[WDI["B9"]]],
        [[WDI["B10"]]],
        [[WDI["B11"]]]
    ],
    [
        [],
        [],
        [],
        []
    ]
    
]

wellsActualStuff = [
    [
        [[WDI["C8"]]],
        [[WDI["C9"]]],
        [[WDI["C10"]]],
        [[WDI["C11"]]]
    ],
    [
        [],
        [],
        [],
        []
    ]
]

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    localXTitle="Time (hrs)",
    localYTitle="Fractional Reporter Displaced",
    colorOffset = 2.5,
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/G1 Activity/DELETE_TEST_1",
    spacer=2
)

AP.generatePlots(timePtsOfInterest, wellsActualStuff,
    localXTitle="Time (hrs)",
    localYTitle="Fractional Reporter Displaced",
    colorOffset = 2.5,
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/G1 Activity/DELETE_TEST_12",
    spacer=2
)