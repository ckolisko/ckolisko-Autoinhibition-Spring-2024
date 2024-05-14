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

fileName = "20231103_AzoR5ControlTests"
filePath = "P:/PhD Work/Assets/Fluorescence Plots/Azobenzene/"
wellPlate1 = AA.WellPlate(filePath+fileName+".txt")

getNormalized = False
hideEarlyData = False
hideLateData = False

invertWellData = False
clrOffset = 2.5 # 1.7 for deep blue, 2.5 for reddish pink

wellCoords = AA.getWellNames("L2", "M5")

timeBreakPts = wellPlate1.getTimeBreaks()
heatCorrectionTimes = []

if getNormalized:
    displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = timeBreakPts[0], upperBound = timeBreakPts[1]) # Find the "low" fluorescence
    scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = timeBreakPts[-1]+5, upperBound = None) # Find the "high" fluorescence
    for wellCoord in wellCoords:
        wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

    if invertWellData:
        AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0

    if len(heatCorrectionTimes)>0:
        AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)
timeBds = [timeBreakPts[-2] if hideEarlyData else 0, timeBreakPts[-1] if hideLateData else None]


timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsGraphSetup = [
    [
        [[
            WDI["L2"],
            WDI["L3"]
        ]],
        [[
            WDI["L4"],
            WDI["L5"]
        ]],
    ],
    [
        [[
            WDI["M2"],
            WDI["M3"]
        ]],
        [[
            WDI["M4"],
            WDI["M5"]
        ]]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsGraphSetup,
    plotXTitles=[], localXTitle="Time (min)", globalXTitle="Azobenzene Controls",
    plotYTitles=[], localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
    colorOffset = clrOffset,
    plotFileName=filePath + fileName + ("" if getNormalized else "_RAW"),
    spacer=2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["G5S1 OFF", "G5S1 2X ON"],
                    ["G5S1 OFF + 2X A5", "G5S1 OFF + 3X AzoR5"]
                ],
                [
                    ["G5S1 2X ON + 3X AzoR5", "G5S1 2X ON + 10X AzoR5"],
                    ["2X A5:3X AzoR5 + G5S1 OFF", "2X A5:10X AzoR5 + G5S1 OFF"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else None
)