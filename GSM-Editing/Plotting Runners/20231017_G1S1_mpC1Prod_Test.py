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

fileName = "20231017_G1S1_mpC1Prod_Test_G1S1"
filePath = "P:/PhD Work/Assets/Fluorescence Plots/G1 Activity/"+fileName+".txt"
wellPlate1 = AA.WellPlate(filePath)

getNormalized = False
hideEarlyData = False
hideLateData = False

invertWellData = False
clrOffset = 1.7 # 1.7 for deep blue, 2.5 for reddish pink

wellCoords = AA.getWellNames("I7", "N10")

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
timeBds = [timeBreakPts[-3] if hideEarlyData else 0, timeBreakPts[-1] if hideLateData else None]


timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsGraphSetup = [
    [
        [[
            WDI["I7"],
            WDI["I8"],
            WDI["I9"],
            WDI["I10"]
        ]],
        [[
            WDI["J7"],
            WDI["J8"],
            WDI["J9"],
            WDI["J10"]
        ]],
        [[
            WDI["K7"],
            WDI["K8"],
            WDI["K9"],
            WDI["K10"]
        ]]
    ],
    [
        [[
            WDI["L7"],
            WDI["L8"],
            WDI["L9"],
            WDI["L10"]
        ]],
        [[
            WDI["M7"],
            WDI["M8"],
            WDI["M9"],
            WDI["M10"]
        ]],[[
            WDI["N7"],
            WDI["N8"],
            WDI["N9"],
            WDI["N10"]
        ]]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsGraphSetup,
    plotXTitles=["G1S1 2X ON", "G1S1 mpBLK 1X Act", "G1S1 mpBLK 2X Act"], localXTitle="Time (min)", globalXTitle="Legend: Conc. of Template, Mult. of 8.92*10^(-3) U/uL RNase H",
    plotYTitles=["4 U/uL RNAP", "1 U/uL RNAP"], localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
    colorOffset = clrOffset,
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/G1 Activity/" + fileName + ("" if getNormalized else "_RAW"),
    spacer=2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["10 nM, 0x", "20 nM, 0x", "10 nM, 1x", "20 nM, 1x"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else None
)