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

fileName = "20231030_ActivityRatesG1G5_S1"
filePath = "P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/"
wellPlate1 = AA.WellPlate(filePath+fileName+".txt")

getNormalized = False
hideEarlyData = False
hideLateData = False

invertWellData = False
clrOffset = 2.5 # 1.7 for deep blue, 2.5 for reddish pink

wellCoords = AA.getWellNames("B12", "G15")

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
            WDI["B12"],
            WDI["B13"],
            WDI["C12"],
            WDI["C13"],
            WDI["D13"],
            WDI["D12"]
        ]],
        [[
            WDI["B14"],
            WDI["B15"],
            WDI["C14"],
            WDI["C15"],
            WDI["D15"],
            WDI["D14"]
        ]]
    ],
    [
        [[
            WDI["E12"],
            WDI["E13"],
            WDI["F12"],
            WDI["F13"],
            WDI["G12"],
            WDI["G13"]
        ]],
        [[
            WDI["E14"],
            WDI["E15"],
            WDI["F14"],
            WDI["F15"],
            WDI["G14"],
            WDI["G15"]
        ]]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsGraphSetup,
    plotXTitles=["G1S1 2X ON", "G5S1 2X ON"], localXTitle="Time (min)", globalXTitle="Legend: Conc. of Template",
    plotYTitles=["4 U/uL RNAP, 1280 nM j.t.", "1 U/uL RNAP, 320 nM j.t."], localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
    colorOffset = clrOffset,
    plotFileName=filePath + fileName + ("" if getNormalized else "_RAW"),
    spacer=2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["10 nM", "20", "40", "75", "100", "150"]
                ],
                [
                    ["5 nM", "10", "15", "25", "40", "60"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else None
)