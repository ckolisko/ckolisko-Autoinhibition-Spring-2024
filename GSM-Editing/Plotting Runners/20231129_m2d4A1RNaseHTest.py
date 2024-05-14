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

fileName = "20231129_m2d4A1RNaseHTest"
filePath = "P:/PhD Work/Assets/Fluorescence Plots/RNase Tests/"
wellPlate1 = AA.WellPlate(filePath+fileName+".txt")
manualAddn = ""

getNormalized = False
hideEarlyData = False
hideLateData = False

invertWellData = False
clrOffset = 2.5 # 1.7 for deep blue, 2.5 for reddish pink

wellCoords = AA.getWellNames("I17", "L18")

timeBreakPts = wellPlate1.getTimeBreaks()
heatCorrectionTimes = []

if getNormalized:
    displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = timeBreakPts[0], upperBound = timeBreakPts[1]) # Find the "low" fluorescence
    scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = timeBreakPts[-2]+5, upperBound = None) # Find the "high" fluorescence
    for wellCoord in wellCoords:
        wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

    if invertWellData:
        AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0

    if len(heatCorrectionTimes)>0:
        AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)
timeBds = [timeBreakPts[-3] if hideEarlyData else 0, timeBreakPts[-2] if hideLateData else None]


timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsGraphSetup = [
    [
        [
            [
                WDI["I17"]
            ],
            [
                WDI["I18"]
            ]
        ],
        [
            [
                WDI["J17"]
            ],
            [
                WDI["J18"]
            ]
        ]
    ],
    [
        [
            [
                WDI["K17"]
            ],
            [
                WDI["K18"]
            ]
        ],
        [
            [
                WDI["L17"]
            ],
            [
                WDI["L18"]
            ]
        ]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsGraphSetup,
    plotXTitles = ["", "+ RNase H"], localXTitle="Time (min)", globalXTitle = "",
    plotYTitles = ["", "+ G4R1 2X ON"], localYTitle = "Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence", globalYTitle = "",
    colorOffset = clrOffset,
    plotFileName = filePath + fileName + (manualAddn if (getNormalized or hideEarlyData or hideLateData) else "_RAW"),
    spacer = 2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["A1", "m2d4A1"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else [-2000, 52000]
)