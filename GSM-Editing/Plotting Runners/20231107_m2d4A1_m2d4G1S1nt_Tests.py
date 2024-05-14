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

fileName = "20231107_m2d4A1_m2d4G1S1nt_Tests"
filePath = "P:/PhD Work/Assets/Fluorescence Plots/Modified Genelets/"
wellPlate1 = AA.WellPlate(filePath+fileName+".txt")

getNormalized = True
hideEarlyData = True
hideLateData = True

invertWellData = False
clrOffset = 2.5 # 1.7 for deep blue, 2.5 for reddish pink

wellCoords = AA.getWellNames("B17", "G19")

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
        [
            [
                WDI["B17"],
                WDI["C17"]
            ],
            [
                WDI["D17"],
                WDI["E17"],
                WDI["F17"],
                WDI["G17"]
            ]
        ],
        [
             [
                WDI["B18"],
                WDI["C18"]
             ],
             [
                WDI["D18"],
                WDI["E18"],
                WDI["F18"],
                WDI["G18"]
            ]
        ]
    ],
    [
        [
            [
                WDI["B19"],
                WDI["C19"]
            ],
            [
                WDI["D19"],
                WDI["E19"],
                WDI["F19"],
                WDI["G19"]
            ]
        ],
        []
    ]
]

AP.generatePlots(timePtsOfInterest, wellsGraphSetup,
    plotXTitles=[], localXTitle="Time (min)", globalXTitle="Legend: State, Modified strands",
    plotYTitles=[], localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
    colorOffset = clrOffset,
    plotFileName=filePath + fileName + ("" if getNormalized else "_RAW"),
    spacer=2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["OFF", "OFF, non-t", "2X ON", "2X ON, non-t", "2X ON, Act", "2X ON, non-t/Act"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else None
)