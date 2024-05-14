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

fileName = "20231212_m2d4A1_3'm5G1S1nt_Tests"
filePath = "P:/PhD Work/Assets/Fluorescence Plots/Modified Genelets/"
wellPlate1 = AA.WellPlate(filePath+fileName+".txt")
manualAddn = ""

getNormalized = False
hideEarlyData = False
hideLateData = False

invertWellData = False
clrOffset = 2.5 # 1.7 for deep blue, 2.5 for reddish pink

wellCoords = AA.getWellNames("B6", "E9")

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
timeBds = [timeBreakPts[-3] if hideEarlyData else 0, timeBreakPts[-2] if hideLateData else timeBreakPts[-2]]


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
                WDI["D6"]
            ],
            [
                WDI["D7"]
            ],
            [],
            [],
            [],
            []
        ],
        [
            [],
            [],
            [
                WDI["B6"]
            ],
            [
                WDI["B7"]
            ],
            [
                WDI["C6"]
            ],
            [
                WDI["C7"]
            ]
        ]
    ],
    [
        [],
        [
            [],
            [],
            [
                WDI["B8"]
            ],
            [
                WDI["B9"]
            ],
            [
                WDI["C8"]
            ],
            [
                WDI["C9"]
            ]
        ]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsGraphSetup,
    plotXTitles = ["Controls (640 nM JT)", "Trials"], localXTitle="Time (min)", globalXTitle = "RNases Multiplier x(0.005 U/uL RNase H, 0.0001 U/uL RNase A/T1)",
    plotYTitles = ["320 nM", "960 nM"], localYTitle = "Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence", globalYTitle = "Trial Junk Template (JT) Concentration",
    colorOffset = clrOffset,
    plotFileName = filePath + fileName + (manualAddn if (getNormalized or hideEarlyData or hideLateData) else "_RAW"),
    spacer = 2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["U-G1S1", "U-G1S1 with 3'm5G1S1-nt"],
                    ["Normal", "m2d4A1", "3'm5G1S1-nt", "Both"]
                ],
                [
                    [],
                    []
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else [-2000, 35000]
)