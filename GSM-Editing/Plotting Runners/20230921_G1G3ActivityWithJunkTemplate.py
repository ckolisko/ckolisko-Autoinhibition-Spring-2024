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

fileName = "P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230921_G1G3ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

getNormalized = True

wellCoords = ["B7", "B8", "B9", "B10", "C7", "C8", "C9", "C10", "D7", "D8", "D9", "D10", "E7", "E8", "E9", "E10", "F7", "F8", "F9", "F10", "G7", "G8", "G9", "G10"]
wellNames = []
heatCorrectionTimes = []

if getNormalized:
    displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2700/60, upperBound = 5000/60) # Find the "low" fluorescence
    scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 46500/60, upperBound = 50000/60) # Find the "high" fluorescence
    for wellCoord in wellCoords:
        wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

    #AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

    if len(heatCorrectionTimes)>0:
        AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [5000/60, 46500/60 if getNormalized else None]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsCursedVersion = [
    [
        [[
            WDI["B7"],
            WDI["B8"],
            WDI["C7"],
            WDI["C8"],
            WDI["D7"],
            WDI["D8"]
        ]],
        [[
            WDI["B9"],
            WDI["B10"],
            WDI["C9"],
            WDI["C10"],
            WDI["D9"],
            WDI["D10"]
        ]]
    ],
    [
        [[
            WDI["E7"],
            WDI["E8"],
            WDI["F7"],
            WDI["F8"],
            WDI["G7"],
            WDI["G8"]
        ]],
        [[
            WDI["E9"],
            WDI["E10"],
            WDI["F9"],
            WDI["F10"],
            WDI["G9"],
            WDI["G10"]
        ]]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    plotXTitles=["G1S1", "G3S1"], globalXTitle="10 nM Template Type", localXTitle="Time (min)",
    plotYTitles=["4 U/uL", "1 U/uL"], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
    colorOffset = 2.5,
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230921_G1G3ActivityWithJunkTemplate" + ("" if getNormalized else "_RAW"),
    spacer=2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["0 nM", "320", "640", "960", "1280", "1500"]
                ],
                [
                    ["0 nM", "40 ", "80", "120", "160", "240"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else None
)