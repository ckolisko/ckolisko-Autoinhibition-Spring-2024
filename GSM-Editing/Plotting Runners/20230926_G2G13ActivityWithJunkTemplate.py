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

fileName = "P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230926_G2G13ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

getNormalized = True

wellCoords = ["I7", "I8", "I9", "I10", "J7", "J8", "J9", "J10", "K7", "K8", "K9", "K10", "L7", "L8", "L9", "L10", "M7", "M8", "M9", "M10", "N7", "N8", "N9", "N10"]
wellNames = []
heatCorrectionTimes = []

if getNormalized:
    displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2500/60, upperBound = 5000/60) # Find the "low" fluorescence
    scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 45000/60, upperBound = None) # Find the "high" fluorescence
    for wellCoord in wellCoords:
        wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

    #AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

    if len(heatCorrectionTimes)>0:
        AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [5000/60, 45000/60 if getNormalized else None]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsCursedVersion = [
    [
        [[
            WDI["I7"],
            WDI["I8"],
            WDI["J7"],
            WDI["M7"],
            WDI["J8"],
            WDI["K7"],
            WDI["K8"]
        ]],
        [[
            WDI["I9"],
            WDI["I10"],
            WDI["J9"],
            WDI["M9"],
            WDI["J10"],
            WDI["K9"],
            WDI["K10"]
        ]]
    ],
    [
        [[
            WDI["L7"],
            WDI["N8"],
            WDI["M8"],
            WDI["N7"],
            WDI["L8"]
        ]],
        [[
            WDI["L9"],
            WDI["N10"],
            WDI["M10"],
            WDI["N9"],
            WDI["L10"]
        ]]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    plotXTitles=["G2S1", "G13S1"], globalXTitle="10 nM Template Type", localXTitle="Time (min)",
    plotYTitles=["4 U/uL", "1 U/uL"], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
    colorOffset = 2.5,
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230926_G2G13ActivityWithJunkTemplate" + ("" if getNormalized else "_RAW"),
    spacer=2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["0 nM", "320", "640", "800", "960", "1280", "1500"]
                ],
                [
                    ["0 nM", "40", "120", "160", "400"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else None
)