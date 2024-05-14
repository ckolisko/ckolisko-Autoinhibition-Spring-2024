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

fileName = "P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230929_G5_10_20ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

getNormalized = False

wellCoords = ["B2", "B3", "B4", "B5", "C2", "C3", "C4", "C5", "D2", "D3", "D4", "D5", "E2", "E3", "E4", "E5", "F2", "F3", "F4", "F5", "G2", "G3", "G4", "G5"]
wellNames = []
heatCorrectionTimes = []

if getNormalized:
    displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 3000/60, upperBound = 5000/60) # Find the "low" fluorescence
    scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 39500/60, upperBound = 48000/60) # Find the "high" fluorescence
    for wellCoord in wellCoords:
        wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

    #AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

    if len(heatCorrectionTimes)>0:
        AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [5000/60, 39500/60 if getNormalized else None]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsCursedVersion = [
    [
        [[
            WDI["B2"],
            WDI["B3"],
            WDI["C2"],
            WDI["C3"],
            WDI["D2"],
            WDI["D3"]
        ]],
        [[
            WDI["B4"],
            WDI["B5"],
            WDI["C4"],
            WDI["C5"],
            WDI["D4"],
            WDI["D5"]
        ]]
    ],
    [
        [[
            WDI["E2"],
            WDI["E3"],
            WDI["F2"],
            WDI["F3"],
            WDI["G2"],
            WDI["G3"]
        ]],
        [[
            WDI["E4"],
            WDI["E5"],
            WDI["F4"],
            WDI["F5"],
            WDI["G4"],
            WDI["G5"]
        ]]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    plotXTitles=["10 nM", "20 nM"], globalXTitle="G5S1 Concentration", localXTitle="Time (min)",
    plotYTitles=["4 U/uL", "1 U/uL"], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
    colorOffset = 2.5,
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230929_G5_10_20ActivityWithJunkTemplate" + ("" if getNormalized else "_RAW"),
    spacer=2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["0 nM", "160", "320", "480", "640", "960"]
                ],
                [
                    ["0 nM", "20 ", "40", "80", "120", "160"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05] if getNormalized else None
)