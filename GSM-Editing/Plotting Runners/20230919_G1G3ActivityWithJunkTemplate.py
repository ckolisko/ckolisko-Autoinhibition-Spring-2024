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

fileName = "C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230919_G1G3ActivityWithJunkTemplate_Optimize_Files/20230919_G1G3ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

getNormalized = False

wellCoords = ["H2", "H3", "H4", "H5", "I2", "I3", "I4", "I5", "J2", "J3", "J4", "J5", "K2", "K3", "K4", "K5", "L2", "L3", "L4", "L5", "M2", "M3", "M4", "M5"]
wellNames = []
heatCorrectionTimes = []

if getNormalized:
    displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 3000/60, upperBound = 5400/60) # Find the "low" fluorescence
    scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 39700/60, upperBound = None) # Find the "high" fluorescence
    for wellCoord in wellCoords:
        wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

    #AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

    if len(heatCorrectionTimes)>0:
        AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [5500/60, 40000/60 if getNormalized else None]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsCursedVersion = [
    [
        [[
            WDI["H2"],
            WDI["H3"],
            WDI["I2"],
            WDI["I3"],
            WDI["J2"],
            WDI["J3"]
        ]],
        [[
            WDI["H4"],
            WDI["H5"],
            WDI["I4"],
            WDI["I5"],
            WDI["J4"],
            []
        ]]
    ],
    [
        [[
            WDI["K2"],
            WDI["K3"],
            WDI["L2"],
            WDI["L3"],
            WDI["M2"],
            WDI["M3"]
        ]],
        [[
            WDI["K4"],
            WDI["K5"],
            WDI["L4"],
            WDI["L5"],
            WDI["M4"],
            []
        ]]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    plotXTitles=["G1S1", "G3S1"], globalXTitle="10 nM Template Type", localXTitle="Time (min)",
    plotYTitles=["4 U/uL", "1 U/uL"], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
    colorOffset = 2.5,
    plotFileName="C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230919_G1G3ActivityWithJunkTemplate_Optimize_Files/20230919_G1G3ActivityWithJunkTemplate3" + ("" if getNormalized else "_RAW"),
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