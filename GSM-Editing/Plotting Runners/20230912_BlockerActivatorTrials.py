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

fileName = "P:/PhD Work/Assets/Fluorescence Plots/G3 Activity/20230912_BlockerActivatorTrials.txt"
wellPlate1 = AA.WellPlate(fileName)

wellCoords = ["B2", "B3", "B4", "C2", "C3", "C4", "C5", "D2", "D3", "D4", "D5", "E2", "E3", "E4", "E5", "F4", "F5"]
wellNames = []
heatCorrectionTimes = []

displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 3800/60, upperBound = 5000/60) # Find the "low" fluorescence
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 35000/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [5000/60, 35000/60]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsCursedVersion = [
    [
        [
            [
                # WDI["C2"],
                # WDI["C3"],
                WDI["C4"],
                WDI["C5"],
            ],
            [
                WDI["B2"],
                # WDI["B3"],
                WDI["B4"],
            ]
        ],
        [
            [
                np.array([]),
                # WDI["C3"],
                WDI["C5"]
            ],
            [
                # WDI["D2"],
                # WDI["D3"],
                WDI["D4"],
                WDI["D5"],
            ],
            []
        ]
    ],
    [
        [
            [
                np.array([]),
                # WDI["C3"],
                WDI["C5"]
            ],
            [
                # WDI["E2"],
                # WDI["E3"],
                WDI["E4"],
                WDI["E5"],
            ],
            []
        ],
        [
            [
                np.array([]),
                WDI["C5"]
            ],
            [
                WDI["F4"],
                WDI["F5"],
            ],
            []
        ]
    ]
]

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    globalXTitle="G3S1 Measurements", localXTitle="Time (min)",
    localYTitle="Fractional Reporter Displaced",
    colorOffset = 2.5,
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/G3 Activity/20230912_BlockerActivatorTrials_G3_NEW",
    spacer=2,
    specialOpts = [
        [
            "Legend",
            [
                [
                    ["5 nM ON", "10 nM ON", "0 nM", "20 nM OFF"]
                ],
                [
                    ["10 nM ON", "20 nM BLK", "20 nM BLK + Act"]
                ]
            ],
            [
                [
                    ["10 nM ON", "20 nM BLK", "20 nM BLK + Act"]
                ],
                [
                    ["10 nM ON", "20 nM BLK", "20 nM BLK + Act"]
                ]
            ]
        ]
    ],
    fixedY=[-0.05, 1.05]
)