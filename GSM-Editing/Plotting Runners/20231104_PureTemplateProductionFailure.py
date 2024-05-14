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

fileName = "20231104_PureTemplateProductionFailure"
filePath = "P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/"
wellPlate1 = AA.WellPlate(filePath+fileName+".txt")

getNormalized = True
hideEarlyData = True
hideLateData = True

invertWellData = False
clrOffset = 2.5 # 1.7 for deep blue, 2.5 for reddish pink

wellCoords = AA.getWellNames("I12", "N15")

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
timeBds = [timeBreakPts[-3] if hideEarlyData else 0, timeBreakPts[-2] if hideLateData else None]


timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsGraphSetup = [[
    [
        [
            [
                WDI["I12"],
                WDI["I13"],
                WDI["J12"],
                WDI["J13"]
            ],
            [
                WDI["I14"],
                WDI["I15"],
                WDI["J14"],
                WDI["J15"]
            ]
        ],
    ]
], [
    [
        [
            [
                WDI["K12"],
                WDI["K13"],
                WDI["L12"],
                WDI["L13"]
            ],
            [
                WDI["K14"],
                WDI["K15"],
                WDI["L14"],
                WDI["L15"]
            ]
        ]
    ]
], [
    [
        [
            [
                WDI["M12"],
                WDI["M13"],
                WDI["N12"],
                WDI["N13"]
            ],
            [
                WDI["M14"],
                WDI["M15"]
            ]
        ]
    ]
]]

modName = ["1 mM NTP + 5 mM MgCl2 Added", "4 U/uL T7 RNAP Added", "5 nM Pure Template Added"]
fileNameAdd = ["NTPMix", "RNAP", "PureTemp"]

for i in np.arange(len(modName)):
    AP.generatePlots(timePtsOfInterest, wellsGraphSetup[i],
        plotXTitles=[modName[i]], localXTitle="Time (min)", globalXTitle="Legend: Pure Template Type and Junk Template Conc.",
        plotYTitles=[], localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
        colorOffset = clrOffset,
        plotFileName=filePath + fileName + "_" + fileNameAdd[i] + ("" if getNormalized else "_RAW"),
        spacer=2,
        specialOpts = [
            [
                "Legend",
                [
                    [
                        ["G1S1, 640 nM j.t.", "960 nM", "1280 nM", "1600 nM", "G5S1, 640 nM j.t.", "960 nM", "1280 nM", "1600 nM"]
                    ]
                ]
            ]
        ],
        fixedY=[-0.05, 1.05] if getNormalized else None
    )