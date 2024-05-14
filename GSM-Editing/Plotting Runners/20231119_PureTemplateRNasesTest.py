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

fileName = "20231119_PureTemplateRNasesTest"
filePath = "P:/PhD Work/Assets/Fluorescence Plots/RNase Tests/"
wellPlate1 = AA.WellPlate(filePath+fileName+".txt")
manualAddn = ""

getNormalized = True
hideEarlyData = False
hideLateData = False

invertWellData = False
clrOffset = 2.5 # 1.7 for deep blue, 2.5 for reddish pink

wellCoords = AA.getWellNames("B21", "J23")

timeBreakPts = wellPlate1.getTimeBreaks()
heatCorrectionTimes = []

if getNormalized:
    displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = timeBreakPts[0], upperBound = timeBreakPts[1]) # Find the "low" fluorescence
    scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = timeBreakPts[-2]+5, upperBound = timeBreakPts[-1]) # Find the "high" fluorescence
    for wellCoord in wellCoords:
        wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

    if invertWellData:
        AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0

    if len(heatCorrectionTimes)>0:
        AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)
timeBds = [timeBreakPts[-3] if hideEarlyData else 0, timeBreakPts[-2] if hideLateData else timeBreakPts[-1]]


timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

# wellsGraphSetup = [
#     [
#         [
#             [
#                 WDI["B21"]
#             ],
#             [
#                 WDI["C21"],
#                 WDI["D21"],
#                 WDI["E21"],
#                 WDI["F21"]
#             ]
#         ],
#         [
#             [
#                 WDI["B22"]
#             ],
#             [
#                 WDI["C22"],
#                 WDI["D22"],
#                 WDI["E22"],
#                 WDI["F22"]
#             ]
#         ],
#         [
#             [
#                 WDI["B23"]
#             ],
#             [
#                 WDI["C23"],
#                 WDI["D23"],
#                 WDI["E23"],
#                 WDI["F23"]
#             ]
#         ]
#     ],
#     [
#         [
#             [
#                 WDI["B21"]
#             ],
#             [
#                 WDI["G21"],
#                 WDI["H21"],
#                 WDI["I21"],
#                 WDI["J21"]
#             ]
#         ],
#         [
#             [
#                 WDI["B22"]
#             ],
#             [
#                 WDI["G22"],
#                 WDI["H22"],
#                 WDI["I22"],
#                 WDI["J22"]
#             ]
#         ],
#         [
#             [
#                 WDI["B23"]
#             ],
#             [
#                 WDI["G23"],
#                 WDI["H23"],
#                 WDI["I23"],
#                 WDI["J23"]
#             ]
#         ]
#     ]
# ]

# AP.generatePlots(timePtsOfInterest, wellsGraphSetup,
#     plotXTitles = ["1 U/uL", "2 U/uL", "4 U/uL"], localXTitle="Time (min)", globalXTitle = "RNAP Concentration",
#     plotYTitles = ["RNase H", "RNase A/T1"], localYTitle = "Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence", globalYTitle = "Type of RNase (Qty in Legend)",
#     colorOffset = clrOffset,
#     plotFileName = filePath + fileName + (manualAddn if (getNormalized or hideEarlyData or hideLateData) else "_RAW"),
#     spacer = 2,
#     specialOpts = [
#         [
#             "Legend",
#             [
#                 [
#                     ["0x", "1x", "2x", "4x", "8x"]
#                 ]
#             ]
#         ]
#     ],
#     fixedY=[-0.05, 1.05] if getNormalized else None
# )

rnaseAConcs = [0, 0.000025, 0.00005, 0.0001, 0.0002]
gNs = []

for i, x in enumerate(rnaseAConcs):
    gNs.append(GSM.GeneletNetwork(
        #                        5
        ["G2", "G1", "G3"],
        ["C1",   "", ""]
    ))

    gNs[-1].setInitialConditions(
        #                            5
        ["ON", "BLK", "ON"],
        [   5,   500, 100],
        GSM.createGeneralClassProperties(
            [0,  20, 100],
            [500, 0, 0]
        ),
        RNAP=1,
        RNaseH=0,
        RNaseA=x,
        standardBLK=False
    )



import GSMRateEvaluation as rE


timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[0]+300)
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
timePtsOfInterest *= 60
WDI = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[0]+300)

rateEvalObj = rE.GSMRateEvaluation(
    gNs,
    [timePtsOfInterest]*len(rnaseAConcs),
    [
        WDI["B21"],
        WDI["G21"],
        WDI["H21"],
        WDI["I21"],
        WDI["J21"]
    ],
    "Blk-C1",
    500,
    [
        {
            "RateType"   : "EnzymeProduction",
            "Target"     : 0,
            "Low Bound"  : -3,
            "High Bound" : 3
        },
        {
            "RateType"   : "RNaseA",
            "Target"     : None,
            "Low Bound"  : -1,
            "High Bound" : 5
        }
    ]
)
rateEvalObj.addTimeDelay()

rateEvalObj.optimizeRates(maxiter=1000)

print("Optimized Rates:\n", rateEvalObj.OptimizedRates)
print("\nOptimized Time Delays:\n", rateEvalObj.OptimizedTimeDelays)