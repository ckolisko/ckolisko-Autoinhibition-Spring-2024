from pathlib import Path
import sys
sys.path.insert(1, str(Path('../')))

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import AssayAnalysis as AA
import ArrayPlots as AP
import GeneralUse as GU
import GeneletSystemModelOriginal as GSM

fileName = "P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230828_G1S1ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

wellCoords = ["J2", "J3", "J4", "K2", "K3", "K4", "L2", "L3", "L4", "M2", "M3", "M4", "N2", "N3", "N4", "O2", "O3", "O4"]
wellNames = []
heatCorrectionTimes = []

displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2000/60, upperBound = 4000/60) # Find the "low" fluorescence
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 34800/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [4000/60, 33130/60]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest
wellsCursedVersion = [
    [
        [[
            WDI["J2"],
            WDI["J3"],
            WDI["J4"],
            WDI["K2"],
            WDI["K3"],
            WDI["K4"],
            WDI["L2"],
            WDI["L3"],
            WDI["L4"]
        ]],
        [[]],
    ],
    [
        [[
            WDI["M2"],
            WDI["M3"],
            WDI["M4"],
            WDI["N2"],
            WDI["N3"],
            WDI["N4"],
            WDI["O2"],
            WDI["O3"],
            WDI["O4"]
        ]],
        [[]]
    ]
    
]

# wellsActualStuff = [
#     [
#         [[WDI["L8"]]],
#         [[WDI["L9"]]],
#         [[WDI["L10"]]],
#         [[WDI["L11"]]]
#     ],
#     [
#         [],
#         [],
#         [],
#         []
#     ]
# ]

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    plotXTitles=["G1S1"], globalXTitle="10 nM Template Type", localXTitle="Time (min)",
    plotYTitles=["4 U/uL", "1 U/uL"], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced",
    colorOffset = 2.5,
    plotFileName="P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230828_G1S1ActivityWithJunkTemplate_Group",
    spacer=2,
    specialOpts = [
        ["Legend", [[["0 nM", "10", "20", "40", "80", "160", "320", "640", "1280"]]]]
    ],
    fixedY=[-0.05, 1.05]
)

# AP.generatePlots(timePtsOfInterest, wellsActualStuff,
#     localXTitle="Time (hrs)",
#     localYTitle="Fractional Reporter Displaced",
#     colorOffset = 2.5,
#     plotFileName="P:/PhD Work/Assets/Fluorescence Plots/G1 Activity/DELETE_TEST_12",
#     spacer=2
# )

sys.exit()
junkTempConc = [0, 10, 20, 40, 80, 160, 320, 640, 1280]
gNs = []

for y in [4,1]:
    for x in junkTempConc:
        gNs.append(GSM.GeneletNetwork(
            #                        5
            ["G2", "G1", "G3"],
            ["C1",   "", ""]
        ))
    
        gNs[-1].setInitialConditions(
            #                            5
            ["ON", "BLK", "ON"],
            [  10,   500, x],
            GSM.createGeneralClassProperties(
                [0,  20, x],
                [500, 0, 0]
            ),
            RNAP=y,
            RNaseH=0,
            standardBLK=False
        )

# Rate constant for S1 -> binding is fine since I am using Coact -> Block-Genelet as the proxy reaction, which has a rate constant of 5 * 10^(-6) 1/nM/s vs.
# 5.71 * 10^(-6) empirically fitted for 15 mM MgCl, 2 mM NTPs

import GSMRateEvaluation as rE


timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[0]+150)
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
timePtsOfInterest *= 60
WDI = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[0]+150)

rateEvalObj = rE.GSMRateEvaluation(
    gNs,
    [timePtsOfInterest]*len(junkTempConc)*2,
    [
        WDI["J2"], WDI["J3"], WDI["J4"],
        WDI["K2"], WDI["K3"], WDI["K4"],
        WDI["L2"], WDI["L3"], WDI["L4"],
        WDI["M2"], WDI["M3"], WDI["M4"],
        WDI["N2"], WDI["N3"], WDI["N4"],
        WDI["O2"], WDI["O3"], WDI["O4"]
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
            "RateType"   : "EnzymeK_M",
            "Target"     : 2,
            "Low Bound"  : -1,
            "High Bound" : 4
        },
        {
            "RateType"   : "EnzymeDecay",
            "Target"     : None,
            "Low Bound"  : 4,
            "High Bound" : 8
        }
    ]
)
rateEvalObj.addTimeDelay()

rateEvalObj.optimizeRates(maxiter=1000)

print("Optimized Rates:\n", rateEvalObj.OptimizedRates)
print("\nOptimized Time Delays:\n", rateEvalObj.OptimizedTimeDelays)