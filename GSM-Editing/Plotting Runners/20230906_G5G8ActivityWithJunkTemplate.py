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

fileName = "P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230906_G5G8ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

getNormalized = False

wellCoords = ["J10", "J11", "J12", "J13", "K10", "K11", "K12", "K13", "L10", "L11", "L12", "L13", "M10", "M11", "N10", "N11", "O10", "O11", "M12", "M13", "N12", "N13", "O12", "O13"]
wellNames = []
heatCorrectionTimes = []

if getNormalized:
    displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2000/60, upperBound = 4500/60) # Find the "low" fluorescence
    scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 31000/60, upperBound = None) # Find the "high" fluorescence
    for wellCoord in wellCoords:
        wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

    #AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

    if len(heatCorrectionTimes)>0:
        AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [4500/60, 31000/60 if getNormalized else None]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest

wellsCursedVersion = [
    [
        [[
            WDI["J10"],
            WDI["J11"],
            WDI["K10"],
            WDI["K11"],
            WDI["L10"],
            WDI["L11"]
        ]],
        [[
            WDI["J12"],
            WDI["J13"],
            WDI["K12"],
            WDI["K13"],
            WDI["L12"],
            WDI["L13"]
        ]]
    ],
    [
        [[
            WDI["M10"],
            WDI["M11"],
            WDI["N10"],
            WDI["N11"],
            WDI["O10"],
            WDI["O11"]
        ]],
        [[
            WDI["M12"],
            WDI["M13"],
            WDI["N12"],
            WDI["N13"],
            WDI["O12"],
            WDI["O13"]
        ]]
    ]
]

# AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
#     plotXTitles=["G5S1", "G8S1"], globalXTitle="10 nM Template Type", localXTitle="Time (min)",
#     plotYTitles=["4 U/uL", "1 U/uL"], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced" if getNormalized else "Raw Fluoresence",
#     colorOffset = 2.5,
#     plotFileName="P:/PhD Work/Assets/Fluorescence Plots/T7 Tests/20230906_G5G8ActivityWithJunkTemplate" + ("" if getNormalized else "_RAW"),
#     spacer=2,
#     specialOpts = [
#         [
#             "Legend",
#             [
#                 [
#                     ["0 nM", "320", "640", "960", "1280", "1500"]
#                 ],
#                 [
#                     ["0 nM", "40 ", "80", "120", "160", "240"]
#                 ]
#             ]
#         ]
#     ],
#     fixedY=[-0.05, 1.05] if getNormalized else None
# )



junkTempConc = [[0,0], [320, 40], [640, 80], [960, 120], [1280, 160], [1500, 240]]
gNs = []

for i, y in enumerate([4,1]):
    for x in junkTempConc:
        gNs.append(GSM.GeneletNetwork(
            #                        5
            ["G2", "G1", "G3"],
            ["C1",   "", ""]
        ))
    
        gNs[-1].setInitialConditions(
            #                            5
            ["ON", "BLK", "ON"],
            [  10,   500, x[i]],
            GSM.createGeneralClassProperties(
                [0,  20, x[i]],
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
        WDI["J10"], WDI["J11"],
        WDI["K10"], WDI["K11"],
        WDI["L10"], WDI["L11"],
        WDI["M10"], WDI["M11"],
        WDI["N10"], WDI["N11"],
        WDI["O10"], WDI["O11"]
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
            "Target"     : 1,
            "Low Bound"  : -1,
            "High Bound" : 4
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