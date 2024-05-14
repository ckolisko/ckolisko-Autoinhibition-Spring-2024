from pathlib import Path
import sys
sys.path.insert(1, str(Path('../../')))

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import AssayAnalysis as AA
import ArrayPlots as AP
import GeneralUse as GU
import GeneletSystemModel as GSM
fileName = "C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230921_G1G3ActivityWithJunkTemplate_Optimize_Files/20230921_G1G3ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

# All well coords for RNAP, except need to ignore M5.
wellCoords = ["B7", "B8", "B9", "B10", "C7", "C8", "C9", "C10", "D7", "D8", "D9", "D10", "E7", "E8", "E9", "E10", "F7", "F8", "F9", "F10", "G7", "G8", "G9", "G10"]
wellNames = []
heatCorrectionTimes = []

displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2000/60, upperBound = 4200/60) # Find the "low" fluorescence
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 46000/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)


timeBds = [4000/60, 45090/60]

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
            WDI["B9"],
            WDI["B10"],
            WDI["C7"],
            WDI["C8"],
            WDI["C9"],
            WDI["C10"],
            WDI["D7"],
            WDI["D8"],
            WDI["D9"],
            WDI["D10"]
        ]],
        [[]],
    ],
    [
        [[
            WDI["E7"],
            WDI["E8"],
            WDI["E9"],
            WDI["E10"],
            WDI["F7"],
            WDI["F8"],
            WDI["F9"],
            WDI["F10"],
            WDI["G7"],
            WDI["G8"],
            WDI["G9"],
            WDI["G10"]
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
    plotFileName="C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230921_G1G3ActivityWithJunkTemplate_Optimize_Files/20230921_G1G3ActivityWithJunkTemplateGraph3",
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

# sys.exit() #Temporary exit.

# From the recipe protocol, Junk Template sect, in nM.
junkTempConc = [0, 320, 640, 960, 1280, 1500]
junkTempConc2 = [0, 40, 80, 120, 160, 240]
gNs = []
# Make 2 copies of each of these, for the different genes, g1s1 and g3s1.
# Wells for RNAP = 4,
for i in range(2):
    for x in junkTempConc:
        gNs.append(GSM.GeneletNetwork(
            #                        5
            # g2 - genelet, g1 - blocker (reporter), g3 - junk
            ["G2", "G1", "G3"],
            ["C1",   "", ""]
        ))
        
        gNs[-1].setInitialConditions(
            #                            5
            # Junk template initial concentration is changing.
            # Initial genelet of 0.0 nM, which goes in bulk, which goes in well.
            # Units in uM
            # S1 F-Q pair is the reporter
            # Junk template is always on.
            ["ON", "BLK", "ON"],
            [  10,   500, x],
                GSM.createGeneralClassProperties(
                # Activator conc vec.
                # 2x on genelet = double activator of genelet amount.
                # goes g1, g2, g3.
                [0,  20, x],
                # Blocker conc vec. 0, except for reporter or specifically blocked already.
                [500, 0, 0]
            ),
            RNAP=4, # T7 RNAP section. careful to iterate over all wells.
            RNaseH=0,
            standardBLK=False # Pretty much always off.
        )
        gNs[-1].modifyRate("EnzymeK_M", 10**1.4593927494173307 , 0) # genelet K_M
        gNs[-1].modifyRate("EnzymeK_M", 10**1.5273483461728665 , 2) # junk K_M
        gNs[-1].modifyRate("RNAPDeath", 10**-3.6033040110165597 , None)
        gNs[-1].modifyRate("DTTDeath", 10**-5.1828082321551445 , None)
        gNs[-1].modifyRate("RNAPRevival", 10**-6.651737126790722e-05 , None)

# Wells for RNAP = 1,
for i in range(2):
    for x in junkTempConc2:
        gNs.append(GSM.GeneletNetwork(
            #                        5
            # g2 - genelet, g1 - blocker (reporter), g3 - junk
            ["G2", "G1", "G3"],
            ["C1",   "", ""]
        ))
        
        gNs[-1].setInitialConditions(
            #                            5
            # Junk template initial concentration is changing.
            # Initial genelet of 0.0 nM, which goes in bulk, which goes in well.
            # Units in uM
            # S1 F-Q pair is the reporter
            # Junk template is always on.
            ["ON", "BLK", "ON"],
            [  10,   500, x],
                GSM.createGeneralClassProperties(
                # Activator conc vec.
                # 2x on genelet = double activator of genelet amount.
                # goes g1, g2, g3.
                [0,  20, x],
                # Blocker conc vec. 0, except for reporter or specifically blocked already.
                [500, 0, 0]
            ),
            RNAP=1, # T7 RNAP section. careful to iterate over all wells.
            RNaseH=0,
            standardBLK=False # Pretty much always off.
        )
        # TODO: Proper locking? Add to actual genelet model? "Lock in" everything except enzyme prod.
        gNs[-1].modifyRate("EnzymeK_M", 10**1.4593927494173307 , 0) # genelet K_M
        gNs[-1].modifyRate("EnzymeK_M", 10**1.5273483461728665 , 2) # junk K_M
        gNs[-1].modifyRate("RNAPDeath", 10**-3.6033040110165597 , None)
        gNs[-1].modifyRate("DTTDeath", 10**-5.1828082321551445 , None)
        gNs[-1].modifyRate("RNAPRevival", 10**-6.651737126790722e-05 , None)

# Rate constant for S1 -> binding is fine since I am using Coact -> Block-Genelet as the proxy reaction, which has a rate constant of 5 * 10^(-6) 1/nM/s vs.
# 5.71 * 10^(-6) empirically fitted for 15 mM MgCl, 2 mM NTPs

import GSMRateEvaluation as rE

numTimePoints = 90
timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[0]+ numTimePoints) #first 90 timepoints.
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
timePtsOfInterest *= 60 
WDI = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[0]+ numTimePoints)

rateEvalObj = rE.GSMRateEvaluation(
    gNs,
    [timePtsOfInterest]*len(junkTempConc)*4,
    # Fitting data
    [
        WDI["B7"], WDI["B8"], WDI["C7"], WDI["C8"], WDI["D7"], WDI["D8"],
        WDI["B9"], WDI["B10"], WDI["C9"], WDI["C10"], WDI["D9"], WDI["D10"],
        WDI["E7"], WDI["E8"], WDI["F7"], WDI["F8"], WDI["G7"], WDI["G8"],
        WDI["E9"], WDI["E10"], WDI["F9"], WDI["F10"], WDI["G9"], WDI["G10"],
    ],
    # Blocker C1 complex
    "Blk-C1",
    # Output scaler scales up fraction on.
    500,
    # Rates of interest. Lock in everything except enzyme production.
    [
        {               
            "RateType"   : "EnzymeProduction",
            "Target"     : 0, # Maps to order genelets were put in. This is g2-c1, aka genelet.
            "Low Bound"  : -3, # = 10^-3 (mL/unit) * (nM/s)
            "High Bound" : 3
        },
    ]
)
rateEvalObj.addTimeDelay()

rateEvalObj.optimizeRates(maxiter=1000)

print("Optimized Rates:\n", rateEvalObj.OptimizedRates)
print("\nOptimized Time Delays:\n", rateEvalObj.OptimizedTimeDelays)