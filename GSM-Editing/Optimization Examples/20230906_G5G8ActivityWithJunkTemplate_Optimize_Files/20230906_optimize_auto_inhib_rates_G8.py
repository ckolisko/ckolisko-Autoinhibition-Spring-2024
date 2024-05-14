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
import GSMRateEvaluation as rE

sys.path.insert(0, '..')
import GSMOptimizationHelperFunctions as GOP


fileName = "C:/Users/kolis/Desktop/Research/LOCALFILES/GSM-Editing/Optimization Examples/20230906_G5G8ActivityWithJunkTemplate_Optimize_Files/20230906_G5G8ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

# All well coords for G5,
wellCoords = ["J12", "J13", "K12", "K13", "L12", "L13", "M12", "M13", "N12", "N13", "O12", "O13"]
wellNames = []
heatCorrectionTimes = []

# Lower and upper are time in minutes.
displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2000/60, upperBound = 4200/60) # Find the "low" fluorescence
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 30900/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)

# Regular time bounds
timeBds = [4000/60, 30900/60]

#Time bound for first hour of activity.
#timeBds = [4000/60, (4000 + 3600)/60]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

# G1 is 7,8 wells.
WDI = wellDataOfInterest
wellsCursedVersion = [
    [
        [[
            WDI["J12"],
            WDI["J13"]
        ]],
        [[
            WDI["K12"],
            WDI["K13"]
        ]],
        [[
            WDI["L12"],
            WDI["L13"]
        ]]
    ],
    [
        [[
            WDI["M12"],
            WDI["M13"]
        ]],
        [[
            WDI["N12"],
            WDI["N13"]
        ]],
        [[
            WDI["O12"],
            WDI["O13"]

        ]]
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

'''
AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    plotXTitles=["B","C","D","E","F","G"], globalXTitle="10 nM Template Type", localXTitle="Time (min)",
    plotYTitles=["", ""], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced",
    colorOffset = 2.5,
    plotFileName="C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230921_G1G3ActivityWithJunkTemplate_Optimize_Files/20230921_G1G3ActivityWithJunkTemplateGraph130",
    spacer=2,
    specialOpts = [
        ["Legend", [[["Well 7","Well 8","Well 9","Well 10"]]]]
    ],
    fixedY=[-0.05, 1.05]
)
'''

# AP.generatePlots(timePtsOfInterest, wellsActualStuff,
#     localXTitle="Time (hrs)",
#     localYTitle="Fractional Reporter Displaced",
#     colorOffset = 2.5,
#     plotFileName="P:/PhD Work/Assets/Fluorescence Plots/G1 Activity/DELETE_TEST_12",
#     spacer=2
# )

# TODO: remove exit.
# sys.exit() #Temporary exit.
# From the recipe protocol, Junk Template sect, in nM.

#J12,J13
junkTempConc = [0,320,640,960,1280,1500,0,40,80,120,160,240]
RNAPArr = [4,4,4,4,4,4,1,1,1,1,1,1]

# Rip these from the list produced by slope optimizer.
# TODO hook up into the function, this stuff is wrong.

#g3 10i 90d tb-120
RNAPRates = [0.342478421701603, 2.436587833164253, 0.9628542888671645, 0.5954340119161966, 0.5479569594643415, 0.5712974308458666, 1.3116318535580918, 1.0117713004484306, 0.5513237950375506, 0.4574437391806116, 0.35702075702075664, 0.19609666222085428]

numGeneletNetworks = len(RNAPRates) # This line is for this specific setup, modify for different traversals.

gNs = []

for i in range(numGeneletNetworks):
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
        [  10,   500, junkTempConc[i]],
        GSM.createGeneralClassProperties(
            # Activator conc vec.
            # 2x on genelet = double activator of genelet amount.
            # goes g1, g2, g3.
            [0,  20, junkTempConc[i]],
            # Blocker conc vec. 0, except for reporter or specifically blocked already.
            [500, 0, 0]
        ),
        RNAP = RNAPArr[i], # T7 RNAP section. careful to iterate over all wells.
        RNaseH=0,
        standardBLK=False # Pretty much always off.
    )
    # Lock in RNAP Prod for each one.
    gNs[-1].modifyRate("RNAPProdRate", RNAPRates[i] , 0) # genelet K_M



# Rate constant for S1 -> binding is fine since I am using Coact -> Block-Genelet as the proxy reaction, which has a rate constant of 5 * 10^(-6) 1/nM/s vs.
# 5.71 * 10^(-6) empirically fitted for 15 mM MgCl, 2 mM NTPs

#TODO this controls time bounds for optimizations.
OptimizationTime = 480 #Optimize for first 480 minutes.
timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[0]+OptimizationTime) #first OptimizationTime min.
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
timePtsOfInterest *= 60 
WDI = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[0]+OptimizationTime)

iterationParameter = 1000

lbound1 = -7
hbound1 = 1
lbound2 = -7
hbound2 = 1
lbound3 = -7
hbound3 = 1
#GOP.OptimizePermutations(gNs, timePtsOfInterest, numGeneletNetworks, WDI, iterationParameter
#                         ,lbound1,hbound1,lbound2,hbound2,lbound3,hbound3)

rateEvalObj = rE.GSMRateEvaluation(
    gNs,
    [timePtsOfInterest]*numGeneletNetworks,
    # Fitting data
    [
        WDI["J12"], WDI["J13"], 
        WDI["K12"], WDI["K13"], 
        WDI["L12"], WDI["L13"], 
        WDI["M12"], WDI["M13"], 
        WDI["N12"], WDI["N13"], 
        WDI["O12"], WDI["O13"], 
    ],
    # Blocker C1 complex
    "Blk-C1",
    # Output scaler scales up fraction on.
    500,
    # Rates of interest.
    [   
     #TODO: We want AutoInhib-Act-Repressor, AutoInhib-Free-Act, and AutoInhib-Act-Genelet
     # Checklist: 
         # 1. RAG
         # 2. RA-
         # 3. R-G
         # 4. -AG
         # 5. R--
         # 6. --G
         # 7. -A-
     
        #{
        #    "RateType"   : "AutoInhib-Act-Repressor",
        #    "Target"     : None,
        #    "Low Bound"  : lbound1,
        #    "High Bound" : hbound1
        #},
        {
            "RateType"   : "AutoInhib-Free-Act",
            "Target"     : None,
            "Low Bound"  : lbound2,
            "High Bound" : hbound2
        },
        #{
        #    "RateType"   : "AutoInhib-Act-Genelet",
        #    "Target"     : 1,
        #    "Low Bound"  : lbound3,
        #    "High Bound" : hbound3
        #},
    ]
)
rateEvalObj.addTimeDelay()

rateEvalObj.optimizeRates(maxiter=iterationParameter) #TODO: Check time for 1 simulation. import time, time . time, have a start time, and then initialize, run 100 simulations, then stop time = time.time, then calculate total time.

print("Optimized Rates:\n", rateEvalObj.OptimizedRates)
print("\nOptimized Time Delays:\n", rateEvalObj.OptimizedTimeDelays)
print("\nOptimization Score: ", rateEvalObj.calcPredCost(rateEvalObj.OptimizedRatesTable))