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

sys.path.insert(0, '..')
import GSMOptimizationHelperFunctions as GOP

fileName = "C:/Users/kolis/Desktop/Research/LOCALFILES/GSM-Editing/Optimization Examples/20230926_G2G13ActivityWithJunkTemplate_Optimize_Files/20230926_G2G13ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

# All well coords for RNAP.
# TODO automize this along with the junk conc and RNAP conc lists in a well object.
wellCoords = ["I9", "I10", "J9", "J10", "K9", "K10", "L9", "L10", "M9", "M10", "N9", "N10"]
wellNames = []
heatCorrectionTimes = []

timeBreakPts = wellPlate1.getTimeBreaks()

# Lower and upper are time in minutes.
displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = timeBreakPts[0], upperBound = timeBreakPts[1]) # Find the "low" fluorescence
# Plus 5 truncates some of the higher data points which were bad anyway (five of them).
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = timeBreakPts[-1] + 5, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)

# Regular time bounds (for graphing)
#timeBds = [4000/60, 40082/60]

# Time bounds 2 hour (for optimization)
timeBds = [5000/60, (3600 + 3600 + 5000) /60]


timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest
# TODO Actually make this a function
# GOP.MakeWellNames(wellDataOfInterest)
wellsCursedVersion = [
    [
        [[
            WDI["I9"],
            WDI["I10"]
        ]],
        [[
            WDI["J9"],
            WDI["J10"]
        ]],
        [[
            WDI["K9"],
            WDI["K10"]
        ]]
    ],
    [
        [[
            WDI["L9"],
            WDI["L10"]
        ]],
        [[
            WDI["M9"],
            WDI["M10"]
        ]],
        [[
            WDI["N9"],
            WDI["N10"]

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

# For every wellDataOFInterest.
# List that holds approximate max slope, trying to account for noise.

#TODO this could easily be a function
maxSlopes = GOP.FindPercentHighestSlope(wellDataOfInterest,timePtsOfInterest, 90, 10)


#Once we have all of the slope, feed the slope through the equation: slope / (RNAP * normalizedByK_Ms / (1 + sum(normalizedByK_Ms))) to get RNAP rates.
'''for i in maxSlopes:
    i = i*RNAP*normalizedByK_Ms/(1+sum(normalizedByK_Ms))'''
    # RNAP is always 4 here.
    # normalizedByK_Ms = ActNGComp/enzymeProdRates["K_M"]
    # ActNGComp is set of initial amounts for being used for [Genelet,reporter,junk] = [10,0,x] TODO: 0 for reporter?
    # enzymeProdRates["K_M"] = [45,45,45]



# More generalizable version of code above. Lists are in order of wells in order we wish to traverse.
# TODO: populate lists with function based on recipe protocols, and make a "well" object to properly link values like these.
# Then, would only need to call a convert slope on each well object.

RNAPArr = [4,4,4,4,4,4,1,1,4,1,1,1]
JunkCon = np.array([0,320,640,960,1280,1500,0,400,800,120,160,40])
K_M = np.array([45,45,45])


# Pass the exact lists to the function.
RNAPProdRates = GOP.ConvertSlopesToRNAPRates(maxSlopes, RNAPArr, K_M, JunkCon)
print("RNAPProdRates:" , RNAPProdRates)

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    plotXTitles=["I","J","K","L","M","N"], globalXTitle="10 nM Template Type", localXTitle="Time (min)",
    plotYTitles=["", ""], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced",
    colorOffset = 2.5,
    plotFileName="C:/Users/kolis/Desktop/Research/LOCALFILES/GSM-Editing/Optimization Examples/20230921_G1G3ActivityWithJunkTemplate_Optimize_Files/20230921_G1G3ActivityWithJunkTemplateActualG3",
    spacer=2,
    specialOpts = [
        ["Legend", [[["Well 9","Well 10"]]]]
    ],
    fixedY=[-0.05, 1.05]
)
