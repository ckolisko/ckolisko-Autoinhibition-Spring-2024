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
import GSMOptimizationHelperFunctions as GOP

    

fileName = "C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230921_G1G3ActivityWithJunkTemplate_Optimize_Files/20230921_G1G3ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

# All well coords for RNAP.
# TODO automize this along with the junk conc and RNAP conc lists in a well object.
wellCoords = ["B7", "B8", "B9", "B10", "C7", "C8", "C9", "C10", "D7", "D8", "D9", "D10", "E7", "E8", "E9", "E10", "F7", "F8", "F9", "F10", "G7", "G8", "G9", "G10"]
wellNames = []
heatCorrectionTimes = []

# Lower and upper are time in minutes.
displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2000/60, upperBound = 4200/60) # Find the "low" fluorescence
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 46000/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)

# Regular time bounds
timeBds = [4000/60, 45090/60]

# Time bounds 1 hour (for graphing only)
# timeBds = [(4000+1200)/60, (3600 + 3600 + 4000 + 1200)/60]


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
            WDI["B7"],
            WDI["B8"],
            WDI["B9"],
            WDI["B10"]
        ]],
        [[
            WDI["C7"],
            WDI["C8"],
            WDI["C9"],
            WDI["C10"]
        ]],
        [[
            WDI["D7"],
            WDI["D8"],
            WDI["D9"],
            WDI["D10"]
        ]]
    ],
    [
        [[
            WDI["E7"],
            WDI["E8"],
            WDI["E9"],
            WDI["E10"]
        ]],
        [[
            WDI["F7"],
            WDI["F8"],
            WDI["F9"],
            WDI["F10"]
        ]],
        [[
            WDI["G7"],
            WDI["G8"],
            WDI["G9"],
            WDI["G10"]

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
# List that holds the fifth max slope.
MaxSlopes5th = GOP.Find5thMax(wellDataOfInterest,timePtsOfInterest)


#Once we have all of the slope, feed the slope through the equation: slope / (RNAP * normalizedByK_Ms / (1 + sum(normalizedByK_Ms))) to get RNAP rates.
'''for i in maxSlopes:
    i = i*RNAP*normalizedByK_Ms/(1+sum(normalizedByK_Ms))'''
    # RNAP is always 4 here.
    # normalizedByK_Ms = ActNGComp/enzymeProdRates["K_M"]
    # ActNGComp is set of initial amounts for being used for [Genelet,reporter,junk] = [10,0,x] TODO: 0 for reporter?
    # enzymeProdRates["K_M"] = [45,45,45]


'''
#TODO: Fix the 4 or 1.
RNAPProdRates = []
RNAP = 4 # For a B7, B8, B9, traversal starts at 4 for this experiment.
K_M = np.array([45,45,45])

# For a B7, B8, B9 traversal. Takes the exact lists from recipe protocol.
JunkCon = np.array([0,320,640,960,1280,1500,0,40,80,120,160,240])
JunkHalf = len(JunkCon)

counter = 0

junk_index = 0
for i in MaxSlopes5th:
    ActNGComp = np.array([10,0,JunkCon[junk_index]])
    # Every other junkindex increments.
    if (counter % 2 == 1):
        junk_index += 1
    normalizedByK_Ms = ActNGComp/K_M
    if (RNAP == 4 and counter == JunkHalf): # Change this for diff setups.
        RNAP = 1
    RNAPProd = ((i * 500 / 60) *  (1 + sum(normalizedByK_Ms))) / (RNAP * normalizedByK_Ms[0])
    print(RNAPProd)
    RNAPProdRates.append(RNAPProd)
    counter += 1
print(RNAPProdRates)
'''
# More generalizable version of code above. Lists are in order of wells in order we wish to traverse.
# TODO: populate lists with function based on recipe protocols, and make a "well" object to properly link values like these.
# Then, would only need to call a convert slope on each well object.
RNAP = [4,4,4,4,4,4,4,4,4,4,4,4,1,1,1,1,1,1,1,1,1,1,1,1] # For a B7, B8, B9, traversal starts at 4 for this experiment.

# TODO Messed up junk array below.
# JunkCon = np.array([0, 0, 320,320,640,640,960,960,1280,1280,1500,1500,0,0,40,40,80,80,120,120,160,160,240,240])
JunkCon = np.array([0,320,0,320, 640,960,640,960, 1280,1500,1280,1500, 0,40,0,40, 80,120,80,120,160,240,160,240])
K_M = np.array([45,45,45])


# Pass the exact lists to the function.
RNAPProdRates = GOP.ConvertSlopesToRNAPRates(MaxSlopes5th, RNAP, K_M, JunkCon)
print(RNAPProdRates)

AP.generatePlots(timePtsOfInterest, wellsCursedVersion,
    plotXTitles=["B","C","D","E","F","G"], globalXTitle="10 nM Template Type", localXTitle="Time (min)",
    plotYTitles=["", ""], globalYTitle="RNAP Concentration", localYTitle="Fractional Reporter Displaced",
    colorOffset = 2.5,
    plotFileName="C:/Users/kolis/Desktop/Research/LOCALFILES/GSM-Editing/Optimization Examples/20230921_G1G3ActivityWithJunkTemplate_Optimize_Files/20230921_G1G3ActivityWithJunkTemplateGraph120MINActual",
    spacer=2,
    specialOpts = [
        ["Legend", [[["Well 7","Well 8","Well 9","Well 10"]]]]
    ],
    fixedY=[-0.05, 1.05]
)
