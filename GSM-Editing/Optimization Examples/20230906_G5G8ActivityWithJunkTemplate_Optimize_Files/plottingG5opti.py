# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:27:57 2024

@author: kolis
"""

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
    
# Generate genelet networks with found RNAP rates to compare to actual data.

# Idea: What if you could just take lists of well data and place them here, in B7, B8, B9 order.
#B7,B8, B9,B10
junkTempConc = [0,320,640,960,1280,1500,0,40,80,120,160,240]
RNAPArr = [4,4,4,4,4,4,1,1,1,1,1,1]

# Rip these from the list produced by slope optimizer.

RNAPRates = [0.3271724351329866, 2.235364757104659, 1.4995707734144543, 1.1884351739690229, 0.9817036037542094, 0.6969121798749746, 1.6846193302557344, 1.598750405712431, 1.3100632172992932, 0.6747656999483433, 0.4021771662460364, 0.20514242456251544]
numGeneletNetworks = len(RNAPRates) # This line is for this specific setup, modify for different traversals.

# Actual graphing part
fileName = "C:/Users/kolis/Desktop/Research/LOCALFILES/GSM-Editing/Optimization Examples/20230906_G5G8ActivityWithJunkTemplate_Optimize_Files/20230906_G5G8ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

# All well coords for RNAP.
# TODO automize this along with the junk conc and RNAP conc lists in a well object.
wellCoords = ["J10", "J11", "K10", "K11", "L10", "L11", "M10", "M11", "N10", "N11", "O10", "O11"]
wellNames = []
heatCorrectionTimes = []

# Lower and upper are time in minutes.
displacer = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 2000/60, upperBound = 5200/60) # Find the "low" fluorescence
scaler = AA.getSteadyStates(wellPlate1, wellCoords, lowerBound = 30900/60, upperBound = None) # Find the "high" fluorescence
for wellCoord in wellCoords:
    wellPlate1.transformWellData(wellCoord, displaceValue = displacer[wellCoord], scalingValue = scaler[wellCoord]) # Perform the normalization

#AA.invertMultiWellData(wellPlate1, wellCoords, reflectVal = 0.5) # Invert around 0.0 - 1.0, comment out if not needed

if len(heatCorrectionTimes)>0:
    AA.imposeMultiHeatCorrections(wellPlate1, wellCoords, heatCorrectionTimes)

wellData = AA.getMultiWellData(wellPlate1, wellCoords, lowerBound=0, upperBound=None)

# Regular time bounds
timeBds = [5000/60, 30900/60]

# Time bounds 1 hour (for graphing only)
#timeBds = [(4000+1200)/60, (3600 + 3600 + 4000 + 1200)/60]



timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

WDI = wellDataOfInterest
#(wellPlate, wellCoords, displaceValue=0, scalingValue=1)


# End of actual graphing copy paste

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

    # Optimized rate lock in section:
    # R does nothing here.     -2.539171784388269 
    gNs[-1].modifyRate("AutoInhib-Act-Repressor", 10**-6.427817427132339 , 1) # genelet K_M
    gNs[-1].modifyRate("AutoInhib-Free-Act", 10**-3.175447035224141, 1) # genelet K_M
    gNs[-1].modifyRate("AutoInhib-Act-Genelet", 10**-5.199293611462969 , 1) # genelet K_M

# Find score of simulated genelet network.

#Simulation Graphs
targFolder = Path("C:/Users/kolis/Desktop/Research/LOCALFILES/GSM-Editing/Optimization Examples/20230906_G5G8ActivityWithJunkTemplate_Optimize_Files/OptimizedGraphs")


simTimeVals = GSM.getTimeValues(timePtsOfInterest[-1] / 60, dataPointNum = len(timePtsOfInterest)) # Returns the time values in seconds (since that's what GSM runs on)

numOfSims = numGeneletNetworks

G2C1Outputs = [] # Initialize groupings of outputs
BlkC1Outputs = []

# Simulate all the genelet networks.
# TODO: Is this correct at all.
for i in np.arange(numOfSims):
    gNs[i].simulate(simTimeVals)
    #G2C1Outputs.append(gNs[i].OutputConcentrations["Act: G2->C1"]) # Gives me the data collected at getTimeValues' points of ON G2C1.
    BlkC1Outputs.append(gNs[i].OutputConcentrations["Blk-C1"]) # Gives me the data collected at getTimeValues' points of ON G1__.


# Place the outputs into the multi-layered vector to prepare them for AP.generatePlots()
plot4DMat = [ # 1st layer, separated by row
    [ # 2nd layer, separated by column
        [[ # 3rd layer, separated by hues (red, green blue)
            BlkC1Outputs[0],
            BlkC1Outputs[1],
            # 4th layer (this vector itself), separated by shade (dark blue, light blue)
        ],[
            WDI["J10"]*500,
            WDI["J11"]*500
        ]],
        [[
            BlkC1Outputs[2],
            BlkC1Outputs[3]
        ],[
           WDI["K10"]*500,
           WDI["K11"]*500 
        ]],
        [[
            BlkC1Outputs[4],
            BlkC1Outputs[5]
        ],[
           WDI["L10"]*500,
           WDI["L11"]*500 
        ]]
    ],
    [
        [[
            BlkC1Outputs[6],
            BlkC1Outputs[7]
        ],[
           WDI["M10"]*500,
           WDI["M11"]*500 
        ]],
        [[
            BlkC1Outputs[8],
            BlkC1Outputs[9]
        ],[
           WDI["N10"]*500,
           WDI["N11"]*500 
        ]],
        [[
            BlkC1Outputs[10],
            BlkC1Outputs[11]

        ],[
           WDI["O10"]*500,
           WDI["O11"]*500 
        ]]
    ]
    
]

# TODO, plot real data on top of this.

AP.generatePlots(simTimeVals/60, plot4DMat,
    plotXTitles=["J", "K","L","M","N","O"], globalXTitle="Genelet Observed", localXTitle="Time (min)",
    localYTitle="Concentration (nM)",
    plotFileName=targFolder / "20230921_opti_graph__Combine_A_", # Name of the file + the folder you want it in. Usually I indicate this elsewhere (top of script) but sometimes I write it here.
                                                      # Default fileSuffix is ".png", can change that though with fileSuffix = "[Whatever you want]"
    specialOpts = [
        ["Legend",
            [ # Same principle as the plot4DMat variable, except only 2 layers (row, column) since the names are just in order of what was plotted to each subplot.
                [
                    ["Simulated Well 10", "Simulated Well 11","Real Well 10","Real Well 11",]
                ]
            ]
        ] 
    ],
    fixedY=(-20, 550) # Locks the y axes, useful in cases where you want to plot between 0 and 1 and don't want confusing axes changing it based on the actual maximum values.
)
