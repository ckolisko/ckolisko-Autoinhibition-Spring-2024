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

# Rip these from the list produced by slope optimizer.

RNAPRates= [0.3119543920229689, 0.012922613765222109, 2.217006299266806, 0.20570012159058532, 0.10939502054522143, 0.08579132016632016, 0.21162042557419777, 0.1657993606048584, 0.12816820276497698, 0.142545136866628, 0.19376311647429198, 0.18177537006359235, 0.3915061071894448, 0.12187450608503238, 0.7622074613718066, 0.21329618723514362, 0.10192375573596896, 0.11393114325440462, 0.13141695028877876, 0.07505635460511861, 0.09489514705249126, 0.05343519699163618, 0.11139515455304927, 0.07225374535657429]
RNAPRatesHalfIndex = len(RNAPRates) // 2 # This line is for this specific setup, modify for different traversals.

gNs = []
RNAPRatesIndex = 0
for x in junkTempConc:
    # Depends on specific experiment, but here we have pairs of 4 until switch to 1 at well 12.
    T7 = 4 if RNAPRatesIndex < RNAPRatesHalfIndex else 1
    for y in [T7, T7]:
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
            RNAP = T7, # T7 RNAP section. careful to iterate over all wells.
            RNaseH=0,
            standardBLK=False # Pretty much always off.
        )
        # Lock in RNAP Prod for each one.
        gNs[-1].modifyRate("RNAPProdRate", RNAPRates[RNAPRatesIndex] , 0) # genelet K_M
        
        # temp auto repress test.
        gNs[-1].modifyRate("AutoInhib-Free-Act", 0.001111111111111 , 1) # genelet K_M
        gNs[-1].modifyRate("AutoInhib-Act-Genelet", 0.001111111111111 , 1) # genelet K_M
        gNs[-1].modifyRate("AutoInhib-Act-Repressor", 0.001111111111111 , 1) # genelet K_M



        RNAPRatesIndex += 1

#Simulation Graphs
targFolder = Path("C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230921_G1G3ActivityWithJunkTemplate_Optimize_Files")

simTime = 2 # Time in hours that the simulations will run
simTimeVals = GSM.getTimeValues(simTime) # Returns the time values in seconds (since that's what GSM runs on)

numOfSims = 24

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
            BlkC1Outputs[2],
            BlkC1Outputs[3]# 4th layer (this vector itself), separated by shade (dark blue, light blue)
        ]],
        [[
            BlkC1Outputs[4],
            BlkC1Outputs[5],
            BlkC1Outputs[6],
            BlkC1Outputs[7]
        ]],
        [[
            BlkC1Outputs[8],
            BlkC1Outputs[9],
            BlkC1Outputs[10],
            BlkC1Outputs[11]
        ]]
    ],
    [
        [[
            BlkC1Outputs[12],
            BlkC1Outputs[13],
            BlkC1Outputs[14],
            BlkC1Outputs[15]
        ]],
        [[
            BlkC1Outputs[16],
            BlkC1Outputs[17],
            BlkC1Outputs[18],
            BlkC1Outputs[19]
        ]],
        [[
            BlkC1Outputs[20],
            BlkC1Outputs[21],
            BlkC1Outputs[22],
            BlkC1Outputs[23]

        ]]
    ]
    
]



AP.generatePlots(simTimeVals/60, plot4DMat,
    plotXTitles=["B", "C","D","E","F","G"], globalXTitle="Genelet Observed", localXTitle="Time (min)",
    localYTitle="Concentration (nM)",
    plotFileName=targFolder / "20230921_opti_graph120MINPredicted", # Name of the file + the folder you want it in. Usually I indicate this elsewhere (top of script) but sometimes I write it here.
                                                      # Default fileSuffix is ".png", can change that though with fileSuffix = "[Whatever you want]"
    specialOpts = [
        ["Legend",
            [ # Same principle as the plot4DMat variable, except only 2 layers (row, column) since the names are just in order of what was plotted to each subplot.
                [
                    ["Well 7", "Well 8","Well 9","Well 10",]
                ]
            ]
        ] 
    ],
    fixedY=(-20, 520) # Locks the y axes, useful in cases where you want to plot between 0 and 1 and don't want confusing axes changing it based on the actual maximum values.
)
