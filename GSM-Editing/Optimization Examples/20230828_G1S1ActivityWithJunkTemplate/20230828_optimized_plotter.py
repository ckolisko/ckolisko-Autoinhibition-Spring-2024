from pathlib import Path
import sys
# The reason I include this is because I have my "GSM Runner" scripts in a folder, so to import GSM I need to add the parent Path.
sys.path.insert(1, str(Path('../../')))
import numpy as np

import GeneletSystemModel as GSM
import ArrayPlots as AP

# Target folder to place images
targFolder = Path("C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230828_G1S1ActivityWithJunkTemplate")
junkTempConc = [0, 10, 20, 40, 80, 160, 320, 640, 1280]
gNs = []

for y in [4,1]:
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
            RNAP=y, # T7 RNAP section. careful to iterate over all wells.
            RNaseH=0,
            standardBLK=False # Pretty much always off.
        )
        # Modify Rates with what was found in optimized run. TODO modify the rates for EP, EK_M's, deaths, rev.
        # def modifyRate(self, rateType, newValue, target=None):
        gNs[-1].modifyRate("EnzymeProduction", 10**8.860303436386519e-05 , 0)
        gNs[-1].modifyRate("EnzymeK_M", 10**1.4593927494173307 , 0) # genelet K_M
        gNs[-1].modifyRate("EnzymeK_M", 10**1.5273483461728665 , 2) # junk K_M
        gNs[-1].modifyRate("RNAPDeath", 10**-3.6033040110165597 , None)
        gNs[-1].modifyRate("DTTDeath", 10**-5.1828082321551445 , None)
        gNs[-1].modifyRate("RNAPRevival", 10**-6.651737126790722e-05 , None)



simTime = 4 # Time in hours that the simulations will run
simTimeVals = GSM.getTimeValues(simTime) # Returns the time values in seconds (since that's what GSM runs on)

numOfSims = 18

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
        [ # 3rd layer, separated by hues (red, green blue)
            G2C1Outputs # 4th layer (this vector itself), separated by shade (dark blue, light blue)
        ],
        [
            BlkC1Outputs
        ]
    ]
]

AP.generatePlots(simTimeVals/60, plot4DMat,
    plotXTitles=["G2C1 ON", "G1__ ON"], globalXTitle="Genelet Observed", localXTitle="Time (min)",
    localYTitle="Concentration (nM)",
    plotFileName=targFolder / "20230828_opti_graph4", # Name of the file + the folder you want it in. Usually I indicate this elsewhere (top of script) but sometimes I write it here.
                                                      # Default fileSuffix is ".png", can change that though with fileSuffix = "[Whatever you want]"
    specialOpts = [
        ["Legend",
            [ # Same principle as the plot4DMat variable, except only 2 layers (row, column) since the names are just in order of what was plotted to each subplot.
                [
                    [str(int(x))+" nM Act 2" for x in junkTempConc]
                ]
            ]
        ] 
    ],
    fixedY=(-20, 520) # Locks the y axes, useful in cases where you want to plot between 0 and 1 and don't want confusing axes changing it based on the actual maximum values.
)
