##### Figure 6 Components #####
#                             #
#      By: Colin Yancey       #
#      Created Sep. 2023      #
#   Last Edited: 09/18/2023   #
#        Schulman Lab         #
#  Johns Hopkins University   #
#                             #
###############################

from pathlib import Path
import sys
# The reason I include this is because I have my "GSM Runner" scripts in a folder, so to import GSM I need to add the parent Path.
sys.path.insert(1, str(Path('../')))
import numpy as np

import GeneletSystemModel as GSM
import ArrayPlots as AP

# Target folder to place images
targFolder = Path("P:/Test Directory")

# Makes a bandpass filter that captures low-strength outputs
newGN = GSM.GeneletNetwork(
    #                        5
    ["G2", "G1"],
    ["C1", ""]
    # When there is something empty in the 2nd line, that probably just means the output doesn't affect anything downstream in the actual circuit.
    # In some cases, that might mean I have intentions to link that node's activity to another genelet circuit and just am using the genelet placeholder as a signalling
    # molecule, similar to what I do in a lot of experiments.
)
newGN.setInitialConditions(
    #                            5
    ["OFF", "BLK"], # Initial annealed conditions of each genelet
    [   50,   60], # Initial concentrations of each genelet
    GSM.createGeneralClassProperties(
        #                   5                         10
        [70, 120], # Activator concentrations (lowest to highest genelet index, ex: the 2 in "G2C1")
        [140,  0]  # Blocker concentrations (lowest to highest genelet index). Counterintuitively, since I named the first genelet G2C1, the activator and blocker
                  # for that node (of which there could be multiple G2 nodes) is actually represented by the second listed activator and blocker number, (120 nM and 0 nM)
    ),
    standardBLK = False, # If this were true, the numbers in blocker list above would have been modified by the script to include additional blocker, 1.5x the
                         # concentration of the total of each genelet index (75 nM for G2__ domains, 90 + 140 nM for G1__ domains). Since this is false, it's just 0 nM and 140 nM.
                         # Regardless of what this is set to, there the initial annealed conditions must be satisfied. I.e., if a G1__ node is BLK, there must be
                         # enough total blocker in solution to exceed the concentration of G1__ BLK or the script will throw an error. This is satisfied as 140 nM > 60 nM.
    
    RNAP = 3.57 # Nowadays, I use 4 U/uL of RNAP. The script's default is 3.57 U/uL but this may change soon. Also, the way we calculate production, as discussed,
                # is in the process of getting overhauled anyway.
)


simTime = 2 # Time in hours that the simulations will run
simTimeVals = GSM.getTimeValues(simTime) # Returns the time values in seconds (since that's what GSM runs on)

numOfSims = 8
maxActivator = 140

simConcList = [maxActivator*i/(numOfSims-1) for i in np.arange(numOfSims)] # I want to plot a bunch of different starting activator values and see what happens.

G2C1Outputs = [] # Initialize groupings of outputs
G1Outputs = []

for i in np.arange(numOfSims):
    newGN.modifyValue("Activator", 1, simConcList[i]) # I'm working on this but basically that 1 indicates the index in the list, so G2's activator (2nd index)...sorry about that.
    newGN.simulate(simTimeVals)
    G2C1Outputs.append(newGN.OutputConcentrations["Act: G2->C1"]) # Gives me the data collected at getTimeValues' points of ON G2C1.
    G1Outputs.append(newGN.OutputConcentrations["Act: G1->"]) # Gives me the data collected at getTimeValues' points of ON G1__.

# Place the outputs into the multi-layered vector to prepare them for AP.generatePlots()
plot4DMat = [ # 1st layer, separated by row
    [ # 2nd layer, separated by column
        [ # 3rd layer, separated by hues (red, green blue)
            G2C1Outputs # 4th layer (this vector itself), separated by shade (dark blue, light blue)
        ],
        [
            G1Outputs
        ]
    ]
]

AP.generatePlots(simTimeVals/60, plot4DMat,
    plotXTitles=["G2C1 ON", "G1__ ON"], globalXTitle="Genelet Observed", localXTitle="Time (min)",
    localYTitle="Concentration (nM)",
    plotFileName=targFolder / "Blocked_To_Activated", # Name of the file + the folder you want it in. Usually I indicate this elsewhere (top of script) but sometimes I write it here.
                                                      # Default fileSuffix is ".png", can change that though with fileSuffix = "[Whatever you want]"
    specialOpts = [
        ["Legend",
            [ # Same principle as the plot4DMat variable, except only 2 layers (row, column) since the names are just in order of what was plotted to each subplot.
                [
                    [str(int(x))+" nM Act 2" for x in simConcList]
                ]
            ]
        ] 
    ],
    fixedY=(-2, 62) # Locks the y axes, useful in cases where you want to plot between 0 and 1 and don't want confusing axes changing it based on the actual maximum values.
)
