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


fileName = "C:/Users/kolis/Desktop/Research/GSM-Editing/Optimization Examples/20230921_G1G3ActivityWithJunkTemplate_Optimize_Files/20230921_G1G3ActivityWithJunkTemplate.txt"
wellPlate1 = AA.WellPlate(fileName)

# All well coords for RNAP, except need to ignore M5.
wellCoords = ["B9", "B10", "C9", "C10", "D9", "D10", "E9", "E10", "F9", "F10", "G9", "G10"]
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

#Time bound for first hour of activity.
#timeBds = [4000/60, (4000 + 3600)/60]

timePtsOfInterest = wellPlate1.getTimeData(lowerBound=timeBds[0], upperBound=timeBds[1])
timePtsOfInterest -= min(timePtsOfInterest) # Set the time to 0 at the start of the run of interest
wellCoordsOfInterest = wellCoords # Whichever wells you actually want to fit

# Make sure this is normalized properly (FRACTION ON or analogous)
wellDataOfInterest = AA.getMultiWellData(wellPlate1, wellCoordsOfInterest, lowerBound=timeBds[0], upperBound=timeBds[1])

# G3 is 9,10 wells.
WDI = wellDataOfInterest
wellsCursedVersion = [
    [
        [[
            WDI["B9"],
            WDI["B10"]
        ]],
        [[
            WDI["C9"],
            WDI["C10"]
        ]],
        [[
            WDI["D9"],
            WDI["D10"]
        ]]
    ],
    [
        [[
            WDI["E9"],
            WDI["E10"]
        ]],
        [[
            WDI["F9"],
            WDI["F10"]
        ]],
        [[
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

#B9,B10
junkTempConc = [0,320,640,960,1280,1500,0,40,80,120,160,240]
RNAPArr = [4,4,4,4,4,4,1,1,1,1,1,1]

# Rip these from the list produced by slope optimizer.
# TODO hook up into the function, this stuff is wrong.

# g3 5i 5b
#RNAPRates = [0.3251609238924649, 0.2056939461715125, 0.1449026559350418, 0.1657993606048584, 0.1663453214763515, 0.18177537006359193, 0.4412780039520986, 0.21329618723514362, 0.10137486947441714, 0.0750518630814648, 0.08117978842131289, 0.07225374535657426]
# g3 10i 5b
#RNAPRates = [0.3238396125205099, 0.15585362284033058, 0.10235676972432792, 0.10263769942205542, 0.09849394034783968, 0.10040925203512709, 0.43669499797633443, 0.1891468152510595, 0.06678558997563547, 0.04218031064655974, 0.04130829941792476, 0.03913744540147772]
#g3 5i 90d tb-120
#RNAPRates = [0.3184163195759182, 0.11820951980666167, 0.07152641739772306, 0.08191277934644799, 0.0758768133050025, 0.09175328203209891, 0.435549246482393, 0.19570762729297814, 0.056996171249565, 0.041307614564217124, 0.03663866557068106, 0.04214801812466831]
#g3 5i 95d tb-120
#RNAPRates = [0.3246481761958855, 0.12688759963373827, 0.07460945263038359, 0.0898079869942983, 0.08755016919807976, 0.10906522203815536, 0.4398048948884603, 0.19989537966016005, 0.0648277062304212, 0.045961993670044415, 0.04094909681429049, 0.04616211508892245]
#g3 10i 90d tb-120
RNAPRates = [0.3132099583491102, 0.11211140965790538, 0.06042749056014538, 0.06266821070481249, 0.05613067923567275, 0.06318858102210588, 0.43055704354450874, 0.1828651867002863, 0.05264531848242257, 0.03403514721136201, 0.024066574443486578, 0.021074009062334163]
#g3 10i 95d tb-120
#RNAPRates = [0.3238396125205099, 0.11750589171257449, 0.06721016807199844, 0.07500447265457887, 0.061285118438655875, 0.07184455102513405, 0.43669499797633443, 0.1891468152510595, 0.05590845805777934, 0.037235032846618235, 0.03196903172343739, 0.02458634390605653]

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
        WDI["B9"], WDI["B10"], 
        WDI["C9"], WDI["C10"], 
        WDI["D9"], WDI["D10"], 
        WDI["E9"], WDI["E10"], 
        WDI["F9"], WDI["F10"], 
        WDI["G9"], WDI["G10"], 
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
        #{
        #    "RateType"   : "AutoInhib-Free-Act",
        #    "Target"     : None,
        #    "Low Bound"  : lbound2,
        #    "High Bound" : hbound2
        #},
        {
            "RateType"   : "AutoInhib-Act-Genelet",
            "Target"     : 1,
            "Low Bound"  : lbound3,
            "High Bound" : hbound3
        },
    ]
)
rateEvalObj.addTimeDelay()

rateEvalObj.optimizeRates(maxiter=iterationParameter) #TODO: Check time for 1 simulation. import time, time . time, have a start time, and then initialize, run 100 simulations, then stop time = time.time, then calculate total time.

print("Optimized Rates:\n", rateEvalObj.OptimizedRates)
print("\nOptimized Time Delays:\n", rateEvalObj.OptimizedTimeDelays)
print("\nOptimization Score: ", rateEvalObj.calcPredCost(rateEvalObj.OptimizedRatesTable))