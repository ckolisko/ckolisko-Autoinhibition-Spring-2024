# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:53:01 2024

@author: kolis
"""
import numpy as np
import math
import GSMRateEvaluation as rE
import itertools


# Finds the 5th max slopes of a list of wells, using timePtsOfInterest.
# wellDataOfInterest is the 2D array of wells and their corresponding data.
# timePtsOfInterest is the timepoints that correspond to the time well data was taken, used to find slope.
def Find5thMax(wellDataOfInterest, timePtsOfInterest):
    numWells = len(wellDataOfInterest)
    if numWells <= 0:
        raise ValueError("Cannot find slope of less then 1 wells.")
    MaxSlopes5th = [0]*numWells
    count = 0

    for i in wellDataOfInterest:
        timeIndex = 0
        #First, second, third, and fourth highest slopes.
        slope1 = 0
        slope2 = 0
        slope3 = 0
        slope4 = 0
        slope5 = 0
        # Iterate through all points of interest, calculating average slope between point and 5 points ahead. 
        # Then, take fifth largest slope.
        # y is fraction florecence, x is minutes.
        for j in range(len(wellDataOfInterest[i]) - 5): # minus 1 avoids out of bounds.
            yCoord1 = wellDataOfInterest[i][j]
            yCoord2 = wellDataOfInterest[i][j + 5]        
            xCoord1 = timePtsOfInterest[timeIndex]
            xCoord2 = timePtsOfInterest[timeIndex + 5]

            timeIndex += 1
            slope = (yCoord2-yCoord1) / (xCoord2-xCoord1)
            # Swap downs if necessary.
            if(slope > slope1):
                temp = slope1
                slope1 = slope
                slope = temp
            if(slope > slope2):
                temp = slope2
                slope2 = slope
                slope = temp
            if(slope > slope3):
                temp = slope3
                slope3 = slope
                slope = temp
            if(slope > slope4):
                temp = slope4
                slope4 = slope
                slope = temp
            if (slope > slope5):
                slope5 = slope
        MaxSlopes5th[count] = slope5
        count += 1
    return MaxSlopes5th

# Finds the percent heighest number, as an alternative to the old 
# version of finding slops, which just found the 5th steepest slope.
# Percent is whole integer
def FindPercentHighestSlope(wellDataOfInterest, timePtsOfInterest, percent, interval):
    numWells = len(wellDataOfInterest)
    if numWells <= 0:
        raise ValueError("Cannot find slope of less then 1 wells.")
        
    bestSlopes = [0]*numWells
    count = 0

    # for each well in well data of interest
    for i in wellDataOfInterest: #TODO make number indexing
        curWellSlopes = [0] * (len(wellDataOfInterest[i]) - interval)
        # for this well, iterate through all points in well, calc average 
        # slopes between current point and 5 points ahead, storing in slopes. 
        # y is fraction of florecence, x is minutes.
        for j in range(len(wellDataOfInterest[i]) - interval): # minus 1 avoids out of bounds.
            yCoord1 = wellDataOfInterest[i][j]
            yCoord2 = wellDataOfInterest[i][j + interval]        
            xCoord1 = timePtsOfInterest[j]
            xCoord2 = timePtsOfInterest[j + interval]

            curWellSlopes[j] = (yCoord2-yCoord1) / (xCoord2-xCoord1)
        # Sort curWellSlopes, take the bottom of 10th decile.
        curWellSlopes.sort()
        decileIndex = int(math.floor(len(curWellSlopes) * (percent / 100)))
        bestSlopes[count] = curWellSlopes[decileIndex]
        count += 1
    return bestSlopes
    

# Converts slopes in units of proportion per minute to units of nM*uL/U RNAP/s
# SlopeList is the list of slopes for each well.
# RNAP is the RNAP concentration in each well.
# K_M is a vector of K_M's for g1, g2, and g3 in the genelet network.
# JunkCon is the list of JunkTemplate concentrations for each well.
def ConvertSlopesToRNAPRates(SlopeList, RNAP, K_M, JunkCon):
    RNAPProdRates = []
    for i in range(len(SlopeList)):
        ActNGComp = np.array([10,0,JunkCon[i]])
        normalizedByK_Ms = ActNGComp/K_M
        RNAPProd = ((SlopeList[i] * 500 / 60) *  (1 + sum(normalizedByK_Ms))) / (RNAP[i] * normalizedByK_Ms[0])
        RNAPProdRates.append(RNAPProd)
        
    return RNAPProdRates

# Takes all of the parameters from original optimization function, but cycles 
# through all permutations to get all optimized rates.
def OptimizePermutations(gNs, timePtsOfInterest, numGeneletNetworks, WDI, iterat,lbound1,hbound1,lbound2,hbound2,lbound3,hbound3):
    
    # call all permutation types, from code 0 to code 6, as defined below. 
    # Codes:
        # 1. RAG = 0
        # 2. R-G = 1
        # 3. RA- = 2
        # 4. -AG = 3
        # 5. R-- = 4
        # 6. --G = 5
        # 7. -A- = 6
    optimizedRatesList = [0] * 7
    optimizedTimeDelayList = [0] * 7
    optimizedRatesListScores = [0] * 7
    rateEvalObjList = getRateEvalObjects(gNs, timePtsOfInterest, numGeneletNetworks, WDI,lbound1,hbound1,lbound2,hbound2,lbound3,hbound3)

    for i in range(7):
        rateEvalObj = rateEvalObjList[i]
        rateEvalObj.addTimeDelay()
        # Actually run simulation.
        rateEvalObj.optimizeRates(maxiter=iterat)
        optimizedRatesList[i] = rateEvalObj.OptimizedRates
        optimizedTimeDelayList[i] = rateEvalObj.OptimizedTimeDelays
        optimizedRatesListScores[i] = rateEvalObj.calcPredCost(rateEvalObj.OptimizedRatesTable)

        
    
    # Print lists out, with the best scores as well.
    for i in range(7):
        print("Permutation:", i , "\nOptimized Rates:\n", optimizedRatesList[i])
        print("\nOptimized Time Delays:\n", optimizedTimeDelayList[i], "\n")
        print("\nOptimization Score: ", optimizedRatesListScores[i], "\n")
    
    
# Creates a rate evaluation object for the given permutation code.
# TODO This does not work
def getRateEvalObjects(gNs, timePtsOfInterest, numGeneletNetworks, WDI,lbound1,hbound1,lbound2,hbound2,lbound3,hbound3):
    ratesDictionary = [
        {
            "RateType"   : "AutoInhib-Act-Repressor",
            "Target"     : None,
            "Low Bound"  : lbound1,
            "High Bound" : hbound1
        },
        {
            "RateType"   : "AutoInhib-Free-Act",
            "Target"     : None,
            "Low Bound"  : lbound2,
            "High Bound" : hbound2
        },
        {
            "RateType"   : "AutoInhib-Act-Genelet",
            "Target"     : 1,
            "Low Bound"  : lbound3,
            "High Bound" : hbound3
        },
    ]
    # Initialize list of dictionaries for rates to be optimized based on permutation code.
    RateEvalList = []
    for i in range(3,0,-1):
        combos = itertools.combinations(ratesDictionary, i)
        for j in combos:
            RateEvalList.append(j)

    
    for i in range(7):
        RateEvalList[i] = rE.GSMRateEvaluation(
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
            RateEvalList[i]    
        )
    return RateEvalList

