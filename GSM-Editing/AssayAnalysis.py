#### Assay Analysis Tools #####
#                             #
#      By: Colin Yancey       #
#      Created Jan. 2022      #
#   Last Edited: 7/20/2023    #
#        Schulman Lab         #
#  Johns Hopkins University   #
#                             #
###############################

import pandas as pd
import numpy as np
from datetime import timedelta
import re
import sys

import ArrayPlots as AP
import GeneralUse as GU

unitChoices = {
    "sec": 1, # Number of seconds in a second
    "min": 60, # Number of seconds in a minute
    "hr": 3600 # Number of seconds in an hour
}

heatCorrectionValues = {
    "TYE665":{
        "Cyt5":{ # Determined from 500 nM of TYE665 fluorophore on CYT5 machine (06/21/23)
            "AmplitudeOvershoot": 0.35,
            "InvTau": 0.09
        }
    }
}

# Returns a list of all well names in the specified rectangle on a well plate
def getWellNames(topLeftName="A1", bottomRightName="P24"):
    letr1 = ord(topLeftName[0])
    letr2 = ord(bottomRightName[0])

    num1 = int(topLeftName[1:])
    num2 = int(bottomRightName[1:])

    return [chr(i)+str(j) for i in np.arange(letr1, letr2+1) for j in np.arange(num1, num2+1)]


# Performs getSteadyState of object wellPlate for a list of wellCoords. See wellPlate.getSteadyState for additional information.
def getSteadyStates(wellPlate, wellCoords, lowerBound=0, upperBound=None):
    ssTab = {}
    for x in wellCoords:
        ssTab[x] = wellPlate.getSteadyState(x, lowerBound, upperBound)
    return ssTab


# Performs transformWellData of object wellPlate for a list of wellCoords. See wellPlate.transformWellData for additional information.
def transformMultiWellData(wellPlate, wellCoords, displaceValue=0, scalingValue=1):
    for x in wellCoords:
        wellPlate.transformWellData(x, displaceValue, scalingValue)


# Performs invertWellData of object wellPlate for a list of wellCoords. See wellplate.invertWellData for additional information.
def invertMultiWellData(wellPlate, wellCoords, reflectVal=None):
    for x in wellCoords:
        wellPlate.invertWellData(x, reflectVal)


# Performs imposeHeatCorrections of object wellPlate for a list of wellCoords. See wellPlate.imposeHeatCorrections for additional information.
def imposeMultiHeatCorrections(wellPlate, wellCoords, timeStamps):
    for x in wellCoords:
        wellPlate.imposeHeatCorrections(x, timeStamps)


# Performs getWellData of object wellPlate for a list of wellCoords. See wellPlate.getWellData for additional information.
def getMultiWellData(wellPlate, wellCoords, lowerBound=0, upperBound=None):
    wDTab = {}
    for x in wellCoords:
        wDTab[x] = wellPlate.getWellData(x, lowerBound, upperBound)
    return wDTab


####### ----- WELL PLATE DATA TRANSFORMATION CLASS ----- #######

class WellPlate:
    
    # Initializes the well plate from the output txt file of the plate reader.
    # fileName: String that gives the file path for the data of interest.
    # timeFormat: The time (first column) will be assumed to be in these units when making calculations.
    # timeUnits: Converts the time in the collected raw data from the timeFormat units to the timeUnits.
        # All time values given hereafter are assumed to be in these units.
    # fluorophore: Dye of interest being observed.
    # plateReader: Plate reader making the measurement.
        # You don't need to give the fluorophore or plateReader if you don't plan on using heat correction.
    # header: Index of the row which contains the title information (i.e., the time label and well names).
        # All files are assumed to only have a time column and all subsequent columns should be well data (make sure this is set in the plate reader's txt output).
    def __init__(self, fileName, timeFormat="sec", timeUnits="min", fluorophore="TYE665", plateReader="Cyt5", header=0):

        tempRaw = pd.read_table(fileName, header=header)

        self.DataColumns = tempRaw.columns[[re.match("^[A-Z][0-9]+$", x)!=None for x in tempRaw.columns]]
        tempRaw = tempRaw.dropna(how='all', subset=self.DataColumns).dropna(axis=1, how='all')

        if timeFormat=="sec":
            self.TimeData = np.array([t/unitChoices[timeUnits] for t in tempRaw["Time"].to_numpy()])
        elif timeFormat=="hrMinSec":
            hrMinSec = np.array([t.split(':') for t in tempRaw["Time"].to_numpy()])
            self.TimeData = timedelta(hours=int(hrMinSec[0]), minutes=int(hrMinSec[1]), seconds=int(hrMinSec[2])).total_seconds()/unitChoices[timeUnits]
        else:
            sys.exit("ERROR: Indicated timeFormat cannot be reconciled with the file's time format.")

        self.RawData = tempRaw.drop("Time", axis=1)
        self.OverflowData = {}

        for coord in self.DataColumns:
            for curValInd, val in enumerate(self.RawData[coord]):
                if val!="OVRFLW":
                    self.RawData[coord][curValInd] = float(val)
                else:
                    print("WARNING: Coordinate "+coord+" contains OVRFLW values.")
                    self.OverflowData[coord] = curValInd-1
                    break
            if not coord in self.OverflowData:
                self.OverflowData[coord] = len(self.RawData[coord])
        
        self.HeatCorrection = heatCorrectionValues[fluorophore][plateReader].copy()


    # Obtains the midpoints in each time break that exceeds the difference threshold
    # timeDiffThreshold: The necessary time difference between two data recordings to be considered and recorded as a time break (in minutes)
    # distInBetween: Linearly interpolated point in between the two times, from 0 (earlier time) to 1 (later time), defaults to 0.5 (midpoint)
    def getTimeBreaks(self, timeDiffThreshold=2, distInBetween=0.5):
        if len(self.TimeData)<=1: # Only one or less data points, thus no breaks.
            return []
        return [self.TimeData[i] + distInBetween*(x - self.TimeData[i]) for i, x in enumerate(self.TimeData[1:]) if (x - self.TimeData[i])>=timeDiffThreshold]


    # Obtains the steady state value (assumed as the "settled" fluorescence value) for a given well in a given time duration.
    # wellCoord: The coordinate of the well (ex: "B2").
    # lowerBound: The lower bound of the recording times to be used in the calculation.
    # upperBound: The upper bound of the recording times to be used in the calculation.
    def getSteadyState(self, wellCoord, lowerBound=0, upperBound=None):
        lowerInd, upperInd = self.getTimeIndexes(lowerBound, upperBound)
        thresholdData = list(self.RawData[wellCoord][lowerInd:upperInd])
        thresholdData.sort()
        return thresholdData[len(thresholdData)//10] # Get the 10th percentile value (lowest to reduce thermal issues, but not too low that it could be a noisy point)


    # Performs an offset and applies a scaling term to a given well's data.
    # wellCoord: The coordinate of the well (ex: "B2").
    # displaceValue: The number to be subtracted from each fluoroscence recording (offsets the data).
    # scalingValue: The number to be divided from each fluorescence recording (can renormalize the data).
    def transformWellData(self, wellCoord, displaceValue=0, scalingValue=1):
        self.RawData.loc[:self.OverflowData[wellCoord], (wellCoord)] -= displaceValue
        self.RawData.loc[:self.OverflowData[wellCoord], (wellCoord)] /= scalingValue


    # Performs a data inversion around a particular value. Useful if your fluorescence values increase when your species of interest decreases.
    # wellCoord: The coordinate of the well (ex: "B2").
    # reflectVal: The value which all data is reflected over. If none, assumes the midpoint of the min and max of the current data.
        # For example, if you wish to invert data which you have normalized to be between 0 and 1, you would choose 0.5.
    def invertWellData(self, wellCoord, reflectVal=None):
        curData = self.RawData.loc[:self.OverflowData[wellCoord], (wellCoord)]
        if reflectVal==None:
            reflectVal = (curData.max() + curData.min)
        self.RawData.loc[:self.OverflowData[wellCoord], (wellCoord)] = 2*reflectVal - curData


    # Attempts to undo the change in fluorescence caused by the plate warming up when placed in the plate reader by dividing out an exponential decay.
    # wellCoord: The coordinate of the well (ex: "B2").
    # timeStamps: A list of the times shortly after which the plate was placed back into the reader.
        # Make sure to at least perform an offset such that the fluorescence is correctly normalized without noise.
    def imposeHeatCorrections(self, wellCoord, timeStamps):
        timeStamps.sort()
        for i in range(len(timeStamps)):
            lowBd, highBd = self.getTimeIndexes(timeStamps[i], None if (i+1)>=len(timeStamps) else timeStamps[i+1]-0.0000001)
            ptsFilter = np.arange(lowBd, highBd)
            timePts = GU.filterList(self.TimeData, ptsFilter)
            self.RawData.loc[ptsFilter, (wellCoord)] /= 1+self.HeatCorrection["AmplitudeOvershoot"]*np.exp(-self.HeatCorrection["InvTau"]*(timePts-timePts[0]))


    # Obtains the time values which happened between the two given bounds.
    # The lower bound of the recording times to be obtained.
    # The upper bound of the recording times to be obtained.
    def getTimeData(self, lowerBound=0, upperBound=None):
        lowBd, highBd = self.getTimeIndexes(lowerBound, upperBound)
        ptsFilter = np.arange(lowBd, highBd)
        return np.array(GU.filterList(self.TimeData, ptsFilter))


    # Obtains the time values which happened between the two given bounds.
    # wellCoord: The coordinate of the well (ex: "B2").
    # The lower time bound of the recordings to be obtained.
    # The upper time bound of the recordings to be obtained.
    def getWellData(self, wellCoord, lowerBound=0, upperBound=None):
        lowBd, highBd = self.getTimeIndexes(lowerBound, upperBound)
        ptsFilter = np.arange(lowBd, highBd)
        return np.array(GU.filterList(list(self.RawData[wellCoord]), ptsFilter))


    # Obtains the indexes of the first time that happened at or after the lower and at or before the upper bound of times.
    # The lower bound of the recording time indexes to be obtained. If none, assumes the minimum time (0).
    # The upper bound of the recording time indexes to be obtained. If none, assumes the maximum time.
    def getTimeIndexes(self, lowerBound=0, upperBound=None):
        if upperBound==None:
            upperBound = self.TimeData.max()
        return (
            np.where(self.TimeData==np.min([i for i in self.TimeData if i>=lowerBound]))[0][0],
            np.where(self.TimeData==np.max([i for i in self.TimeData if i<=upperBound]))[0][0]+1
        )


    # Adjusts the heat correction modifier by manually assigning an expected change in fluorescence and time constant for the exponential decay.
    # ampCorrection: The percent expected change from the starting fluorescence to the fully heated state's fluorescence.
    # invTau: The time it takes to reach (1 - 1/e) or roughly ~63% of the change in fluorescence (with 100% in this case being the fully heated state).
    def adjustHeatCorrection(self, ampCorrection=None, tau=None):
        if ampCorrection != None:
            self.HeatCorrection["AmplitudeOvershoot"] = ampCorrection
        if tau != None:
            self.HeatCorrection["InvTau"] = 1/tau