from scipy.optimize import minimize
import sys
import numpy as np
import copy
import matplotlib.pyplot as plt
import GeneletSystemModel as GSM


class GSMRateEvaluation:

    # Initializes the GSM rate fitting instance.
    # geneletNetworks: A GSM GeneletNetwork or list of GeneletNetwork instances.
    # timePoints: A numpy array, or a list of multiple numpy arrays of time points (in sec).
    # fittingData: A numpy array, or a list of multiple numpy arrays of corresponding data to the time points.
    # outputValue: The molecule/complex name in question to be measured when performing the fit.
    # outputScalar: The normalization to be applied so that the output values match in expected scale with the fitting data.
    # ratesOfInterest: A list containing multiple dictionary entrieos detailing each rate's name and upper/lower legal values (as powers of 10).
    def __init__(self, geneletNetworks, timePoints, fittingData, outputValue, outputScaler, ratesOfInterest):
        
        if isinstance(geneletNetworks, list):
            self.DataCount = len(geneletNetworks)
            if self.DataCount == sum([1 for x in geneletNetworks if isinstance(x, object)]):
                if self.DataCount == len(timePoints) and self.DataCount == len(fittingData):
                    self.GeneletNetworks = geneletNetworks # List of genelet networks created with matching conditions to the data
                    self.TimePoints = timePoints # List of time point sets recorded
                    self.FittingData = fittingData # List of empirical data sets recorded
                else:
                    sys.exit("Multiple genelet networks should be paired with the same number of time and data instances.")
            else:
                sys.exit("All networks must be GeneletNetwork instances.")
        elif isinstance(geneletNetworks, object):
            self.DataCount = 1
            if isinstance(timePoints, list) and isinstance(fittingData, list): # Reformat the data for general use in this class
                self.GeneletNetworks = [geneletNetworks]
                self.TimePoints = [timePoints]
                self.FittingData = [fittingData]
            else:
                sys.exit("A single genelet network should be paired with only one instance of time and data.")
        else:
            sys.exit("Incorrect genelet network format. Networks must be provided as a GeneletNetwork instance. Multiple networks must be provided in a list.")

        self.OutputValue = outputValue
        self.OutputScaler = outputScaler

        self.RatesOfInterest = ratesOfInterest


    # Adds a time delay variable to all fits to account for time misalignment with when the experimental data was taken (not immediately put in plate reader)
    # timeDelayLB: Lower bound of the time delay. Default is 0 seconds.
    # timeDelayUB: Upper bound of the time delay. Default is 600 seconds. (ex: Experiment was spiked 0 to 600 seconds before being placed into the reader.)
    def addTimeDelay(self, timeDelayLB = 0, timeDelayUB = 600):
        self.TimeDelayBounds = (timeDelayLB, timeDelayUB)


    def updateRates(self, newRateVals):
        for i, x in enumerate(self.RatesOfInterest):
            for gN in self.GeneletNetworks:
                # Nelder-Mead fit has a difficult time adjusting the learning rate to explore large magnitudes effectively, so I account for this when defining rates.
                gN.modifyRate(x["RateType"], 10**newRateVals[i], x["Target"])



    def calcPredCost(self, y):
        self.updateRates(y)
        totCost = 0
        for i, gN in enumerate(self.GeneletNetworks):
            reqTimes = copy.copy(self.TimePoints[i])
            if hasattr(self, "TimeDelayBounds"):
                reqTimes += y[i-self.DataCount] # Tagged to the end of the optimization parameters so pick from the end in ascending order
            gN.simulate(reqTimes)
            totCost += np.sum((gN.OutputConcentrations[self.OutputValue]/self.OutputScaler - self.FittingData[i])**2)
        #    plt.plot(self.TimePoints[i], gN.OutputConcentrations[self.OutputValue]/self.OutputScaler)
        #    plt.plot(self.TimePoints[i], self.FittingData[i])
        #    plt.title(i)
        #    plt.ylim(-.05, 1.05)
        #    plt.show()
        
        #print("totcost:", totCost)
        return totCost
    

    def optimizeRates(self, maxiter=2000):
        
        bds = [(x["Low Bound"], x["High Bound"]) for x in self.RatesOfInterest]
        if hasattr(self, "TimeDelayBounds"):
            bds = bds + [self.TimeDelayBounds]*self.DataCount
        res = minimize(self.calcPredCost, [(x[0]+x[1])/2 for x in bds], method="Nelder-Mead", options={'maxiter': maxiter}, bounds=bds)
        self.OptimizedRatesTable = res["x"]
        
        self.OptimizedRates = {}
        for i, itm in enumerate(self.RatesOfInterest):
            self.OptimizedRates[itm["RateType"]+" -> "+str(itm["Target"])] = res["x"][i]
            
        if hasattr(self, "TimeDelayBounds"):
            self.OptimizedTimeDelays = res["x"][len(self.RatesOfInterest):]