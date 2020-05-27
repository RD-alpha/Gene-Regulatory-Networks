import numpy as np
import matplotlib.pyplot as plt
from config import *
from GeneClass import *

def simulate(stepSize, simLength, concs, species):
    # initialize simulation variables
    totalSteps = int(simLength / stepSize)
    concHistory = np.empty([totalSteps, len(concs)])

    for m in range(totalSteps):
        concHistory[m] = concs
        dconcs = np.zeros(len(concs))

        # select all species which are gene products and vary their concentration
        for gene in species:
            dconcs[gene.speciesId] = gene.getDerivative(concs)
        concs = concs + dconcs * stepSize
        concs[concs < 0] = 0

        concs = event(m * stepSize, concs)

    return concHistory.T


def graphSim():
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    for gene in speciesList:
        ax2.plot(times, concentrationHist[gene.speciesId])
    ax2.legend(speciesNames[~signals])

    for signalHistory in concentrationHist[signals]:
        ax1.plot(times, signalHistory)
    ax1.legend(speciesNames[signals])

    plt.show()


step, endTime, speciesNames, speciesList, concentrations, signals = setup()
concentrationHist = simulate(step, endTime, concentrations, speciesList)

times = step*np.array(range(0, int(endTime / step)))
graphSim()