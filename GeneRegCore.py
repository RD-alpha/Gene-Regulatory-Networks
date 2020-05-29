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


def graphSim(saveLoc, plotInfo):
    fig, axes = plt.subplots(len(plotInfo), 1, sharex=True)

    if len(plotInfo) > 1:
        for i in range(len(plotInfo)):
            for speciesHistory in concentrationHist[plotInfo[i]]:
                axes[i].plot(times, speciesHistory)
            axes[i].legend(speciesNames[plotInfo[i]])
    else:
        for speciesHistory in concentrationHist[plotInfo[0]]:
            axes.plot(times, speciesHistory)
        axes.legend(speciesNames[tuple(plotInfo)])

    # save plot to saveLoc, path specified by user
    if saveLoc != -1:
        plt.savefig(saveLoc)

    plt.show()


step, endTime, path, speciesNames, speciesList, concentrations, plotInfo = setup()
concentrationHist = simulate(step, endTime, concentrations, speciesList)

times = step * np.array(range(0, int(endTime / step)))
graphSim(path, plotInfo)
