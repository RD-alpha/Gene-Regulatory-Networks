import numpy as np
import matplotlib.pyplot as plt
from GeneClass import *


def setup():
    species = np.ndarray(concCount, dtype=np.object)
    # args: speciesId, operatorCount, alpha, betas, possible operator states, thresholds, step behaviour (True)
    # or michaelis behaviour (False), species which binds to operators, hill coefficient
    species[1] = Gene(1, speciesNames[1], 1, [0], 0.05, [0.5], [[True]], [1], [True], [2])
    species[2] = Gene(2, speciesNames[2], 1, [1], 0.03, [0.5], [[True]], [10], [False], [2])
    species[3] = Gene(3, speciesNames[3], 2, [1, 2], 0.05, 0.3 * np.array([1]), [[True, False]], [10, 5],
                           [False, False], [2, 8])

    return species


def event(t, systemState):
    if t == 0:
        systemState[0] = 0
    elif t == 5:
        systemState[0] = 1
    elif t == 25:
        systemState[0] = 0
    elif t == 200:
        systemState[0] = 1
    elif t == 500:
        systemState[0] = 0
    return systemState


# def simulate(stepSize, simLength, concs, signalConcs):
def simulate(stepSize, simLength, concs, species):
    # initialize simulation variables
    totalSteps = int(simLength / stepSize)
    concHistory = np.empty([totalSteps, len(concs)])

    for m in range(totalSteps):
        concHistory[m] = concs
        dconcs = np.zeros(concCount)

        # select all species which are gene products and vary their concentration
        for gene in species[~signals]:
            dconcs[gene.speciesId] = gene.getDerivative(concs)
        concs = concs + dconcs * stepSize
        concs[concs < 0] = 0

        concs = event(m * stepSize, concs)

    return concHistory


# geneNames holds the labels of all genes
# concentrations holds the concentrations of all substances
speciesNames = np.array(["S", "X", "Y", "Z"])
concCount = len(speciesNames)
concentrations = np.array([0, 0, 0, 0])

# signals holds boolean of whether a species is a signal (True) or a gene product (False)
signals = np.full(concCount, False)
signals[0] = True

# step is the time advanced every timestep
# endTime is the time at which the simulation is finished
step = 0.1
endTime = 600

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

speciesList = setup()
concentrationHist = simulate(step, endTime, concentrations, speciesList).T
steps = np.array(range(0, int(endTime / step)))
times = steps * step

for gene in speciesList[~signals]:
    ax2.plot(times, concentrationHist[gene.speciesId])
ax2.legend(speciesNames[~signals])

for signalHistory in concentrationHist[signals]:
    ax1.plot(times, signalHistory)
ax1.legend(speciesNames[signals])

plt.show()
