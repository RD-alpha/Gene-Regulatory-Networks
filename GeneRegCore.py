import numpy as np
import matplotlib.pyplot as plt

# geneNames holds the labels of all genes
# concentrations holds the concentrations of all gene products
# signalNames holds the labels of all signals
# signals records the signal concentration
geneNames = np.array(["S", "X", "Y", "Z"])
geneCount = len(geneNames)
concentrations = np.zeros(geneCount)
concentrations[0] = 1
# step is the time advanced every timestep
# endTime is the time at which the simulation is finished
step = 0.1
endTime = 400


class Gene:
    def __init__(self, id, name, opAmount, alpha, betas, states, thresholds, logic, ligands, coeff):
        self.geneId = id
        self.geneName = name
        self.operatorCount = opAmount

        self.alpha = alpha
        self.betas = betas
        self.operatorStates = states

        self.thresholds = thresholds
        self.logic = logic
        self.operatorLigands = ligands
        self.hillsCoeff = coeff

    def logicFunc(self, opID, systemState):
        if systemState[self.operatorLigands[opID]] >= self.thresholds[opID]:
            return 1
        else:
            return 0

    def hill(self, opID, systemState):
        c = (systemState[self.operatorLigands[opID]] / self.thresholds[opID]) ** self.hillsCoeff[opID]
        return c / (1 + c)

    def getOccupancy(self, opID, state, systemState):
        if self.logic[opID]:
            p = self.logicFunc(opID, systemState)
        else:
            p = self.hill(opID, systemState)
        if state:
            return p
        else:
            return 1 - p

    def getDerivative(self, systemState):
        rate = -self.alpha * systemState[self.geneId]

        for h in range(len(self.operatorStates)):
            p = 1
            for k in range(len(self.operatorStates[h])):
                p *= self.getOccupancy(k, self.operatorStates[h][k], systemState)
            rate += self.betas[h] * p
        return rate


def setup():
    genes = np.ndarray(geneCount, dtype=np.object)
    genes[0] = Gene(0, "S", 0, 0, [0], [], [], [], [], [])
    genes[1] = Gene(1, "X", 1, 0.05, [0.5], [[True]], [0], [True], [0], [0])
    genes[2] = Gene(2, "Y", 1, 0.04, [0.4], [[True]], [8], [False], [1], [10])
    genes[3] = Gene(3, "Z", 2, 0.02, [0.3], [[True, False]], [8, 5], [False, False], [1, 2], [10, 10])

    return genes


#    def getEvent(t, signalChange):
#        if t == 5:
#            signalChange[0] = 1
#        if t == 20:
#            signalChange[0] = 0
#        if t == 200:
#            signalChange[0] = 1
#        return signalChange


# def simulate(stepSize, simLength, concs, signalConcs):
def simulate(stepSize, simLength, concs, genes):
    # initialize simulation variables
    totalSteps = int(simLength / stepSize)
    concHistory = np.empty([totalSteps, len(concs)])
    #    signalLog = np.empty([totalSteps, len(signalConcs)])

    for m in range(totalSteps):
        concHistory[m] = concs
        #       signalLog[m] = signalConcs

        # dconcs = getDerivatives(concs, signalConcs)
        dconcs = np.zeros(geneCount)
        for n in range(geneCount):
            dconcs[n] = genes[n].getDerivative(concs)

        concs = concs + dconcs * stepSize
        concs[concs < 0] = 0

    #        signalConcs = getEvent(m * step, signalConcs)

    #    return concHistory, signalLog
    return concHistory


# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

# concentrationHist, signalHistory = simulate(step, endTime, concentrations, signals)
geneList = setup()
concentrationHist = simulate(step, endTime, concentrations, geneList)
times = np.array(range(0, int(endTime / step))) * step

# for i in range(0, geneCount):
#    ax2.plot(times, concentrationHist.T[i])
# ax2.legend(geneNames)

# for i in range(0, signalCount):
#    ax1.plot(times, signalHistory.T[i])
# ax1.legend(signalNames)

for i in range(0, geneCount):
    plt.plot(times, concentrationHist.T[i])
plt.legend(geneNames)
plt.show()
