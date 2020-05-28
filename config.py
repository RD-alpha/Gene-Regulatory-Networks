import numpy as np
from GeneClass import *


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


def setup():
    # --------------
    # set up simulation parameters
    # --------------

    # step is the time advanced every timestep
    # endTime is the time at which the simulation is finished
    step = 0.1
    endTime = 600

    # -------------
    # set up gene network
    # -------------

    # speciesNames holds the labels of all species
    # concentrations holds the concentrations of all substances at time 0
    speciesNames = np.array(["S", "X", "Y", "Z"])
    concentrations = np.array([0, 0, 0, 0])

    # signals holds boolean of whether a species is a signal (True) or a gene product (False)
    signals = np.array([True, False, False, False])

    #species is initialized here, it holds an array of all genes
    species = np.ndarray(3, dtype=np.object)

    # genes are assigned to species
    # species[n] = Gene() initializes an empty gene
    # species[n] = Gene(filename="geneX.txt") loads from file /Genes/geneX.txt
    # species[n].save("geneX.txt") saves gene to file /Genes/geneX.txt for later use
    species[0] = Gene(filename="geneX.txt")
    species[1] = Gene(filename="geneY.txt")
    species[2] = Gene(filename="geneZ.txt")

    return step, endTime, speciesNames, species, concentrations, signals
