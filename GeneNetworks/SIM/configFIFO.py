import numpy as np
from GeneClass import *


def event(t, systemState):
    if t == 0:
        systemState[0] = 0
    elif t == 100:
        systemState[0] = 1
    elif t == 700:
        systemState[0] = 0

    return systemState


def setup():
    # --------------
    # set up simulation parameters
    # --------------

    # step is the time advanced every timestep
    # endTime is the time at which the simulation is finished
    # save is the folder in which to save output image e.g. "Genes/NegFFL/NegFFL.png" set to -1 if no saving wanted
    step = 0.1
    endTime = 1400
    saveFile = "GeneNetworks/SIM/FIFO.png"
    plotInfo = [[0], [1, 2], [3, 4, 5]]

    # -------------
    # set up gene network
    # -------------

    # speciesNames holds the labels of all species
    # concentrations holds the concentrations of all substances at time 0
    speciesNames = np.array(["S", "X", "Y", "Z1", "Z2", "Z3"])
    concentrations = np.array([0, 0, 0, 0, 0, 0])

    # signals holds boolean of whether a species is a signal (True) or a gene product (False)
    # signals = np.array([True, False, False, False, False])

    # species is initialized here, it holds an array of all genes
    # species = np.full(len(speciesNames[~signals]), None)

    # genes are assigned to species
    # species[n] = Gene() initializes an empty gene
    # species[n].attribute = value will set the value of an attribute. The following must be set (refer to readme):
    # speciesId, alpha, operatorLigands, thresholds, logic, hillsCoeff, operatorStates, betas
    # species[n] = Gene(filename="geneX") loads from file /GeneNetworks/geneX.txt
    # species[n].save("geneX") saves gene to file GeneNetworks/geneX for later use
    # saveNetwork(species,"LIFO") will save any changes done to the network after loading to GeneNetworks/LIFO.txt
    dir_path="SIM/FIFO"
    species = loadNetwork(dir_path)

    return step, endTime, saveFile, dir_path, speciesNames, species, concentrations, plotInfo
