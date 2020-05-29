import numpy as np
from GeneClass import *


def event(t, systemState):
    return systemState


def setup():
    # --------------
    # set up simulation parameters
    # --------------

    # step is the time advanced every timestep
    # endTime is the time at which the simulation is finished
    # save is the folder in which to save output image e.g. "GeneNetworks/NegFFL/NegFFL.png" set to -1 if no saving wanted
    step = 0.1
    endTime = 1000
    saveFile = "GeneNetworks/Oscillator/clock.png"
    plotInfo = [[0, 1, 2]]  # plot info is a list of lists. every list contains the species to be graphed on a subplot

    # speciesNames holds the labels of all species
    # concentrations holds the concentrations of all substances at time 0
    speciesNames = np.array(["X", "Y", "Z"])
    concentrations = np.array([0, 5, 10])

    # genes are assigned to species
    # species = loadNetwork("LIFO/LIFO") will load network from file "GeneNetworks/LIFO/LIFO.txt"
    # saveNetwork(species,"LIFO/LIFO") will save any changes done to the network after loading to GeneNetworks/LIFO/LIFO.txt
    species = loadNetwork("Oscillator/clock")

    return step, endTime, saveFile, speciesNames, species, concentrations, plotInfo
