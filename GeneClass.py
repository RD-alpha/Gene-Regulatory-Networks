import json


def saveNetwork(geneList, path):
    data = [None] * len(geneList)
    for i in range(len(geneList)):
        data[i] = geneList[i].__dict__
    with open("GeneNetworks/" + path + ".txt", "w") as json_file:
        json.dump(data, json_file, indent=4)


def loadNetwork(path):
    with open(path) as json_file:
        data = json.load(json_file)
        output = [None]*len(data)
        for i in range(len(data)):
            output[i] = Gene(obj=data[i])
    return output


class Gene:
    def __init__(self, **kwargs):
        self.speciesName = None
        self.operatorLigands = None
        self.thresholds = None
        self.hillsCoeff = None
        self.alpha = None
        self.logic = None
        self.speciesId = None
        self.operatorStates = None
        self.betas = None

        for key, value in kwargs.items():
            if key == "filename":
                with open("GeneNetworks/" + value + ".txt") as json_file:
                    data = json.load(json_file)
                    self.__dict__ = data
            elif key == "obj":
                self.__dict__ = value

    def save(self, path):
        with open("GeneNetworks/" + path + ".txt", 'w') as outfile:
            json.dump(self.__dict__, outfile, indent=4)

    def logicFunc(self, opID, systemState):
        return 1 if systemState[self.operatorLigands[opID]] >= self.thresholds[opID] else 0

    def hill(self, opID, systemState):
        c = (systemState[self.operatorLigands[opID]] / self.thresholds[opID]) ** self.hillsCoeff[opID]
        return c / (1 + c)

    def getOccupancy(self, opID, state, systemState):
        if self.logic[opID]:
            p = self.logicFunc(opID, systemState)
        else:
            p = self.hill(opID, systemState)
        return p if state else 1 - p

    def getDerivative(self, systemState):
        rate = -self.alpha * systemState[(self.speciesId)]

        for h in range(len(self.operatorStates)):
            p = 1
            for k in range(len(self.operatorStates[h])):
                p *= self.getOccupancy(k, self.operatorStates[h][k], systemState)
            rate += self.betas[h] * p
        return rate
