import sys
from PySide2 import QtCore, QtWidgets, QtGui
from GeneClass import *
import json
import numpy as np
from config import *
import matplotlib.pyplot as plt
from GeneRegCore import *
import os

class MyWidget(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Gene Regulatory Networks")
        self.comboBox=QtWidgets.QComboBox(self)
        
        for x in speciesList:
            self.comboBox.addItem("Gene "+x.speciesName)
            
        self.comboBox.addItem("Add new gene...")
        
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.comboBox)
        
        self.boundbox_gene=QtWidgets.QHBoxLayout()
        self.layout.addLayout(self.boundbox_gene)
        
        self.gene_properties=QtWidgets.QFormLayout()
        self.operatorspecifics=QtWidgets.QFormLayout()
        self.boundbox_gene.addLayout(self.gene_properties)
        self.boundbox_gene.addLayout(self.operatorspecifics)
        self.setLayout(self.layout)
        
        self.nameLineEdit=QtWidgets.QLineEdit()
        self.alphaLineEdit=QtWidgets.QLineEdit()
        self.opCountSpin=QtWidgets.QSpinBox()
        self.spIDSpin=QtWidgets.QSpinBox()
        self.commitButtonL=QtWidgets.QPushButton("Commit")
        self.gene_properties.addRow(self.tr("&Name:"), self.nameLineEdit)
        self.gene_properties.addRow(self.tr("&Species ID:"), self.spIDSpin)
        self.gene_properties.addRow(self.tr("&Alpha value:"), self.alphaLineEdit)
        self.gene_properties.addRow(self.tr("&Operator count:"), self.opCountSpin)
        self.gene_properties.addWidget(self.commitButtonL)
        self.commitButtonL.clicked.connect(self.commitEditGene)
        
        self.ligandCombo=QtWidgets.QComboBox()
        #self.betaLineEdit=QtWidgets.QLineEdit()
        self.thresholdLineEdit=QtWidgets.QLineEdit()
        self.logicCheck=QtWidgets.QCheckBox()
        self.hillsLineEdit=QtWidgets.QLineEdit()
        self.buttonEdit=QtWidgets.QPushButton("Edit operator states")
        self.opIdSpin=QtWidgets.QSpinBox()
        self.commitButtonR=QtWidgets.QPushButton("Commit")
        self.operatorspecifics.addRow(self.tr("&Operator:"), self.ligandCombo)
        self.operatorspecifics.addRow(self.tr("&Species Id"),self.opIdSpin)
        #self.operatorspecifics.addRow(self.tr("&Beta value:"),self.betaLineEdit)
        self.operatorspecifics.addRow(self.tr("&Thresholds:"),self.thresholdLineEdit)
        self.operatorspecifics.addRow(self.tr("&Logic:"), self.logicCheck)
        self.operatorspecifics.addRow(self.tr("&Hills coefficient:"),self.hillsLineEdit)
        self.operatorspecifics.addWidget(self.buttonEdit)
        self.operatorspecifics.addWidget(self.commitButtonR)
        self.logicCheck.stateChanged.connect(self.hillsUpdate)
        self.commitButtonR.clicked.connect(self.commitEditOperator)
        
        
        self.plotGraphButton=QtWidgets.QPushButton("Draw Graph")
        self.browseButton=QtWidgets.QPushButton("Browse")
        self.graphButtonLayout=QtWidgets.QHBoxLayout()
        self.layout.addLayout(self.graphButtonLayout)
        self.graphButtonLayout.addWidget(self.plotGraphButton)
        self.graphButtonLayout.addWidget(self.browseButton)
        self.plotGraphButton.clicked.connect(self.execPlot)
        self.browseButton.clicked.connect(self.browseDialog)
        
        self.geneSelect()
        self.operatorSelect()
        self.ligandCombo.activated.connect(self.operatorSelect)
        self.buttonEdit.clicked.connect(self.findStates)
        self.comboBox.activated.connect(self.geneSelect)
    def geneSelect(self):
        comboBoxCount=self.comboBox.count()
        if (self.comboBox.currentIndex()+1)==comboBoxCount:
            self.nameGene=""
            askNameGene=dialogNewGene()
            askNameGene.selectType(True)
            askNameGene.exec_()
            print(askNameGene.DialogCode)
            self.comboBox.insertItem(comboBoxCount-1,self.nameGene)
            self.comboBox.setCurrentIndex(comboBoxCount-1)
            speciesList.append(Gene())
        self.currentGene=speciesList[self.comboBox.currentIndex()];
        self.nameLineEdit.setText(self.currentGene.speciesName)
        self.alphaLineEdit.setText(str(self.currentGene.alpha))
        self.spIDSpin.setValue(self.currentGene.speciesId)
        self.opCountSpin.setValue(self.currentGene.operatorCount)
        self.ligandCombo.clear()
        for x in self.currentGene.operatorLigands:
            self.ligandCombo.addItem(speciesNames[x])
        self.ligandCombo.addItem("Add new operator...")
    def operatorSelect(self):
        ligandBoxCount=self.ligandCombo.count()
        if (self.ligandCombo.currentIndex()+1)==ligandBoxCount:
            self.nameLigand=""
            askNameLigand=dialogNewGene()
            askNameLigand.selectType(False)
            askNameLigand.exec_()
            self.ligandCombo.insertItem(ligandBoxCount-1,self.nameLigand)
            self.ligandCombo.setCurrentIndex(ligandBoxCount-1)
        self.selectedIndex=self.ligandCombo.currentIndex()
        #self.betaLineEdit.setText(str(self.currentGene.betas[self.selectedIndex]))
        self.thresholdLineEdit.setText(str(self.currentGene.thresholds[self.selectedIndex]))
        self.logicCheck.setChecked(self.currentGene.logic[self.selectedIndex])
        self.hillsLineEdit.setText(str(self.currentGene.hillsCoeff[self.selectedIndex]))
    def hillsUpdate(self):
        if self.logicCheck.checkState():
            self.hillsLineEdit.setReadOnly(True)
        else:
            self.hillsLineEdit.setReadOnly(False)
    def findStates(self):
        self.boolCombos=np.full((2**(self.ligandCombo.count()-1),self.ligandCombo.count()-1),1)
        self.boolPerm=np.full((self.ligandCombo.count()-1),0)
        for i in range(0, 2**(self.ligandCombo.count()-1)):
            self.bin_num=str(bin(i)[2:])
            for j in range(1,len(self.bin_num)+1):
                self.boolPerm[-j]=self.bin_num[-j]
            self.boolCombos[i]=self.boolPerm
        popupTable=dialogOperatorTable(self.ligandCombo.count()-1,self.boolCombos,self.currentGene)
        popupTable.exec_()
        
            
    def commitEditGene(self):
        self.currentGene.speciesName=self.nameLineEdit.text()
        self.currentGene.alpha=float(self.alphaLineEdit.text())
        self.currentGene.speciesId=self.spIDSpin.value()
        self.currentGene.operatorCount=self.opCountSpin.value()
        speciesList[self.comboBox.currentIndex()]=self.currentGene
        saveNetwork(speciesList,dir_path)
        
    def commitEditOperator(self):
        self.opIndex=self.ligandCombo.currentIndex()
        self.currentGene.thresholds[self.opIndex]=float(self.thresholdLineEdit.text())
        self.currentGene.hillsCoeff[self.opIndex]=float(self.hillsLineEdit.text())
        self.currentGene.logic[self.opIndex]=self.logicCheck.isChecked()
        speciesList[self.comboBox.currentIndex()]=self.currentGene
        saveNetwork(speciesList,dir_path)
    def browseDialog(self):
        self.path = QtWidgets.QFileDialog.getOpenFileName(self, self.tr("Find file"),os.path.dirname(os.path.abspath(__file__))+"\\GeneNetworks")[0]
        step, endTime, saveFilePath, dir_path ,speciesNames, speciesList, concentrations, plotInfo = setup(self.path)
        self.comboBox.clear()
        for x in speciesList:
            self.comboBox.addItem("Gene "+x.speciesName) 
        self.comboBox.addItem("Add new gene...")
        self.geneSelect()
        self.operatorSelect()
    def execPlot(self):
        step, endTime, saveFilePath, dir_path ,speciesNames, speciesList, concentrations, plotInfo = setup(self.path)
        concentrationHist = simulate(step, endTime, concentrations, speciesList)
        times = step * np.array(range(0, int(endTime / step)))
        graphSim(saveFilePath, plotInfo,concentrationHist,times,speciesNames)

                

class dialogNewGene(QtWidgets.QDialog):
    def __init__(self,parent=None):
        super().__init__(parent)
        self.layout=QtWidgets.QVBoxLayout()
        self.edit=QtWidgets.QLineEdit()
        self.button=QtWidgets.QPushButton("Done")

        self.layout.addWidget(self.edit)
        self.layout.addWidget(self.button)

        self.setLayout(self.layout)
    def selectType(self,is_gene):
        self.isGene=is_gene
        self.cont_init()
    def cont_init(self):
        if self.isGene:
            self.setWindowTitle("New gene")
        else:
            self.setWindowTitle("New operator")
        self.button.setDefault(False)
        self.button.clicked.connect(self.giveName)

    def giveName(self):
        if self.isGene:
            widget.nameGene=self.edit.text()

        else:
            widget.nameLigand=self.edit.text()
        self.accept()

class dialogOperatorTable(QtWidgets.QDialog):
    def __init__(self,no_operators,boolCombos,selectedGene,parent=None):
        super().__init__(parent)
        
        self.selectedGene=selectedGene
        self.no_operators=no_operators
        
        self.layout=QtWidgets.QHBoxLayout()
        self.setLayout(self.layout)
        self.setWindowTitle("Edit operator state table")
        self.button=QtWidgets.QPushButton("Done")
        
        self.layoutsArray=np.full((2**(self.no_operators)),QtWidgets.QFormLayout())
        for k in range(2**(self.no_operators)):
            self.layoutsArray[k]=QtWidgets.QFormLayout()
        self.textWidgets=np.full(2**(self.no_operators),QtWidgets.QLabel())
        for k in range(2**(self.no_operators)):
            self.textWidgets[k]=QtWidgets.QLabel()
            self.textWidgets[k].setText(str(boolCombos[k]))
            
        self.checkBoxWidgets=np.full(2**(self.no_operators),QtWidgets.QCheckBox())
        for k in range(2**(self.no_operators)):
            self.checkBoxWidgets[k]=QtWidgets.QCheckBox()
     
        self.betaLineEditWidget=np.full(2**(self.no_operators),QtWidgets.QLineEdit())
        for k in range(2**(self.no_operators)):
            self.betaLineEditWidget[k]=QtWidgets.QLineEdit()
            
        self.no_states=len(selectedGene.operatorStates)
        self.active_k=np.full(self.no_states,0)
        for j in range(self.no_states):
            self.running_str=""
            self.savedState=selectedGene.operatorStates[j]
            for i in self.savedState:
                self.running_str=self.running_str+ str(bin(i))[2:]
                
            self.active_k[j]=int(self.running_str,2) 
        for k in range(2**(self.no_operators)):
            self.textWidgets[k].setText(str(boolCombos[k]))
            self.layout.addLayout(self.layoutsArray[k])
            self.layoutsArray[k].addRow(self.tr("&Operator state"),self.textWidgets[k])
            self.checkBoxWidgets[k].setChecked(False)
            self.layoutsArray[k].addRow(self.tr("&Active"),self.checkBoxWidgets[k])
            self.layoutsArray[k].addRow(self.tr("Beta value"),self.betaLineEditWidget[k])
            
        for z in range(0,len(self.active_k)):
            k=self.active_k[z]
            self.checkBoxWidgets[k].setChecked(True)
            self.betaLineEditWidget[k].setText(str(selectedGene.betas[z]))
        self.layout.addWidget(self.button)
        self.button.clicked.connect(self.commitEditTable)
        
    def commitEditTable(self):
        self.array_length=0
        for j in range(2**(self.no_operators)):
            if self.checkBoxWidgets[j].isChecked():
                self.array_length=self.array_length+1
        self.newStates=np.full((self.array_length,self.no_operators),False)
        self.counter=0
        for j in range(2**(self.no_operators)):
            if self.checkBoxWidgets[j].isChecked():
                self.newStates[self.counter]=[bool(int(char)) for char in str(bin(j))[2:]]
                self.counter=self.counter+1
        self.newBetas=np.full(self.array_length,0.0)
        self.counter=0
        for j in range(2**(self.no_operators)):
            if self.checkBoxWidgets[j].isChecked():
                self.newBetas[self.counter]=float(self.betaLineEditWidget[j].text())
                self.counter=self.counter+1
                
        self.selectedGene.operatorStates=self.newStates.tolist()
        self.selectedGene.betas=self.newBetas.tolist()
        speciesList[widget.comboBox.currentIndex()]=self.selectedGene
        saveNetwork(speciesList,dir_path)
        self.accept()


        
if __name__ == "__main__":
    init_path=os.path.dirname(os.path.abspath(__file__))+"\\GeneNetworks\\SIM\\FIFO.txt"
    step, endTime, saveFilePath, dir_path ,speciesNames, speciesList, concentrations, plotInfo = setup(init_path)
#UI stuff, all very tedious
    app = QtWidgets.QApplication([])
    widget = MyWidget()
    widget.resize(800, 600)
    widget.show()
    sys.exit(app.exec_())
