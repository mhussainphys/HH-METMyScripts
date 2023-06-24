import ROOT as r
import os
import collections
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import cm
r.gSystem.Load('libGenVector')

MassPoints = [500]
BaseDirectory = '/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_HToGG_fastSimJmeCorrection/'
FilenameList = [BaseDirectory + 'SMS-TChiHH_mChi-' + str(mass) + '_mLSP-1_HToGG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUSummer16v3Fast_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1.root' for mass in MassPoints]
InFileList = [r.TFile.Open(Filename ,'READ') for Filename in FilenameList]
TreeList = [InFile.Get('Events') for InFile in InFileList]
eventCounter = 10000

def indices(lst, item):
    return [i for i, x in enumerate(lst) if x == item]

#################################
###   Function to be tested  ####
#################################


def TruthPhotonList(event):
    TruthPhotonList = []
    MotherIndexList = []
    TruthPhotonHGGList = [] 
    NumTruthPhotons = 0
    
    for i in range(len(event.GenPart_pdgId)):
        MotherIDindex = event.GenPart_genPartIdxMother[i]
        if event.GenPart_pdgId[i] == 22 and MotherIDindex != -1 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] == 25:
            TruthPhoton = r.Math.XYZTVector()
            TruthPhoton.SetCoordinates(event.GenPart_pt[i],event.GenPart_eta[i],event.GenPart_phi[i],event.GenPart_mass[i])
            TruthPhotonList.append(TruthPhoton)
            MotherIndexList.append(event.GenPart_genPartIdxMother[i])
            NumTruthPhotons = NumTruthPhotons + 1
    
    if NumTruthPhotons == 2 and len(Counter(MotherIndexList).keys()) == 1: 
        TruthPhotonHGGList = TruthPhotonList
    elif NumTruthPhotons == 3:
        for items, counts in Counter(MotherIndexList).items():
            if counts == 2:
                MotherIndexRepeatedTwice = items
                TruthPhotonsIndicesToGrab = indices(MotherIndexList, MotherIndexRepeatedTwice)
                for index in TruthPhotonsIndicesToGrab: TruthPhotonHGGList.append(TruthPhotonList[index])
    elif NumTruthPhotons == 4: 
        TruthPhotonHGGList = TruthPhotonList
    
    return TruthPhotonHGGList



for tree in TreeList:
    NumEvent = 0
    for event in tree:
        TruthPhotonHGGList = TruthPhotonList(event)
        NumEvent += 1
        if NumEvent == eventCounter: break
    print 'Filename Processed: ', FilenameList[TreeList.index(tree)]

