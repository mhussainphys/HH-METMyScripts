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
eventCounter = 1


for tree in TreeList:
    NumEvent = 0

    for event in tree:
		
		for i in range(len(event.Photon_pt)):
			RecoPhoton1 = r.Math.PtEtaPhiMVector()
			RecoPhoton1.SetCoordinates(event.Photon_pt[i],event.Photon_eta[i],event.Photon_phi[i],0)
			RecoPhoton2 = r.Math.XYZTVector()
			RecoPhoton2.SetCoordinates(event.Photon_pt[i],event.Photon_eta[i],event.Photon_phi[i],0)
			RecoPhoton3 = r.Math.XYZTVector()
			RecoPhoton3.SetCoordinates(227.874938965, 0.857788085938, 0.0795593261719, 0)
			print np.sqrt((RecoPhoton2.Eta() - 0.857788085938)**2 + (RecoPhoton2.Phi() - 0.0795593261719)**2)
			print r.Math.VectorUtil.DeltaR(RecoPhoton1, RecoPhoton3)
			print event.Photon_pt[i],event.Photon_eta[i],event.Photon_phi[i],0
			print RecoPhoton1.Pt(), RecoPhoton1.Eta(), RecoPhoton1.Phi(), RecoPhoton1.M()
			print RecoPhoton2.Pt(), RecoPhoton2.Eta(), RecoPhoton2.Phi(), RecoPhoton2.M()

		NumEvent += 1
		if NumEvent == eventCounter: break
    print 'Filename Processed: ', FilenameList[TreeList.index(tree)]
