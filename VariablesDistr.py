import ROOT as r
import os
r.gSystem.Load('libGenVector')

#Looping over the entire dataset and chaining them
BaseDirectory = '/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_HToGG_fastSimJmeCorrection/'
chain = r.TChain('Events')
for filename in os.listdir(BaseDirectory):
    f = os.path.join(BaseDirectory, filename)
    if os.path.isfile(f):
        chain.Add(f)
        break #Debug mode: For quick processing, just one file

#Define Histograms
Gamma_PtHist = r.TH1D('Photon Pt',   'Photon Pt distribution for nJet=2 and nPhoton=2 cut;   Photon_Pt[GeV];   Entries',500,0,1500) 
MET_PtHist = r.TH1D('MET Pt',   'Pt_miss distribution for nJet=2 and nPhoton=2 cut;   Pt_miss[GeV];   Entries',500,0,1500) 

#Define Lorentz Vectors
PhotonArray = [r.TLorentzVector(), r.TLorentzVector()]
BJetArray = [r.TLorentzVector(), r.TLorentzVector()]
ReconstructedHiggsPhoton = r.TLorentzVector()
ReconstructedHiggsBJet = r.TLorentzVector()

#Initializing 4 vectors
PhotonArray[0].SetPxPyPzE(0,0,0,1) 
PhotonArray[1].SetPxPyPzE(0,0,0,1) 
BJetArray[0].SetPxPyPzE(0,0,0,1)
BJetArray[1].SetPxPyPzE(0,0,0,1)

#Looping over the events
for event in chain:
    if (event.nPhoton>2):
        PhotonArray[0].SetPtEtaPhiM(event.Photon_pt[0],event.Photon_eta[0],event.Photon_phi[0],0)
        PhotonArray[1].SetPtEtaPhiM(event.Photon_pt[1],event.Photon_eta[1],event.Photon_phi[1],0)
        
        #Filling Histograms
        Gamma_PtHist.Fill(PhotonArray[0].Pt())
        Gamma_PtHist.Fill(PhotonArray[1].Pt())
        MET_PtHist.Fill(event.MET_pt) 

#Writing Final Histograms
out1 = r.TFile('Gamma_PtHist.root','RECREATE')
Gamma_PtHist.Write()
out1.Close()

out2 = r.TFile('MET_PtHist.root','RECREATE')
MET_PtHist.Write()
out2.Close()

#Drawing the Histograms
Gamma_PtHist.Draw()
#MET_PtHist.Draw()
r.gPad.Update()
raw_input('Please press enter to continue.')


