import ROOT
import os

directory = '/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_HToGG_fastSimJmeCorrection/'

chain = ROOT.TChain('Events')
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        chain.Add(f)

m_gg    = ROOT.TH1D('m_gg',   'Mass of 2 photon system;   m_{2g} [GeV];   Entries',70,60,340)
photons = [ROOT.TLorentzVector(), ROOT.TLorentzVector()]
H = ROOT.TLorentzVector()

for event in chain:
    photons[0].SetPxPyPzE(0,0,0,1) 
    photons[1].SetPxPyPzE(0,0,0,1) 
    if (event.nPhoton>=2):
        photons[0].SetPtEtaPhiM(event.Photon_pt[0],event.Photon_eta[0],event.Photon_phi[0],0);
        photons[1].SetPtEtaPhiM(event.Photon_pt[1],event.Photon_eta[1],event.Photon_phi[1],0);
        m_gg.Fill((photons[0]+photons[1]).M())

out = ROOT.TFile('higgs_peak.root','RECREATE')
m_gg.Write()
out.Close()
m_gg.Draw()
ROOT.gPad.Update()
raw_input('Please press enter to continue.')