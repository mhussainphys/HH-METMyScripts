import ROOT as r
import os
r.gSystem.Load('libGenVector')

directory = '/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_HToGG_fastSimJmeCorrection/'

chain = r.TChain('Events')
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        chain.Add(f)

DeltaR_pt = r.TH2D('DeltaR_pt','Separation between two photons from Higgs decay vs Reco Higgs Pt;   Reco Higgs Pt [GeV];   Delta R [Radians]',500,0,1500,500,0,10)

photons = [r.TLorentzVector(), r.TLorentzVector()]
H = r.TLorentzVector()

for event in chain:
    photons[0].SetPxPyPzE(0,0,0,1) 
    photons[1].SetPxPyPzE(0,0,0,1) 
    if (event.nJet==2 and event.nPhoton==2):
        photons[0].SetPtEtaPhiM(event.Photon_pt[0],event.Photon_eta[0],event.Photon_phi[0],0)
        photons[1].SetPtEtaPhiM(event.Photon_pt[1],event.Photon_eta[1],event.Photon_phi[1],0)
        HiggsRecoMass = (photons[0]+photons[1]).M()
        if HiggsRecoMass > 100 and HiggsRecoMass < 150:
            H = photons[0]+photons[1]
            R = r.Math.VectorUtil.DeltaR(photons[0],photons[1])
            HiggsPt = H.Pt() 
            DeltaR_pt.Fill(HiggsPt,R)

out = r.TFile('DeltaR_pt_wholedataset.root','RECREATE')
DeltaR_pt.Write()
out.Close()
DeltaR_pt.Draw()



