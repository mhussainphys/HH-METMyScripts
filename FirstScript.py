import ROOT
infile = ROOT.TFile('/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_HToGG_fastSimJmeCorrection/SMS-TChiHH_mChi-200_mLSP-1_HToGG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUSummer16v3Fast_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1.root','READ')
infile.Events.Draw('Photon_mass','nJet==2&&nPhoton==2&&MET_pt>150','')




chain = ROOT.TChain(’tree’)
chain.Add('/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_HToGG_fastSimJmeCorrection/SMS-TChiHH_mChi-200_mLSP-1_HToGG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISummer16NanoAODv7__PUSummer16v3Fast_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1.root')

# Setup Histograms and TLorentzVectors
m_4l    = ROOT.TH1D(’m_4l’,   ’Mass of 4l system;   m_{4l} [GeV];   Entries’,70,60,340)
m_4e    = ROOT.TH1D(’m_4e’,   ’Mass of 4e system;   m_{4e} [GeV];   Entries’,70,60,340)
m_4mu   = ROOT.TH1D(’m_4mu’,  ’Mass of 4mu system;  m_{4mu} [GeV];  Entries’,70,60,340)
m_2e2mu = ROOT.TH1D(’m_2e2mu’,’Mass of 2e2mu system;m_{2e2mu} [GeV];Entries’,70,60,340)
leptons = [ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()]
H = ROOT.TLorentzVector()
for event in chain:
nels, nmus, nleps = 0, 0, 0
charge = 0
leptons[0].SetPxPyPzE(0,0,0,1) leptons[1].SetPxPyPzE(0,0,0,1) leptons[2].SetPxPyPzE(0,0,0,1) leptons[3].SetPxPyPzE(0,0,0,1) H.SetPxPyPzE(0,0,0,1)
j=0
while j < event.els_pt.size(): #loop over electrons
        #only count electrons that pass some quality criteria
        if (event.els_pt[j]>7 and abs(event.els_sceta[j])<2.5 and event.els_reliso[j]<0.35 and
     event.els_sigid[j]):
            #write your code to process electrons here
nels += 1
            nleps += 1
        j += 1
j=0
while j < event.mus_pt.size(): #loop over muons
        #only count muons that pass some quality criteria
        if (event.mus_pt[j]>5 and abs(event.mus_eta[j])<2.4 and event.mus_reliso[j]<0.35 and
     event.mus_sigid[j]):
            #write your code to process muons here
nmus += 1
            nleps += 1
        j += 1
    if (nleps == 4):
        #write your code to reconstruct Higgs mass, check charge == 0,
        #and fill relevant histograms here
# Write final histograms to file
out = ROOT.TFile(’higgs_peak.root’,’RECREATE’)
m_4l.Write()
m_4e.Write()
m_4mu.Write()
m_2e2mu.Write()
out.Close()



TLorentzVector ph1, ph2;
ph1.SetPtEtaPhiM(Photon_pt[0],Photon_eta[0],Photon_phi[0],0);
ph2.SetPtEtaPhiM(Photon_pt[1],Photon_eta[1],Photon_phi[1],0);
(ph1+ph2).M()

