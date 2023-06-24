def DeltaRPhotonNearestJet(photon, event):
    jet = r.TLorentzVector()
    MinDeltaR = []
    for i in range(len(photon)):
        DeltaR = []
        for j in range(len(event.Jet_pt)):
            jet.SetPtEtaPhiM(event.Jet_pt[j],event.Jet_eta[j],event.Jet_phi[j],event.Jet_mass[j])
            DeltaR.append(r.Math.VectorUtil.DeltaR(photon[i],jet))
        MinDeltaR.append(min(DeltaR))
    return MinDeltaR


def ChainMultipleFiles(NumFiles):
    directory = '/net/cms17/cms17r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_HToGG_fastSimJmeCorrection/'
    chain = r.TChain('Events')
    NumFilesCounter = 0 
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            chain.Add(f)
        NumFilesCounter = NumFilesCounter + 1
        if NumFilesCounter == NumFiles: break 
    return chain


'''
If you find a higgs: record it's index in a list of list.
Loop through GenParMotherid, if it matches the index in the higgs list, assign it to the list (daughter particle) of list (higgs) of list (event)
'''
def HiggsPieChart(tree, HiggsDecayDict):
    HiggsTreeList = []
    DaughterParticleTreeList =[]
    PieChart = []
    explode =[]
    colors = []
    LegendLabel = []
    eventCounter = 0
    TotalDecays = 0
    DebugMode = False
    for event in tree:
        HiggsEventList = []
        DaughterParticleEventList = [[],[]]
        for i in range(len(event.GenPart_pdgId)):
            if event.GenPart_pdgId[i] == 25:
                HiggsEventList.append(i)
        for j in range(len(event.GenPart_genPartIdxMother)):
            if event.GenPart_genPartIdxMother[j] in HiggsEventList:
                if HiggsEventList.index(event.GenPart_genPartIdxMother[j]) == 0:
                    DaughterParticleEventList[0].append(event.GenPart_pdgId[j])
                elif HiggsEventList.index(event.GenPart_genPartIdxMother[j]) == 1:
                    DaughterParticleEventList[1].append(event.GenPart_pdgId[j])
        DaughterParticleTreeList.append([DaughterParticleEventList[0],DaughterParticleEventList[1]])
        HiggsTreeList.append(HiggsEventList)
        eventCounter = eventCounter + 1
        if DebugMode:
            if eventCounter == 100: break;
    for event in DaughterParticleTreeList:
        for daughterpair in event:
            for value in HiggsDecayDict.values():
                if collections.Counter(daughterpair) == collections.Counter(value[0]):
                    value[1] = value[1] + 1

    for key in HiggsDecayDict.keys(): 
        TotalDecays = TotalDecays + HiggsDecayDict.get(key)[1]
        print 'Number of decays in channel %s is %d' %(key, HiggsDecayDict.get(key)[1])
    print 'Number of Total Decays: ', TotalDecays
    print 'Total Number of Events: ', eventCounter
    for value in HiggsDecayDict.values(): PieChart.append(value[1]) 
    mylabels = [key for key in HiggsDecayDict.keys()]
    for value in HiggsDecayDict.values(): explode.append(value[2]) 
    #colors=sns.color_palette('muted', 11)
    cs=cm.Set1(np.arange(11)/11.)
    fig, ax = plt.subplots()
    ax.pie(PieChart, colors=cs, startangle=90)
    ax.set_title("Higgs decay into different channels")
    for i in range(len(PieChart)): 
        Perc = (float(PieChart[i])/float(TotalDecays))*100
        LegendLabel.append((mylabels[i] + " - " + "%.2f%%") % Perc) 
    ax.legend(LegendLabel, loc="best")
    plt.tight_layout()   
    ax.axis('equal')   
    plt.show()

def RoottreeInspect(tree):
    TotalNumTruthPhotons = 0
    TotalNumHiggs = 0
    TotalNumTruthB = 0
    TotalNumEventsProcessesed = 0
    TotalNumTwoTruthPhotonOneMother = 0
    TotalOneTwoTruthPhotonOneMother = 0
    TotalNumZeroTruthPhotons = 0
    NormalEvents = 0
    GGDecay = 0 
    BBDecay = 0
    WWDecay = 0 
    tauDecay = 0   
    ZZDecay = 0
    GlDecay = 0 
    CCDecay = 0
    SSDecay = 0
    MuMuDecay = 0
    Notrace = 0
    print 'Inspecting Root file...'
    for event in tree:
        TruthPhotonMotherIndexList = []
        TruthBMotherIndexList = []
        NumHiggsEvent = 0
        NumBEvent = 0
        NumGEvent = 0
        HiggsIndex = []
        OtherHiggsDecayList = []
        for i in range(len(event.GenPart_pdgId)):
            MotherIDindex = event.GenPart_genPartIdxMother[i]
            if event.GenPart_pdgId[i] == 22 and MotherIDindex != -1 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] == 25:
                TotalNumTruthPhotons = TotalNumTruthPhotons + 1
                TruthPhotonMotherIndexList.append(event.GenPart_genPartIdxMother[i])
                NumGEvent = NumGEvent + 1
            elif event.GenPart_pdgId[i] == 25:
                HiggsIndex.append(i)
                TotalNumHiggs = TotalNumHiggs + 1
                NumHiggsEvent = NumHiggsEvent + 1
            elif event.GenPart_pdgId[i] == 5 or event.GenPart_pdgId[i] == -5:
                if MotherIDindex != -1 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] == 25:
                    TotalNumTruthB = TotalNumTruthB + 1
                    NumBEvent = NumBEvent + 1
        '''elif len(TruthPhotonMotherIndexList)==2 and len(Counter(TruthPhotonMotherIndexList).keys()) == 1:
            TotalNumTwoTruthPhotonOneMother = TotalNumTwoTruthPhotonOneMother + 1
        elif len(TruthPhotonMotherIndexList)==1:
            TotalNumOneTruthPhotonOneMother = TotalNumOneTruthPhotonOneMother + 1
        elif len(TruthPhotonMotherIndexList)==0:
            TotalNumZeroTruthPhotons = TotalNumZeroTruthPhotons + 1'''

        if NumHiggsEvent == 2 and NumGEvent == 2 and len(Counter(TruthPhotonMotherIndexList).keys()) == 1: 
            for index in HiggsIndex: 
                if index != Counter(TruthPhotonMotherIndexList).keys():
                    OtherHiggsIndex = index
            for i in range(len(event.GenPart_pdgId)):
                if event.GenPart_genPartIdxMother[i] == OtherHiggsIndex:
                    OtherHiggsDecayList.append(event.GenPart_pdgId[i])
        if OtherHiggsDecayList == [22,22]:
            GGDecay = GGDecay + 1
        elif OtherHiggsDecayList == [5,-5]:
            BBDecay = BBDecay + 1
        elif OtherHiggsDecayList == [24,-24]:
            WWDecay = WWDecay + 1
        elif OtherHiggsDecayList == [15,-15]:
            tauDecay = tauDecay + 1
        elif OtherHiggsDecayList == [21,21]:
            GlDecay = GlDecay + 1   
        elif OtherHiggsDecayList == [23,23]:
            ZZDecay = ZZDecay + 1   
        elif OtherHiggsDecayList == [4,-4]:
            CCDecay = CCDecay + 1 
        elif OtherHiggsDecayList == [3,-3]:
            SSDecay = SSDecay + 1 
        elif OtherHiggsDecayList == [13,-13]:
            MuMuDecay = MuMuDecay + 1 
        elif OtherHiggsDecayList == []:
            Notrace = Notrace + 1   
        else:
            print OtherHiggsDecayList

        TotalNumEventsProcessesed = TotalNumEventsProcessesed + 1
        if TotalNumEventsProcessesed == 2000: break
    print 'Total Num of Events Processed: ', TotalNumEventsProcessesed
    print 'Total Num of Gen Higgs: ', TotalNumHiggs
    print 'Total Num of Truth Photons: ', TotalNumTruthPhotons
    print 'Total Number of Truth B: ', TotalNumTruthB
    print 'Total Number of Normal Events: ', NormalEvents
    print GGDecay
    y = np.array([GGDecay, BBDecay, WWDecay, tauDecay, GlDecay, ZZDecay, CCDecay, Notrace])
    mylabels = ["gg", "bb", "WW", "tau tau", "glgl","ZZ", "cc", "ss","Mu Mu", "No trace"]
    plt.pie(y, labels = mylabels, autopct='%1.1f%%')
    plt.show()











#HistogramsList = TaggedHistogram(TreeList)

#chain = ChainMultipleFiles(10)
#HiggsPieChart(chain, HiggsDecayDict)

#for tree in TreeList:
#    HiggsPieChart(tree, HiggsDecayDict)

'''
#photon = r.TLorentzVector()
#jet = r.TLorentzVector()
for tree in TreeList:
    EventCounter = 0
    for event in tree:
        for i in range(len(event.Photon_pt)):
            DeltaR = []
            photon = r.TLorentzVector()
            photon.SetPtEtaPhiM(event.Photon_pt[i],event.Photon_eta[i],event.Photon_phi[i],0)
            for j in range(len(event.Jet_pt)):
                jet = r.TLorentzVector()
                jet.SetPtEtaPhiM(event.Jet_pt[j],event.Jet_eta[j],event.Jet_phi[j],event.Jet_mass[j])
                DeltaR.append(r.Math.VectorUtil.DeltaR(photon,jet))
            HistogramsList[TreeList.index(tree)].Fill(min(DeltaR))
        EventCounter+=1
#        if EventCounter == 0000: break
'''


'''

for tree in TreeList:
    
    for event in tree:
        
        if (event.nPhoton>=2):
            NewPhotonList = []
            NewPhotonPtList = []
            photonList = []
            ph = [r.Math.XYZTVector(), r.Math.XYZTVector()]
            for i in range(len(event.Photon_pt)):            
            
                photons = r.Math.XYZTVector()
                photons.SetCoordinates(event.Photon_pt[i],event.Photon_eta[i],event.Photon_phi[i],0)
                photonList.append(photons)
            
            MinDeltaR = DeltaRPhotonNearestJet(photonList, event)

            # Get index in MinDeltaR > 0.3
            indexList = [idx for idx, val in enumerate(MinDeltaR) if val > 0.3]

            if len(indexList) >= 2:
                
                # Get photons from the photonList with the index in indexList and make a new photon list, make a photon pt list from this new list
                for index in indexList: NewPhotonList.append(photonList[index]) 
                for photon in NewPhotonList: NewPhotonPtList.append(photon.Pt()) 

                # Get two highest pt photon from pt list, and get corresponding photons index, get those photons from photon list
                NewPhotonPtList.sort(reverse=True)
                ph[0] = NewPhotonList[NewPhotonPtList.index(NewPhotonPtList[0])]
                ph[1] = NewPhotonList[NewPhotonPtList.index(NewPhotonPtList[1])]
                
                Mass= (ph[0]+ph[1]).M()
                print ph[0].Pt(), ph[1].Pt(), Mass

                HistogramsList[TreeList.index(tree)].Fill((ph[0]+ph[1]).M())

'''
'''

for tree in TreeList:
    
    for event in tree:
        
        if (event.nPhoton>=2):
            NewPhotonList = []
            NewPhotonPtList = []
            photonList = []
            ph = [r.Math.XYZTVector(), r.Math.XYZTVector()]
            for i in range(len(event.Photon_pt)):            
            
                photons = r.Math.XYZTVector()
                photons.SetCoordinates(event.Photon_pt[i],event.Photon_eta[i],event.Photon_phi[i],0)
                photonList.append(photons)
            
            MinDeltaR = DeltaRPhotonNearestJet(photonList, event)

            # Get index in MinDeltaR > 0.3
            indexList = [idx for idx, val in enumerate(MinDeltaR) if val > 0.3]

            if len(indexList) == 2:
                
                Mass= (photonList[0]+photonList[1]).M()
                print photonList[0].Pt(), photonList[1].Pt(), Mass

                HistogramsList[TreeList.index(tree)].Fill((photonList[0]+photonList[1]).M())
'''
'''
for tree in TreeList:
    
    for event in tree:
        
        if (event.nPhoton>=2):
            NewPhotonList = []
            NewPhotonPtList = []
            photonList = [r.Math.XYZTVector(), r.Math.XYZTVector()]
            ph = [r.Math.XYZTVector(), r.Math.XYZTVector()]      
            photonList[0].SetCoordinates(event.Photon_pt[0],event.Photon_eta[0],event.Photon_phi[0],0)
            photonList[1].SetCoordinates(event.Photon_pt[1],event.Photon_eta[1],event.Photon_phi[1],0)
            
            MinDeltaR = DeltaRPhotonNearestJet(photonList, event)

            # Get index in MinDeltaR > 0.3
            indexList = [idx for idx, val in enumerate(MinDeltaR) if val > 0.3]

            if indexList == [0,1]:
                mass = (photonList[0]+photonList[1]).M()
                print photonList[0].Pt(), photonList[1].Pt(), mass
                HistogramsList[TreeList.index(tree)].Fill((photonList[0]+photonList[1]).M())
'''


