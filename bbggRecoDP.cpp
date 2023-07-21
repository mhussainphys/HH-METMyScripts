#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>
#include <regex>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "TLegend.h" // Controls error level reporting
#include "TVector2.h"
#include "Math/Vector4D.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ROOT::Math;

void GetOptions(int argc, char *argv[]);

namespace{
  float lumi = 1;
}

NamedFunc add_cut(NamedFunc & current_cut, NamedFunc additional_cut) {
  current_cut = current_cut && additional_cut;
  return current_cut;
}

/*

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////// Defining do_reconstruction function to choose veto objects, and choose objects to reconstruct Higg //////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  
  Tuple information: 
  
  All the following objects pass the object selection/definition criteria of the 2b analysis and ZY analysis. These can be tweaked further.
  
    Output:
      0: b1_idx : Leading pt ak4 b jet id 
      1: b2_idx : Sub-leading pt ak4 b jet id
      2: g1_idx : Leading pt photon id
      3: g2_idx : Sub-leading pt photon id 
      4: b_idxs : All ak4 b jet IDs
      5: g_idxs : All photon IDs
      6: mu_veto_idxs : All veto muon IDs
      7: e_veto_idxs : All veto electron IDs
      8: ak8bjet_idx : Leading pt AK8 jet id
      9: ak8bjet_idxs : All ak8 b jet IDs

    Input: 
      muons: pt, eta, dz, dxy, pfreliso, mediumId 
      electrons: pt, eta, phi, dz, dxy, pfreliso, cutBased
      photons: pt, eta, isEb, isEE, electronVeto, photon_id
      Jets: pt, btagDeepFlavB, eta, phi
      Fatjets: pt, btagDDBvL, eta, phi, mass, msoftdrop
*/

tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> > do_reconstruction (vector<float> * mu_pt, vector<float> * mu_eta, vector<float> * mu_dz, vector<float> * mu_dxy, vector<float> * mu_pfreliso, vector<bool> * mu_mediumId,
                               vector<float> * e_pt, vector<float> * e_eta, vector<float> * e_phi, vector<float> * e_dz, vector<float> * e_dxy, vector<float> * e_pfreliso, vector<int> * e_cutBased, 
                               vector<float> * photon_pt, vector<float> * photon_eta, vector<bool> * photon_isEb, vector<bool> * photon_isEE, vector<bool> * photon_electronVeto, vector<float> * photon_id,
                               vector<float> * jet_pt, vector<float> * jet_btagDeepFlavB, vector<float> * jet_eta, vector<float> * jet_phi, 
                               vector<float> * FatJet_pt, vector<float> * FatJet_btagDDBvL, vector<float> * FatJet_eta, vector<float> * FatJet_phi, vector<float> * FatJet_mass, vector<float> * FatJet_msoftdrop
                               ) {

  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> > result = {-1, -1, -1, -1, {}, {}, {}, {}, -1, {}};
  
  //////////// Select veto objects ///////// 

  // Select veto muons
  vector<int> & mu_veto_idxs = get<6>(result);
  mu_veto_idxs.reserve(mu_pt->size());
  for (unsigned iPart = 0; iPart<mu_pt->size(); iPart++) {
    if (mu_pt->at(iPart) <= 10) continue;
    if (fabs(mu_eta->at(iPart)) > 2.4) continue;
    if (fabs(mu_dz->at(iPart)) > 0.5) continue;
    if (fabs(mu_dxy->at(iPart)) > 0.2) continue;
    if (mu_pfreliso->at(iPart) >= 0.2) continue;
    if (!mu_mediumId -> at(iPart)) continue;
    mu_veto_idxs.push_back(iPart);
  }
  get<6>(result) = mu_veto_idxs;
  
  // Select veto electrons
  vector<int> & e_veto_idxs = get<7>(result);
  e_veto_idxs.reserve(e_pt->size());
  for (unsigned iPart = 0; iPart<e_pt->size(); iPart++) {
    if (e_pt->at(iPart) <= 10) continue;
    if (fabs(e_eta->at(iPart)) > 2.5) continue;
    if (fabs(e_dz->at(iPart)) > 0.5) continue;
    if (fabs(e_dxy->at(iPart)) > 0.2) continue;
    if (e_pfreliso->at(iPart) >= 0.1) continue;
    bool isBarrel = fabs(e_eta -> at(iPart)) <= 1.479;
    if ((isBarrel && fabs(e_dxy -> at(iPart))>0.05) || (!isBarrel && fabs(e_dxy -> at(iPart))>0.10)) continue; 
    if ((isBarrel && fabs(e_dz -> at(iPart))>0.10) || (!isBarrel && fabs(e_dz -> at(iPart))>0.20)) continue;
    if (e_cutBased -> at(iPart) != 1) continue;
    e_veto_idxs.push_back(iPart);
  }
  get<7>(result) = e_veto_idxs;

  ////////// Stop reconstruction if there are veto objects in the event /////////
  if ((e_veto_idxs.size()!=0)||(mu_veto_idxs.size()!=0)) return result;


  ///////// Select objects for Higgs reconstruction /////////

  // Select photons 
  vector<int> photon_idxs = get<5>(result);
  photon_idxs.reserve(photon_pt->size());
  for (unsigned iPart = 0; iPart<photon_pt->size(); iPart++) {
    // Selection
    if (photon_pt->at(iPart) <= 15) continue;
    //if (photon_lep_drmin->at(iPart) <= 0.3) continue;
    if (fabs(photon_eta->at(iPart)) >= 2.5) continue;
    if (!(photon_isEb->at(iPart) || photon_isEE->at(iPart))) continue;
    if (!photon_electronVeto->at(iPart)) continue;
    if (fabs(photon_eta->at(iPart)) < 1.4442) {
      if (photon_id->at(iPart) <= -0.4) continue;
    } else if (fabs(photon_eta->at(iPart)) > 1.566 && fabs(photon_eta->at(iPart)) < 2.5) {
      if (photon_id->at(iPart) <= -0.58) continue;
    }
    // Save photon
    photon_idxs.push_back(iPart);
  }
  get<5>(result) = photon_idxs;

  // Select b jets
  vector<int> & ak4bjet_idxs = get<4>(result);
  ak4bjet_idxs.reserve(jet_pt->size());
  bool jet_islep = false

  for (unsigned iPart = 0; iPart<jet_pt->size(); iPart++) {
    if (jet_pt->at(iPart) <= 30) continue;
    if (jet_btagDeepFlavB->at(iPart) < 0.3093) continue;
    if (fabs(jet_eta->at(iPart)) > 2.4) continue;

    for(int iel = 0; iel<e_pt->size(); iel++){
      if (dR(e_eta->at(iel), jet_eta->at(iPart), e_phi->at(iel), jet_phi->at(iPart))<0.4 &&
          fabs(jet_pt>at(iPart) - e_pt->at(iel))/e_pt->at(iel) < 1) jet_islep = true;
    }
    if (jet_islep == true){
      jet_islep = false; 
      continue;
    }

    // save ak4 b jet
    ak4bjet_idxs.push_back(iPart);
  }

  vector<int> & ak8bjet_idxs = get<9>(result);
  ak8bjet_idxs.reserve(jet_pt->size());

  for (unsigned iPart = 0; iPart<jet_pt->size(); iPart++) {
    if (FatJet_pt->at(iPart) <= 300) continue;
    //if (FatJet_btagDDBvL->at(iPart) < 0.632) continue;
    if (fabs(FatJet_eta->at(iPart)) > 2.4) continue;
    if (FatJet_msoftdrop->at(iPart) <= 60 && FatJet_msoftdrop->at(iPart) >= 260) continue;

    // save ak8 fat jet
    ak8bjet_idxs.push_back(iPart);
  }

  get<4>(result) = ak4bjet_idxs;
  get<9>(result) = ak8bjet_idxs


  // Choose leading and subleading photon pt

  if (photon_idxs.size() >= 2) {
    float max_photon_pt = 0;
    float submax_photon_pt = 0;
    for (unsigned iIdx=0; iIdx < photon_idxs.size(); ++iIdx) {
      int photon_idx = photon_idxs[iIdx];
      if (photon_pt->at(photon_idx) > max_photon_pt) {
        maxpt_photon_id = photon_idx
        max_photon_pt = photon_pt->at(photon_idx);
      }
      else if (photon_pt->at(photon_idx) > submax_photon_pt) {
        submax_photon_id = photon_idx
        submax_photon_pt = photon_pt->at(photon_idx);
      }
    }
    get<2>(result) = maxpt_photon_id;
    get<3>(result) = submax_photon_id;
  } 
  
  // Choose leading and subleading b jet

  if (ak4bjet_idxs.size() >= 2) {
    float max_ak4bjet_pt = 0;
    float submax_ak4bjet_pt = 0;
    for (unsigned iIdx=0; iIdx < ak4bjet_idxs.size(); ++iIdx) {
      int ak4bjet_idx = ak4bjet_idxs[iIdx];
      if (jet_pt->at(ak4bjet_idx) > max_ak4bjet_pt) {
        maxpt_ak4bjet_id = ak4bjet_idx
        max_ak4bjet_pt = jet_pt->at(ak4bjet_idx);
      }
      else if (jet_pt->at(ak4bjet_idx) > submax_ak4bjet_pt) {
        submax_ak4bjet_id = ak4bjet_idx
        submax_ak4bjet_pt = jet_pt->at(ak4bjet_idx);
      }
    }
    get<0>(result) = maxpt_ak4bjet_id;
    get<1>(result) = submax_ak4bjet_id;
  }

  // Choose leading Ak8 b jet

  if (ak8bjet_idxs.size() >= 1) {
    float max_ak8bjet_pt = 0;
    for (unsigned iIdx=0; iIdx < ak8bjet_idxs.size(); ++iIdx) {
      int ak8bjet_idx = ak8bjet_idxs[iIdx];
      if (FatJet_pt->at(ak8bjet_idx) > max_ak8bjet_pt) {
        maxpt_ak8bjet_id = ak8bjet_idx
        max_ak8bjet_pt = FatJet_pt->at(ak8bjet_idx);
      }
    }
    get<8>(result) = maxpt_ak8bjet_id
  }

  return result;
}


/*

  /////////////////////////////////////////////////////////////////////////////////
  ///////////////// Defining NamedFuncs to reconstruct Higgs mass /////////////////
  /////////////////////////////////////////////////////////////////////////////////  
  
  Tuple information: 
  
  All the following objects pass the object selection/definition criteria of the 2b analysis and ZY analysis. These can be tweaked further.
  
    Output:
      0: b1_idx : Leading pt ak4 b jet id 
      1: b2_idx : Sub-leading pt ak4 b jet id
      2: g1_idx : Leading pt photon id
      3: g2_idx : Sub-leading pt photon id 
      4: b_idxs : All ak4 b jet IDs
      5: g_idxs : All photon IDs
      6: mu_veto_idxs : All veto muon IDs
      7: e_veto_idxs : All veto electron IDs
      8: ak8bjet_idx : Leading pt AK8 jet id
      9: ak8bjet_idxs : All ak8 b jet IDs

    Input: 
      muons: pt, eta, dz, dxy, pfreliso, mediumId 
      electrons: pt, eta, phi, dz, dxy, pfreliso, cutBased
      Jets: pt, btagDeepFlavB, eta, phi
      Fatjets: pt, btagDDBvL, eta, phi, mass, msoftdrop
*/


const NamedFunc is_gg("is_gg",[](const Baby &b) -> NamedFunc::ScalarType{
  float is_gg_recon = 0;
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  if ((get<2>(recon_info) != -1) && (get<3>(recon_info) != -1)) is_gg_recon = 1;
  return is_gg_recon;
});

const NamedFunc is_ak4bb("is_ak4bb",[](const Baby &b) -> NamedFunc::ScalarType{
  float is_ak4bb_recon = 0;
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  if ((get<0>(recon_info) != -1) && (get<1>(recon_info) != -1)) is_ak4bb_recon = 1;
  return is_ak4bb_recon;
});

const NamedFunc is_ak8b("is_ak8b",[](const Baby &b) -> NamedFunc::ScalarType{
  float is_ak8b_recon = 0;
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  if (get<8>(recon_info) != -1) is_ak8b_recon = 1;
  return is_ak8b_recon;
});

const NamedFunc get_gg_m("get_gg_m",[](const Baby &b) -> NamedFunc::ScalarType{
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  int photon_idx1 = get<2>(recon_info);
  int photon_idx2 = get<3>(recon_info);

  LorentzVector<PtEtaPhiM4D<double> > photon_1;
  LorentzVector<PtEtaPhiM4D<double> > photon_2;
  
  photon_1.SetCoordinates(b.Photon_pt()->at(photon_idx1), b.Photon_eta()->at(photon_idx1), b.Photon_phi()->at(photon_idx1), 0.0);
  photon_2.SetCoordinates(b.Photon_pt()->at(photon_idx2), b.Photon_eta()->at(photon_idx2), b.Photon_phi()->at(photon_idx2), 0.0);

  LorentzVector<PtEtaPhiM4D<double> > gg = photon_1 + photon_2;
  return gg.mass();
});

const NamedFunc get_ak4bb_m("get_ak4bb_m",[](const Baby &b) -> NamedFunc::ScalarType{
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  int ak4b_idx1 = get<0>(recon_info);
  int ak4b_idx2 = get<1>(recon_info);

  LorentzVector<PtEtaPhiM4D<double> > ak4b_1;
  LorentzVector<PtEtaPhiM4D<double> > ak4b_2;
  
  ak4b_1.SetCoordinates(b.Jet_pt()->at(ak4b_idx1), b.Jet_eta()->at(ak4b_idx1), b.Jet_phi()->at(ak4b_idx1), b.Jet_mass()->at(ak4b_idx1));
  ak4b_2.SetCoordinates(b.Jet_pt()->at(ak4b_idx2), b.Jet_eta()->at(ak4b_idx2), b.Jet_phi()->at(ak4b_idx2), b.Jet_mass()->at(ak4b_idx2));

  LorentzVector<PtEtaPhiM4D<double> > ak4bb = ak4b_1 + ak4b_2;
  return ak4bb.mass();
});


const NamedFunc get_ak8b_m("get_ak8b_m",[](const Baby &b) -> NamedFunc::ScalarType{
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  int ak8b_idx1 = get<8>(recon_info);
  LorentzVector<PtEtaPhiM4D<double> > ak8b;  
  ak8b.SetCoordinates(b.FatJet_pt()->at(ak8b_idx1), b.FatJet_eta()->at(ak8b_idx1), b.FatJet_phi()->at(ak8b_idx1), b.FatJet_mass()->at(ak8b_idx1));
  return ak8b.mass();

});


/*

  /////////////////////////////////////////////////////////////////////////////////
  ///////////////// Defining NamedFuncs to get photon, b jets pts /////////////////
  /////////////////////////////////////////////////////////////////////////////////  
  
  Tuple information: 
  
  All the following objects pass the object selection/definition criteria of the 2b analysis and ZY analysis. These can be tweaked further.
  
    Output:
      0: b1_idx : Leading pt ak4 b jet id 
      1: b2_idx : Sub-leading pt ak4 b jet id
      2: g1_idx : Leading pt photon id
      3: g2_idx : Sub-leading pt photon id 
      4: b_idxs : All ak4 b jet IDs
      5: g_idxs : All photon IDs
      6: mu_veto_idxs : All veto muon IDs
      7: e_veto_idxs : All veto electron IDs
      8: ak8bjet_idx : Leading pt AK8 jet id
      9: ak8bjet_idxs : All ak8 b jet IDs

    Input: 
      muons: pt, eta, dz, dxy, pfreliso, mediumId 
      electrons: pt, eta, phi, dz, dxy, pfreliso, cutBased
      Jets: pt, btagDeepFlavB, eta, phi
      Fatjets: pt, btagDDBvL, eta, phi, mass, msoftdrop

*/

const NamedFunc get_lead_photon_pt("get_lead_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  int photon_idx = get<2>(recon_info);
  LorentzVector<PtEtaPhiM4D<double> > photon (b.photon_pt()->at(photon_idx), b.photon_eta()->at(photon_idx), b.photon_phi()->at(photon_idx), 0.0);
  return photon.pt();
});

const NamedFunc get_sublead_photon_pt("get_sublead_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  int photon_idx = get<3>(recon_info);
  LorentzVector<PtEtaPhiM4D<double> > photon (b.photon_pt()->at(photon_idx), b.photon_eta()->at(photon_idx), b.photon_phi()->at(photon_idx), 0.0);
  return photon.pt();
});

const NamedFunc get_lead_ak4b_pt("get_lead_ak4b_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  int ak4b_idx = get<0>(recon_info);
  LorentzVector<PtEtaPhiM4D<double> > ak4b (b.Jet_pt()->at(ak4b_idx), b.Jet_eta()->at(ak4b_idx), b.Jet_phi()->at(ak4b_idx), b.Jet_mass()->at(ak4b_idx));
  return ak4b.pt();
});

const NamedFunc get_sublead_ak4b_pt("get_sublead_ak4b_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  int ak4b_idx = get<1>(recon_info);
  LorentzVector<PtEtaPhiM4D<double> > ak4b (b.Jet_pt()->at(ak4b_idx), b.Jet_eta()->at(ak4b_idx), b.Jet_phi()->at(ak4b_idx), b.Jet_mass()->at(ak4b_idx));
  return ak4b.pt();
});

const NamedFunc get_lead_ak8b_pt("get_lead_ak8b_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  tuple<int, int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, int, vector<int> >  recon_info = do_reconstruction(b.Muon_pt(), b.Muon_eta(), b.Muon_dz(), b.Muon_dxy(), b.Muon_miniPFRelIso_all(), b.Muon_mediumId(),
                     b.Electron_pt(), b.Electron_eta(), b.Electron_phi(), b.Electron_dz(), b.Electron_dxy(), b.Electron_miniPFRelIso_all(), b.Electron_cutBased(),
                     b.Photon_pt(), b.Photon_eta(), b.Photon_isScEtaEB(), b.Photon_isScEtaEE(), b.Photon_electronVeto(), b.Photon_mvaID_WP90(),
                     b.Jet_pt(), b.Jet_btagDeepFlavB(), b.Jet_eta(), b.Jet_phi(), 
                     b.FatJet_pt(), b.FatJet_btagDDBvL(), b.FatJet_eta(), b.FatJet_phi(), b.FatJet_mass(), b.FatJet_msoftdrop()
                     );
  int ak8b_idx = get<8>(recon_info);
  LorentzVector<PtEtaPhiM4D<double> > ak8b (b.FatJet_pt()->at(ak8b_idx), b.FatJet_eta()->at(ak8b_idx), b.FatJet_phi()->at(ak8b_idx), b.FatJet_mass()->at(ak8b_idx));
  return ak8b.pt();
});


/*

  ///////////////////////////////////////////////// 
  ///////////////// Main function ///////////////// 
  ///////////////////////////////////////////////// 
  
  Tuple information: 
  
  All the following objects pass the object selection/definition criteria of the 2b analysis and ZY analysis. These can be tweaked further.
  
    Output:
      0: b1_idx : Leading pt ak4 b jet id 
      1: b2_idx : Sub-leading pt ak4 b jet id
      2: g1_idx : Leading pt photon id
      3: g2_idx : Sub-leading pt photon id 
      4: b_idxs : All ak4 b jet IDs
      5: g_idxs : All photon IDs
      6: mu_veto_idxs : All veto muon IDs
      7: e_veto_idxs : All veto electron IDs
      8: ak8bjet_idx : Leading pt AK8 jet id
      9: ak8bjet_idxs : All ak8 b jet IDs

    Input: 
      muons: pt, eta, dz, dxy, pfreliso, mediumId 
      electrons: pt, eta, phi, dz, dxy, pfreliso, cutBased
      Jets: pt, btagDeepFlavB, eta, phi
      Fatjets: pt, btagDDBvL, eta, phi, mass, msoftdrop

*/


int main(int argc, char *argv[]){
  
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);
  time_t begtime, endtime;
  time(&begtime);


  /*
      
      Defining all the data sample paths
      
      Signal:

      0: SMS-TChiHH_mChi-150_mLSP-1_HToGG
      1: SMS-TChiHH_mChi-500_mLSP-1_HToGG
      2: SMS-TChiHH_mChi-1000_mLSP-1_HToGG
      
      Background:

      0: DiPhotonJetsBox_MGG-80toInf
      1: DiPhotonJetsBox1BJet_MGG-80toInf
      2: DiPhotonJetsBox2BJets_MGG-80toInf
      3: GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf
      4: GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf
      5: TTGG_0Jets
      6: TTGJets
      7: TTTo2L2Nu
      8: QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf
      9: QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf
  
  */
  
  Palette colors("txt/colors.txt","default");
  Process::Type sig  = Process::Type::signal;
  Process::Type bkg  = Process::Type::background;

  std::vector<std::string> BkgSampleNames = {"DiPhotonJetsBox_MGG-80toInf", 
                                            "DiPhotonJetsBox1BJet_MGG-80toInf", 
                                            "DiPhotonJetsBox2BJets_MGG-80toInf",
                                            "GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf",
                                            "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf",
                                            "TTGG_0Jets",
                                            "TTGJets",
                                            "TTTo2L2Nu",
                                            "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf",
                                            "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf"};

  std::vector<std::string> SigSampleNames = {"SMS-TChiHH_mChi-150_mLSP-1_HToGG", 
                                            "SMS-TChiHH_mChi-500_mLSP-1_HToGG", 
                                            "SMS-TChiHH_mChi-1000_mLSP-1_HToGG"};

  string bfolderbkg("/net/cms37/data1/mhussain/HH-MET/DataSample/Background/");
  string bfoldersig("/net/cms37/data1/mhussain/HH-MET/DataSample/Signal/");
  
  string mchi_150(bfoldersig + SigSampleNames[0] + "*");
  string mchi_500(bfoldersig + SigSampleNames[1] + "*");
  string mchi_1000(bfoldersig + SigSampleNames[2] + "*");
  
  string Diphoton(bfolderbkg + BkgSampleNames[0] + "/*");
  string Diphoton1B(bfolderbkg + BkgSampleNames[1] + "/*");
  string Diphoton2B(bfolderbkg + BkgSampleNames[2] + "/*");
  string GJet1(bfolderbkg + BkgSampleNames[3] + "/*");
  string GJet2(bfolderbkg + BkgSampleNames[4] + "/*");
  string TTGG(bfolderbkg + BkgSampleNames[5] + "/*");
  string TTGJets(bfolderbkg + BkgSampleNames[6] + "/*");
  string TTTo2L2Nu(bfolderbkg + BkgSampleNames[7] + "/*");
  string QCD_Pt2(bfolderbkg + BkgSampleNames[8] + "/*");
  string QCD_Pt1(bfolderbkg + BkgSampleNames[9] + "/*");


  // Weighting the samples by cross section and luminosity 

  float Luminosity = 36.31;

  std::vector<std::float> xsecSig = {58.94448,
                                    0.486923,
                                    0.01354};
  
  std::vector<std::float> xsecBkg = {2995.9381,
                                    29.719735,
                                    17.697494,
                                    7959.152,
                                    31313.744,
                                    0.6285261,
                                    31452.1689761,
                                    24949.2295261,
                                    32232,
                                    113100.0
                                    };
  
  const NamedFunc wgt("wgt", [](const Baby &b) -> NamedFunc::ScalarType{

    int signal_scale_factor = 1000;

    if (b.SampleTypeString() == SigSampleNames[0]) return  xsecSig[0] * signal_scale_factor /b.GetEntries();
    if (b.SampleTypeString() == SigSampleNames[1]) return  xsecSig[1] * signal_scale_factor/b.GetEntries();
    if (b.SampleTypeString() == SigSampleNames[2]) return  xsecSig[2] * signal_scale_factor/b.GetEntries();

    if (b.SampleTypeString() == BkgSampleNames[0]) return  xsecBkg[0]/b.GetEntries();
    if (b.SampleTypeString() == BkgSampleNames[1]) return  xsecBkg[1]/b.GetEntries();
    if (b.SampleTypeString() == BkgSampleNames[2]) return  xsecBkg[2]/b.GetEntries();
    if (b.SampleTypeString() == BkgSampleNames[3]) return  xsecBkg[3]/b.GetEntries();
    if (b.SampleTypeString() == BkgSampleNames[4]) return  xsecBkg[4]/b.GetEntries();
    if (b.SampleTypeString() == BkgSampleNames[5]) return  xsecBkg[5]/b.GetEntries();
    if (b.SampleTypeString() == BkgSampleNames[6]) return  xsecBkg[6]/b.GetEntries();
    if (b.SampleTypeString() == BkgSampleNames[7]) return  xsecBkg[7]/b.GetEntries();   
    if (b.SampleTypeString() == BkgSampleNames[8]) return  xsecBkg[8]/b.GetEntries();
    if (b.SampleTypeString() == BkgSampleNames[9]) return  xsecBkg[9]/b.GetEntries(); 
    
    else return -1.;

  });


  // Setting up Plotmaker

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::overflow)
          .YTitleOffset(1.75)
          .AutoYAxis(true)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .FileExtensions({"pdf"});
  
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {log_stack};

  auto proc_mchi_150 = Process::MakeShared<Baby_nano>(SigSampleNames[0], sig, TColor::GetColor("#ff0000"), {mchi_150}, "1");
  auto proc_mchi_500 = Process::MakeShared<Baby_nano>(SigSampleNames[1], sig, TColor::GetColor("#00ff00"), {mchi_500}, "1");
  auto proc_mchi_1000 = Process::MakeShared<Baby_nano>(SigSampleNames[2], sig, TColor::GetColor("#0000ff"), {mchi_1000}, "1");

  auto proc_Diphoton = Process::MakeShared<Baby_nano>(BkgSampleNames[0], bkg, TColor::GetColor("#006600"), {Diphoton}, "1");
  auto proc_DiPhoton1B = Process::MakeShared<Baby_nano>(BkgSampleNames[1], bkg, TColor::GetColor("#00cc00"), {Diphoton1B}, "1");
  auto proc_DiPhoton2B = Process::MakeShared<Baby_nano>(BkgSampleNames[2], bkg, TColor::GetColor("#ffcc00"), {Diphoton2B}, "1");
  auto proc_GJet1 = Process::MakeShared<Baby_nano>(BkgSampleNames[3], bkg, TColor::GetColor("#99dcff"), {GJet1}, "1");
  auto proc_GJet2 = Process::MakeShared<Baby_nano>(BkgSampleNames[4], bkg, TColor::GetColor("#9999ff"), {GJet2}, "1");
  auto proc_TTGG = Process::MakeShared<Baby_nano>(BkgSampleNames[5], bkg, TColor::GetColor("#0194da"), {TTGG}, "1");
  auto proc_TTGJets = Process::MakeShared<Baby_nano>(BkgSampleNames[6], bkg, TColor::GetColor("#09ba01"), {TTGJets}, "1");
  auto proc_TTTo2L2Nu = Process::MakeShared<Baby_nano>(BkgSampleNames[7], bkg, TColor::GetColor("#deba87"), {TTTo2L2Nu}, "1");
  auto proc_QCDPt1 = Process::MakeShared<Baby_nano>(BkgSampleNames[8], bkg, TColor::GetColor("#fe6200"), {QCD_Pt1}, "1");
  auto proc_QCDPt2 = Process::MakeShared<Baby_nano>(BkgSampleNames[9], bkg, TColor::GetColor("#ff6600"), {QCD_Pt2}, "1");
  
  
  proc_mchi_150->SetLineWidth(3);
  proc_mchi_500->SetLineWidth(3);
  proc_mchi_1000->SetLineWidth(3);

  vector<shared_ptr<Process>> procs_sig = {proc_mchi_150, proc_mchi_500, proc_mchi_1000};
  vector<shared_ptr<Process>> procs_bkg = {proc_QCDPt1, proc_QCDPt2, proc_TTGJets, proc_GJet2, proc_TTTo2L2Nu, proc_GJet1, proc_DiPhoton0B, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGG};
  vector<shared_ptr<Process>> procs = {proc_mchi_150, proc_mchi_500, proc_mchi_1000, proc_QCDPt1, proc_QCDPt2, proc_TTGJets, proc_GJet2, proc_TTTo2L2Nu, proc_GJet1, proc_DiPhoton0B, proc_DiPhoton1B, proc_DiPhoton2B, proc_TTGG};

  PlotMaker pm;

  NamedFunc baseline = "1";

  pm.Push<Hist1D>(Axis(160, 100, 180, "llphoton_m[0]", "llg mass [GeV]", {}, {}), baseline, procs_sig, ops).Weight(weight).Tag("FixName:pico_llg_m");
  pm.Push<Hist1D>(Axis(80, 100, 180, "llphoton_m[0]", "llg mass [GeV]", {}, {}), baseline, procs_bkg, ops).Weight(weight).Tag("FixName:bkg_pico_llg_m");
  pm.Push<Hist1D>(Axis(5, 0, 5, "nJet", "nJet", {}), cutflow.at(1), procs, ops).Weight(wgt).Tag("nJet");

  pm.min_print_ = true;
  pm.MakePlots(Luminosity);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;

}
























