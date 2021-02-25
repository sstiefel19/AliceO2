// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "Analysis/Jet.h"
#include "Analysis/JetFinder.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include <TH1F.h>
#include <TMath.h>
#include <TString.h>
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;
//using namespace o2::framework::expressions;


struct ProcessGammaJet {

  //general event histograms
  OutputObj<TH1I> hist_events_N{TH1I("events_N", "events_N", 10, 0, 10)};
  OutputObj<TH1F> hist_events_z{TH1F("events_z", "events_z", 100, -10., 10.)};
  OutputObj<TH2F> hist_events_xy{TH2F("events_xy", "events_xy", 100, -1., 1., 100, -1., 1.)};
  OutputObj<TH1I> hist_tracks_NGood{TH1I("tracks_NGood", "tracks_NGood", 1000, 0, 10000)};
  OutputObj<TH2F> hist_tracks_etaphi{TH2F("tracks_etaphi", "tracks_etaphi", 102, -2.01, 2.01, 100, 0., 2. * M_PI)};

  //photon histograms
  OutputObj<TH1F> hist_calo_pt{TH1F("calo_pt", "calo_pt", 500, 0., 100.)};
  OutputObj<TH1F> hist_calo_E{TH1F("calo_E", "calo_E", 500, 0., 100.)};
  OutputObj<TH2F> hist_calo_etaphi{TH2F("calo_etaphi", "calo_etaphi", 102, -2.01, 2.01, 100, 0., 2. * M_PI)};

  //jet histograms
  OutputObj<TH1F> hist_jet_pt{TH1F("jet_pt", "jet_pt", 500, 0., 100.)};
  OutputObj<TH2F> hist_jet_etaphi{TH2F("jet_etaphi", "jet_etaphi", 102, -2.01, 2.01, 100, 0., 2. * M_PI)};

  //gamma-jet histograms
  OutputObj<TH2F> hist_gammajet_detadphi{TH2F("gammajet_detadphi", "gammajet_detadphi", 102, -2.01, 2.01, 100, 0., 2. * M_PI)};
  OutputObj<TH1F> hist_gammajet_dE{TH1F("gammajet_dE", "gammajet_dE", 100,0.,100.)};
  OutputObj<TH1F> hist_gammajet_xj{TH1F("gammajet_xj", "gammajet_xj", 100,0.,10.)};
  OutputObj<TH1F> hist_gammajet_ratioE{TH1F("gammajet_ratioE", "gammajet_ratioE", 100,0.,10.)};

  //meson candidate histograms
  OutputObj<TH2F> hist_meson_masspt{TH2F("meson_masspt", "meson_masspt", 100, 0., 1., 250, 0., 100.)};
  

  //the following doesn't work, why?
  //hist_NEvents->GetXaxis()->SetBinLabel(0,"accepted");

 //use a filter for tracks
  Filter trackCuts = aod::track::pt >= 0.15f && aod::track::eta > -0.9f && aod::track::eta < 0.9f;

  //initialize the jet finder
  std::vector<fastjet::PseudoJet> jets;
  std::vector<fastjet::PseudoJet> inputParticles;
  JetFinder jetFinder;
  void init(InitContext const&)
  {
    jetFinder.jetR = 0.4;
  }

  //function for calculating deltaphi in the correct interval
  template <typename T> T relativePhi(T phi1, T phi2) {
    if (phi1 < -TMath::Pi()) {
      phi1 += (2 * TMath::Pi());
    } else if (phi1 > TMath::Pi()) {
      phi1 -= (2 * TMath::Pi());
    }
    if (phi2 < -TMath::Pi()) {
      phi2 += (2 * TMath::Pi());
    } else if (phi2 > TMath::Pi()) {
      phi2 -= (2 * TMath::Pi());
    }
    T deltaPhi = phi2 - phi1;
    if (deltaPhi < -TMath::Pi()) {
      deltaPhi += (2 * TMath::Pi());
    } else if (deltaPhi > TMath::Pi()) {
      deltaPhi -= (2 * TMath::Pi());
    }
    return deltaPhi;
  }

  void process(o2::aod::Collision const& collision, soa::Filtered<aod::Tracks> const& tracks) //o2::aod::Calos const& calos
  {
    jets.clear();
    inputParticles.clear();
    //up the events histogram
    hist_events_N->Fill(0);
    hist_events_z->Fill(collision.posZ());
    hist_events_xy->Fill(collision.posX(),collision.posY());
    //loop for charged tracks

    //std::vector<aod::Tracks> photonvec;
    //photonvec.clear();
    auto nGoodTracks = 0;
    for (auto& track : tracks) {
      nGoodTracks++;
      hist_tracks_etaphi->Fill(track.eta(),track.phi());

      auto energy = std::sqrt(track.p() * track.p() + JetFinder::mPion * JetFinder::mPion);
      inputParticles.emplace_back(track.px(), track.py(), track.pz(), energy);
      inputParticles.back().set_user_index(track.globalIndex());
      //photonvec.emplace_back(track);
    }
    hist_tracks_NGood->Fill(nGoodTracks);
    //need to find out how to get the energy
    //for (auto& calo : calos) {
      //hist_calo_pt->Fill(calo.CaloType());
    //}

    //jet finder
    fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
    for (const auto& jet : jets) {
      hist_jet_pt->Fill(jet.pt());
      hist_jet_etaphi->Fill(jet.eta(),jet.phi());
      //here starts the gamma-jet correlation calculations
      //for now the track is assumed to be the photon..
      if(jet.pt()>5.){
        for (auto& track : tracks) {
          if(track.pt()>1.){ //this should be replaced by the photon candidates
            auto trackphi = 0.0;
            auto deta = track.eta()-jet.eta();
            auto dphi = TMath::Abs(relativePhi(jet.phi(), trackphi));
            hist_gammajet_detadphi->Fill(deta,dphi);
            if(dphi> TMath::Pi() * 3./4.){
              //energy of the assumed photon (wrong)
              auto trackE = std::sqrt(track.p() * track.p());
              auto jetE = jet.E();
              hist_gammajet_dE->Fill(abs(trackE-jetE));
              hist_gammajet_xj->Fill(jet.pt()/track.pt());
              hist_gammajet_ratioE->Fill(jetE/trackE);
            }
          }
        }
      }
    }
    //loop for photons
    //todo still have to remove double counting
    // for (int i = 0; i<photonvec.size() ;i++) {
    //   for (int j = i+1; j<photonvec.size() ;j++) {
    //     // get track momenta
    //     //array<float, 3> photon1 = {photonvec.at(i).px(), photonvec.at(i).py(), photonvec.at(i).pz()};
    //     // for now assuming massless particles
    //     //auto theta = ;
    //     //auto invariantmass = sqrt(2.* e1 e2 * (1-costheta));
    //     hist_meson_masspt->Fill();
    //   }
    // }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<ProcessGammaJet>("ga-gammajet"),
  };
}
