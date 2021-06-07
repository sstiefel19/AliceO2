
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include "AnalysisDataModel/StrangenessTables.h"
#include "AnalysisDataModel/HFSecondaryVertex.h" // for BigTracks

#include "AnalysisDataModel/PID/PIDResponse.h"
#include "AnalysisDataModel/PID/PIDTPC.h"

#include <TH1F.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidRespTPCEl, aod::pidRespTPCPi, aod::McTrackLabels>;


template <typename T>
Double_t computePsiPair(T const& theV0)
{
  TVector3 lPos(theV0.pxpos(), theV0.pypos(), theV0.pzpos());
  TVector3 lNeg(theV0.pxneg(), theV0.pyneg(), theV0.pzneg());

  Double_t lXsi = lPos.Angle(lNeg);

  if (lXsi > 0.0) {
    Double_t lDeltaTheta = lNeg.Theta() - lPos.Theta();
    Double_t lArg = lDeltaTheta / lXsi;

    if (lArg >= -1. && lArg <= 1.) {
      return asin(lArg);
    } else {
      return 100.;
    }
  }
  return -100.;
}

struct GammaConversionsMC {

  Configurable<float> fSinglePtCut{"fSinglePtCut", 0.04, "minimum daughter track pt"};

  Configurable<float> fEtaCut{"fEtaCut", 0.8, "accepted eta range"};

  Configurable<float> fMinR{"fMinR", 5., "minimum conversion radius of the V0s"};
  Configurable<float> fMaxR{"fMaxR", 180., "maximum conversion radius of the V0s"};

  Configurable<float> fPIDnSigmaBelowElectronLine{"fPIDnSigmaBelowElectronLine", -3., "minimum sigma electron PID for V0 daughter tracks"};

  Configurable<float> fPIDnSigmaAboveElectronLine{"fPIDnSigmaAboveElectronLine", 3., "maximum sigma electron PID for V0 daughter tracks"};

  Configurable<float> fPIDnSigmaAbovePionLine{"fPIDnSigmaAbovePionLine", 3., "minimum sigma to be over the pion line for low momentum tracks"}; //case 4: 3.0sigma, 1.0 sigma at high momentum
  
  Configurable<float> fPIDnSigmaAbovePionLineHighP{"fPIDnSigmaAbovePionLineHighP", 1., "minimum sigma to be over the pion line for high momentum tracks"};

  Configurable<float> fPIDMinPnSigmaAbovePionLine{"fPIDMinPnSigmaAbovePionLine", 0.4, "minimum track momentum to apply any pion rejection"}; //case 7:  // 0.4 GeV

  Configurable<float> fPIDMaxPnSigmaAbovePionLine{"fPIDMaxPnSigmaAbovePionLine", 8., "border between low and high momentum pion rejection"}; //case 7:  // 8. GeV

  Configurable<float> fQtPtMax{"fQtPtMax", 0.11, "up to fQtMax, multiply the pt of the V0s by this value to get the maximum qt "};

  Configurable<float> fQtMax{"fQtMax", 0.040, "maximum qt"};
  Configurable<float> fMaxPhotonAsymmetry{"fMaxPhotonAsymmetry", 0.95, "maximum photon asymetry"};
  Configurable<float> fPsiPairCut{"fPsiPairCut", 0.1, "maximum psi angle of the track pair"};

  Configurable<float> fCosPAngleCut{"fCosPAngleCut", 0.85, "mimimum cosinus of the pointing angle"}; // case 4

  HistogramRegistry registry{
    "registry",
    {
      {"IsPhotonSelected", "IsPhotonSelected", {HistType::kTH1F, {{12, -0.5f, 11.5f}}}},

      {"beforeCuts/hPtRec_before", "hPtRec_before", {HistType::kTH1F, {{100, 0.0f, 25.0f}}}},
      {"beforeCuts/hEtaRec_before", "hEtaRec_before", {HistType::kTH1F, {{1000, -2.f, 2.f}}}},
      {"beforeCuts/hPhiRec_before", "hEtaRec_before", {HistType::kTH1F, {{1000, 0.f, 2.f * M_PI}}}},
      {"beforeCuts/hConvPointR_before", "hConvPointR_before", {HistType::kTH1F, {{800, 0.f, 200.f}}}},

      {"hPtRec", "hPtRec", {HistType::kTH1F, {{100, 0.0f, 25.0f}}}},
      {"hEtaRec", "hEtaRec", {HistType::kTH1F, {{1000, -2.f, 2.f}}}},
      {"hPhiRec", "hEtaRec", {HistType::kTH1F, {{1000, 0.f, 2.f * M_PI}}}},
      {"hConvPointR", "hConvPointR", {HistType::kTH1F, {{800, 0.f, 200.f}}}},

      {"hTPCdEdxSig", "hTPCdEdxSig", {HistType::kTH2F, {{150, 0.03f, 20.f}, {400, -10.f, 10.f}}}},
      {"hTPCdEdx", "hTPCdEdx", {HistType::kTH2F, {{150, 0.03f, 20.f}, {800, 0.f, 200.f}}}},

      {"hArmenteros", "hArmenteros", {HistType::kTH2F, {{200, -1.f, 1.f}, {250, 0.f, 25.f}}}},
      {"hPsiPtRec", "hPsiPtRec", {HistType::kTH2F, {{500, -2.f, 2.f}, {100, 0.f, 25.f}}}},

      {"hCosPAngle", "hCosPAngle", {HistType::kTH1F, {{1000, -1.f, 1.f}}}},

      // resolution histos
      {"resolutions/hPtRes", "hPtRes_Rec-MC", {HistType::kTH1F, {{100, -0.5f, 0.5f}}}},
      {"resolutions/hEtaRes", "hEtaRes_Rec-MC", {HistType::kTH1F, {{100, -2.f, 2.f}}}},
      {"resolutions/hPhiRes", "hPhiRes_Rec-MC", {HistType::kTH1F, {{100, -M_PI, M_PI}}}},

      {"resolutions/hConvPointRRes", "hConvPointRRes_Rec-MC", {HistType::kTH1F, {{100, -200.f, 200.f}}}},
      {"resolutions/hConvPointAbsoluteDistanceRes", "hConvPointAbsoluteDistanceRes", {HistType::kTH1F, {{100, -0.0f, 200.f}}}},
    },
  };

  enum photonCuts {
    kPhotonIn = 0,
    kTrackEta,
    kTrackPt,
    kElectronPID,
    kPionRejLowMom,
    kPionRejHighMom,
    kV0Radius,
    kArmenteros,
    kPsiPair,
    kCosinePA,
    kPhotonOut
  };

  std::vector<TString> fPhotCutsLabels{
    "kPhotonIn",
    "kTrackEta",
    "kTrackPt",
    "kElectronPID",
    "kPionRejLowMom",
    "kPionRejHighMom",
    "kV0Radius",
    "kArmenteros",
    "kPsiPair",
    "kCosinePA",
    "kPhotonOut"};
    
  void init(InitContext const&)
  {
    TAxis* lXaxis = registry.get<TH1>(HIST("IsPhotonSelected"))->GetXaxis();
    for (size_t i = 0; i < fPhotCutsLabels.size(); ++i) {
      lXaxis->SetBinLabel(i + 1, fPhotCutsLabels[i]);
    }
  }

  void process(aod::Collision const& theCollision,
               aod::V0Datas const& theV0s,
               tracksAndTPCInfo const& theTracks,
               aod::McParticles const& theMcParticles)
  {
    for (auto& lV0 : theV0s) {
      float lPhiV0Rec = static_cast<float>(M_PI) + std::atan2(-lV0.py(), -lV0.px());

      // fill some QA histograms before any cuts
      registry.fill(HIST("beforeCuts/hPtRec_before"), lV0.pt());
      registry.fill(HIST("beforeCuts/hEtaRec_before"), lV0.eta());
      registry.fill(HIST("beforeCuts/hPhiRec_before"), lPhiV0Rec);
      registry.fill(HIST("beforeCuts/hConvPointR_before"), lV0.v0radius());
      registry.fill(HIST("IsPhotonSelected"), kPhotonIn);

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfo>(); //positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfo>(); //negative daughter

      // single track eta cut
      if (TMath::Abs(lTrackPos.eta()) > fEtaCut || TMath::Abs(lTrackNeg.eta()) > fEtaCut) {
        registry.fill(HIST("IsPhotonSelected"), kTrackEta);
        continue;
      }

      // single track pt cut
      if (lTrackPos.pt() < fSinglePtCut || lTrackNeg.pt() < fSinglePtCut) {
        registry.fill(HIST("IsPhotonSelected"), kTrackPt);
        continue;
      }

      if (!(selectionPIDTPC_track(lTrackPos) && selectionPIDTPC_track(lTrackNeg))) {
        continue;
      }

      if (lV0.v0radius() < fMinR || lV0.v0radius() > fMaxR) {
        registry.fill(HIST("IsPhotonSelected"), kV0Radius);
        continue;
      }

      if (!ArmenterosQtCut(lV0.alpha(), lV0.qtarm(), lV0.pt())) {
        registry.fill(HIST("IsPhotonSelected"), kArmenteros);
        continue;
      }

      Double_t lPsiPair = computePsiPair(lV0);
      if (TMath::Abs(lPsiPair) > fPsiPairCut) {
        // todo: rethink this
        if (TMath::Abs(lPsiPair) > 90.) {
          registry.fill(HIST("hPsiPtRec"), lPsiPair, lV0.pt());
        }
        registry.fill(HIST("IsPhotonSelected"), kPsiPair);
        continue;
      }

      auto lV0CosinePA = RecoDecay::CPA(array{theCollision.posX(), theCollision.posY(), theCollision.posZ()}, array{lV0.x(), lV0.y(), lV0.z()}, array{lV0.px(), lV0.py(), lV0.pz()});
      if (lV0CosinePA < fCosPAngleCut) {
        registry.fill(HIST("IsPhotonSelected"), kCosinePA);
        continue;
      }

      registry.fill(HIST("IsPhotonSelected"), kPhotonOut);

      registry.fill(HIST("hPtRec"), lV0.pt());
      registry.fill(HIST("hEtaRec"), lV0.eta());
      registry.fill(HIST("hPhiRec"), lPhiV0Rec);
      registry.fill(HIST("hConvPointR"), lV0.v0radius());

      registry.fill(HIST("hTPCdEdxSig"), lTrackNeg.p(), lTrackNeg.tpcNSigmaEl());
      registry.fill(HIST("hTPCdEdxSig"), lTrackPos.p(), lTrackPos.tpcNSigmaEl());

      registry.fill(HIST("hTPCdEdx"), lTrackNeg.p(), lTrackNeg.tpcSignal());
      registry.fill(HIST("hTPCdEdx"), lTrackPos.p(), lTrackPos.tpcSignal());

      registry.fill(HIST("hArmenteros"), lV0.alpha(), lV0.qtarm());
      registry.fill(HIST("hPsiPtRec"), lPsiPair, lV0.pt());

      registry.fill(HIST("hCosPAngle"), lV0CosinePA);

      // now check if we have an actual photon
      // todo: verify it is enough to check only mother0 being equal
      if (lTrackPos.mcParticle().mother0() > -1 &&
          lTrackPos.mcParticle().mother0() == lTrackNeg.mcParticle().mother0()) {
        auto mother = theMcParticles.iteratorAt(lTrackPos.mcParticle().mother0());

        if (mother.pdgCode() == 22) {

          registry.fill(HIST("resolutions/hPtRes"), lV0.pt() - mother.pt());
          registry.fill(HIST("resolutions/hEtaRes"), lV0.eta() - mother.eta());
          registry.fill(HIST("resolutions/hPhiRes"), lPhiV0Rec - mother.phi());

          TVector3 lConvPointRec(lV0.x(), lV0.y(), lV0.z());
          TVector3 lPosTrackVtxMC(lTrackPos.mcParticle().vx(), lTrackPos.mcParticle().vy(), lTrackPos.mcParticle().vz());
          // take the origin of the positive mc track as conversion point (should be identical with negative, verified this on a few photons)
          TVector3 lConvPointMC(lPosTrackVtxMC);

          registry.fill(HIST("resolutions/hConvPointRRes"), lConvPointRec.Perp() - lConvPointMC.Perp());
          registry.fill(HIST("resolutions/hConvPointAbsoluteDistanceRes"), TVector3(lConvPointRec - lConvPointMC).Mag());
        }
      } 
    } 
  }

  Bool_t ArmenterosQtCut(Double_t theAlpha, Double_t theQt, Double_t thePt)
  {
    // in AliPhysics this is the cut for if fDo2DQt && fDoQtGammaSelection == 2
    Float_t lQtMaxPtDep = fQtPtMax * thePt;
    if (lQtMaxPtDep > fQtMax) {
      lQtMaxPtDep = fQtMax;
    }
    if (!(TMath::Power(theAlpha / fMaxPhotonAsymmetry, 2) + TMath::Power(theQt / lQtMaxPtDep, 2) < 1)) {
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename T>
  bool selectionPIDTPC_track(const T& theTrack)
  {
    // TPC Electron Line
    if (theTrack.tpcNSigmaEl() < fPIDnSigmaBelowElectronLine || theTrack.tpcNSigmaEl() > fPIDnSigmaAboveElectronLine) {
      registry.fill(HIST("IsPhotonSelected"), kElectronPID);
      return kFALSE;
    }

    // TPC Pion Line
    if (theTrack.p() > fPIDMinPnSigmaAbovePionLine) {
      // low pt Pion rej
      if (theTrack.p() < fPIDMaxPnSigmaAbovePionLine) {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLine) {
          registry.fill(HIST("IsPhotonSelected"), kPionRejLowMom);
          return kFALSE;
        }
      }
      // High Pt Pion rej
      else {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineHighP) {
          registry.fill(HIST("IsPhotonSelected"), kPionRejHighMom);
          return kFALSE;
        }
      }
    }
    return kTRUE;
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GammaConversionsMC>(cfgc)};
}
