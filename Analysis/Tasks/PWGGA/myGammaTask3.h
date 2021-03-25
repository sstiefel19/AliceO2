
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

using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidRespTPC, aod::McTrackLabels>;


using std::cout;
using std::endl;


struct ATask
{
    // aodConversionCutnumber  06000008d00000001100000000
    // V0ReaderPhotonCutnumber 00000070000000000500004000
    // from trainconfig 681    0dm00009a4770c00amd0404000
    /* cuts to be implemented:
     * 
     * from V0ReaderPhotonCutnumber
     * SinglePtCut 6 = case 6:  // 0.04 GeV  + min gamma pT cut of 10 MeV
        fSinglePtCut = 0.040;
        fPtCut       = 0.01;
     * 
     * from photon cutnumber from trainconfig 681 0dm00009a4770c00amd0404000
     * eta cut d -> actually not used?! AcceptanceCut completely commented out in AliConversionPhotonCuts
     * MinRCut m (also maxR , 5-180cm)
     * ClsTPCCut 9  [left out for the moment]
     * ededxSigmaCut a
     * pidedxSigmaCut 4
     * piMomdedxSigmaCut 7
     * piMaxMomdedxSigmaCut 7
     * TOFelectronPID c do later
     * QtMaxCut a
     * Chi2GammaCut m
     * PsiPair d
     * CosinePointingAngle 4
     * RejectToCloseV0s 4
     *  
     * 
     * 
    */


    // make a pt cut
    Configurable<float> fSinglePtSq{"fSinglePtSq", TMath::Power(0.04,2), "fSinglePtSq"};

    Configurable<float> fEtaCut{"fEtaCut", 0.8, "fEtaCut"};

    Configurable<float> fMinR{"fMinR", 5., "fMinR"};
    Configurable<float> fMaxR{"fMaxR", 180., "fMaxR"};

    // do later, sth like fraction of findable clusters 
    Configurable<float> fMinClsTPCToF{"fMinClsTPCToF", 0.6, "fMinClsTPCToF"};

    Configurable<int> fUseCorrectedTPCClsInfo{"fUseCorrectedTPCClsInfo", 1, "fUseCorrectedTPCClsInfo"};

    Configurable<bool> fDoElecDeDxPostCalibration{"fDoElecDeDxPostCalibration", false, "fDoElecDeDxPostCalibration"};
    
    Configurable<float> fPIDnSigmaBelowElectronLine{"fPIDnSigmaBelowElectronLine", -3., "fPIDnSigmaBelowElectronLine"};

    Configurable<float> fPIDnSigmaAboveElectronLine{"fPIDnSigmaAboveElectronLine", 3., "fPIDnSigmaAboveElectronLine"};

    //case 4:  // 3.0sigma, 1.0 sigma at high pt
    Configurable<float> fPIDnSigmaAbovePionLine{"fPIDnSigmaAbovePionLine", 3., "fPIDnSigmaAbovePionLine"};

    Configurable<float> fPIDnSigmaAbovePionLineHighPt{"fPIDnSigmaAbovePionLineHighPt", 1., "fPIDnSigmaAbovePionLineHighPt"};

     //case 7:  // 0.4 GeV
    Configurable<float> fPIDMinPnSigmaAbovePionLine{"fPIDMinPnSigmaAbovePionLine", 0.4, "fPIDMinPnSigmaAbovePionLine"};

    //case 7:  // 8. GeV
    Configurable<float> fPIDMaxPnSigmaAbovePionLine{"fPIDMaxPnSigmaAbovePionLine", 8., "fPIDMaxPnSigmaAbovePionLine"};

    
    //~ // case 10:  //a
    Configurable<float> fQtPtMax{"fQtPtMax", 0.11, "fQtPtMax"};

    Configurable<float> fQtMax{"fQtMax", 0.040, "fQtMax"};
    
    Configurable<int> fDoQtGammaSelection{"fDoQtGammaSelection", 2, "fDoQtGammaSelection"};

    Configurable<bool> fDo2DQt{"fDo2DQt", true, "fDo2DQt"};
    
    Configurable<float> fMaxPhotonAsymmetry{"fMaxPhotonAsymmetry", 0.95, "fMaxPhotonAsymmetry"};
    
    Configurable<float> fPsiPairCut{"fPsiPairCut", 0.1, "fPsiPairCut"};


    //case 22: // m for exp cut (fDo2DPsiPairChi2 = 2)
     Float_t   fChi2CutConversion = 30.;
     Float_t   fChi2CutConversionExpFunc = -0.15;
          
    // case 4:
    Configurable<float> fCosPAngleCut{"fCosPAngleCut", 0.85, "fCosPAngleCut"};

    //case 4:
    Configurable<bool> fDoToCloseV0sCut{"fDoToCloseV0sCut", true, "fDoToCloseV0sCut"};

    Configurable<bool> fDoDoubleCountingCut{"fDoDoubleCountingCut", true, "fDoDoubleCountingCut"};

    Configurable<float> fDeltaR{"fDeltaR", 6., "fDeltaR"};

    Configurable<float> fOpenAngle{"fOpenAngle", 0.02, "fOpenAngle"};


    // some variables for counting    
    Size_t i = 0;
    Size_t nPassed =0;
    Size_t nConfirmedGammas=0;

    void init(InitContext const&);
        
    void process(aod::Collision const& collision, 
                 aod::V0Datas const& v0s, 
                 tracksAndTPCInfo const& tracks,
                 aod::McParticles const& mcParticles);
    
    // Armenteros Qt Cut
    Bool_t ArmenterosQtCut(Double_t theAlpha, Double_t theQt, Double_t thePt);
    
    
    template <typename T>
    bool selectionPIDTPC_track(const T& theTrack);
    
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
        
        
        {"hPsiPtRec", "hPsiPtRec", {HistType::kTH2F, {{500, -2.f, 2.f}, {100, 0.f, 25.f}}}},
        
        {"hTPCdEdxSigafter", "hTPCdEdxSigafter", {HistType::kTH2F, {{150,0.03f,20.f}, {400, -10.f,10.f}}}},
        
        {"hTPCdEdxafter", "hTPCdEdxafter", {HistType::kTH2F, {{150,0.03f,20.f}, {800,0.f,200.f}}}},
        
        {"hArmenterosafter", "hArmenterosafter", {HistType::kTH2F, {{200,-1.f,1.f}, {250,0.f,25.f}}}},
        
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
            kPhotonIn=0,
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
};
