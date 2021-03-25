
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

using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::pidRespTPC, aod::McTrackLabels>;


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

    Configurable<float> fLineCutZRSlope{"fLineCutZRSlope", tan(2*atan(exp(-fEtaCut))), "fLineCutZRSlope"};


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

    Configurable<float> fMinPhotonAsymmetry{"fMinPhotonAsymmetry", 0., "fMinPhotonAsymmetry"};

    Configurable<float> fMaxPhotonAsymmetry{"fMaxPhotonAsymmetry", 0.95, "fMaxPhotonAsymmetry"};

    //~ // case 10:  //a
    Configurable<float> fQtPtMax{"fQtPtMax", 0.11, "fQtPtMax"};

    Configurable<float> fQtMax{"fQtMax", 0.040, "fQtMax"};

    Configurable<int> fDoQtGammaSelection{"fDoQtGammaSelection", 2, "fDoQtGammaSelection"};

    Configurable<bool> fDo2DQt{"fDo2DQt", true, "fDo2DQt"};


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
    Size_t nPID = 0;
    Size_t nArmenteros = 0;
    Size_t nCosPAngleCut =0;
    Size_t nPassed =0;
    Size_t nConfirmedGammas=0;

    void init(InitContext const&);
        
    void process(aod::Collision const& collision, 
                 aod::V0Datas const& v0s, 
                 tracksAndTPCInfo const& tracks,
                 aod::McParticles const& mcParticles);
    
    // Armenteros Qt Cut
    template <typename T>
    Bool_t ArmenterosQtCut(T const& theV0);
    
    template <typename T>
    bool selectionPIDTPC_V0(const T& theTrack);
    
    template <typename T>
    bool selectionPIDTPC_track(const T& theTrack);
    
    HistogramRegistry registry{
    "registry",
    {
        {"IsPhotonSelected", "IsPhotonSelected", {HistType::kTH1F, {{12, -0.5f, 11.5f}}}},
        
        {"hPtRec", "hPtRec", {HistType::kTH1F, {{100, 0.0f, 25.0f}}}},
        {"hPtRes", "hPtRes_Rec-MC", {HistType::kTH1F, {{100, -0.5f, 0.5f}}}},
        
        {"hConvPointR", "hConvPointR", {HistType::kTH1F, {{100, 0.f, 300.f}}}},
        {"hConvPointRRes", "hConvPointRRes_Rec-MC", {HistType::kTH1F, {{100, -200.f, 200.f}}}},
        {"hConvPointAbsoluteDistanceRes", "hConvPointAbsoluteDistanceRes", {HistType::kTH1F, {{100, -0.0f, 200.f}}}},
        
        {"hEtaRes", "hEtaRes_Rec-MC", {HistType::kTH1F, {{100, -3.14f, 3.14f}}}},
        {"hPhiRes", "hPhiRes_Rec-MC", {HistType::kTH1F, {{100, -3.14f, 3.14f}}}},
    },
  };
  
    enum photonCuts {
            kPhotonIn=0,
            kOnFly,
            kNoTracks,
            kNoV0,
            kTrackCuts,
            kdEdxCuts,
            kConvPointFail,
            kPhotonCuts,
            kEventPlane,
            kPhotonOut
    };
    
    std::vector<TString> fPhotCutsLabels{
        "kPhotonIn",
        "kOnFly",
        "kNoTracks",
        "kNoV0",
        "kTrackCuts",
        "kdEdxCuts",
        "kConvPointFail",
        "kPhotonCuts",
        "kEventPlane",
        "kPhotonOut"};
};

// just the histos from ATask to fill them for all confirmed photons without additional cuts
struct BTask
{
    void process(aod::Collision const& collision, 
                 aod::V0Datas const& v0s, 
                 tracksAndTPCInfo const& tracks,
                 aod::McParticles const& mcParticles);
        
    HistogramRegistry registry{
    "registry",
    {
        {"hPtRec", "hPtRec", {HistType::kTH1F, {{100, 0.0f, 25.0f}}}},
        {"hPtRes", "hPtRes_Rec-MC", {HistType::kTH1F, {{100, -0.5f, 0.5f}}}},
        
        {"hConvPointR", "hConvPointR", {HistType::kTH1F, {{100, 0.f, 300.f}}}},
        {"hConvPointRRes", "hConvPointRRes_Rec-MC", {HistType::kTH1F, {{100, -200.f, 200.f}}}},
        {"hConvPointAbsoluteDistanceRes", "hConvPointAbsoluteDistanceRes", {HistType::kTH1F, {{100, -0.0f, 200.f}}}},
        
        {"hEtaRes", "hEtaRes_Rec-MC", {HistType::kTH1F, {{100, -3.14f, 3.14f}}}},
        {"hPhiRes", "hPhiRes_Rec-MC", {HistType::kTH1F, {{100, -3.14f, 3.14f}}}},
    },
  };
};
