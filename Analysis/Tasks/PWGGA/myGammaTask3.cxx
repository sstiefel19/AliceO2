// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "myGammaTask3.h"


// This is a very simple example showing how to create an histogram
// FIXME: this should really inherit from AnalysisTask but
//        we need GCC 7.4+ for that

  
template <typename T>
TVector3 getMomPos(T const& theV0)
{
    return TVector3(theV0.pxpos(),theV0.pypos(),theV0.pzpos());
}
 
template <typename T>
TVector3 getMomNeg(T const& theV0)
{
    return TVector3(theV0.pxneg(),theV0.pyneg(),theV0.pzneg());
}

// todo: rethink the handling of lXsi<=0 and lArg with +-100
template <typename T>
Double_t computePsiPair(T const& theV0)
{
    TVector3 lPos(getMomPos(theV0));
    TVector3 lNeg(getMomNeg(theV0));
    
    Double_t lXsi = lPos.Angle(lNeg);
    
    if (lXsi > 0.0 )
    {
        Double_t lDeltaTheta = lNeg.Theta() - lPos.Theta();
        Double_t lArg        = lDeltaTheta / lXsi;
        
        if (lArg >= -1. && lArg <= 1.)
        {
            return asin(lArg);
        }
        else
        {
            return 100.;
        }
    }
    return -100.; 
}
 
// probably can be improved
// checks for division by 0?
template <typename T>
void computeArmenteros(T const& theV0, Double_t &theAlpha, Double_t &theQt)
{
    TVector3 lPos(getMomPos(theV0));
    TVector3 lNeg(getMomNeg(theV0));
    
    TVector3 lMother = lPos + lNeg; // SFS verify this is legit
    
    // SFS is this the most efficient way?
    Double_t thetaP = acos((lPos*lMother) / (lPos.Mag() * lMother.Mag()));
    Double_t thetaN = acos((lNeg*lMother) / (lNeg.Mag() * lMother.Mag()));
    
    // could also calculate this directly from scalar products
    theAlpha = ((lPos.Mag())*cos(thetaP) - (lNeg.Mag())*cos(thetaN)) / 
               ((lPos.Mag())*cos(thetaP) + (lNeg.Mag())*cos(thetaN))  ;
    
    // and this with cross product 
    theQt = lPos.Mag() * sin(thetaP);
}


// Armenteros Qt Cut
Bool_t ATask::ArmenterosQtCut(Double_t theAlpha, Double_t theQt, Double_t thePt){  
    
    if(fDo2DQt){
       if(fDoQtGammaSelection==2){
        Float_t qtMaxPtDep = fQtPtMax*thePt;
        if (qtMaxPtDep > fQtMax)
          qtMaxPtDep      = fQtMax;
        if ( !(TMath::Power(theAlpha/fMaxPhotonAsymmetry,2)+TMath::Power(theQt/qtMaxPtDep,2) < 1) ){
          return kFALSE;
        }
      }
    }     
    return kTRUE;
}


  
// returns false if track not selected
// in original from run2 there is an if(fDoElecDeDxPostCalibration)
template <typename T>
bool ATask::selectionPIDTPC_track(const T& theTrack)
{
    // TPC Electron Line
    if( theTrack.tpcNSigmaEl() < fPIDnSigmaBelowElectronLine || theTrack.tpcNSigmaEl() > fPIDnSigmaAboveElectronLine)
    {
        registry.fill(HIST("IsPhotonSelected"), kElectronPID);
        return kFALSE;
    }
    
    // TPC Pion Line
    if( theTrack.p()>fPIDMinPnSigmaAbovePionLine)
    {
        // low pt Pion rej
        if (theTrack.p()<fPIDMaxPnSigmaAbovePionLine )
        {
            if( theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi()<fPIDnSigmaAbovePionLine)
            {
                registry.fill(HIST("IsPhotonSelected"), kPionRejLowMom);
                return kFALSE;
            }
        }
        // High Pt Pion rej
        else
        {
            if( theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi()<fPIDnSigmaAbovePionLineHighPt)
            {
                registry.fill(HIST("IsPhotonSelected"), kPionRejHighMom);
                return kFALSE;
            }
        }
    }
    return kTRUE;  
}



void ATask::init(InitContext const&)
{
    TAxis *lXaxis = registry.get<TH1>(HIST("IsPhotonSelected"))->GetXaxis();
    for (size_t i=0; i< fPhotCutsLabels.size(); ++i)
    {
        lXaxis->SetBinLabel(i+1, fPhotCutsLabels[i]);
    }
}
  
void ATask::process(aod::Collision const& collision, 
                    aod::V0Datas const& v0s, 
                    tracksAndTPCInfo const& tracks,
                    aod::McParticles const& mcParticles)
{
    for (auto& v0 : v0s) {
        ++i;
        float lPhiV0Rec = static_cast<float>(M_PI) + std::atan2(-v0.py(), -v0.px());
        
        // fill some QA histograms before any cuts
        registry.fill(HIST("beforeCuts/hPtRec_before"), v0.pt());
        registry.fill(HIST("beforeCuts/hEtaRec_before"), v0.eta());
        registry.fill(HIST("beforeCuts/hPhiRec_before"), lPhiV0Rec);
        registry.fill(HIST("beforeCuts/hConvPointR_before"), v0.v0radius());

        registry.fill(HIST("IsPhotonSelected"), kPhotonIn);
        
        auto trackPos = v0.template posTrack_as<tracksAndTPCInfo>(); //positive daughter
        auto trackNeg = v0.template negTrack_as<tracksAndTPCInfo>(); //positive daughter
        
        // single track eta cut
        if (TMath::Abs(trackPos.eta()) > fEtaCut || TMath::Abs(trackNeg.eta()) > fEtaCut)
        {
            registry.fill(HIST("IsPhotonSelected"), kTrackEta);
            continue;
        }
        
        // single track pt cut    
        if ( (TMath::Power(v0.pxpos(),2) + TMath::Power(v0.pypos(),2) < fSinglePtSq) ||
             (TMath::Power(v0.pxneg(),2) + TMath::Power(v0.pyneg(),2) < fSinglePtSq)
           )
        {
            registry.fill(HIST("IsPhotonSelected"), kTrackPt);
            continue;
        } 
        
        if (!(selectionPIDTPC_track(trackPos) && selectionPIDTPC_track(trackNeg)))
        {
            continue;
        }
        
        if (v0.v0radius() < fMinR || v0.v0radius() > fMaxR)
        {
                registry.fill(HIST("IsPhotonSelected"), kV0Radius);
                continue;
        }
        
        Double_t lAlpha = -1.;
        Double_t lQt = -1.;
        computeArmenteros(v0, lAlpha, lQt);
        //todo: check if values of lAlpha and Qt need to be checked
        if (!ArmenterosQtCut(lAlpha, lQt, v0.pt()))
        {
            registry.fill(HIST("IsPhotonSelected"), kArmenteros);
            continue;
        }
        
        Double_t lPsiPair = computePsiPair(v0);
        if (TMath::Abs(lPsiPair) > fPsiPairCut)
        {   
            // todo: rethink this
            if (TMath::Abs(lPsiPair) > 90.)
            {
                registry.fill(HIST("hPsiPtRec"), lPsiPair, v0.pt());
            }
            registry.fill(HIST("IsPhotonSelected"), kPsiPair);
            continue;
        }
        
        auto lV0CosinePA = RecoDecay::CPA(array{collision.posX(), collision.posY(), collision.posZ()}, array{v0.x(), v0.y(), v0.z()}, array{v0.px(), v0.py(), v0.pz()});
        if (lV0CosinePA < fCosPAngleCut)
        {   
            registry.fill(HIST("IsPhotonSelected"), kCosinePA);
            continue;
        } 
        
        ++nPassed;  
        registry.fill(HIST("IsPhotonSelected"), kPhotonOut);
        
        registry.fill(HIST("hPtRec"), v0.pt());
        registry.fill(HIST("hEtaRec"), v0.eta());
        registry.fill(HIST("hPhiRec"), lPhiV0Rec);
        registry.fill(HIST("hConvPointR"), v0.v0radius());

        registry.fill(HIST("hPsiPtRec"), lPsiPair, v0.pt());
        
        registry.fill(HIST("hTPCdEdxSigafter"), trackNeg.p(), trackNeg.tpcNSigmaEl());
        registry.fill(HIST("hTPCdEdxSigafter"), trackPos.p(), trackPos.tpcNSigmaEl());

        registry.fill(HIST("hTPCdEdxafter"), trackNeg.p(), trackNeg.tpcSignal());
        registry.fill(HIST("hTPCdEdxafter"), trackPos.p(), trackPos.tpcSignal());
        
        registry.fill(HIST("hArmenterosafter"), lAlpha, lQt);
        
        registry.fill(HIST("hCosPAngle"), lV0CosinePA);

        cout << "SFS signals: " << trackNeg.tpcSignal() << " " << trackPos.tpcSignal() << endl;

        // now check if we have an actual photon
        // todo: verify it is enough to check only mother0 being equal
        if (trackPos.mcParticle().mother0() > -1 && 
            trackPos.mcParticle().mother0() == trackNeg.mcParticle().mother0())
        {
            auto mother = mcParticles.iteratorAt(trackPos.mcParticle().mother0());
                        
            LOGF(info, "SFS mother: %d", mother.pdgCode());
            if (mother.pdgCode() == 22) 
            {
                ++nConfirmedGammas;
                
                registry.fill(HIST("resolutions/hPtRes"), v0.pt() - mother.pt());
                registry.fill(HIST("resolutions/hEtaRes"), v0.eta() - mother.eta()); 
                registry.fill(HIST("resolutions/hPhiRes"), lPhiV0Rec - mother.phi());
                
                TVector3 lConvPointRec(v0.x(), v0.y(), v0.z());
                TVector3 lPosTrackVtxMC(trackPos.mcParticle().vx(), trackPos.mcParticle().vy(), trackPos.mcParticle().vz());
                // take the origin of the positive mc track as conversion point (should be identical with negative, verified this on a few photons)
                TVector3 lConvPointMC(lPosTrackVtxMC);
                
                registry.fill(HIST("resolutions/hConvPointRRes"), lConvPointRec.Perp() - lConvPointMC.Perp());
                registry.fill(HIST("resolutions/hConvPointAbsoluteDistanceRes"), TVector3(lConvPointRec - lConvPointMC).Mag());
            }
        }
    }
    cout << "SFS total v0s: " << i << " passed all cuts " << nPassed << " ratio " << nPassed / (Float_t)i << " nConfirmedGammas: " << nConfirmedGammas << " purity: " << nConfirmedGammas / (Float_t)nPassed <<  endl;
}




WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{ 
    return WorkflowSpec{
        adaptAnalysisTask<ATask>(cfgc, TaskName{"myGammaTask3"})
    };
}
