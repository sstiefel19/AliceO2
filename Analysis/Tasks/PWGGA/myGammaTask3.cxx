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

  

// can I use move constructor?
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
 
// probably can be improved
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
template <typename T>
Bool_t ATask::ArmenterosQtCut(T const& theV0){  
    
    Double_t lAlpha = -1.;
    Double_t lQt = -1;
    computeArmenteros(theV0, lAlpha, lQt);
    //todo: check if values of lAlpha and Qt need to be checked
    
    if(fDo2DQt){
       if(fDoQtGammaSelection==2){
        Float_t qtMaxPtDep = fQtPtMax*theV0.pt();
        if (qtMaxPtDep > fQtMax)
          qtMaxPtDep      = fQtMax;
        if ( !(TMath::Power(lAlpha/fMaxPhotonAsymmetry,2)+TMath::Power(lQt/qtMaxPtDep,2) < 1) ){
          return kFALSE;
        }
      }
    }     
    return kTRUE;
}

template <typename T>
bool ATask::selectionPIDTPC_V0(const T &theV0)
{
    auto trackPos = theV0.template posTrack_as<tracksAndTPCInfo>(); //positive daughter
    auto trackNeg = theV0.template negTrack_as<tracksAndTPCInfo>(); //negative daughter
    return selectionPIDTPC_track(trackPos) && selectionPIDTPC_track(trackNeg);
}
  
// returns false if track not selected
// in original from run2 there is an if(fDoElecDeDxPostCalibration)
template <typename T>
bool ATask::selectionPIDTPC_track(const T& theTrack)
{
    // TPC Electron Line
    if( theTrack.tpcNSigmaEl() < fPIDnSigmaBelowElectronLine || theTrack.tpcNSigmaEl() > fPIDnSigmaAboveElectronLine)
    {
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
                return kFALSE;
            }
        }
        // High Pt Pion rej
        else
        {
            if( theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi()<fPIDnSigmaAbovePionLineHighPt)
            {
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
        
        registry.fill(HIST("IsPhotonSelected"), kPhotonIn);
        
        auto trackPos = v0.template posTrack_as<tracksAndTPCInfo>(); //positive daughter
        auto trackNeg = v0.template negTrack_as<tracksAndTPCInfo>(); //positive daughter
        
        // single track pt cut    
        if ( (TMath::Power(v0.pxpos(),2) + TMath::Power(v0.pypos(),2) < fSinglePtSq) ||
             (TMath::Power(v0.pxneg(),2) + TMath::Power(v0.pyneg(),2) < fSinglePtSq)
           )
        {
            registry.fill(HIST("IsPhotonSelected"), kTrackCuts);
            continue;
        } 
        
        // single track eta cut
        if (TMath::Abs(trackPos.eta()) > fEtaCut || TMath::Abs(trackNeg.eta()) > fEtaCut)
        {
            registry.fill(HIST("IsPhotonSelected"), kTrackCuts);
            continue;
        }
        
        if (!selectionPIDTPC_V0(v0))
        {
            registry.fill(HIST("IsPhotonSelected"), kdEdxCuts);
            ++nPID;
            continue;
        }
        
        if (v0.v0radius() < fMinR || v0.v0radius() > fMaxR)
        {
                registry.fill(HIST("IsPhotonSelected"), kPhotonCuts);
                continue;
        }
        
        if (!ArmenterosQtCut(v0))
        {
            registry.fill(HIST("IsPhotonSelected"), kPhotonCuts);
            ++nArmenteros;
            continue;
        }
        
        auto V0CosinePA = RecoDecay::CPA(array{collision.posX(), collision.posY(), collision.posZ()}, array{v0.x(), v0.y(), v0.z()}, array{v0.px(), v0.py(), v0.pz()});
        if (V0CosinePA < fCosPAngleCut)
        {   
            registry.fill(HIST("IsPhotonSelected"), kPhotonCuts);
            ++nCosPAngleCut;
            continue;
        } 
        
        ++nPassed;  
        registry.fill(HIST("IsPhotonSelected"), kPhotonOut);
        
        registry.fill(HIST("hPtRec"), v0.pt());
        
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
                                
                registry.fill(HIST("hPtRes"), v0.pt() - mother.pt());
                
                TVector3 lConvPointRec(v0.x(), v0.y(), v0.z());
                TVector3 lPosTrackVtxMC(trackPos.mcParticle().vx(), trackPos.mcParticle().vy(), trackPos.mcParticle().vz());
                // take the origin of the positive mc track as conversion point (should be identical with negative, verified this on a few photons)
                TVector3 lConvPointMC(lPosTrackVtxMC);
                
                registry.fill(HIST("hConvPointR"), lConvPointRec.Perp());
                registry.fill(HIST("hConvPointRRes"), lConvPointRec.Perp() - lConvPointMC.Perp());
                registry.fill(HIST("hConvPointAbsoluteDistanceRes"), TVector3(lConvPointRec - lConvPointMC).Mag());
                
                float lPhiV0Rec = static_cast<float>(M_PI) + std::atan2(-v0.py(), -v0.px());
                registry.fill(HIST("hEtaRes"), v0.eta() - mother.eta()); 
                registry.fill(HIST("hPhiRes"), lPhiV0Rec - mother.phi());
            }
        }
    }
    cout << "SFS total v0s: " << i << " nPID: " << nPID << " nArmenteros: " << nArmenteros << " nCosPAngleCut: " << nCosPAngleCut << " passed all cuts " << nPassed << " ratio " << nPassed / (Float_t)i << " nConfirmedGammas: " << nConfirmedGammas << " purity: " << nConfirmedGammas / (Float_t)nPassed <<  endl;
}

// the same but fill the histos for all V0s that are confirmed gammas
void BTask::process(aod::Collision const& collision, 
                    aod::V0Datas const& v0s, 
                    tracksAndTPCInfo const& tracks,
                    aod::McParticles const& mcParticles)
{    
    for (auto& v0 : v0s) {
        
        auto trackPos = v0.template posTrack_as<tracksAndTPCInfo>(); //positive daughter
        auto trackNeg = v0.template negTrack_as<tracksAndTPCInfo>(); //positive daughter
        
        // now check if we have an actual photon
        // todo: verify it is enough to check only mother0 being equal
        if (trackPos.mcParticle().mother0() > -1 && 
            trackPos.mcParticle().mother0() == trackNeg.mcParticle().mother0())
        {
            auto mother = mcParticles.iteratorAt(trackPos.mcParticle().mother0());
                        
            if (mother.pdgCode() == 22) 
            {
                registry.fill(HIST("hPtRec"), v0.pt());
                registry.fill(HIST("hPtRes"), v0.pt() - mother.pt());
                
                TVector3 lConvPointRec(v0.x(), v0.y(), v0.z());
                TVector3 lPosTrackVtxMC(trackPos.mcParticle().vx(), trackPos.mcParticle().vy(), trackPos.mcParticle().vz());
                // take the origin of the positive mc track as conversion point (should be identical with negative, verified this on a few photons)
                TVector3 lConvPointMC(lPosTrackVtxMC);
                
                registry.fill(HIST("hConvPointR"), lConvPointRec.Perp());
                registry.fill(HIST("hConvPointRRes"), lConvPointRec.Perp() - lConvPointMC.Perp());
                registry.fill(HIST("hConvPointAbsoluteDistanceRes"), TVector3(lConvPointRec - lConvPointMC).Mag());
                
                float lPhiV0Rec = static_cast<float>(M_PI) + std::atan2(-v0.py(), -v0.px());
                registry.fill(HIST("hEtaRes"), v0.eta() - mother.eta()); 
                registry.fill(HIST("hPhiRes"), lPhiV0Rec - mother.phi()); 
            }
        }
    }
}


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{ 
    return WorkflowSpec{
        adaptAnalysisTask<ATask>(cfgc, TaskName{"myGammaTask3"}),
        adaptAnalysisTask<BTask>(cfgc, TaskName{"myGammaTask3_resolutionsFromAllGammas"})
    };
}
