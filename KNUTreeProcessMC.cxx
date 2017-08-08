#include "KNUTree.h"

#include <iostream>
#include <cstring>

#include "FileManager.hh"
#include "AMSRootSupport.hh"
#include "AMSRootParticleFactory.hh"
#include "AMSGeometry.h"
#include "Utilities.hh"
#include "SlowControlLookup.hh"
#include "ParticleId.hh"
#include "TrdKCluster.h"
#include "amschain.h"
#include "bcorr.h"
#include "TrCharge.h"
#include "TrSim.h"
#include "TrExtAlignDB.h"
#include "TVector3.h"
#include "SplineTrack.hh"

/**
 * @brief Process events
 * @return
 */
void KNUTree::ProcessMC()
{
  const bool debugMode = false;
  KNUOUT << "Begin to loop over " << nEvents << " events. [MC Mode]" << endl;
  KNUOUT << "Initializing variables ... " << endl;
  KNUTree::InitLoop();

  // Disabling reading settings from file
  TRMCFFKEY.ReadFromFile  = 0;
  TRFITFFKEY.ReadFromFile = 0;
  TRFITFFKEY.magtemp      = 0;

  const bool setAMSRootDefaults = true;
  AMSRootSupport amsRootSupport(AC::MCRun, setAMSRootDefaults);

  Analysis::EventFactory&              eventFactory = amsRootSupport.EventFactory();
  Analysis::AMSRootParticleFactory& particleFactory = amsRootSupport.ParticleFactory();

  //for( unsigned long int e = 0; e < 10; e++ )
  for( unsigned long int e = 0; e < nEvents; e++ )
  {
    // Print out the working status briefly.
    if( (int)e % 100 == 0 || e == nEvents - 1 )
    {
      KNUOUT << "Processing " << e << " out of " << nEvents << " (" << (float)e/nEvents*100. << "%), "
        << "nSelected = " << nSelected << ", currently, selection efficiency is (" << (float)nSelected/e*100. << "%) " << endl;
        //<< "nSelectedParticles = " << nSelectedParticles << ", out of nExaminedParticles = " << nExaminedParticles<< endl;
    }

    // Get event
    pAMSEvent = pChain->GetEvent(e);
    pHeader   = &( pAMSEvent->fHeader );

    if( !pAMSEvent ){ if( debugMode ) { KNUOUT << "pAMSEvent is NULL!" << endl; } continue; }
    hEvtCounter->Fill(0);
    if( pAMSEvent->nParticle() != 1 ) continue;
    hEvtCounter->Fill(1);
    pParticle = pAMSEvent->pParticle(0);
    if( !IsGoodParticle(pParticle) ) continue;
    hEvtCounter->Fill(2);

    nSelected++;
    // End of Cut

    // Save variables
    nRun            = pAMSEvent->Run();
    nEvent          = pAMSEvent->Event();
    nLevel1         = pAMSEvent->nLevel1();
    nParticle       = pAMSEvent->nParticle();
    nCharge         = pAMSEvent->nCharge();
    nTrTrack        = pAMSEvent->nTrTrack();
    nTrdTrack       = pAMSEvent->nTrdTrack();
    nAntiCluster    = pAMSEvent->nAntiCluster();
    nRichRing       = pAMSEvent->nRichRing();
    nRichRingB      = pAMSEvent->nRichRingB();
    nBeta           = pAMSEvent->nBeta();
    nBetaB          = pAMSEvent->nBetaB();
    nBetaH          = pAMSEvent->nBetaH();
    nEcalShower     = pAMSEvent->nEcalShower();
    nVertex         = pAMSEvent->nVertex();

    nAntiMCCluster  = pAMSEvent->nAntiMCCluster();
    nTrMCCluster    = pAMSEvent->nTrMCCluster();
    nTofMCCluster   = pAMSEvent->nTofMCCluster();
    nTofMCPmtHit    = pAMSEvent->nTofMCPmtHit();
    nEcalMCHit      = pAMSEvent->nEcalMCHit();
    nTrdMCCluster   = pAMSEvent->nTrdMCCluster();
    nRichMCCluster  = pAMSEvent->nRichMCCluster();
    nMCTrack        = pAMSEvent->nMCTrack();
    nMCEventg       = pAMSEvent->nMCEventg();

    // Object preparation
    pTrTrack = pParticle->pTrTrack();
    int id_inner = pTrTrack->iTrTrackPar(1, 3, 23);
    pBetaH = pParticle->pBetaH();
    pEcalShower = pParticle->pEcalShower();
    MCEventgR* pMCEventg = pAMSEvent->GetPrimaryMC();
    particleID = pMCEventg->Particle;
    trkID = pMCEventg->trkID;
    parentID = pMCEventg->parentID;
    for( int j = 0; j < 3; ++j )
    {
      genPosition[j] = pMCEventg->Coo[j];
      genDirCosine[j] = pMCEventg->Dir[j];
    }
    genMass = pMCEventg->Mass;
    genCharge = pMCEventg->Charge;
    genMomentum = pMCEventg->Momentum;

    // ACSoft related lines
    particleFactory.SetAMSTrTrackR(pTrTrack);
    particleFactory.SetAMSBetaHR(pBetaH);
    particleFactory.SetAMSEcalShowerR(pEcalShower);

    if( !amsRootSupport.SwitchToSpecificTrackFitById(id_inner) ) continue;
    Analysis::Event& event = amsRootSupport.BuildEvent(pChain, pAMSEvent);

    // Only do this if you need access to TRD segments/tracks and vertices
    eventFactory.PerformTrdTracking(event);
    eventFactory.PerformTrdVertexFinding(event);

    // If you want to use TrdQt, you shold always add Analysis::CreateSplineTrack, so it can use the
    // Tracker track extrapolation to pick up the TRD hits in the tubes and calculate the path length
    // and also Analysis::FillTrdQt so that the likelihoods are calculated.
    // Additionally you should pass Analysis::CreateTrdTrack if you want to examine the TRD tracks
    // found by the previous PerformTrdTracking() call, via the Analysis::Particle interface.

    int productionSteps = 0;
    productionSteps |= Analysis::CreateSplineTrack;
    productionSteps |= Analysis::CreateTrdTrack;
    //productionSteps |= Analysis::FillTrdQt;
    eventFactory.FillParticles(event, productionSteps);

    const Analysis::Particle* acParticle = event.PrimaryParticle();
    assert(acParticle);

    // Tracker
    if( id_inner < 0 ) cout << "ERROR!" << endl;
    // Inner tracker
    trkFitCodeInner              = id_inner;
    trkRigidityInner             = pTrTrack->GetRigidity(id_inner);
    trkRigidityInverseErrorInner = pTrTrack->GetErrRinv(id_inner);
    trkReducedChisquareInnerX    = pTrTrack->GetNormChisqX(id_inner);
    trkReducedChisquareInnerY    = pTrTrack->GetNormChisqY(id_inner);

    // Upper half of inner tracker
    trkFitCodeUpperInner = pTrTrack->iTrTrackPar(1, 1, 23);
    if( trkFitCodeUpperInner > 0 )
    {
      trkRigidityUpperInner             = pTrTrack->GetRigidity(trkFitCodeUpperInner);
      trkRigidityInverseErrorUpperInner = pTrTrack->GetErrRinv(trkFitCodeUpperInner);
      trkReducedChisquareUpperInnerX    = pTrTrack->GetNormChisqX(trkFitCodeUpperInner);
      trkReducedChisquareUpperInnerY    = pTrTrack->GetNormChisqY(trkFitCodeUpperInner);
    }

    // Lower half of inner tracker
    trkFitCodeLowerInner = pTrTrack->iTrTrackPar(1, 2, 23);
    if( trkFitCodeLowerInner > 0 )
    {
      trkRigidityLowerInner             = pTrTrack->GetRigidity(trkFitCodeLowerInner);
      trkRigidityInverseErrorLowerInner = pTrTrack->GetErrRinv(trkFitCodeLowerInner);
      trkReducedChisquareLowerInnerX    = pTrTrack->GetNormChisqX(trkFitCodeLowerInner);
      trkReducedChisquareLowerInnerY    = pTrTrack->GetNormChisqY(trkFitCodeLowerInner);
    }

    // L1 + inner tracker
    trkFitCodeL1Inner = pTrTrack->iTrTrackPar(1, 5, 20);
    if( trkFitCodeL1Inner > 0 )
    {
      trkRigidityL1Inner             = pTrTrack->GetRigidity(trkFitCodeL1Inner);
      trkRigidityInverseErrorL1Inner = pTrTrack->GetErrRinv(trkFitCodeL1Inner);
      trkReducedChisquareL1InnerX    = pTrTrack->GetNormChisqX(trkFitCodeL1Inner);
      trkReducedChisquareL1InnerY    = pTrTrack->GetNormChisqY(trkFitCodeL1Inner);
    }

    // L9 + inner tracker
    trkFitCodeL9Inner = pTrTrack->iTrTrackPar(1, 6, 20);
    if( trkFitCodeL9Inner > 0 )
    {
      trkRigidityL9Inner             = pTrTrack->GetRigidity(trkFitCodeL9Inner);
      trkRigidityInverseErrorL9Inner = pTrTrack->GetErrRinv(trkFitCodeL9Inner);
      trkReducedChisquareL9InnerX    = pTrTrack->GetNormChisqX(trkFitCodeL9Inner);
      trkReducedChisquareL9InnerY    = pTrTrack->GetNormChisqY(trkFitCodeL9Inner);
    }

    // Full span tracker
    trkFitCodeFS = pTrTrack->iTrTrackPar(1, 7, 20);
    if( trkFitCodeFS > 0 )
    {
      trkRigidityFS             = pTrTrack->GetRigidity(trkFitCodeFS);
      trkRigidityInverseErrorFS = pTrTrack->GetErrRinv(trkFitCodeFS);
      trkReducedChisquareFSX    = pTrTrack->GetNormChisqX(trkFitCodeFS);
      trkReducedChisquareFSY    = pTrTrack->GetNormChisqY(trkFitCodeFS);
    }

    // Maximum span tracker
    trkFitCodeMS = pTrTrack->iTrTrackPar(1, 0, 20);
    if( trkFitCodeMS > 0 )
    {
      trkRigidityMS             = pTrTrack->GetRigidity(trkFitCodeMS);
      trkRigidityInverseErrorMS = pTrTrack->GetErrRinv(trkFitCodeMS);
      trkReducedChisquareMSX    = pTrTrack->GetNormChisqX(trkFitCodeMS);
      trkReducedChisquareMSY    = pTrTrack->GetNormChisqY(trkFitCodeMS);
    }

    for(int ii = 0; ii < 9; ii++)
    {
      trkLayerJQ[ii] = 0.;
      trkEdepLayerJ[ii] = 0.;
    }

    AMSPoint amspnt;
    AMSDir amsdir;
    for( int ilayer = 0; ilayer < 9; ilayer++)
    {
      trkLayerJQ[ilayer] = pTrTrack->GetLayerJQ(ilayer+1);
      pTrRecHit = pTrTrack->GetHitLJ(ilayer+1);
      if( !pTrRecHit )
      {
        trkEdepLayerJXSideOK[ilayer] = -1;
        trkEdepLayerJYSideOK[ilayer] = -1;
      }
      else
      {
        if( pTrRecHit->GetEdep(0) != 0 ) trkEdepLayerJXSideOK[ilayer] = 1;
        else trkEdepLayerJXSideOK[ilayer] = 0;
        if( pTrRecHit->GetEdep(1) != 0 ) trkEdepLayerJYSideOK[ilayer] = 1;
        else trkEdepLayerJYSideOK[ilayer] = 0;
        trkEdepLayerJ[ilayer] += pTrRecHit->GetEdep(0) + pTrRecHit->GetEdep(1);
      }
      pTrTrack->InterpolateLayerJ(ilayer+1, amspnt, amsdir);
      trkPosXLJ[ilayer]     = amspnt.x();
      trkPosYLJ[ilayer]     = amspnt.y();
      trkPosZLJ[ilayer]     = amspnt.z();
      trkDirThetaLJ[ilayer] = amsdir.gettheta();
      trkDirPhiLJ[ilayer]   = amsdir.getphi();
    }

    trkCharge = pTrTrack->GetQ();
    std::vector<like_t> like;
    if( fabs(trkRigidityInner) < 3 && tofBeta < 1)
    {
      trkZ = pTrTrack->GetZ(like, tofBeta);
      trkInnerZ = pTrTrack->GetInnerZ(like, tofBeta);
      for(int k = 0; k < 9; k++) trkLayerJZ[k] = pTrTrack->GetInnerZ(like, k, tofBeta);
    }
    else
    {
      trkZ = pTrTrack->GetZ(like);
      trkInnerZ = pTrTrack->GetInnerZ(like);
      for(int k = 0; k < 9; k++) trkLayerJZ[k] = pTrTrack->GetInnerZ(like, k);
    }
    trkInnerCharge = pTrTrack->GetInnerQ();
    trkHasExtLayers = pTrTrack->HasExtLayers();

    // 
    //
    // TOF Section
    //
    //
    tofNCluster         = pAMSEvent->nTofCluster();
    tofNClusterH        = pAMSEvent->nTofClusterH();
    if(!pBetaH) { KNUERR << "Event: " << e << " NULL pBetaH pointer! This should not be happen! Exit the program!" << std::endl; exit(1); }
    tofBeta             = pBetaH->GetBeta();
    tofInvBetaErr       = pBetaH->GetEBetaV();
    tofNormEBetaV       = pBetaH->GetNormEBetaV();
    tofBetaS            = pBetaH->GetBetaS();
    tofBetaC            = pBetaH->GetBetaC();
    tofEBetaCV          = pBetaH->GetEBetaCV();
    tofMass             = pBetaH->GetMass();
    tofMassError        = pBetaH->GetEMass();
    tofNUsedHits        = pBetaH->NTofClusterH();
    tofNUnusedHits      = tofNCluster - tofNUsedHits;
    if( pBetaH->IsGoodBeta()   == true ) isGoodBeta = 1;
    else isGoodBeta     = 0;
    if( pBetaH->IsTkTofMatch() == true ) isTkTofMatch = 1;
    else isTkTofMatch   = 0;
    tofChisqT           = pBetaH->GetChi2T();
    tofChisqC           = pBetaH->GetChi2C();
    tofReducedChisqT    = pBetaH->GetNormChi2T();
    tofReducedChisqC    = pBetaH->GetNormChi2C();
    tofSumHit           = pBetaH->GetSumHit();
    tofNUsedClusterH    = pBetaH->GetUseHit();
    for( int k = 0; k < 4; k++ )
    {
      tofPatternOnLayer[k] = pBetaH->GetPattern(k);
      tofTimeOnLayer[k] = pBetaH->GetTime(k);
      tofETimeOnLayer[k] = pBetaH->GetETime(k);
      tofETCooOnLayer[k] = pBetaH->GetETCoo(k);
      tofTResidualOnLayer[k] = pBetaH->GetTResidual(k);
      tofTkTFLenOnLayer[k] = pBetaH->GetTkTFLen(k);
      tofCResidualXOnLayer[k] = pBetaH->GetCResidual(k, 0);
      tofCResidualYOnLayer[k] = pBetaH->GetCResidual(k, 1);
    }
    tofCharge = pBetaH->GetQ(tofNUsedLayersForQ, tofChargeRMS);
    tofChargeOnLayer[0] = pBetaH->GetQL(0);
    tofChargeOnLayer[1] = pBetaH->GetQL(1);
    tofChargeOnLayer[2] = pBetaH->GetQL(2);
    tofChargeOnLayer[3] = pBetaH->GetQL(3);
    tofZ = pBetaH->GetZ(nTofLUsedForZ, probTOFZ);

    for(int i = 0; i < 4; ++i)
    {
      TofClusterHR* pTofClusterH = pBetaH->GetClusterHL(i);
      if( pTofClusterH != 0 )
      {
        tofDepositedEnergyOnLayer[i] = pTofClusterH->GetEdep();
      }
      else
        tofDepositedEnergyOnLayer[i] = -1;
    }

    int ncls[4] = {0, 0, 0, 0};
    nTofClustersInTime  = pAMSEvent->GetNTofClustersInTime(pBetaH, ncls);

    tofAcCharge      = acParticle->TofCharge();
    tofAcUpperCharge = acParticle->UpperTofCharge();
    tofAcLowerCharge = acParticle->LowerTofCharge();

    // TOF MC section
    vector<TofMCClusterR>* pvTofMCCluster;
    vector<int> score;
    score.clear();
    bool survivedAtL[4] = { false, false, false, false };

    /*
    cout << "Event : " << e << " | Particle : " << Form("%+3d", particleID) <<
      " | GtrkID : " << Form("%+5d", trkID) <<
      " | ParentNo : " << Form("%+5d", parentID) << endl;
      */
    if( nTofMCCluster != 0 )
    {
      //for( std::vector<TofMCClusterR>::iterator it = pvTofMCCluster->begin(); it !- pvTofMCCluster->end(); ++it )
      for( Int_t i = 0; i < nTofMCCluster; ++i )
      {
        TofMCClusterR* pTofMCCluster = pAMSEvent->pTofMCCluster(i);
        switch(pTofMCCluster->GetLayer())
        {
          case 0:
            if( 
              particleID == pTofMCCluster->Particle &&
              trkID == pTofMCCluster->GtrkID &&
              parentID == pTofMCCluster->ParentNo )
              survivedAtL[0] = true;
            break;
          case 1:
            if(
              particleID == pTofMCCluster->Particle &&
              trkID == pTofMCCluster->GtrkID &&
              parentID == pTofMCCluster->ParentNo )
              survivedAtL[1] = true;
            break;
          case 2:
            if(
              particleID == pTofMCCluster->Particle &&
              trkID == pTofMCCluster->GtrkID &&
              parentID == pTofMCCluster->ParentNo )
              survivedAtL[2] = true;
            break;
          case 3:
            if(
              particleID == pTofMCCluster->Particle &&
              trkID == pTofMCCluster->GtrkID &&
              parentID == pTofMCCluster->ParentNo )
              survivedAtL[3] = true;
            break;
          default:
            break;
        }

        /*    Information print-out for debugging
        cout << "TofMCCluster---| i : " << Form("%2d",i) <<
          " | Layer : "                 << pTofMCCluster->GetLayer() <<
          " | Bar : "                   << Form("%3d", pTofMCCluster->GetBar()) <<
          " | Idsoft : "                << pTofMCCluster->Idsoft <<
          " | ParentNo : "              << Form("%+5d", pTofMCCluster->ParentNo) <<
          " | GtrkID : "                << Form("%+5d", pTofMCCluster->GtrkID) <<
          " | Particle : "              << Form("%+3d", pTofMCCluster->Particle) <<
          " | Coo_X(cm) : "             << Form("%+7.4f", pTofMCCluster->Coo[0]) <<
          " | Coo_Y(cm) : "             << Form("%+7.4f", pTofMCCluster->Coo[1]) <<
          " | Coo_Z(cm) : "             << Form("%+7.4f", pTofMCCluster->Coo[2]) <<
          " | TOF(sec) : "              << Form("%7.4f", pTofMCCluster->TOF) <<
          " | Beta : "                  << Form("%7.4f", pTofMCCluster->Beta) <<
          " | Edep(GeV) : "             << Form("%7.4f", pTofMCCluster->Edep) <<
          " | EdepR(GeV) : "            << Form("%7.4f", pTofMCCluster->EdepR) <<
          " | Step(cm) : "              << Form("%7.4f", pTofMCCluster->Step) << endl;
          */
      }
      for( int i = 0; i < pBetaH->NTofClusterH(); ++i )
      {
        TofClusterHR* pTofClusterH = pBetaH->GetClusterHL(i);
        /*    Information print-out for debugging
        cout << "TofClusterH----| i : " << Form("%2d",i) <<
          " | Layer : "                 << pTofClusterH->Layer <<
          " | Bar : "                   << Form("%3d", pTofClusterH->Bar) <<
          " | Coo_X(cm) : "             << Form("%+7.4f", pTofClusterH->Coo[0]) <<
          " | Coo_Y(cm) : "             << Form("%+7.4f", pTofClusterH->Coo[1]) <<
          " | Coo_Z(cm) : "             << Form("%+7.4f", pTofClusterH->Coo[2]) <<
          " | AEdep : "                 << Form("%+7.4f", pTofClusterH->AEdep) <<
          " | DEdep : "                 << Form("%+7.4f", pTofClusterH->DEdep) << endl;
          */
      }
    }
    if( survivedAtL[0] && survivedAtL[1] && survivedAtL[2] && survivedAtL[3] )
    {
      //cout << "==============Event: " << e << ", SURVIVED!!================" << endl;
      survived_tof = 1;
    }
    else survived_tof = 0;

    //
    //
    // ECAL shower Section
    //
    //

    // Save ECAL shower related variables in case of EcalShowerR object exists.
    if( pEcalShower )
    {
      isEcalAvailable = 1;
      showerEnergyCorrected = pEcalShower->GetCorrectedEnergy(2, 2);
      if( showerEnergyCorrected < 0.5 ) continue;
      showerEnergyD = (pEcalShower->EnergyD)/1000.;
      for(int i = 0; i < pEcalShower->NEcal2DCluster(); ++i)
      {
        for(int j = 0; j < pEcalShower->pEcal2DCluster(i)->NEcalCluster(); ++j)
        {
          if (pEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Edep > 0)
            showerEnergyDL[pEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Plane] += pEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Edep;
          else
            showerEnergyDL[pEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Plane] = -1;
        }
      }
      showerEnergyE = pEcalShower->EnergyE;
      showerBDT     = pEcalShower->GetEcalBDT();
      showerCofG[0] = pEcalShower->CofG[0];
      showerCofG[1] = pEcalShower->CofG[1];
      showerCofG[2] = pEcalShower->CofG[2];
      AMSPoint trackPoint;
      AMSDir   trackDirection;
      pTrTrack->Interpolate(pEcalShower->CofG[2], trackPoint, trackDirection, id_inner);
      pEcalShower->NormaliseVariableLAPP();
      showerCofGDist = std::sqrt(std::pow(trackPoint[0] - showerCofG[0], 2) + std::pow( trackPoint[1] - showerCofG[1], 2));
      showerCofGdX   = std::fabs(trackPoint[0] - showerCofG[0]);
      showerCofGdY   = std::fabs(trackPoint[1] - showerCofG[1]);
    }

    trdPElectronToProtonLogLikelihoodRatio = acParticle->CalculateElectronProtonLikelihood();
    trdPHeliumToProtonLogLikelihoodRatio   = acParticle->CalculateHeliumProtonLikelihood();
    trdPHeliumToElectronLogLikelihoodRatio = acParticle->CalculateHeliumElectronLikelihood();

    // query Xenon partial pressure
    bool pXeOk = false;
    acEventTime = event.EventTime();
    double pXe = Utilities::SlowControlLookup::Self()->QueryXenonPressure(acEventTime, pXeOk);

    // prepare for likelihood calculation
    const Analysis::Particle::TrdHitsVector& trdHybridHits = acParticle->TrdHitsFromTrdAndTrackerTracks();
    Analysis::TrdLikelihoodCalculation lh(trdHybridHits, pXe);

    // calculate electron/proton likelihood and check that pdfs could be retrieved as expected
    trdPElectronLikelihood = lh.ComputeTrdLikelihood(ParticleId::Electron, acParticle->Rigidity());
    trdPProtonLikelihood   = lh.ComputeTrdLikelihood(ParticleId::Proton, acParticle->Rigidity());
    trdPHeliumLikelihood   = lh.ComputeTrdLikelihood(ParticleId::Helium, acParticle->Rigidity());
    bool TrdPdfQueriesOk   = lh.PdfQueriesOk();

    // number of total TRD hits available
    trdPNHybridHits = acParticle->NumberOfTrdHitsFromTrdAndTrackerTracks();

    // number of active layers used for likelihood generation
    trdPNActiveLayers = lh.NumberOfActiveLayers();

    // number of unsuitable hits filtered out before the likelihood calculation
    trdPNLowAdcHits = lh.NumberOfLowAdcHits();
    trdPNLowDxHits = lh.NumberOfLowDxHits();

    acTrackerTrdTrackDistanceAtTrdCenterInSigmas = acParticle->TrackerTrdTrackDistanceAtTrdCenterInSigmas();
    acTrackerTrdTrackAngularDistanceInSigmas     = acParticle->TrackerTrdTrackAngularDistanceInSigmas();

    pRichRing = pParticle->pRichRing();
    if( pRichRing )
    {
      isRichAvailable            = 1;
      richRebuild                =   (int)pRichRing->Rebuild();
      richIsGood                 =   (int)pRichRing->IsGood();
      richIsClean                =   (int)pRichRing->IsClean();
      richIsNaF                  =   (int)pRichRing->IsNaF();
      richTileIndex              =   (int)pRichRing->getTileIndex();
      richUsedHits               =   (int)pRichRing->getUsedHits();
      richRingWidth              = (float)pRichRing->RingWidth();
      richDistanceTileBorder     = (float)pRichRing->DistanceTileBorder();
      richChargeCorrections      =   (int)pRichRing->buildChargeCorrections();
      richPMTCorrectionsFailed   =   (int)pRichRing->PmtCorrectionsFailed();
      richWidth                  =        pRichRing->getWidth();
      richNHits                  =        pRichRing->getHits();
      richNPMTsOnRing            =        pRichRing->getPMTs();
      richBeta                   =        pRichRing->getBeta();
      richBetaError              =        pRichRing->getBetaError();
      richChargeSquared          =        pRichRing->getCharge2Estimate();
      richKolmogorovProbability  =        pRichRing->getProb();
      richPhotoelectrons         =        pRichRing->getPhotoElectrons();
      richExpectedPhotoelectrons =        pRichRing->getExpectedPhotoelectrons();
      const float* tmpRichTrackEmissionPoint = pRichRing->getTrackEmissionPoint();
      richTrackEmissionPoint[0]  = tmpRichTrackEmissionPoint[0];  // X
      richTrackEmissionPoint[1]  = tmpRichTrackEmissionPoint[1];  // Y
      richTrackEmissionPoint[2]  = tmpRichTrackEmissionPoint[2];  // Z
      richTrackEmissionPoint[3]  = tmpRichTrackEmissionPoint[3];  // Theta
      richTrackEmissionPoint[4]  = tmpRichTrackEmissionPoint[4];  // Phi

      AMSPoint trackPointRICH;
      AMSDir   trackDirAtRICH;
      float zrichradiatortop = AC::AMSGeometry::ZRICHradiator;
      pTrTrack->Interpolate(zrichradiatortop, trackPointRICH, trackDirAtRICH, id_inner);
      richTrackerTrackOnRadiatorX= trackPointRICH[0];
      richTrackerTrackOnRadiatorY= trackPointRICH[1];
      /*
      cout << "RICH: " << richTrackEmissionPoint[0] << endl;
      cout << "RICH: " << richTrackEmissionPoint[1] << endl;
      cout << "RICH: " << richTrackEmissionPoint[2] << endl;
      cout << "RICH: " << richTrackEmissionPoint[3] << endl;
      cout << "RICH: " << richTrackEmissionPoint[4] << endl;
      cout << "RICH: " << sizeof(richTrackEmissionPoint) << endl;
      cout << "RICH: ====================================" << endl;
      */
      richTheta                  =        pRichRing->getTrackTheta();
      richPhi                    =        pRichRing->getTrackPhi();
      richBetaConsistency        =        pRichRing->getBetaConsistency();
      richReflectedHits          =        pRichRing->getReflectedHits();
      richPMTChargeConsistency   =        pRichRing->getPMTChargeConsistency();
    }

    // TRD vertex cut
    trdNVertex = 0;
    const std::vector<Analysis::TrdVertex>& verticesXZ = event.TrdVerticesXZ();
    const std::vector<Analysis::TrdVertex>& verticesYZ = event.TrdVerticesYZ();

    for( std::vector<Analysis::TrdVertex>::const_iterator xzIter = verticesXZ.begin(); xzIter != verticesXZ.end(); ++xzIter)
    {
      const Analysis::TrdVertex& xzVertex = *xzIter;
      for( std::vector<Analysis::TrdVertex>::const_iterator yzIter = verticesYZ.begin(); yzIter != verticesYZ.end(); ++yzIter)
      {
        const Analysis::TrdVertex& yzVertex = *yzIter;
        if( std::max( xzVertex.NumberOfSegments(), yzVertex.NumberOfSegments() ) < 3 )
          continue;
        if( fabs( xzVertex.Z() - yzVertex.Z() ) < fabs( xzVertex.ErrorZ() + yzVertex.ErrorZ() ) )
        {
          trdNVertex++;
        }
      }
    }
    pTrdTrack = pParticle->pTrdTrack();
    int trdNUsedHits = 0;
    int trdNUsedSegment = 0;
    if(!pTrdTrack)
    {
      trdNClusters      = -1;
      trdNUsedHits      = -1;
      trdNUsedSegment   = -1;
      trdNUsedHits      = -1;
      trdNTracks        = -1;
      trdTrackTheta     = -9.;
      trdTrackPhi       = -9.;
      trdTrackPattern   = -9;
      trdTrackCharge    = -9;
      for(int i = 0; i < 20; ++i) trdTrackEdepL[i] = -9.;
    }
    else
    {
      trdNClusters = pAMSEvent->nTrdCluster();
      trdNUsedSegment = pTrdTrack->NTrdSegment();
      for(int i = 0; i < trdNUsedSegment; i++)
      {
        pTrdSegment = pTrdTrack->pTrdSegment(i);
        trdNUsedHits += pTrdSegment->NTrdCluster();
      }
      trdNUnusedHits = trdNClusters - trdNUsedHits;

      trdNTracks  = pAMSEvent->nTrdTrack();
      trdTrackTheta   = pTrdTrack->Theta;
      trdTrackPhi     = pTrdTrack->Phi;
      trdTrackPattern = pTrdTrack->Pattern;
      trdTrackCharge  = pTrdTrack->Q;
      trdTrackChi2    = pTrdTrack->Chi2;
      trdTrackMeanDepositedEnergy = 0.;
      int ntrdlayers = 0;
      for(int i = 0; i < pTrdTrack->NTrdSegment(); i++)
      {
        for(int j = 0; j < pTrdTrack->pTrdSegment(i)->NTrdCluster(); j++)
        {
          int trdLayer  = pTrdTrack->pTrdSegment(i)->pTrdCluster(j)->Layer;
          float trdEdep = pTrdTrack->pTrdSegment(i)->pTrdCluster(j)->EDep;
          trdTrackEdepL[trdLayer] = trdEdep;
          trdTrackTotalDepositedEnergy += trdEdep;
          ntrdlayers++;
        }
      }
      trdTrackMeanDepositedEnergy = trdTrackTotalDepositedEnergy/(float)ntrdlayers;

      // Check whether TrkInner meet with TrdTrack at UTOF
      AMSPoint trackPointAtUTOF;
      AMSDir   trackDirAtUTOF;
      float ztofupper = AC::AMSGeometry::ZTOFUpper;
      pTrTrack->Interpolate(ztofupper, trackPointAtUTOF, trackDirAtUTOF, id_inner);
      AMSPoint trdTrackPointAtUTOF;
      AMSDir   trdTrackDirAtUTOF;
      pTrdTrack->Interpolate(ztofupper, trdTrackPointAtUTOF, trdTrackDirAtUTOF);
      trdTrackDeviationXWithInnerTrk = trackPointAtUTOF[0] - trdTrackPointAtUTOF[0];
      trdTrackDeviationYWithInnerTrk = trackPointAtUTOF[1] - trdTrackPointAtUTOF[1];
      trdTrackDistanceBetweenInnerTrk = std::sqrt( std::pow( trdTrackDeviationXWithInnerTrk, 2) + std::pow( trdTrackDeviationYWithInnerTrk, 2) );

      const Analysis::SplineTrack* pSplTrack = acParticle->GetSplineTrack();
      TVector3 v = pSplTrack->InterpolateToZ(AC::AMSGeometry::ZTOFUpper);
      ResidualXBetweenInnerTrackAndSplineTrack = v.X() - trackPointAtUTOF[0];
      ResidualYBetweenInnerTrackAndSplineTrack = v.Y() - trackPointAtUTOF[1];
    }

    // Below of this for TrdK

    trdKNRawHits = pAMSEvent->NTrdRawHit();
    if(trdKNRawHits > 0)
    {
      double trdKLikelihoodRatio[3] = {-1., -1., -1.};
      double trdKLikelihood[3]      = {-1., -1., -1.};

      pTrdKCluster = new TrdKCluster(pAMSEvent, pTrTrack, id_inner);
      if(!pTrdKCluster) continue;

      trdKIsReadAlignmentOK = pTrdKCluster->IsReadAlignmentOK;
      trdKIsReadCalibOK     = pTrdKCluster->IsReadCalibOK;
      trdKIsValid           = pTrdKCluster->GetLikelihoodRatio_TrTrack(15, trdKLikelihoodRatio, trdKLikelihood, trdKNHits, trdKTotalPathLength, trdKTotalAmp, -1, 0);
      if(trdKIsValid != 0 && trdKNHits != 0)
      {
        pTrdKCluster->CalculateTRDCharge();
        trdKCharge                              = pTrdKCluster->GetTRDCharge();
        trdKChargeError                         = pTrdKCluster->GetTRDChargeError();
        trdKNUsedHitsForCharge                  = pTrdKCluster->GetQTRDHitCollectionNuclei().size();
        trdKElectronToProtonLogLikelihoodRatio  = trdKLikelihoodRatio[0];
        trdKHeliumToElectronLogLikelihoodRatio  = trdKLikelihoodRatio[1];
        trdKHeliumToProtonLogLikelihoodRatio    = trdKLikelihoodRatio[2];
        trdKElectronLikelihood                  = trdKLikelihood[0];
        trdKProtonLikelihood                    = trdKLikelihood[1];
        trdKHeliumLikelihood                    = trdKLikelihood[2];

        AMSPoint trExtraP0;
        AMSDir   trExtraDir;
        pTrdKCluster->GetTrTrackExtrapolation(trExtraP0, trExtraDir);

        for(int l = 0; l < pTrdKCluster->NHits(); l++)
        {
          TrdKHit* pTrdKHit = pTrdKCluster->GetHit(l);
          int tmpL = 0;
          tmpL = pTrdKHit->TRDHit_Layer;
          float tmpAmp = 0;
          tmpAmp = pTrdKHit->TRDHit_Amp;
          float tmpPathLength = 0;
          tmpPathLength = pTrdKHit->Tube_Track_3DLength(&trExtraP0, &trExtraDir);

          if( tmpAmp < 15 || tmpPathLength <= 0 ) continue;
          if( pTrdKHit->IsAligned == 0 || pTrdKHit->IsCalibrated == 0 ) continue;
          trdKAmpLayer[tmpL] += tmpAmp;
        }
      }
      delete pTrdKCluster;
    }


    /*
     *
     *  MC particle tracking related code piece
     *
     *
     */

    conversionIdContainer.clear();

    // ParticleR class coordinate information
    //
    //       %% OLD SCHEME to J SCHEME
    //        J SCHEME(index)
    //        Layer           1[0] 2[1] 3[2] 4[3] 5[4] 6[5] 7[6] 8[7] 9[8]
    //        OLD SCHEME(index)
    //        Layer           8[7] 1[0] 2[1] 3[2] 4[3] 5[4] 6[5] 7[6] 9[8]

    const bool debug_mc = true;
    if( debug_mc == true ) {
      KNUOUT << "Event: " << e << " ========== " << endl;
      KNUOUT << "Coordinate information for each detector" << endl;
      KNUOUT << "TRK L1:     X: " << pParticle->TrCoo[7][0] << ", Y: " << pParticle->TrCoo[7][1] << ", Z: " << pParticle->TrCoo[7][2] << endl;
      KNUOUT << "TRD Top:    X: " << pParticle->TRDCoo[1][0] << ", Y: " << pParticle->TRDCoo[1][1] << ", Z: " << pParticle->TRDCoo[1][2] << endl;
      KNUOUT << "TRD Center: X: " << pParticle->TRDCoo[0][0] << ", Y: " << pParticle->TRDCoo[0][1] << ", Z: " << pParticle->TRDCoo[0][2] << endl;
      KNUOUT << "TOF L1:     X: " << pParticle->TOFCoo[0][0] << ", Y: " << pParticle->TOFCoo[0][1] << ", Z: " << pParticle->TOFCoo[0][2] << endl;
      KNUOUT << "TOF L2:     X: " << pParticle->TOFCoo[1][0] << ", Y: " << pParticle->TOFCoo[1][1] << ", Z: " << pParticle->TOFCoo[1][2] << endl;
      KNUOUT << "TRK L2:     X: " << pParticle->TrCoo[0][0] << ", Y: " << pParticle->TrCoo[0][1] << ", Z: " << pParticle->TrCoo[0][2] << endl;
      KNUOUT << "TRK L3:     X: " << pParticle->TrCoo[1][0] << ", Y: " << pParticle->TrCoo[1][1] << ", Z: " << pParticle->TrCoo[1][2] << endl;
      KNUOUT << "TRK L4:     X: " << pParticle->TrCoo[2][0] << ", Y: " << pParticle->TrCoo[2][1] << ", Z: " << pParticle->TrCoo[2][2] << endl;
      KNUOUT << "TRK L5:     X: " << pParticle->TrCoo[3][0] << ", Y: " << pParticle->TrCoo[3][1] << ", Z: " << pParticle->TrCoo[3][2] << endl;
      KNUOUT << "TRK L6:     X: " << pParticle->TrCoo[4][0] << ", Y: " << pParticle->TrCoo[4][1] << ", Z: " << pParticle->TrCoo[4][2] << endl;
      KNUOUT << "TRK L7:     X: " << pParticle->TrCoo[5][0] << ", Y: " << pParticle->TrCoo[5][1] << ", Z: " << pParticle->TrCoo[5][2] << endl;
      KNUOUT << "TRK L8:     X: " << pParticle->TrCoo[6][0] << ", Y: " << pParticle->TrCoo[6][1] << ", Z: " << pParticle->TrCoo[6][2] << endl;
      KNUOUT << "TOF L3:     X: " << pParticle->TOFCoo[2][0] << ", Y: " << pParticle->TOFCoo[2][1] << ", Z: " << pParticle->TOFCoo[2][2] << endl;
      KNUOUT << "TOF L4:     X: " << pParticle->TOFCoo[3][0] << ", Y: " << pParticle->TOFCoo[3][1] << ", Z: " << pParticle->TOFCoo[3][2] << endl;
      KNUOUT << "RICH Rad.:  X: " << pParticle->RichCoo[0][0] << ", Y: " << pParticle->RichCoo[0][1] << ", Z: " << pParticle->RichCoo[0][2] << endl;
      KNUOUT << "RICH PMT:   X: " << pParticle->RichCoo[1][0] << ", Y: " << pParticle->RichCoo[1][1] << ", Z: " << pParticle->RichCoo[1][2] << endl;
      KNUOUT << "TRK L9:     X: " << pParticle->TrCoo[8][0] << ", Y: " << pParticle->TrCoo[8][1] << ", Z: " << pParticle->TrCoo[8][2] << endl;
      KNUOUT << "ECAL Enter: X: " << pParticle->EcalCoo[0][0] << ", Y: " << pParticle->EcalCoo[0][1] << ", Z: " << pParticle->EcalCoo[0][2] << endl;
      KNUOUT << "ECAL COFG:  X: " << pParticle->EcalCoo[1][0] << ", Y: " << pParticle->EcalCoo[1][1] << ", Z: " << pParticle->EcalCoo[1][2] << endl;
      KNUOUT << "ECAL Exit:  X: " << pParticle->EcalCoo[2][0] << ", Y: " << pParticle->EcalCoo[2][1] << ", Z: " << pParticle->EcalCoo[2][2] << endl << endl;
    }

    double tolerance_distance = 0.1;
    KNUOUT << "MCCluster Info" << endl;
    TrMCClusterR* pTrMCCluster = NULL;

    for( int j = 0; j < nTrMCCluster; j++ )
    {
      pTrMCCluster = pAMSEvent->pTrMCCluster(j);
      if( debug_mc == true ) { KNUOUT << "TrMCCluster(" << j << "): X: " << pTrMCCluster->GetXgl().x() << ", Y: " << pTrMCCluster->GetXgl().y() << ", Z: " << pTrMCCluster->GetXgl().z(); }

      /*
      int iTrCooMap[9] = { 7, 0, 1, 2, 3, 4, 5, 6, 8 };
      for( int iLayer = 0; iLayer < 9; iLayer++ )
      {
        if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[iTrCooMap[iLayer]][2] ) < tolerance_distance )
        {
          double dist = 0;
          double deltax = 0;
          double deltay = 0;
          deltax = pParticle->TrCoo[iTrCooMap[iLayer]][0] - pTrMCCluster->GetXgl().x();
          deltay = pParticle->TrCoo[iTrCooMap[iLayer]][1] - pTrMCCluster->GetXgl().y();
          dist = sqrt( deltax * deltax + deltay * deltay );
          cout << Form(" -> L%d hit, Distance = %f", iLayer, dist);
          if( dist < 1 )
          {
            cout << " --> NEAR hit";
            if( !pTrMCCluster->IsPrimary() ) { cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId); }
            else cout << " survived." << endl;
          }
          else
            cout << " --> FAR hit." << endl;
        }
      }
      */
      // case : L1
      if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[7][2] ) < tolerance_distance )
      {
        if( debug_mc == true ) cout << " -> L1 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[7][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[7][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true ) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true ) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true ) {
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else
            if( debug_mc == true ) cout << " survived. " << endl;
        }
        else
          if( debug_mc == true) cout << " --> FAR hit" << endl;
      }
      // case : L2
      else if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[0][2] ) < tolerance_distance )
      {
        if( debug_mc == true) cout << " -> L2 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[0][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[0][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true) {
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else if( debug_mc == true ) cout << " survived. " << endl;
        }
        else if( debug_mc == true ) cout << " --> FAR hit" << endl;
      }
      // case : L3
      else if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[1][2] ) < tolerance_distance )
      {
        if( debug_mc == true ) cout << " -> L3 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[1][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[1][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true ) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true ) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true ) {
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else
            if( debug_mc == true ) cout << " survived. " << endl;
        }
        else
          if( debug_mc == true ) cout << " --> FAR hit" << endl;
      }
      // case : L4
      else if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[2][2] ) < tolerance_distance )
      {
        if( debug_mc == true ) cout << " -> L4 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[2][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[2][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true ) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true ) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true ) {
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else
            if( debug_mc == true ) cout << " survived. " << endl;
        }
        else
          if( debug_mc == true ) cout << " --> FAR hit" << endl;
      }
      // case : L5
      else if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[3][2] ) < tolerance_distance )
      {
        if( debug_mc == true ) cout << " -> L5 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[3][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[3][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true ) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true ) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true ) { 
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else
            if( debug_mc == true ) cout << " survived. " << endl;
        }
        else
          if( debug_mc == true ) cout << " --> FAR hit" << endl;
      }
      // case : L6
      else if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[4][2] ) < tolerance_distance )
      {
        if( debug_mc == true) cout << " -> L6 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[4][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[4][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true ) {
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else
            if( debug_mc == true ) cout << " survived. " << endl;
        }
        else
          if( debug_mc == true ) cout << " --> FAR hit" << endl;
      }
      // case : L7
      else if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[5][2] ) < tolerance_distance )
      {
        if( debug_mc == true) cout << " -> L7 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[5][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[5][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true ) {
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else
            if( debug_mc == true ) cout << " survived. " << endl;
        }
        else
          if( debug_mc == true ) cout << " --> FAR hit" << endl;
      }
      // case : L8
      else if( fabs( pTrMCCluster->GetXgl().z() - pParticle->TrCoo[6][2] ) < tolerance_distance )
      {
        if( debug_mc == true) cout << " -> L8 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[6][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[6][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true ) {
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else
            if( debug_mc == true ) cout << " survived. " << endl;
        }
        else
          if( debug_mc == true ) cout << " --> FAR hit" << endl;
      }
      // case : L9
      else if( fabs( pTrMCCluster->GetXgl().z() == pParticle->TrCoo[8][2] ) < tolerance_distance )
      {
        if( debug_mc == true) cout << " -> L9 hit, Distance = ";
        double dist = 0;
        double deltax = 0;
        double deltay = 0;
        deltax = pParticle->TrCoo[8][0] - pTrMCCluster->GetXgl().x();
        deltay = pParticle->TrCoo[8][1] - pTrMCCluster->GetXgl().y();
        dist = sqrt( deltax * deltax + deltay * deltay );
        if( debug_mc == true) cout << dist;
        if( dist < 1 )
        {
          if( debug_mc == true) cout << " --> NEAR hit";
          if( !pTrMCCluster->IsPrimary() ){
            if( debug_mc == true ) {
              cout << "  [CONVERTED!] PID: " << pTrMCCluster->GetPart() << endl; conversionId = pTrMCCluster->GetPart(); conversionIdContainer.push_back(conversionId);
            }
          }
          else
            if( debug_mc == true ) cout << " survived. " << endl;
        }
        else
          if( debug_mc == true ) cout << " --> FAR hit" << endl;
      }
      else
      {
        if( debug_mc == true) cout << " can't find layer. " << endl;
      }
    }

    conversionIdSet.clear();
    pair<std::set<int>::iterator, bool> pr;
    for( std::vector<int>::const_iterator it = conversionIdContainer.begin(); it != conversionIdContainer.end(); ++it )
    {
      conversionIdSet.insert(*it);
    }
    for( std::set<int>::const_iterator it = conversionIdSet.begin(); it != conversionIdSet.end(); ++it )
      if( debug_mc == true ) cout << *it << endl;

    RichMCClusterR* pRichMCCluster = NULL;
    bool isPrimaryExistInRICH = false;
    for( int j = 0 ; j < nRichMCCluster; j++ )
    {
      pRichMCCluster = pAMSEvent->pRichMCCluster(j);
      if( pRichMCCluster->Id == pAMSEvent->GetPrimaryMC()->Particle )
      {
        isPrimaryExistInRICH = true;
      }
      if( pRichMCCluster->Id != 50 && pRichMCCluster->Id != pAMSEvent->GetPrimaryMC()->Particle )
      {
        if( debug_mc == true )
        {
          KNUOUT << "Event: " << e << " / iRichMCCluster: " << j << endl;
          KNUOUT << "ID: " << pRichMCCluster->Id << endl;
          KNUOUT << "Primary ID: " << pAMSEvent->GetPrimaryMC()->Particle << endl;
          KNUOUT << "isPrimaryExistInRICH: " << isPrimaryExistInRICH << endl;
        }
      }
    }
    outTree->Fill();
  }
}
