#include "KNUTree.h"

#include <iostream>
#include <cstring>

#include "FileManager.hh"
#include "AMSRootSupport.hh"
#include "AMSRootParticleFactory.hh"
#include "Utilities.hh"
#include "SlowControlLookup.hh"
#include "TrdKCluster.h"
#include "amschain.h"
#include "bcorr.h"
#include "TrSim.h"
#include "TrExtAlignDB.h"

#define KNUOUT std::cout << "[KNUTree::Loop] "
#define KNUERR std::cerr << "[KNUTree::Loop] Error: "

void KNUTree::Loop()
{/*{{{*/
  KNUOUT << "Begin to loop over " << nevt << " events." << endl;
  KNUTree::InitializeLoop();

  // Latest RTI database (pass6)
  AMSSetupR::RTI::UseLatest(6);
  // Latest alignment
  TkDBc::UseFinal();
  // Disabling reading settings from file
  TRMCFFKEY.ReadFromFile = 0;
  TRFITFFKEY.ReadFromFile = 0;
  TRFITFFKEY.magtemp = 0;

  const bool setAMSRootDefaults = true;
  AMSRootSupport amsRootSupport(AC::ISSRun, setAMSRootDefaults);
  Analysis::EventFactory& eventFactory = amsRootSupport.EventFactory();
  Analysis::AMSRootParticleFactory& particleFactory = amsRootSupport.ParticleFactory();

  for( unsigned long int e = 0; e < nevt; e++ )
  //for( unsigned long int e = 0; e < 5000; e++ )
  {
    pAMSEvent = pChain->GetEvent(e);
    pHeader = &(pAMSEvent->fHeader);

    if( !pAMSEvent ) continue;
    hEvtCounter->Fill(0);
    int ptlIdx;
    if( isMC == false )
    {
      if( this->IsBadRun(pAMSEvent) ) this->End(1);
      hEvtCounter->Fill(1);
      if( !this->IsScienceRun(pAMSEvent) ) continue;
      hEvtCounter->Fill(2);
      //if( !this->IsHardwareStatusGood(pAMSEvent) ) continue;
      //hEvtCounter->Fill(3);
      if( this->IsUnbiasedPhysicsTriggerEvent(pAMSEvent) ) continue;
      hEvtCounter->Fill(3);
      if( pAMSEvent->nParticle() != 1 ) continue;
      hEvtCounter->Fill(4);
      if( !this->IsTrkAlignmentGood(pAMSEvent) ) continue;
      hEvtCounter->Fill(5);
      ptlIdx = this->GetGoodParticleIndex(pAMSEvent);
      if( ptlIdx == -1 ) continue;
      hEvtCounter->Fill(6);
    }
    else
    {
      if( this->IsBadRun(pAMSEvent) ) this->End(1);
      hEvtCounter->Fill(1);
      if( this->IsUnbiasedPhysicsTriggerEvent(pAMSEvent) ) continue;
      hEvtCounter->Fill(2);
      if( pAMSEvent->nParticle() != 1 ) continue;
      hEvtCounter->Fill(3);
      if( !this->IsTrkAlignmentGood(pAMSEvent) ) continue;
      hEvtCounter->Fill(4);
      ptlIdx = this->GetGoodParticleIndex(pAMSEvent);
      if( ptlIdx == -1 ) continue;
      hEvtCounter->Fill(5);
    }

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

    if( isMC == false )
    {
      liveTime        = pAMSEvent->LiveTime();
      utcTime         = (float)pHeader->UTCTime(0);
      orbitAltitude   = pHeader->RadS;
      orbitLatitude   = pHeader->ThetaS;
      orbitLongitude  = pHeader->PhiS;
      orbitLatitudeM  = pHeader->ThetaM;
      orbitLongitudeM = pHeader->PhiM;
      velR            = pHeader->VelocityS;
      velTheta        = pHeader->VelTheta;
      velPhi          = pHeader->VelPhi;
      yaw             = pHeader->Yaw;
      pitch           = pHeader->Pitch;
      roll            = pHeader->Roll;
      zenithAngle     = pHeader->Zenith();
      isInSAA         = pAMSEvent->IsInSAA();

      double tmpglong;
      double tmpglat;
      gCoordCalcResult  = pAMSEvent->GetGalCoo(
        gCoordCalcResultBit,
        tmpglong,
        tmpglat);
      gLongitude  = tmpglong;
      gLatitude   = tmpglat;

      double tmpSunPosAzimuth;
      double tmpSunPosElevation;
      sunPosCalcResult   = pHeader->getSunAMS(tmpSunPosAzimuth, tmpSunPosElevation);
      sunPosAzimuth      = tmpSunPosAzimuth;
      sunPosElevation    = tmpSunPosElevation;

      AMSPoint solarArray;
      isInShadow         = pAMSEvent->isInShadow(solarArray, 0);
      solarArrayCoord[0] = solarArray.x();
      solarArrayCoord[1] = solarArray.y();
      solarArrayCoord[2] = solarArray.z();
    }

    pParticle = pAMSEvent->pParticle(ptlIdx);
    if(!pParticle) continue;
    //ptlCharge           = (unsigned int)pParticle->Charge;
    // Unfold below to see the particle coordinates
    /*
    ptlMomentum         = pParticle->Momentum;
    ptlTheta            = pParticle->Theta;
    ptlPhi              = pParticle->Phi;
    ptlCoo[0]           = pParticle->Coo[0];
    ptlCoo[1]           = pParticle->Coo[1];
    ptlCoo[2]           = pParticle->Coo[2];
    ptlTofCoo[0][0]     = pParticle->TOFCoo[0][0];
    ptlTofCoo[0][1]     = pParticle->TOFCoo[0][1];
    ptlTofCoo[0][2]     = pParticle->TOFCoo[0][2];
    ptlTofCoo[1][0]     = pParticle->TOFCoo[1][0];
    ptlTofCoo[1][1]     = pParticle->TOFCoo[1][1];
    ptlTofCoo[1][2]     = pParticle->TOFCoo[1][2];
    ptlTofCoo[2][0]     = pParticle->TOFCoo[2][0];
    ptlTofCoo[2][1]     = pParticle->TOFCoo[2][1];
    ptlTofCoo[2][2]     = pParticle->TOFCoo[2][2];
    ptlTofCoo[3][0]     = pParticle->TOFCoo[3][0];
    ptlTofCoo[3][1]     = pParticle->TOFCoo[3][1];
    ptlTofCoo[3][2]     = pParticle->TOFCoo[3][2];
    ptlTofTrLength[0]   = pParticle->TOFTLength[0];
    ptlTofTrLength[1]   = pParticle->TOFTLength[1];
    ptlTofTrLength[2]   = pParticle->TOFTLength[2];
    ptlTofTrLength[3]   = pParticle->TOFTLength[3];
    ptlAntiCoo[0][0]    = pParticle->AntiCoo[0][0];
    ptlAntiCoo[0][1]    = pParticle->AntiCoo[0][1];
    ptlAntiCoo[0][2]    = pParticle->AntiCoo[0][2];
    ptlAntiCoo[0][3]    = pParticle->AntiCoo[0][3];
    ptlAntiCoo[0][4]    = pParticle->AntiCoo[0][4];
    ptlAntiCoo[1][0]    = pParticle->AntiCoo[1][0];
    ptlAntiCoo[1][1]    = pParticle->AntiCoo[1][1];
    ptlAntiCoo[1][2]    = pParticle->AntiCoo[1][2];
    ptlAntiCoo[1][3]    = pParticle->AntiCoo[1][3];
    ptlAntiCoo[1][4]    = pParticle->AntiCoo[1][4];
    ptlEcalCoo[0][0]    = pParticle->EcalCoo[0][0];
    ptlEcalCoo[0][1]    = pParticle->EcalCoo[0][1];
    ptlEcalCoo[0][2]    = pParticle->EcalCoo[0][2];
    ptlEcalCoo[1][0]    = pParticle->EcalCoo[1][0];
    ptlEcalCoo[1][1]    = pParticle->EcalCoo[1][1];
    ptlEcalCoo[1][2]    = pParticle->EcalCoo[1][2];
    ptlEcalCoo[2][0]    = pParticle->EcalCoo[2][0];
    ptlEcalCoo[2][1]    = pParticle->EcalCoo[2][1];
    ptlEcalCoo[2][2]    = pParticle->EcalCoo[2][2];
    ptlTrCoo[0][0]      = pParticle->TrCoo[0][0];
    ptlTrCoo[0][1]      = pParticle->TrCoo[0][1];
    ptlTrCoo[0][2]      = pParticle->TrCoo[0][2];
    ptlTrCoo[1][0]      = pParticle->TrCoo[1][0];
    ptlTrCoo[1][1]      = pParticle->TrCoo[1][1];
    ptlTrCoo[1][2]      = pParticle->TrCoo[1][2];
    ptlTrCoo[2][0]      = pParticle->TrCoo[2][0];
    ptlTrCoo[2][1]      = pParticle->TrCoo[2][1];
    ptlTrCoo[2][2]      = pParticle->TrCoo[2][2];
    ptlTrCoo[3][0]      = pParticle->TrCoo[3][0];
    ptlTrCoo[3][1]      = pParticle->TrCoo[3][1];
    ptlTrCoo[3][2]      = pParticle->TrCoo[3][2];
    ptlTrCoo[4][0]      = pParticle->TrCoo[4][0];
    ptlTrCoo[4][1]      = pParticle->TrCoo[4][1];
    ptlTrCoo[4][2]      = pParticle->TrCoo[4][2];
    ptlTrCoo[5][0]      = pParticle->TrCoo[5][0];
    ptlTrCoo[5][1]      = pParticle->TrCoo[5][1];
    ptlTrCoo[5][2]      = pParticle->TrCoo[5][2];
    ptlTrCoo[6][0]      = pParticle->TrCoo[6][0];
    ptlTrCoo[6][1]      = pParticle->TrCoo[6][1];
    ptlTrCoo[6][2]      = pParticle->TrCoo[6][2];
    ptlTrCoo[7][0]      = pParticle->TrCoo[7][0];
    ptlTrCoo[7][1]      = pParticle->TrCoo[7][1];
    ptlTrCoo[7][2]      = pParticle->TrCoo[7][2];
    ptlTrCoo[8][0]      = pParticle->TrCoo[8][0];
    ptlTrCoo[8][1]      = pParticle->TrCoo[8][1];
    ptlTrCoo[8][2]      = pParticle->TrCoo[8][2];
    ptlTrdCoo[0][0]     = pParticle->TRDCoo[0][0];
    ptlTrdCoo[0][1]     = pParticle->TRDCoo[0][1];
    ptlTrdCoo[0][2]     = pParticle->TRDCoo[0][2];
    ptlTrdCoo[1][0]     = pParticle->TRDCoo[1][0];
    ptlTrdCoo[1][1]     = pParticle->TRDCoo[1][1];
    ptlTrdCoo[1][2]     = pParticle->TRDCoo[1][2];
    ptlRichCoo[0][0]    = pParticle->RichCoo[0][0];
    ptlRichCoo[0][1]    = pParticle->RichCoo[0][1];
    ptlRichCoo[0][2]    = pParticle->RichCoo[0][2];
    ptlRichCoo[1][0]    = pParticle->RichCoo[1][0];
    ptlRichCoo[1][1]    = pParticle->RichCoo[1][1];
    ptlRichCoo[1][2]    = pParticle->RichCoo[1][2];
    ptlRichPath[0]      = pParticle->RichPath[0];
    ptlRichPath[1]      = pParticle->RichPath[1];
    ptlRichPathBeta[0]  = pParticle->RichPathBeta[0];
    ptlRichPathBeta[1]  = pParticle->RichPathBeta[1];
    ptlRichLength       = pParticle->RichLength;
    ptlRichParticles    = pParticle->RichParticles;
    */
    ptlCutOffStoermer   = pParticle->CutoffS;
    ptlCutOffDipole     = pParticle->Cutoff;

    // 
    //
    // TOF Section
    //
    //
    tofNCluster         = pAMSEvent->nTofCluster();
    tofNClusterH        = pAMSEvent->nTofClusterH();
    pBetaH              = pParticle->pBetaH();
    if(!pBetaH) { KNUERR << "NULL pBetaH pointer! This should not be happen! Exit the program!" << std::endl; exit(1); }
    tofBeta             = pBetaH->GetBeta();
    tofInvBetaErr       = pBetaH->GetEBetaV();
    tofMass             = pBetaH->GetMass();
    tofMassError        = pBetaH->GetEMass();
    tofNUsedHits        = pBetaH->NTofClusterH();
    tofNUnusedHits      = tofNCluster - tofNUsedHits;
    if( pBetaH->IsGoodBeta()   == true ) isGoodBeta = 1;
    else isGoodBeta     = 0;
    if( pBetaH->IsTkTofMatch() == true ) isTkTofMatch = 1;
    else isTkTofMatch   = 0;
    tofReducedChisqT    = pBetaH->GetNormChi2T();
    tofReducedChisqC    = pBetaH->GetNormChi2C();
    pBetaH->GetQ(tofNUsedLayersForQ, tofCharge);
    tofChargeOnLayer[0] = pBetaH->GetQL(0);
    tofChargeOnLayer[1] = pBetaH->GetQL(1);
    tofChargeOnLayer[2] = pBetaH->GetQL(2);
    tofChargeOnLayer[3] = pBetaH->GetQL(3);

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

    //
    //
    // Tracker Section
    //
    //
    pTrTrack = pParticle->pTrTrack();
    if(!pTrTrack) continue;
    trkCharge       = pTrTrack->GetQ();
    trkInnerCharge  = pTrTrack->GetInnerQ();
    trkHasExtLayers = pTrTrack->HasExtLayers();

    int z_inn = floor(trkCharge + .5);
    if( z_inn < 1 ) z_inn = 1;
    float mass = TrFit::Mproton; // Proton
    if( z_inn >= 2) mass = TrFit::Mhelium; // Helium
    if( z_inn >= 3) mass = TrFit::Mproton*2*z_inn; // Ion Z>=3
    int refit = 3;
    // ISS or MC
    bool isiss = (pAMSEvent->nMCEventg() == 0);
    // For ISS, enable ion linearity (charge Z >= 3) and dz correction
    if( isiss )
    {
      TRCLFFKEY.UseNonLinearity = 0;
      TRFITFFKEY.Zshift = 2;
      if( z_inn >= 3)
        TrClusterR::SetLinearityCorrection();
      else
        TrClusterR::UnsetLinearityCorrection();
      refit = 3;
    }
    // MC Smearing and Refit
    else
    {
      TrExtAlignDB::SmearExtAlign();
      TRCLFFKEY.UseSensorAlign = 0;
      TRFITFFKEY.Zshift = -1;
      refit = 3;
    }
    int mfit = pTrTrack->iTrTrackPar(1, 5, 20+refit, mass, z_inn);
    float rigidity = 0;
    if( mfit >= 0 )
    {
      if( isiss ) rigidity = pTrTrack->GetCorrectedRigidity(mfit);
      else rigidity = pTrTrack->GetRigidity(mfit);
      // Magnet temperature correction for ISS pass4
      float bcor = pTrTrack->GetBcorr(mfit); // positive
      if( isiss && bcor == 1 )
      { // Not applied
        float bcor1 = 1, bcor2 = 1;
        int bret1 = MagnetVarp::btempcor(bcor1, 0, 1);
        int bret2 = MagnetVarp::btempcor(bcor2, 0, 2);
        if( bret1 == 0 && bret2 == 0 ) bcor = ( bcor1 + bcor2 ) / 2;
        else if( bret1 != 0 && bret2 == 0 ) bcor = bcor2;
        rigidity *= bcor;
      }
    }

    trkFitCodeL1Inner              = mfit;
    trkRigidityL1Inner             = rigidity;
    trkRigidityInverseErrorL1Inner = pTrTrack->GetErrRinv(trkFitCodeL1Inner);
    trkReducedChisquareL1InnerX    = pTrTrack->GetNormChisqX(trkFitCodeL1Inner);
    trkReducedChisquareL1InnerY    = pTrTrack->GetNormChisqY(trkFitCodeL1Inner);

    int id_maxspan = pTrTrack->iTrTrackPar(1, 0, 20+refit, mass, z_inn);
    if( id_maxspan >= 0 )
    {
      if( isiss ) trkRigidityMS = pTrTrack->GetCorrectedRigidity(id_maxspan);
      else trkRigidityMS = pTrTrack->GetRigidity(id_maxspan);
      // Magnet temperature correction for ISS pass4
      float bcor = pTrTrack->GetBcorr(id_maxspan); // positive
      if( isiss && bcor == 1 )
      {
        float bcor1 = 1, bcor2 = 1;
        int bret1 = MagnetVarp::btempcor(bcor1, 0, 1);
        int bret2 = MagnetVarp::btempcor(bcor2, 0, 2);
        if( bret1 == 0 && bret2 == 0 ) bcor = ( bcor1 + bcor2 ) / 2;
        else if( bret1 != 0 && bret2 == 0 ) bcor = bcor2;
        trkRigidityMS *= bcor;
      }
      trkFitCodeMS              = id_maxspan;
      trkRigidityInverseErrorMS = pTrTrack->GetErrRinv(id_maxspan);
      trkReducedChisquareMSX    = pTrTrack->GetNormChisqX(id_maxspan);
      trkReducedChisquareMSY    = pTrTrack->GetNormChisqY(id_maxspan);
    }
    else
    {
      trkFitCodeMS              = id_maxspan;
      trkRigidityMS             = -99999.;
      trkRigidityInverseErrorMS = -99999.;
      trkReducedChisquareMSX    = -99999.;
      trkReducedChisquareMSY    = -99999.;
    }


    int id_inner = pTrTrack->iTrTrackPar(1, 3, 20+refit, mass, z_inn);
    if( id_inner >= 0 )
    {
      if( isiss ) trkRigidityInner = pTrTrack->GetCorrectedRigidity(id_inner);
      else trkRigidityInner = pTrTrack->GetRigidity(id_inner);
      // Magnet temperature correction for ISS pass4
      float bcor = pTrTrack->GetBcorr(id_inner); // positive
      if( isiss && bcor == 1 )
      {
        float bcor1 = 1, bcor2 = 1;
        int bret1 = MagnetVarp::btempcor(bcor1, 0, 1);
        int bret2 = MagnetVarp::btempcor(bcor2, 0, 2);
        if( bret1 == 0 && bret2 == 0 ) bcor = ( bcor1 + bcor2 ) / 2;
        else if( bret1 != 0 && bret2 == 0 ) bcor = bcor2;
        trkRigidityInner *= bcor;
      }
      trkFitCodeInner              = id_inner;
      trkRigidityInverseErrorInner = pTrTrack->GetErrRinv(id_inner);
      trkReducedChisquareInnerX    = pTrTrack->GetNormChisqX(id_inner);
      trkReducedChisquareInnerY    = pTrTrack->GetNormChisqY(id_inner);
    }
    else
    {
      trkFitCodeInner              = id_inner;
      trkRigidityInner             = -99999.;
      trkRigidityInverseErrorInner = -99999.;
      trkReducedChisquareInnerX    = -99999.;
      trkReducedChisquareInnerY    = -99999.;
    }


    for(int ii = 0; ii < 9; ii++)
    {
      trkLayerJQ[ii] = 0.;
      trkEdepLayerJ[ii] = 0.;
    }

    for( int ilayer = 0; ilayer <9; ilayer++)
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
    }

    //
    //
    // ECAL shower Section
    //
    //

    // Save ECAL shower related variables in case of EcalShowerR object exists.
    EcalShowerR* pEcalShower = pParticle->pEcalShower();
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
      pTrTrack->Interpolate(pEcalShower->CofG[2], trackPoint, trackDirection, id_maxspan);
      pEcalShower->NormaliseVariableLAPP();
      showerCofGDist = std::sqrt(std::pow(trackPoint[0] - showerCofG[0], 2) + std::pow( trackPoint[1] - showerCofG[1], 2));
      showerCofGdX   = std::fabs(trackPoint[0] - showerCofG[0]);
      showerCofGdY   = std::fabs(trackPoint[1] - showerCofG[1]);
    }

    // ACSoft related lines
    particleFactory.SetAMSTrTrackR(pTrTrack);
    particleFactory.SetAMSEcalShowerR(pEcalShower);
    particleFactory.SetAMSBetaHR(pBetaH);

    if( !amsRootSupport.SwitchToSpecificTrackFitById(id_maxspan) ) continue;
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

    const Analysis::Particle* particle = event.PrimaryParticle();
    assert(particle);

    tofUpperCharge = particle->UpperTofCharge();
    tofLowerCharge = particle->LowerTofCharge();

    trdPElectronToProtonLogLikelihoodRatio = particle->CalculateElectronProtonLikelihood();
    trdPHeliumToProtonLogLikelihoodRatio   = particle->CalculateHeliumProtonLikelihood();
    trdPHeliumToElectronLogLikelihoodRatio = particle->CalculateHeliumElectronLikelihood();

    pRichRing = pParticle->pRichRing();
    if( pRichRing )
    {
      isRichAvailable            = 1;
      richRebuild                = (int)pRichRing->Rebuild();
      richIsGood                 = (int)pRichRing->IsGood();
      richIsClean                = (int)pRichRing->IsClean();
      richIsNaF                  = (int)pRichRing->IsNaF();
      richUsedHits               = (int)pRichRing->getUsedHits();
      richRingWidth              = (float)pRichRing->RingWidth();
      richNHits                  = pRichRing->getHits();
      richNPMTsOnRing            = pRichRing->getPMTs();
      richBeta                   = pRichRing->getBeta();
      richBetaError              = pRichRing->getBetaError();
      richChargeSquared          = pRichRing->getCharge2Estimate();
      richKolmogorovProbability  = pRichRing->getProb();
      richPhotoelectrons         = pRichRing->getPhotoElectrons();
      richExpectedPhotoelectrons = pRichRing->getExpectedPhotoelectrons();
      richTheta                  = pRichRing->getTrackTheta();
      richPhi                    = pRichRing->getTrackPhi();
    }

    pTrdTrack = pParticle->pTrdTrack();
    int trdNUsedHits = 0;
    int trdNUsedSegment = 0;
    if(!pTrdTrack)
    {
      trdNClusters = -1;
      trdNUsedHits = -1;
      trdNUsedSegment = -1;
      trdNUsedHits = -1;
      trdNTracks = -1;
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
    }

    // Below of this for TrdK

    trdKNRawHits = pAMSEvent->NTrdRawHit();
    if(trdKNRawHits > 0)
    {
      double trdKLikelihoodRatio[3] = {-1., -1., -1.};
      double trdKLikelihood[3]      = {-1., -1., -1.};

      pTrdKCluster = new TrdKCluster(pAMSEvent, pTrTrack, id_maxspan);
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

    // The following code lines are about ACC

    /*
    pAMSEvent->RebuildAntiClusters();
    accNHits = pAMSEvent->nAntiCluster();

    if( accNHits != 0 )
    {
      for( int i = 0; i < accNHits; ++i)
      {
        AntiClusterR* anti = pAMSEvent->pAntiCluster(i);
        if( !anti ) continue;
        if( anti->Edep != 0 )
          accEdep.push_back(anti->Edep);
        else
          accEdep.push_back(-1);
        accTime.push_back(anti->time);
      }
    }

    accSector.clear();
    accTime.clear();
    accHitPosZ.clear();
    accChi2.clear();
    accNPairs.clear();
    accUnfoldedHitPosZ.clear();
    accRawCharge.clear();
    accHitPosZfromADC.clear();
    accUnfoldedHitPosZfromADC.clear();
    accEdep.clear();
    accTimePG.clear();
    */

    if( e % 100000 == 0 || e == nevt - 1 )
      KNUOUT << "Processed " << e << " out of " << nevt << " (" << (float)e/nevt*100. << "%)" << endl;

    nSelected++;
    outTree->Fill();
    if( isMC == false )
      hEvtCounter->Fill(7);
    else
      hEvtCounter->Fill(6);
  } // End of the Loop
}/*}}}*/

/**
 * @brief
 * @return
 */
bool KNUTree::IsGoodLiveTime(AMSEventR* thisEvent)
{/*{{{*/
  return true;
}/*}}}*/



/**
 * @brief
 * @return
 */
bool KNUTree::IsInSouthAtlanticAnomaly(AMSEventR* thisEvent)
{/*{{{*/
  return false;
}/*}}}*/



/**
 * @brief
 * @return
 */
bool KNUTree::IsInSolarArrays(AMSEventR* thisEvent)
{/*{{{*/
  return false;
}/*}}}*/

