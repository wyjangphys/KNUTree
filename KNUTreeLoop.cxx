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

#define KNUOUT std::cout << "[KNUTree::Loop] "
#define KNUERR std::cerr << "[KNUTree::Loop] Error: "

void KNUTree::Loop()
{/*{{{*/
  KNUOUT << "Begin to loop over " << nevt << " events." << endl;

  AMSSetupR::RTI::UseLatest();
  TkDBc::UseFinal();

  const bool setAMSRootDefaults = true;
  AMSRootSupport amsRootSupport(AC::ISSRun, setAMSRootDefaults);
  Analysis::EventFactory& eventFactory = amsRootSupport.EventFactory();
  Analysis::AMSRootParticleFactory& particleFactory = amsRootSupport.ParticleFactory();

  AMSEventR* pev = NULL;
  HeaderR* header = NULL;
  for( unsigned long int e = 0; e < nevt; e++ )
  {
    pev = pChain->GetEvent(e);
    header = &(pev->fHeader);

    if( !pev ) continue;
    hEvtCounter->Fill(0);
    if( this->IsBadRun(pev) ) this->End(1);
    hEvtCounter->Fill(1);
    if( !this->IsScienceRun(pev) ) continue;
    hEvtCounter->Fill(2);
    if( !this->IsHardwareStatusGood(pev) ) continue;
    hEvtCounter->Fill(3);
    if( this->IsUnbiasedPhysicsTriggerEvent(pev) ) continue;
    hEvtCounter->Fill(4);
    if( pev->nParticle() != 1 ) continue;
    hEvtCounter->Fill(5);
    if( !this->IsTrkAlignmentGood(pev) ) continue;
    hEvtCounter->Fill(6);
    int ptlIdx = this->GetGoodParticleIndex(pev);
    if( ptlIdx == -1 ) continue;
    hEvtCounter->Fill(7);

    nRun            = pev->Run();
    nEvent          = pev->Event();
    nLevel1         = pev->nLevel1();
    nParticle       = pev->nParticle();
    nCharge         = pev->nCharge();
    nTrTrack        = pev->nTrTrack();
    nTrdTrack       = pev->nTrdTrack();
    nAntiCluster    = pev->nAntiCluster();
    nRichRing       = pev->nRichRing();
    nRichRingB      = pev->nRichRingB();
    nBeta           = pev->nBeta();
    nBetaB          = pev->nBetaB();
    nBetaH          = pev->nBetaH();
    nShower         = pev->nEcalShower();
    nVertex         = pev->nVertex();
    liveTime        = pev->LiveTime();
    utcTime         = (float)header->UTCTime(0);
    orbitAltitude   = header->RadS;
    orbitLatitude   = header->ThetaS;
    orbitLongitude  = header->PhiS;
    orbitLatitudeM  = header->ThetaM;
    orbitLongitudeM = header->PhiM;
    velR            = header->VelocityS;
    velTheta        = header->VelTheta;
    velPhi          = header->VelPhi;
    yaw             = header->Yaw;
    pitch           = header->Pitch;
    roll            = header->Roll;
    zenithAngle     = header->Zenith();
    isInSAA         = pev->IsInSAA();

    double tmpglong;
    double tmpglat;
    gCoordCalcResult  = pev->GetGalCoo(
        gCoordCalcResultBit,
        tmpglong,
        tmpglat);
    gLongitude  = tmpglong;
    gLatitude   = tmpglat;

    double tmpSunPosAzimuth;
    double tmpSunPosElevation;
    sunPosCalcResult   = header->getSunAMS(tmpSunPosAzimuth, tmpSunPosElevation);
    sunPosAzimuth      = tmpSunPosAzimuth;
    sunPosElevation    = tmpSunPosElevation;

    AMSPoint solarArray;
    isInShadow         = pev->isInShadow(solarArray, 0);
    solarArrayCoord[0] = solarArray.x();
    solarArrayCoord[1] = solarArray.y();
    solarArrayCoord[2] = solarArray.z();

    ParticleR* pParticle = pev->pParticle(ptlIdx);
    if(!pParticle) continue;
    ptlCharge           = (unsigned int)pParticle->Charge;
    // Unfold below to see the particle coordinates
    ptlMomentum         = pParticle->Momentum;/*{{{*/
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
    ptlCutOffStoermer   = pParticle->CutoffS;
    ptlCutOffDipole     = pParticle->Cutoff;/*}}}*/

    tofNCluster         = pev->nTofCluster();
    tofNClusterH        = pev->nTofClusterH();
    BetaHR* pBetaH       = pParticle->pBetaH();
    //if(!pBeta) continue;
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
    for(int i = 0; i < pBetaH->NTofClusterH(); ++i)
    {
      if( pBetaH->TestExistHL(i) == true )
        tofDepositedEnergyOnLayer[i] = pBetaH->GetClusterHL(i)->GetEdep();
      else
        tofDepositedEnergyOnLayer[i] = -1;
    }

    int ncls[4] = {0, 0, 0, 0};
    nTofClustersInTime  = pev->GetNTofClustersInTime(pBetaH, ncls);


    TrTrackR* pTrTrack = pParticle->pTrTrack();
    if(!pTrTrack) continue;
    int id_maxspan = pTrTrack->iTrTrackPar(1, 0, 20);

    trkFitCodeMS              = id_maxspan;
    trkRigidityMS             = pTrTrack->GetRigidity(id_maxspan);
    trkRigidityInverseErrorMS = pTrTrack->GetErrRinv(id_maxspan);
    trkReducedChisquareMSX    = pTrTrack->GetNormChisqX(id_maxspan);
    trkReducedChisquareMSY    = pTrTrack->GetNormChisqY(id_maxspan);

    int id_inner = pTrTrack->iTrTrackPar(1, 3, 20);
    if( id_inner >= 0 )
    {
      trkFitCodeInner              = id_inner;
      trkRigidityInner             = pTrTrack->GetRigidity(id_inner);
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


    TrRecHitR* pTrRecHit = NULL;
    for(int ii = 0; ii < 9; ii++) trkEdepLayerJ[ii] = 0.;

    for( int ilayer = 0; ilayer <9; ilayer++)
    {
      pTrRecHit = pTrTrack->GetHitLJ(ilayer);
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
    trkCharge       = pTrTrack->GetQ();
    trkInnerCharge  = pTrTrack->GetInnerQ();
    trkHasExtLayers = pTrTrack->HasExtLayers();

    TrdTrackR* pTrdTrack = pParticle->pTrdTrack();
    //if( !pTrdTrack ) continue;

    // Initialize ECAL shower related variables
    showerEnergyD         = -1;
    for(int i = 0; i < 18; ++i) showerEnergyDL[i] = 0;
    showerEnergyE         = -1;
    showerEnergyCorrected = -1;
    showerBDT             = -1;
    showerCofG[0]         = -99;
    showerCofG[1]         = -99;
    showerCofG[2]         = -99;
    showerCofGDist        = -99;
    showerCofGdX          = -99;
    showerCofGdY          = -99;

    // Save ECAL shower related variables in case of EcalShowerR object exists.
    EcalShowerR* pEcalShower = pParticle->pEcalShower();
    if( pEcalShower )
    {
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
    Analysis::Event& event = amsRootSupport.BuildEvent(pChain, pev);

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

    RichRingR* richRing = pParticle->pRichRing();
    if(!richRing)
    {
      richRebuild               = -1;
      richIsGood                = -1;
      richIsClean               = -1;
      richIsNaF                 = -1;
      richRingWidth             = -1;
      richUsedHits              = -1;
      richNHits                 = -1;
      richNPMTsOnRing           = -1;
      richBeta                  = -1;
      richBetaError             = -1;
      richChargeSquared         = -1;
      richKolmogorovProbability = -1;
      richPhotoelectrons        = -1;
      richExpectedPhotoelectrons = -1;
      richTheta                 = -9;
      richPhi                   = -9;
    }
    else
    {
      richRebuild                = (int)richRing->Rebuild();
      richIsGood                 = (int)richRing->IsGood();
      richIsClean                = (int)richRing->IsClean();
      richIsNaF                  = (int)richRing->IsNaF();
      richUsedHits               = (int)richRing->getUsedHits();
      richRingWidth              = (float)richRing->RingWidth();
      richNHits                  = richRing->getHits();
      richNPMTsOnRing            = richRing->getPMTs();
      richBeta                   = richRing->getBeta();
      richBetaError              = richRing->getBetaError();
      richChargeSquared          = richRing->getCharge2Estimate();
      richKolmogorovProbability  = richRing->getProb();
      richPhotoelectrons         = richRing->getPhotoElectrons();
      richExpectedPhotoelectrons = richRing->getExpectedPhotoelectrons();
      richTheta                  = richRing->getTrackTheta();
      richPhi                    = richRing->getTrackPhi();
    }

    TrdTrackR* trdTrack = pParticle->pTrdTrack();
    if(!trdTrack) continue;
    trdNClusters = pev->nTrdCluster();
    int trdNUsedHits = 0;
    int trdNUsedSegment = trdTrack->NTrdSegment();
    for(int i = 0; i < trdNUsedSegment; i++)
    {
      TrdSegmentR* pTrdSegment = trdTrack->pTrdSegment(i);
      trdNUsedHits += pTrdSegment->NTrdCluster();
    }
    trdNUnusedHits = trdNClusters - trdNUsedHits;

    trdNTracks  = pev->nTrdTrack();
    if(!trdTrack)
    {
      trdTrackTheta     = -9.;
      trdTrackPhi       = -9.;
      trdTrackPattern   = -9;
      trdTrackCharge    = -9;
      for(int i = 0; i < 20; ++i) trdTrackEdepL[i] = -9.;
    }
    else
    {
      trdTrackTheta   = trdTrack->Theta;
      trdTrackPhi     = trdTrack->Phi;
      trdTrackPattern = trdTrack->Pattern;
      trdTrackCharge  = trdTrack->Q;
      trdTrackChi2    = trdTrack->Chi2;
      trdTrackMeanDepositedEnergy = 0.;
      int ntrdlayers = 0;
      for(int i = 0; i < trdTrack->NTrdSegment(); i++)
      {
        for(int j = 0; j < trdTrack->pTrdSegment(i)->NTrdCluster(); j++)
        {
          int trdLayer = trdTrack->pTrdSegment(i)->pTrdCluster(j)->Layer;
          float trdEdep = trdTrack->pTrdSegment(i)->pTrdCluster(j)->EDep;
          trdTrackEdepL[trdLayer] = trdEdep;
          trdTrackTotalDepositedEnergy += trdEdep;
          ntrdlayers++;
        }
      }
      trdTrackMeanDepositedEnergy = trdTrackTotalDepositedEnergy/(float)ntrdlayers;
    }

    // Below of this for TrdK
    trdKNRawHits                                = -1;
    trdKIsReadAlignmentOK                       = -1;
    trdKIsReadCalibOK                           = -1;
    trdKNHits                                   = -1;
    trdKIsValid                                 = -1;
    trdKElectronToProtonLogLikelihoodRatio      = -1;
    trdKHeliumToElectronLogLikelihoodRatio      = -1;
    trdKHeliumToProtonLogLikelihoodRatio        = -1;
    trdKCharge                                  = -1;
    trdKChargeError                             = -1;
    trdKNUsedHitsForCharge                      = -1;
    for(int k = 0; k < 20; k++) trdKAmpLayer[k] = -1;
    trdKTotalPathLength                         = -1;
    trdKElectronLikelihood                      = -1;
    trdKProtonLikelihood                        = -1;
    trdKHeliumLikelihood                        = -1;

    trdKNRawHits = pev->NTrdRawHit();
    if(trdKNRawHits > 0)
    {
      double trdKLikelihoodRatio[3] = {-1., -1., -1.};
      double trdKLikelihood[3]      = {-1., -1., -1.};

      TrdKCluster* trdK = new TrdKCluster(pev, pTrTrack, id_maxspan);
      if(!trdK) continue;

      trdKIsReadAlignmentOK = trdK->IsReadAlignmentOK;
      trdKIsReadCalibOK     = trdK->IsReadCalibOK;
      trdKIsValid           = trdK->GetLikelihoodRatio_TrTrack(15, trdKLikelihoodRatio, trdKLikelihood, trdKNHits, trdKTotalPathLength, trdKTotalAmp, -1, 0);
      if(trdKIsValid != 0 && trdKNHits != 0)
      {
        trdK->CalculateTRDCharge();
        trdKCharge                              = trdK->GetTRDCharge();
        trdKChargeError                         = trdK->GetTRDChargeError();
        trdKNUsedHitsForCharge                  = trdK->GetQTRDHitCollectionNuclei().size();
        trdKElectronToProtonLogLikelihoodRatio  = trdKLikelihoodRatio[0];
        trdKHeliumToElectronLogLikelihoodRatio  = trdKLikelihoodRatio[1];
        trdKHeliumToProtonLogLikelihoodRatio    = trdKLikelihoodRatio[2];
        trdKElectronLikelihood                  = trdKLikelihood[0];
        trdKProtonLikelihood                    = trdKLikelihood[1];
        trdKHeliumLikelihood                    = trdKLikelihood[2];

        AMSPoint trExtraP0;
        AMSDir   trExtraDir;
        trdK->GetTrTrackExtrapolation(trExtraP0, trExtraDir);

        for(int l = 0; l < trdK->NHits(); l++)
        {
          TrdKHit* trdKHit = trdK->GetHit(l);
          int tmpL = 0;
          tmpL = trdKHit->TRDHit_Layer;
          float tmpAmp = 0;
          tmpAmp = trdKHit->TRDHit_Amp;
          float tmpPathLength = 0;
          tmpPathLength = trdKHit->Tube_Track_3DLength(&trExtraP0, &trExtraDir);

          if( tmpAmp < 15 || tmpPathLength <= 0 ) continue;
          if( trdKHit->IsAligned == 0 || trdKHit->IsCalibrated == 0 ) continue;
          trdKAmpLayer[tmpL] += tmpAmp;
        }
      }
    }

    // The following code lines are about ACC

    pev->RebuildAntiClusters();
    accNHits = pev->nAntiCluster();

    if( accNHits != 0 )
    {
      for( int i = 0; i < accNHits; ++i)
      {
        AntiClusterR* anti = pev->pAntiCluster(i);
        if( !anti ) continue;
        if( anti->Edep != 0 )
          accEdep.push_back(anti->Edep);
        else
          accEdep.push_back(-1);
        accTime.push_back(anti->time);
      }
    }


    if( e % 100000 == 0 || e == nevt - 1 )
      KNUOUT << "Processed " << e << " out of " << nevt << " (" << (float)e/nevt*100. << "%)" << endl;

    nProcessedNumber = nProcessed;
    outTree->Fill();

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

    nProcessed++;
    hEvtCounter->Fill(8);
  }
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

