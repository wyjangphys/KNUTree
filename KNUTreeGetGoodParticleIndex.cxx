#include "KNUTree.h"

/**
 * @brief
 * @return
 */
int KNUTree::GetGoodParticleIndex(AMSEventR* thisEvent)
{/*{{{*/
  bool debugMode     = false;
  int  nParticle     = 0;
  int  nGoodParticle = 0;
  int iresult        = -1;
  float maxRigidity  = 0;
  vector<int> iGoodParticle;
  iGoodParticle.clear();

  ParticleR*   pParticle   = NULL;
  BetaHR*      pBetaH      = NULL;
  TrTrackR*    pTrTrack    = NULL;
  EcalShowerR* pEcalShower = NULL;
  TrdTrackR*   pTrdTrack   = NULL;
  bool         isGoodBetaH   = false;
  bool         isTkTofMatch  = false;
  bool         isTkTrdMatch  = false;
  bool         isTkRingMatch = false;
  bool         result        = false;

  nParticle = thisEvent->nParticle();

  if( nParticle < 1 ) return -1; // In case of no reconstructed particle return -1!

  if( debugMode ) KNUERR << "Event: (" << thisEvent->Event() << " / " << nEvents << ")" << "nParticle: " << nParticle;
  for( int i = 0; i < nParticle; ++i )
  {
    // Object pointer check
    pParticle   = thisEvent->pParticle(i);
    pBetaH      = pParticle->pBetaH();
    pTrTrack    = pParticle->pTrTrack();
    pTrdTrack   = pParticle->pTrdTrack();

    if( pParticle == NULL ) { if( debugMode ) { KNUERR << "pParticle is NULL" << endl; } continue; }
    if( pBetaH    == NULL ) { if( debugMode ) { KNUERR << "pBetaH is NULL"    << endl; } continue; }
    if( pTrTrack  == NULL ) { if( debugMode ) { KNUERR << "pTrTrack is NULL"  << endl; } continue; }
    if( pTrdTrack == NULL ) { if( debugMode ) { KNUERR << "pTrdTrack is NULL" << endl; } continue; }

    if( pTrTrack->iTrTrackPar(1, 3, 0) < 0 ) { if( debugMode ) { KNUERR << "TrkFitCode is wrong!"            << endl; } continue; }
    if( !IsGoodBeta( pBetaH ) )              { if( debugMode ) { KNUERR << "IsGoodBeta test failed"          << endl; } continue; }
    if( !pBetaH->IsGoodBeta() )              { if( debugMode ) { KNUERR << "BetaH::IsGoodBeta test failed"   << endl; } continue; }
    if( !pBetaH->IsTkTofMatch() )            { if( debugMode ) { KNUERR << "BetaH::IsTkTofMatch test failed" << endl; } continue; }
    iGoodParticle.push_back(i);
    if( debugMode ) { KNUOUT << "Particle index " << i << " is pushed back in the GoodParticle container." << endl; }
  }

  // Determine which index has highest rigidity
  for( int i = 0; i < iGoodParticle.size(); i++ )
  {
    if( debugMode ) { KNUOUT << "[KNUTree::GetGoodParticleIndex()] " << i << " is good particle index "; }
    pParticle = thisEvent->pParticle(iGoodParticle[i]);
    pTrTrack = pParticle->pTrTrack();
    int id_inner = pTrTrack->iTrTrackPar(1, 3, 23);
    float rigidity = pTrTrack->GetRigidity(id_inner);
    if( debugMode ) { cout << "with rigidity(inner): " << rigidity << endl; }

    if( rigidity > maxRigidity )
    {
      maxRigidity = rigidity;
      iresult = iGoodParticle[i];
    }
  }
  //KNUOUT << "ParticleR index: " << iresult << " is a good particle!!" << endl;
  return iresult;
}/*}}}*/

std::vector<int> KNUTree::GetGoodParticles(AMSEventR* thisEvent)
{/*{{{*/
  bool debugMode = false;
  int nGoodParticle = 0;
  std::vector<int> GoodParticleIndexContainer;

  //if( nParticle < 1 ) return 0;
  if( debugMode ) KNUERR << "Event: " << thisEvent->Event() << " / " << "nParticle: " << nParticle;

  for( int i = 0; i < nParticle; ++i )
  {
    pParticle = thisEvent->pParticle(i);
    if( !pParticle ){ hGoodParticleCounter->Fill(1); if( debugMode ) { cerr << " / CUTTED / Reason: pParticle is NULL.\n"; } continue; } // Exit for no ParticleR pointer!
    pBetaH      = pParticle->pBetaH();
    if( !pBetaH    ){ hGoodParticleCounter->Fill(2); if( debugMode ){ cerr << " / CUTTED / Reason: pBeta is NULL.\n"; } continue; } // Exit for no BetaHR pointer!
    pTrTrack    = pParticle->pTrTrack();
    if( !pTrTrack  ){ hGoodParticleCounter->Fill(3); if( debugMode ){ cerr << " / CUTTED / Reason: pTrTrack is NULL.\n"; } continue; } // Exit for no TrTrackR pointer!

    // Check whether the particle meets preselection criteria.
    //
    // BetaHR::IsGoodBeta() will return true when track
    // passes 4 layers of TOF measured hits and the 
    // track and TOF hit match.
    if( !pBetaH->IsGoodBeta() ) { continue; }

    GoodParticleIndexContainer.push_back(i);
  }

  return GoodParticleIndexContainer;
}/*}}}*/
