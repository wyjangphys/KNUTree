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
  int iGoodParticle;
  //int  iGoodParticle = -1;

  ParticleR*   pParticle           = NULL;
  BetaHR*      pBetaH              = NULL;
  TrTrackR*    pTrTrack            = NULL;
  EcalShowerR* pEcalShower         = NULL;
  TrdTrackR*   pTrdTrack           = NULL;
  RichRingR*   pRichRing           = NULL;
  bool         BetaHRGoodBetaTest  = false;
  bool         KNUTreeGoodBetaTest = false;
  bool         GoodBeta            = false;
  bool         TkTofMatch          = false;
  bool         GoodTrTrack         = false;
  bool         ShowerTkMatch       = false;
  bool         FidVolumeTest       = false;
  bool         GoodRing            = false;
  bool         result              = false;

  nParticle = thisEvent->nParticle();

  if( nParticle < 1 ) return -1; // In case of no reconstructed particle return -1!

  if( debugMode ) KNUERR << "Event: " << thisEvent->Event() << " / " << "nParticle: " << nParticle;
  for( int i = 0; i < nParticle; ++i )
  {
    // Object pointer check
    pParticle   = thisEvent->pParticle(i);
    if( !pParticle ){ hGoodParticleCounter->Fill(1); if( debugMode ){ cerr << " / CUTTED / Reason: pParticle is NULL.\n"; return -2; }} // Exit for no ParticleR pointer!
    pBetaH      = pParticle->pBetaH();
    if( !pBetaH    ){ hGoodParticleCounter->Fill(2); if( debugMode ){ cerr << " / CUTTED / Reason: pBeta is NULL.\n"; return -3; }} // Exit for no BetaHR pointer!
    pTrTrack    = pParticle->pTrTrack();
    if( !pTrTrack  ){ hGoodParticleCounter->Fill(3); if( debugMode ){ cerr << " / CUTTED / Reason: pTrTrack is NULL.\n"; return -4; }} // Exit for no TrTrackR pointer!

    // Check whether the particle meets preselection criteria.
    //
    // BetaHR::IsGoodBeta() will return true when track passes 4 layers of TOF measured hits.
    // KNUTree::IsGoodBeta() is here to reject events with beta less than 0.4 and to check TOF build type.
    if( IsGoodBeta( pBetaH ) )
    {
      GoodBeta = true;
    }
    else
    {
      GoodBeta = false;
      if( debugMode ) cerr << " / CUTTED / Reason: KNUTree::IsGoodBeta() test failed.\n";
      continue;
    }
    if( pBetaH->IsGoodBeta() )
      GoodBeta = true;
    else
    {
      GoodBeta = false;
      if( debugMode ) cerr << " / CUTTED / Reason: BetaHR::IsGoodBeta() test failed.\n";
      continue;
    }

    // BetaHR::IsTkTofMatch() will return true when Tracker track and TOF geometry agrees.
    if( pBetaH->IsTkTofMatch() )
      TkTofMatch = true;
    else
    {
      TkTofMatch = false;
      if( debugMode ) cerr << " / CUTTED / Reason: BetaHR::IsTkTofMatch() test failed.\n";
      continue;
    }

    nGoodParticle++;

    //iGoodParticle.push_back(i);
    iGoodParticle = i;
    cout << " / iParticle: " << i << " is selected as a good particle candidate.\n";
  }
  return iGoodParticle;
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
