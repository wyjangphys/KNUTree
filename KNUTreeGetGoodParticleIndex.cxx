#include "KNUTree.h"

/**
 * @brief
 * @return
 */
int KNUTree::GetGoodParticleIndex(AMSEventR* thisEvent)
{/*{{{*/
  bool debugMode = false;
  int nParticle = 0;
  int nGoodParticle = 0;
  int iGoodParticle = -1;

  ParticleR*   pParticle       = NULL;
  BetaHR*      pBetaH          = NULL;
  TrTrackR*    pTrTrack        = NULL;
  EcalShowerR* pEcalShower     = NULL;
  TrdTrackR*   pTrdTrack       = NULL;
  RichRingR*   pRichRing       = NULL;
  bool         GoodBeta        = false;
  bool         TkTofMatch      = false;
  bool         GoodTrTrack     = false;
  bool         ShowerTkMatch   = false;
  bool         FidVolumeTest   = false;
  bool         GoodRing        = false;
  bool         result          = false;

  nParticle = thisEvent->nParticle();

  if( nParticle < 1 )
    return -1;

  if( debugMode ) cerr << "Event : " << thisEvent->Event() << " / ";
  for( int i = 0; i < nParticle; ++i )
  {
    pParticle   = thisEvent->pParticle(i);
    if( debugMode ) cerr << "ParticleR : ";
    if( !pParticle   )
    {
      if( debugMode ) cerr << "0";
      return -2;
    }
    else
      if( debugMode ) cerr << "1";
    pBetaH      = pParticle->pBetaH();
    if( debugMode ) cerr << " / BetaHR : ";
    if( !pBetaH      )
    {
      if( debugMode ) cerr << "0";
      return -3;
    }
    else
      if( debugMode ) cerr << "1";
    pTrTrack    = pParticle->pTrTrack();
    if( debugMode ) cerr << " / TrTrackR : ";
    if( !pTrTrack    )
    {
      if( debugMode ) cerr << "0";
      return -4;
    }
    else
      if( debugMode ) cerr << "1";
    pEcalShower = pParticle->pEcalShower();
    if( debugMode ) cerr << " / EcalShowerR : ";
    if( !pEcalShower )
    {
      if( debugMode ) cerr << "0";
      return -5;
    }
    else
      if( debugMode ) cerr << "1";
    pTrdTrack = pParticle->pTrdTrack();
    if( !pTrdTrack )
    {
      if( debugMode ) cerr << "0";
      return -6;
    }
    else
      if( debugMode ) cerr << "1";
//    pRichRing   = pParticle->pRichRing();

    // BetaHR::IsGoodBeta() will return true when track passes 4 layers of TOF measured hits.
    // BetaHR::IsTkTofMatch() will return true when Tracker track and TOF geometry are match.
    GoodBeta      = pBetaH->IsGoodBeta();
    TkTofMatch    = pBetaH->IsTkTofMatch();
    GoodTrTrack   = IsGoodTrTrack(pTrTrack);
    ShowerTkMatch = IsShowerTrackMatched(pEcalShower, pTrTrack);
    FidVolumeTest = IsTrackInsideEcalFiducialVolume(pTrTrack);

    if( debugMode ) cerr << " / Result: " << GoodBeta << TkTofMatch << GoodTrTrack << ShowerTkMatch << FidVolumeTest;
    result = GoodBeta && TkTofMatch && GoodTrTrack && ShowerTkMatch && FidVolumeTest;
    if( debugMode ) cerr << " / Final : " << result;
    nGoodParticle++;
    if( nGoodParticle > 1 )
    {
      cout << "More than two good particles" << endl;
      return -1; // return -1 if there are more than one good particle
    }
    else if( result == false ) return -1;
    else
      iGoodParticle = i;
  }
  if( debugMode ) cerr << " / iGoodParticle : " << iGoodParticle << endl;
  return iGoodParticle;
}/*}}}*/

