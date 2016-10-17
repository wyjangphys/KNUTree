#include "KNUTree.h"

/**
 * @brief
 * @return
 */
int KNUTree::GetGoodParticleIndex(AMSEventR* thisEvent)
{/*{{{*/
  bool debugMode = true;
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
    return -1;           // In case of no reconstructed particle return -1 and exit!

  if( debugMode ) cerr << "Event : " << thisEvent->Event() << " / ";
  for( int i = 0; i < nParticle; ++i )
  {
    pParticle   = thisEvent->pParticle(i);
    if( !pParticle ){ hGoodParticleCounter->Fill(1); if( debugMode ){ cerr << "-2"; return -2; }} // Exit for no ParticleR pointer!
    pBetaH      = pParticle->pBetaH();
    if( !pBetaH    ){ hGoodParticleCounter->Fill(2); if( debugMode ){ cerr << "-3"; return -3; }} // Exit for no BetaHR pointer!
    pTrTrack    = pParticle->pTrTrack();
    if( !pTrTrack  ){ hGoodParticleCounter->Fill(3); if( debugMode ){ cerr << "-4"; return -4; }} // Exit for no TrTrackR pointer!
    //pTrdTrack   = pParticle->pTrdTrack();
    //if( !pTrdTrack ){ hGoodParticleCounter->Fill(4); if( debugMode ){ cerr << "-6"; return -6; }} // Exit for no TrdTrackR pointer!

    /*
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
      */
    /* ECAL related cuts are removed since we don't need to require ECAL information for deuteron analysis
     * 2016.6.9
     */
    /*
    pEcalShower = pParticle->pEcalShower();
    if( debugMode ) cerr << " / EcalShowerR : ";
    if( !pEcalShower )
    {
      if( debugMode ) cerr << "0";
      return -5;
    }
    else
      if( debugMode ) cerr << "1";
      */

    /*
    pTrdTrack = pParticle->pTrdTrack();
    if( !pTrdTrack )
    {
      if( debugMode ) cerr << "0";
      return -6;
    }
    else
      if( debugMode ) cerr << "1";
    */
//    pRichRing   = pParticle->pRichRing();

//    GoodBeta      = pBetaH->IsGoodBeta(); // BetaHR::IsGoodBeta() will return true when track passes 4 layers of TOF measured hits.
    GoodBeta      = IsGoodBeta(pBetaH);     // To reject events with beta less than 0.4 and to check build type.
    //TkTofMatch    = pBetaH->IsTkTofMatch(); // BetaHR::IsTkTofMatch() will return true when Tracker track and TOF geometry are match.
    GoodTrTrack   = IsGoodTrTrack(pTrTrack);
    /* ECAL related cuts are removed since we don't need to require ECAL information for deuteron analysis
     * 2016.6.9
     */
//    ShowerTkMatch = IsShowerTrackMatched(pEcalShower, pTrTrack);
//    FidVolumeTest = IsTrackInsideEcalFiducialVolume(pTrTrack);

    //if( debugMode ) cerr << " / Result: " << GoodBeta << TkTofMatch << GoodTrTrack << ShowerTkMatch << FidVolumeTest;
    //result = GoodBeta && TkTofMatch && GoodTrTrack;
    result = GoodBeta && GoodTrTrack;
    //if( debugMode ) cerr << " / Final : " << result;

    /*
    cout << "Result: " << endl;
    cout << "GoodBeta: " << GoodBeta;
    cout << " / TkTofMatch: " << TkTofMatch;
    cout << " / GoodTrTrack: " << GoodTrTrack;
    cout << " / Final: " << result << endl << endl;
    */

    nGoodParticle++;
    if( nGoodParticle > 1 )
    {
      if( debugMode )
      {
        KNUERR << "More than two good particles" << endl;
        return -1; // return -1 if there are more than one good particle
      }
    }
    else if( result == false )
    {
      if( debugMode )
      {
        KNUERR << "Failed to find a good particle!" << std::endl;
        KNUERR << "GoodBeta: " << GoodBeta << std::endl;
        KNUERR << "TkTofMatch: " << TkTofMatch << std::endl;
        KNUERR << "GoodTrTrack: " << GoodTrTrack << std::endl;
        return -1;
      }
    }
    else
      iGoodParticle = i;
  }
  //if( debugMode ) cerr << " / iGoodParticle : " << iGoodParticle << endl;
  return iGoodParticle;
}/*}}}*/

