#include "KNUTree.h"

/**
 * @brief
 * @return
 */
bool KNUTree::IsShowerTrackMatched(EcalShowerR* thisShower, TrTrackR* thisTrack)
{/*{{{*/
  if( !thisShower || !thisTrack ) return false;

  int      id_inner = thisTrack->iTrTrackPar(1, 3, 0);
  AMSPoint tkPoint;
  AMSDir   tkDir;

  thisTrack->Interpolate(thisShower->CofG[2], tkPoint, tkDir, id_inner);
  double dX = TMath::Abs(tkPoint[0] - thisShower->CofG[0] );
  double dY = TMath::Abs(tkPoint[1] - thisShower->CofG[1] );

  if( dX < 3.6 && dY < 7.2 ) return true;
  else return false;
}/*}}}*/

