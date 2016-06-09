#include "KNUTree.h"

/**
 * @brief
 * @return
 */
bool KNUTree::IsTrackInsideEcalFiducialVolume(TrTrackR* thisTrack)
{/*{{{*/
  if( !thisTrack) return false;

  int id_fullspan = thisTrack->iTrTrackPar(1, 7, 0);
  AMSPoint entryPoint, exitPoint;
  AMSDir   entryDir,   exitDir;

  float zEntry = -142.792;
  float zExit  = -158.457;

  thisTrack->Interpolate(zEntry, entryPoint, entryDir, id_fullspan);
  thisTrack->Interpolate(zExit,  exitPoint,  exitDir,  id_fullspan);

  bool Entry_in_32_4 = (TMath::Abs( entryPoint.x() )<32.4)  && (TMath::Abs( entryPoint.y() )<32.4);
  bool Exit_in_32_4  = (TMath::Abs( exitPoint.x() ) <32.4)  && (TMath::Abs( exitPoint.y() ) <32.4);
  bool Entry_in_31_5 = (TMath::Abs( entryPoint.x() )<31.5)  && (TMath::Abs( entryPoint.y() )<31.5);
  bool Exit_in_31_5  = (TMath::Abs( exitPoint.x() ) <31.5)  && (TMath::Abs( exitPoint.y() ) <31.5);

  // Request: Shower axis in ECAL volume (Entry&Exit<32.4), at least Entry||Exit within 1 cell (0.5 Moliere radius) from the border
  bool inacc = (Exit_in_32_4 && Entry_in_32_4) && ( Exit_in_31_5 || Entry_in_31_5 );
  if( inacc ) return true;
  return false;
}/*}}}*/

