#include "KNUTree.h"

/**
 * @brief
 * @return
 */
bool KNUTree::IsGoodTrTrack(TrTrackR* thisTrack)
{/*{{{*/
  // 1. THE EVENT SHOULD HAVE MAXSPAN TRACK
  bool debugMode = false;

  if( !thisTrack )
  {
    if( debugMode ) cerr << "No track pointer. Reject it." << endl;
    return false;
  }
  if( thisTrack->IsFake() )
  {
    if( debugMode ) cerr << "Found a fake track. Reject it." << endl;
    return false;
  }

  int id_maxspan;
  id_maxspan = thisTrack->iTrTrackPar(1, 0, 20);

  /*
  if( id_fullspan < 0 )
  {
    if( debugMode == true )
    {
      switch(id_fullspan)
      {
        case -1: cerr << "[FS]The requested fit cannot be performed on this track." << endl; break;
        case -2: cerr << "[FS]The requested fit it is not available without refitting." << endl; break;
        case -3: cerr << "[FS]The refit failed." << endl; break;
        case -4: cerr << "[FS]Should not happen!! Contact the developers!." << endl; break;
        case -5: cerr << "[FS]The refit failed because there was a problem retrieving dynamic alignment for Ext Layers" << endl; break;
        default: break;
      }
    }

    return false;
  }
  */

  if( id_maxspan < 0 )
  {
    if( debugMode == true )
    {
      switch( id_maxspan )
      {
        case -1: cerr << "[MS]The requested fit cannot be performed on this track." << endl; break;
        case -2: cerr << "[MS]The requested fit it is not available without refitting." << endl; break;
        case -3: cerr << "[MS]The refit failed." << endl; break;
        case -4: cerr << "[MS]Should not happen!! Contact the developers!." << endl; break;
        case -5: cerr << "[MS]The refit failed because there was a problem retrieving dynamic alignment for Ext Layers" << endl; break;
        default: break;
      }
    }

    return false;
  }

  /*
  if( id_inner < 0 )
  {
    if( debugMode )
    {
      switch( id_inner )
      {
        case -1: cerr << "[IN]The requested fit cannot be performed on this track." << endl; break;
        case -2: cerr << "[IN]The requested fit it is not available without refitting." << endl; break;
        case -3: cerr << "[IN]The refit failed." << endl; break;
        case -4: cerr << "[IN]Should not happen!! Contact the developers!." << endl; break;
        case -5: cerr << "[IN]The refit failed because there was a problem retrieving dynamic alignment for Ext Layers" << endl; break;
        default: break;
      }
    }

    return false;
  }
  */

  bool hitOnLayerJ[9];
  for(int i = 0; i < 9; i++) hitOnLayerJ[i] = thisTrack->TestHitBitsJ(i, id_maxspan);

  if( !(hitOnLayerJ[2] || hitOnLayerJ[3]) ) { if( debugMode ){cerr << "No hit on layer 2/3 " << endl;} return false; }
  if( !(hitOnLayerJ[4] || hitOnLayerJ[5]) ) { if( debugMode ){cerr << "No hit on layer 4/5 " << endl;} return false; }
  if( !(hitOnLayerJ[6] || hitOnLayerJ[7]) ) { if( debugMode ){cerr << "No hit on layer 6/7 " << endl;} return false; }

  // Reject events with rigidity less than 0.5 GV
  if( thisTrack->GetRigidity(id_maxspan) < 0.5 ) return false;

  return true;
}/*}}}*/

