#include "KNUTree.h"

/**
 * @brief   This function determines whether the given track object is useful or not.
 * @return  Good track: true / Bad track: false
 */
bool KNUTree::IsGoodTrTrack(TrTrackR* thisTrack)
{/*{{{*/
  bool debugMode = true;
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
  // Tracks should be refitted properly in order to be used for analysis for ions.
  // To refit correctly, mass and charge information is required.
  float q_inn = thisTrack->GetInnerQ();
  int z_inn = floor(q_inn + 0.5);
  float mass;
  if( z_inn < 1 ) z_inn = 1;
  mass = TrFit::Mproton;
  if( z_inn >= 2 ) mass = TrFit::Mhelium;
  if( z_inn >= 3 ) mass = TrFit::Mproton*2*z_inn;
  int refit = 3;

  if( isMC == false ) // For ISS, enable ion linearity (charge Z >= 3) and dZ correction
  {
    TRCLFFKEY.UseNonLinearity = 0; // Off default linearity correction option in datacard (it speed up refit)
    TRFITFFKEY.Zshift = 2; // Ion dZ correction.
    if( z_inn >= 3 )
      TrClusterR::SetLinearityCorrection(); // Linearity correction (important for ions with 2<Z<13)
    else
      TrClusterR::UnsetLinearityCorrection(); // Off Linearity correction (Z<3)
    refit = 3;
  }
  else  // MC smearing and refit
  {
    TrExtAlignDB::SmearExtAlign(); // MC smear ext-layer
    TRCLFFKEY.UseSensorAlign = 0;
    TRFITFFKEY.Zshift = -1;
    refit = 3;
  }

  // 1. THE EVENT SHOULD HAVE L1 + INNER PATTERN TRACK
  int id_l1inner = thisTrack->iTrTrackPar(1, 5, 20 + refit, mass, z_inn);
  if( id_l1inner < 0 )
  {
    if( debugMode )
    {
      switch( id_l1inner )
      {
        case -1: cerr << "[L1Inner]The requested fit cannot be performed on this track." << endl; break;
        case -2: cerr << "[L1Inner]The requested fit it is not available without refitting." << endl; break;
        case -3: cerr << "[L1Inner]The refit failed." << endl; break;
        case -4: cerr << "[L1Inner]Should not happen!! Contact the developers!." << endl; break;
        case -5: cerr << "[L1Inner]The refit failed because there was a problem retrieving dynamic alignment for Ext Layers" << endl; break;
        default: break;
      }
    }
    return false;
  }


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

  /*
  int id_maxspan;
  id_maxspan = thisTrack->iTrTrackPar(1, 0, 23);

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
  */

  int id_inner;
  id_inner = thisTrack->iTrTrackPar(1, 3, 20+refit, mass, z_inn);
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

  bool xside[9];
  bool yside[9];
  for( int ilayer = 0; ilayer < 9; ilayer++)
  {
    pTrRecHit = thisTrack->GetHitLJ(ilayer+1); // argument should be in the range from 1-9
    if( !pTrRecHit ) continue;
    else
    {
      if( pTrRecHit->GetEdep(0) != 0 ){xside[ilayer] = true;}
      else xside[ilayer] = false;
      if( pTrRecHit->GetEdep(1) != 0 ){yside[ilayer] = true; }
      else yside[ilayer] = false;
    }
  }

  if( !yside[1] ) { if( debugMode ) { KNUERR << "No hit on layer 2" << endl; } return false; }
  if( !(yside[2] || yside[3]) ) { if( debugMode ){ KNUERR << "No hit on layer 3/4" << endl; } return false; }
  if( !(yside[4] || yside[5]) ) { if( debugMode ){ KNUERR << "No hit on layer 5/6" << endl; } return false; }
  if( !(yside[6] || yside[7]) ) { if( debugMode ){ KNUERR << "No hit on layer 7/8" << endl; } return false; }
  if( !(xside[0] && yside[0]) ) { if( debugMode ){ KNUERR << "No L1X & L1Y hit" << endl; } return false; }
  float l1inner_chi2y = thisTrack->GetNormChisqY( id_l1inner );
  float inner_chi2y   = thisTrack->GetNormChisqY( id_inner );
  if( !( l1inner_chi2y < 10 && l1inner_chi2y - inner_chi2y < 10 ) ) { if( debugMode ){ KNUERR << "Chisquare in condition does not fitted." << endl; } return false; }

  /*
  bool hitOnLayerJ[9];
  for(int i = 0; i < 9; i++) hitOnLayerJ[i] = thisTrack->TestHitBitsJ(i, id_maxspan);

  if( !(hitOnLayerJ[2] || hitOnLayerJ[3]) ) { if( debugMode ){cerr << "No hit on layer 2/3 " << endl;} return false; }
  if( !(hitOnLayerJ[4] || hitOnLayerJ[5]) ) { if( debugMode ){cerr << "No hit on layer 4/5 " << endl;} return false; }
  if( !(hitOnLayerJ[6] || hitOnLayerJ[7]) ) { if( debugMode ){cerr << "No hit on layer 6/7 " << endl;} return false; }
  */

  // Reject events with rigidity less than 0.5 GV
  if( thisTrack->GetRigidity(id_l1inner) < 0.05 ) return false;

  return true;
}/*}}}*/

