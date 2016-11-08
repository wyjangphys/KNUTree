#include "KNUTree.h"

/**
 * @brief  This function checks two aspects from TOF measured information.
 *
 *         First is hit configuration. BetaR::Pattern variable contains 
 *         information for the configuration of TOF hit pattern.
 *         If BetaR::Pattern is greater than 5, it means there are at least 
 *         three hits on TOF layers out of 4 layers.
 *
 *         Second is it's velocity. We restricted our interest of particle velocity to the range of 0.4 < beta <~ 1.
 *         So we discarded slow particles to ensure the accuracy of measurement.
 *
 * @return true : The BetaR object has good pattern and the speed of the 
 *              particle is enough to be analyzed.
 *         false : The BetaR object doesn't have good hit pattern or the 
 *              speed of the particle is not enough to be analyzed.
 */
bool KNUTree::IsGoodBeta(BetaHR* thisBeta)
{
  bool debugMode = false;
  if( !thisBeta ) return false;

  if( thisBeta->GetBeta() < 0.4 )
  {
    if( debugMode ) KNUERR << "This event has beta less than 0.4." << endl;
    return false;
  }

  /*
  if( thisBeta->GetBuildType() != 1 )
  {
    if( debugMode ) KNUERR << "This event doesn't have BuildType = 1." << endl;
    return false;
  }
  */

  return true;
}

