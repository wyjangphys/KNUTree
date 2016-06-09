#include "KNUTree.h"

/**
 * @brief This function returns hardware status.
 * @return true : Hardware status is good / false : Hardware status is not good
 */
bool KNUTree::IsHardwareStatusGood(AMSEventR* thisEvent)
{/*{{{*/
  int nDaqEvent = thisEvent->nDaqEvent();
  bool goodHW = true;
  bool error = false;

  for(int i = 0; i < nDaqEvent; i++)
  {
    DaqEventR* daqEvent = thisEvent->pDaqEvent(i);
    for( int iJINJ = 0; iJINJ < 4; iJINJ++ ) error |= (bool)( ( daqEvent->JINJStatus[iJINJ]>>8 ) & 0x7F );
    for( int iJErr = 0; iJErr < 24;iJErr++ ) error |= (bool)(   daqEvent->JError[iJErr] & 0x7F );
    if(error) goodHW &= false;
  }

  if( goodHW == false ) KNUERR << "Event: " << thisEvent->Event() << " has bad HW status." << endl;
  return goodHW;
}/*}}}*/

