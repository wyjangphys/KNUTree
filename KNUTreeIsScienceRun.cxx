#include "KNUTree.h"

/**
 * @brief This function checks whether the run containing input event is science-run or not.
 * @return true : The event is a science-run / false : The run containing input event is not a science-run
 */
bool KNUTree::IsScienceRun(AMSEventR* thisEvent)
{/*{{{*/
  HeaderR* header = &(thisEvent->fHeader);
  if( !header ) return false;
  if( ( header->RunType >> 12 ) != 0xf ) return false;
  return true;
}/*}}}*/

bool KNUTree::IsScienceRun(AMSChain* chain)
{
  AMSEventR* event = chain->GetEvent(1);
  HeaderR* header = &(event->fHeader);
  if( !header )
  {
    cout << "[KNUTree::IsScienceRun] NO HEADER" << endl;
    return false;
  }
  if( ( header->RunType >> 12 ) != 0xf )
  {
    cout << "[KNUTree::IsScienceRun]" << (header->RunType>>12) << endl;
    return false;
  }
  return true;
}
