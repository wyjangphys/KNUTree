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
