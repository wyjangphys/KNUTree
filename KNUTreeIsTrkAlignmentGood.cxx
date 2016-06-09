#include "KNUTree.h"

/**
 * @brief
 * @return
 */
bool KNUTree::IsTrkAlignmentGood(AMSEventR* thisEvent)
{/*{{{*/
  AMSPoint pn1, pn9, pd1, pd9;
  thisEvent->GetRTIdL1L9(0, pn1, pd1, thisEvent->UTime(), 60);
  thisEvent->GetRTIdL1L9(1, pn9, pd9, thisEvent->UTime(), 60);
  if(pd1.y() > 35 || pd9.y() > 45)
    return false;
  else
    return true;
}/*}}}*/

