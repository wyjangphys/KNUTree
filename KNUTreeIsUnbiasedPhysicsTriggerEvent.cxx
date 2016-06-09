#include "KNUTree.h"

/**
 * @brief This function checks whether the event triggered by unbiased physics trigger or not.
 * @return true : The trigger is unbiased physics trigger / false : The trigger is physics trigger
 */
bool KNUTree::IsUnbiasedPhysicsTriggerEvent(AMSEventR* thisEvent)
{/*{{{*/
  int nLevel1 = thisEvent->nLevel1();
  bool unbiased = false;

  for(int i = 0; i < nLevel1; i++)
  {
    Level1R* pLevel1 = thisEvent->pLevel1(i);
    bitset<8> physicsBitPattern(pLevel1->PhysBPatt);

    if(physicsBitPattern.test(1) || physicsBitPattern.test(2) || physicsBitPattern.test(3) || physicsBitPattern.test(4) || physicsBitPattern.test(5))
      unbiased |= false;
    else
      unbiased |= true;
  }

  return unbiased;
}/*}}}*/

