#include "KNUTree.h"

/**
 * @brief This function checks whether the number of ACC hits exceed 5 or not.
 * @return true : Number of ACC hits smaller than or equal to 5. false : Number of ACC hits greater than 5.
 */
bool KNUTree::IsACCPatternGood(AMSEventR* thisEvent)
{/*{{{*/
  int nACCHit = 0;
  int nLevel1 = thisEvent->nLevel1();

  for(int i = 0; i < nLevel1; i++)
  {
    Level1R* pLevel1 = thisEvent->pLevel1(i);
    for( int j = 0; j < 8; j++)
    {
      if(((pLevel1->AntiPatt>>i)&1)==1) nACCHit++;
    }
  }
  if( nACCHit > 5 ) return false;
  return true;
}/*}}}*/
