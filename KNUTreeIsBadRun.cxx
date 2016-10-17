#include "KNUTree.h"

/**
 * @brief   This function checks whether current event contained in bad-run.
 * @return  true : The event is contained in bad-run. (Therefore, it should be rejected.) / false : The event is good to analyze further.
 */
bool KNUTree::IsBadRun(AMSEventR* event)
{/*{{{*/
  unsigned int runNumber = (unsigned int)event->fHeader.Run;
  if( runNumber == 1306219312 || runNumber == 1306219522 || runNumber == 1306233745 )
  {
    KNUERR << "The program is aborted since the current run: " << runNumber << " is a bad run(1)." << std::endl;
    KNUERR << "Reason : The run is either 1306219312 or 1306219522 or 1306233745" << std::endl;
    return true;
  }
  else if( runNumber >= 1307125541 && runNumber <= 1307218054 )
  {
    KNUERR << "The program is aborted since the current run: " << runNumber << " is a bad run(2)." << std::endl;
    KNUERR << "Reason: The run is contained in the region from 1307125541 to 1307218054" << std::endl;
    return true;
  }
  else if( runNumber == 1321198167 )
  {
    KNUERR << "The program is aborted since the current run: " << runNumber << " is a bad run(3)." << std::endl;
    KNUERR << "Reason: The run 1321198167 is a bad run." << std::endl;
    return true;
  }
  else if( event->isBadRun( runNumber ) )
  {
    KNUERR << "The program is aborted since the current run: " << runNumber << " is a bad run(4)." << std::endl;
    KNUERR << "REason: The run: " << runNumber << " is determined as a bad run by AMSEventR::isBadRun() function." << std::endl;
    return true;
  }
  else return false;
}/*}}}*/
