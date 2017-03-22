/**
 * @file main.cxx
 * @brief Define main function here.
 * @author Wooyoung Jang, wyjang@knu.ac.kr
 * @date 2016.04.15
 *
 * This source code defining the main function of the KNU standard analysis template.
 */

#include <iostream>
#include <ctime>
#include <root.h>

#include "KNUTree.h"

int main(int argc, char* argv[])
{
  KNUOUT << "This is KNU standard AMS analysis template code." << std::endl;

  TStopwatch clock;
  clock.Start();
  KNUOUT << "TStopwatch is started." << std::endl;

  KNUTree knuTree("KNUTree", argc, argv);
  knuTree.Begin();
  knuTree.Init();
  //knuTree.Loop();
  //knuTree.LoopMultiPtl();
  if( knuTree.isMC == true )
  {
    cout << " !!!!!!!!! THIS IS MC RUN !!!!!!!!!!!!!!!!!! " << endl;
    knuTree.ProcessMC();
  }
  else
  {
    cout << " !!!!!!!! This is ISS run!!!!!!" << endl;
    knuTree.Process();
  }
  knuTree.End();

  clock.Stop();
  KNUOUT << "TStopwatch report : " << std::endl;
  clock.Print();
  return 0;
}
