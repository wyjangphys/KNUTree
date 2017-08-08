/**
 * @file main.cxx
 * @brief Define main function here.
 * @author Wooyoung Jang, wyjang@knu.ac.kr
 * @date 2016.04.15
 *
 * This source code defining the main function of the KNU standard analysis template.
 */

#include <iostream>
#include <string>
#include <cstdio>
#include <ctime>
#include <root.h>

#include "KNUTree.h"

const std::string CurrentDateTime()
{
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}

int main(int argc, char* argv[])
{
  KNUOUT << "This is KNU standard AMS analysis template code." << std::endl;

  TStopwatch clock;
  clock.Start();
  KNUOUT << "TStopwatch is started." << std::endl;
  KNUOUT << "Current time : " << CurrentDateTime() << std::endl;

  KNUTree knuTree("KNUTree", argc, argv);
  knuTree.Begin();
  knuTree.Init();
  //knuTree.Loop();
  //knuTree.LoopMultiPtl();
  if( knuTree.isMC == true )
  {
    KNUOUT << " !!!!!!!!! THIS IS MC RUN !!!!!!!!!!!!!!!!!! " << endl;
    knuTree.ProcessMC();
  }
  else
  {
    KNUOUT << " !!!!!!!! This is ISS run!!!!!!" << endl;
    knuTree.Process();
  }
  knuTree.End();

  clock.Stop();
  KNUOUT << "TStopwatch report : " << std::endl;
  KNUOUT << "Current time : " << CurrentDateTime() << std::endl;
  clock.Print();
  return 0;
}
