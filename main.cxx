/**
 * @file main.cxx
 * @brief main function definition file.
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

  KNUTree knuTree("KNUTree", argc, argv);
  knuTree.Begin();
  knuTree.Init();
  //knuTree.Loop();
  knuTree.LoopMultiPtl();
  knuTree.End();

  return 0;
}
