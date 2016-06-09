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

#include "amschain.h"

#include <TROOT.h>
#include <TUnixSystem.h>
#include <TChain.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include "KNUTree.h"

int main(int argc, char* argv[])
{
  KNUOUT << "This is KNU standard AMS analysis template program." << std::endl;
  KNUERR << "Test for console error out." << std::endl;

  KNUTree knuTree("KNUTree", argc, argv);
  knuTree.Begin();
  knuTree.Init();
  knuTree.Loop();
  knuTree.End();
}
