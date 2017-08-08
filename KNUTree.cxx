#include "KNUTree.h"

#include <iostream>
#include <cstring>

#include "FileManager.hh"
#include "AMSRootSupport.hh"
#include "AMSRootParticleFactory.hh"
#include "Utilities.hh"
#include "SlowControlLookup.hh"
#include "TrdKCluster.h"
#include "amschain.h"

/**
 * @brief Basic constructor
 */
KNUTree::KNUTree(char* _name, int argc, char* argv[])
{/*{{{*/
  SetName(_name);
  if( argc == 1 ){ outputFileName = new char[ strlen("testrun.root")+1]; strncpy(outputFileName, "testrun.root", strlen("testrun.root")+1); }
  else if( argc == 2 ){ outputFileName = new char[ strlen("testrun.root")+1]; strncpy(outputFileName, "testrun.root", strlen("testrun.root")+1); }
  else if( argc == 3 ){ outputFileName = new char[ strlen(argv[2])+1]; strncpy(outputFileName, argv[2], strlen(argv[2])+1); }
  else if( argc == 4 ){ outputFileName = new char[ strlen(argv[2])+1]; strncpy(outputFileName, argv[2], strlen(argv[2])+1); }
  nTreesInChain = 0;
  pChain = new AMSChain();
  RegisterFileToChain(argc, argv, pChain);

  nEvents             = pChain->GetEntries();
  nCut                = 0;
  nExaminedParticles  = 0;
  nSelected           = 0;
  nSelectedParticles  = 0;
  nRun                = 0;
  nEvent              = 0;
  nLevel1             = 0;
  nParticle           = 0;
  nCharge             = 0;
  nTrTrack            = 0;
  nTrdTrack           = 0;
  nAntiCluster        = 0;
  nTofClustersInTime  = 0;
  nRichRing           = 0;
  nRichRingB          = 0;
  nBeta               = 0;
  nBetaB              = 0;
  nBetaH              = 0;
  nEcalShower         = 0;
  nVertex             = 0;
  particleID          = 0;
  trkID               = 0;
  parentID            = 0;
  for(int j = 0; j < 3; j++ )
  {
    genPosition[j] = 0;
    genDirCosine[j] = 0;
  }
  genMass             = 0;
  genCharge           = 0;
  genMomentum         = 0;
  particleType        = 0;
  liveTime            = 0;
  utcTime             = 0;
  orbitAltitude       = 0;
  orbitLatitude       = 0;
  orbitLongitude      = 0;
  orbitLatitudeM      = 0;
  orbitLongitudeM     = 0;
  velR                = 0;
  velTheta            = 0;
  velPhi              = 0;
  yaw                 = 0;
  pitch               = 0;
  roll                = 0;
  gLongitude          = 0;
  gLatitude           = 0;
  gCoordCalcResult    = 0;
  gCoordCalcResultBit = 0;
  sunPosAzimuth       = 0;
  sunPosElevation     = 0;
  sunPosCalcResult    = 0;
  unixTime            = 0;
  acEventTime         = 0.;
  solarArrayCoord[0]  = 0;
  solarArrayCoord[1]  = 0;
  solarArrayCoord[2]  = 0;
  isInShadow          = 0;
  zenithAngle         = 0;
  isInSAA             = 0;
  //ptlCharge = 0;
  /*
  ptlMomentum = 0;
  ptlTheta = 0;
  ptlPhi = 0;
  ptlCoo[0]           = 0;
  ptlCoo[1]           = 0;
  ptlCoo[2]           = 0;
  ptlTofCoo[0][0]     = 0;
  ptlTofCoo[0][1]     = 0;
  ptlTofCoo[0][2]     = 0;
  ptlTofCoo[1][0]     = 0;
  ptlTofCoo[1][1]     = 0;
  ptlTofCoo[1][2]     = 0;
  ptlTofCoo[2][0]     = 0;
  ptlTofCoo[2][1]     = 0;
  ptlTofCoo[2][2]     = 0;
  ptlTofCoo[3][0]     = 0;
  ptlTofCoo[3][1]     = 0;
  ptlTofCoo[3][2]     = 0;
  ptlTofTrLength[0]   = 0;
  ptlTofTrLength[1]   = 0;
  ptlTofTrLength[2]   = 0;
  ptlTofTrLength[3]   = 0;
  ptlAntiCoo[0][0]    = 0;
  ptlAntiCoo[0][1]    = 0;
  ptlAntiCoo[0][2]    = 0;
  ptlAntiCoo[0][3]    = 0;
  ptlAntiCoo[0][4]    = 0;
  ptlAntiCoo[1][0]    = 0;
  ptlAntiCoo[1][1]    = 0;
  ptlAntiCoo[1][2]    = 0;
  ptlAntiCoo[1][3]    = 0;
  ptlAntiCoo[1][4]    = 0;
  ptlEcalCoo[0][0]    = 0;
  ptlEcalCoo[0][1]    = 0;
  ptlEcalCoo[0][2]    = 0;
  ptlEcalCoo[1][0]    = 0;
  ptlEcalCoo[1][1]    = 0;
  ptlEcalCoo[1][2]    = 0;
  ptlEcalCoo[2][0]    = 0;
  ptlEcalCoo[2][1]    = 0;
  ptlEcalCoo[2][2]    = 0;
  ptlTrCoo[0][0]      = 0;
  ptlTrCoo[0][1]      = 0;
  ptlTrCoo[0][2]      = 0;
  ptlTrCoo[1][0]      = 0;
  ptlTrCoo[1][1]      = 0;
  ptlTrCoo[1][2]      = 0;
  ptlTrCoo[2][0]      = 0;
  ptlTrCoo[2][1]      = 0;
  ptlTrCoo[2][2]      = 0;
  ptlTrCoo[3][0]      = 0;
  ptlTrCoo[3][1]      = 0;
  ptlTrCoo[3][2]      = 0;
  ptlTrCoo[4][0]      = 0;
  ptlTrCoo[4][1]      = 0;
  ptlTrCoo[4][2]      = 0;
  ptlTrCoo[5][0]      = 0;
  ptlTrCoo[5][1]      = 0;
  ptlTrCoo[5][2]      = 0;
  ptlTrCoo[6][0]      = 0;
  ptlTrCoo[6][1]      = 0;
  ptlTrCoo[6][2]      = 0;
  ptlTrCoo[7][0]      = 0;
  ptlTrCoo[7][1]      = 0;
  ptlTrCoo[7][2]      = 0;
  ptlTrCoo[8][0]      = 0;
  ptlTrCoo[8][1]      = 0;
  ptlTrCoo[8][2]      = 0;
  ptlTrdCoo[0][0]     = 0;
  ptlTrdCoo[0][1]     = 0;
  ptlTrdCoo[0][2]     = 0;
  ptlTrdCoo[1][0]     = 0;
  ptlTrdCoo[1][1]     = 0;
  ptlTrdCoo[1][2]     = 0;
  ptlRichCoo[0][0]    = 0;
  ptlRichCoo[0][1]    = 0;
  ptlRichCoo[0][2]    = 0;
  ptlRichCoo[1][0]    = 0;
  ptlRichCoo[1][1]    = 0;
  ptlRichCoo[1][2]    = 0;
  ptlRichPath[0]      = 0;
  ptlRichPath[1]      = 0;
  ptlRichPathBeta[0]  = 0;
  ptlRichPathBeta[1]  = 0; 
  ptlRichLength       = 0;
  ptlRichParticles    = 0;
  ptlCutOffStoermer   = 0;
  ptlCutOffDipole     = 0;
  */
  //ptlCutOffMax[0] = 0;
  //ptlCutOffMax[1] = 0;
  isEcalAvailable       = 0;
  showerEnergyD         = 0;
  showerEnergyDL[0]     = 0;
  showerEnergyDL[1]     = 0;
  showerEnergyDL[2]     = 0;
  showerEnergyDL[3]     = 0;
  showerEnergyDL[4]     = 0;
  showerEnergyDL[5]     = 0;
  showerEnergyDL[6]     = 0;
  showerEnergyDL[7]     = 0;
  showerEnergyDL[8]     = 0;
  showerEnergyDL[9]     = 0;
  showerEnergyDL[10]    = 0;
  showerEnergyDL[11]    = 0;
  showerEnergyDL[12]    = 0;
  showerEnergyDL[13]    = 0;
  showerEnergyDL[14]    = 0;
  showerEnergyDL[15]    = 0;
  showerEnergyDL[16]    = 0;
  showerEnergyDL[17]    = 0;
  showerEnergyE         = 0;
  showerEnergyCorrected = 0;
  showerBDT             = 0;
  showerCofG[0]         = 0;
  showerCofG[1]         = 0;
  showerCofG[2]         = 0;
  showerCofGDist        = 0;
  showerCofGdX          = 0;
  showerCofGdY          = 0;
  tofNCluster           = 0;
  tofNClusterH          = 0;
  tofNUsedHits          = 0;
  tofNUnusedHits        = 0;
  tofNUsedLayersForQ    = 0;
  tofBeta               = 0;
  tofInvBetaErr         = 0;
  tofMass               = 0;
  tofMassError          = 0;
  isGoodBeta            = 0;
  isTkTofMatch          = 0;
  tofReducedChisqC      = 0;
  tofReducedChisqT      = 0;
  tofDepositedEnergyOnLayer[0] = 0;
  tofDepositedEnergyOnLayer[1] = 0;
  tofDepositedEnergyOnLayer[2] = 0;
  tofDepositedEnergyOnLayer[3] = 0;
  tofEstimatedChargeOnLayer[0] = 0;
  tofEstimatedChargeOnLayer[1] = 0;
  tofEstimatedChargeOnLayer[2] = 0;
  tofEstimatedChargeOnLayer[3] = 0;
  tofCharge = 0;
  tofChargeRMS = 0;
  tofUpperCharge = 0;
  tofLowerCharge = 0;
  tofAcCharge = 0;
  tofAcUpperCharge = 0;
  tofAcLowerCharge = 0;
  tofChargeOnLayer[0] = 0;
  tofChargeOnLayer[1] = 0;
  tofChargeOnLayer[2] = 0;
  tofChargeOnLayer[3] = 0;
  trkFitCodeL1Inner = -1;
  trkRigidityL1Inner = -9999;
  trkRigidityInverseErrorL1Inner = -9999;
  trkReducedChisquareL1InnerX = -9999;
  trkReducedChisquareL1InnerY = -9999;
  trkFitCodeL9Inner = -1;
  trkRigidityL9Inner = -9999;
  trkRigidityInverseErrorL9Inner = -9999;
  trkReducedChisquareL9InnerX = -9999;
  trkReducedChisquareL9InnerY = -9999;
  trkFitCodeUpperInner = -1;
  trkRigidityUpperInner = -9999;
  trkRigidityInverseErrorUpperInner = -9999;
  trkReducedChisquareUpperInnerX = -9999;
  trkReducedChisquareUpperInnerY = -9999;
  trkFitCodeLowerInner = -1;
  trkRigidityLowerInner = -9999;
  trkRigidityInverseErrorLowerInner = -9999;
  trkReducedChisquareLowerInnerX = -9999;
  trkReducedChisquareLowerInnerY = -9999;
  trkFitCodeFS = -1;
  trkRigidityFS = -9999;
  trkRigidityInverseErrorFS = -9999;
  trkReducedChisquareFSX = -9999;
  trkReducedChisquareFSY = -9999;
  trkFitCodeMS = -1;
  trkRigidityMS = -9999;
  trkRigidityInverseErrorMS = -9999;
  trkReducedChisquareMSY = -9999;
  trkReducedChisquareMSX = -9999;
  trkFitCodeInner = -1;
  trkRigidityInverseErrorInner = -9999;
  trkReducedChisquareInnerX = -9999;
  trkReducedChisquareInnerY = -9999;
  trkLayerJQ[0] = 0;
  trkLayerJQ[1] = 0;
  trkLayerJQ[2] = 0;
  trkLayerJQ[3] = 0;
  trkLayerJQ[4] = 0;
  trkLayerJQ[5] = 0;
  trkLayerJQ[6] = 0;
  trkLayerJQ[7] = 0;
  trkLayerJQ[8] = 0;
  trkEdepLayerJXSideOK[0] = 0;
  trkEdepLayerJXSideOK[1] = 0;
  trkEdepLayerJXSideOK[2] = 0;
  trkEdepLayerJXSideOK[3] = 0;
  trkEdepLayerJXSideOK[4] = 0;
  trkEdepLayerJXSideOK[5] = 0;
  trkEdepLayerJXSideOK[6] = 0;
  trkEdepLayerJXSideOK[7] = 0;
  trkEdepLayerJXSideOK[8] = 0;
  trkEdepLayerJYSideOK[0] = 0;
  trkEdepLayerJYSideOK[1] = 0;
  trkEdepLayerJYSideOK[2] = 0;
  trkEdepLayerJYSideOK[3] = 0;
  trkEdepLayerJYSideOK[4] = 0;
  trkEdepLayerJYSideOK[5] = 0;
  trkEdepLayerJYSideOK[6] = 0;
  trkEdepLayerJYSideOK[7] = 0;
  trkEdepLayerJYSideOK[8] = 0;
  trkEdepLayerJ[0] = 0;
  trkEdepLayerJ[1] = 0;
  trkEdepLayerJ[2] = 0;
  trkEdepLayerJ[3] = 0;
  trkEdepLayerJ[4] = 0;
  trkEdepLayerJ[5] = 0;
  trkEdepLayerJ[6] = 0;
  trkEdepLayerJ[7] = 0;
  trkEdepLayerJ[8] = 0;
  for( int i = 0; i < 9; i++ )
  {
    trkPosXLJ[i] = 0;
    trkPosYLJ[i] = 0;
    trkPosZLJ[i] = 0;
    trkDirThetaLJ[i] = 0;
    trkDirPhiLJ[i] = 0;
  }
  trkCharge = 0;
  trkInnerCharge = 0;
  trkHasExtLayers = 0;
  richRebuild = 0;
  richIsGood = 0;
  richIsClean = 0;
  richIsNaF = 0;
  richUsedHits = 0;
  richRingWidth = 0;
  richNHits = 0;
  richNPMTsOnRing = 0;
  richBeta = 0;
  richBetaError = 0;
  richChargeSquared = 0;
  richKolmogorovProbability = 0;
  richPhotoelectrons = 0;
  richPhotoElectrons = 0;
  richExpectedPhotoelectrons = 0;
  richExpectedPhotoElectrons = 0;
  richNUsedHits = 0;
  richTheta = 0;
  richPhi = 0;
  trdNClusters = 0;
  trdNUnusedHits = 0;
  trdNTracks = 0;
  trdNVertex = 0;
  trdTrackPhi = 0;
  trdTrackTheta = 0;
  trdTrackChi2 = 0;
  trdTrackPattern = 0;
  trdTrackCharge = 0;
  trdTrackEdepL[0] = 0;
  trdTrackEdepL[1] = 0;
  trdTrackEdepL[2] = 0;
  trdTrackEdepL[3] = 0;
  trdTrackEdepL[4] = 0;
  trdTrackEdepL[5] = 0;
  trdTrackEdepL[6] = 0;
  trdTrackEdepL[7] = 0;
  trdTrackEdepL[8] = 0;
  trdTrackEdepL[9] = 0;
  trdTrackEdepL[10] = 0;
  trdTrackEdepL[11] = 0;
  trdTrackEdepL[12] = 0;
  trdTrackEdepL[13] = 0;
  trdTrackEdepL[14] = 0;
  trdTrackEdepL[15] = 0;
  trdTrackEdepL[16] = 0;
  trdTrackEdepL[17] = 0;
  trdTrackEdepL[18] = 0;
  trdTrackEdepL[19] = 0;
  trdTrackMeanDepositedEnergy = 0;
  trdTrackTotalDepositedEnergy = 0;
  trdQtNActiveLayer = 0;
  trdQtIsCalibrationGood = 0;
  trdQtIsInsideTrdGeometricalAcceptance = 0;
  trdQtIsValid = 0;
  trdQtNActiveStraws = 0;
  trdQtNActiveLayers = 0;
  trdQtNTRDVertex = 0;
  trdQtElectronToProtonLogLikelihoodRatio = 0;
  trdQtHeliumToElectronLogLikelihoodRatio = 0;
  trdQtHeliumToProtonLogLikelihoodRatio = 0;
  trdKNRawHits = 0;
  trdKIsReadAlignmentOK = 0;
  trdKIsReadCalibOK = 0;
  trdKNHits = 0;
  trdKIsValid = 0;
  trdKElectronToProtonLogLikelihoodRatio = 0;
  trdKHeliumToElectronLogLikelihoodRatio = 0;
  trdKHeliumToProtonLogLikelihoodRatio = 0;
  trdKCharge = 0;
  trdKChargeError = 0;
  trdKNUsedHitsForCharge = 0;
  trdKAmpLayer[0] = 0;
  trdKAmpLayer[1] = 0;
  trdKAmpLayer[2] = 0;
  trdKAmpLayer[3] = 0;
  trdKAmpLayer[4] = 0;
  trdKAmpLayer[5] = 0;
  trdKAmpLayer[6] = 0;
  trdKAmpLayer[7] = 0;
  trdKAmpLayer[8] = 0;
  trdKAmpLayer[9] = 0;
  trdKAmpLayer[10] = 0;
  trdKAmpLayer[11] = 0;
  trdKAmpLayer[12] = 0;
  trdKAmpLayer[13] = 0;
  trdKAmpLayer[14] = 0;
  trdKAmpLayer[15] = 0;
  trdKAmpLayer[16] = 0;
  trdKAmpLayer[17] = 0;
  trdKAmpLayer[18] = 0;
  trdKAmpLayer[19] = 0;
  trdKTotalPathLength = 0;
  trdKTotalAmp = 0;
  trdKElectronLikelihood = 0;
  trdKHeliumLikelihood = 0;
  trdPElectronToProtonLogLikelihoodRatio = 0;
  trdPHeliumToElectronLogLikelihoodRatio = 0;
  trdPHeliumToProtonLogLikelihoodRatio = 0;
  trdPElectronLikelihood = 0;
  trdPProtonLikelihood = 0;
  trdPHeliumLikelihood = 0;
  trdPNHybridHits = 0;
  trdPNActiveLayers = 0;
  trdPNLowAdcHits = 0;
  trdPNLowDxHits = 0;
  acTrackerTrdTrackDistanceAtTrdCenterInSigmas = 0;
  acTrackerTrdTrackAngularDistanceInSigmas = 0;
}/*}}}*/

/**
 * @brief Destructor
 */
KNUTree::~KNUTree()
{/*{{{*/
  //delete[] name;
  //delete pChain;
  //delete hEvtCounter;     // <<< This line causes segmentation violation. Don't know why.
  //delete pOutputTFile;
  //delete[] outputFileName;// <<< This line causes seg. viol?
  //delete outTree;         // <<< This line cuases segmentation violation. Don't know why.
}/*}}}*/

/**
 * @brief Set the name of the object
 * @return true: when success to set the name / false: when failed to set the name
 */
bool KNUTree::SetName(char* _name)
{/*{{{*/
  this->name = new char[strlen(_name)+1];
  if( strncpy(this->name, _name, strlen(_name)+1) )
    return true;
  else
    return false;
}/*}}}*/

/**
 * @brief Return the name of the object
 * @return pointer to the name object
 */
const char* KNUTree::GetName()
{/*{{{*/
  return name;
}/*}}}*/

/**
 * @brief 
 * @return
 */
bool KNUTree::RegisterFileToChain(int argc, char* argv[], AMSChain* amsChain)
{/*{{{*/
  /*************************************************************************************************************
   *
   * FILE OPENING
   *
   *************************************************************************************************************/

  /*
   * 1. Test mode
   *   1-1. Skirmish test
   *     Use the designated run to test the code. No arguments are required. ( argc == 1)
   *   1-2. Skirmish test within given number of events
   *     User should pass the number of events to be analyzed as the program argument.
   *     ./a.out <number of events to be analyzed>  ( argc == 2 )
   *   1-3. Test using user-given run file within given number of events
   *     ./a.out <run file to be analyzed> <output file> <number of events to be analyzed> ( argc == 4 )
   *
   * 2. Job submission mode
   *   2-1. Standard job submission mode
   *     ./a.out <list file which contains list of files to be analyzed> <output file> ( argc == 3 )
   */

  //char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1366261348.00000001.root";  // Bad run example
  //char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1362134524.00000001.root";
  //char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1308865423.00000001.root";
  //char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1307005643.00000001.root";   // Path of test run
  //char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1323051106.00000001.root";   // Path of test run
  //char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/MC/AMS02/2014/d.B1030/d.pl1.0_520_GG_Blic/872728226.00000001.root";
  char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1319745123.00000001.root";
  char inputFileName[256];      // File path for single run. (Cat 3.)

  if( argc == 1 )
  {
    KNUOUT << "RUN MODE : Single Test Run (Cat. 1)" << std::endl;
    if( skirmishRunPath[31] == 'M' )
    {
        isMC = true;
        KNUOUT << "==== THE PROGRAM IS RUNNING AS MONTE CARLO MODE ====" << std::endl;
    }
    else
    {
        isMC = false;
        KNUOUT << "==== THE PROGRAM IS RUNNING AS ISS DATA PROCESSING MODE ====" << std::endl;
    }

    if( amsChain->Add(skirmishRunPath) != 1 )
    {
      KNUERR << " ERROR    : File open error, [" << skirmishRunPath << "] can not be found!" << std::endl;
      return false;
    }
    else
    {
      KNUOUT << "File [" << skirmishRunPath << "] opened successfully." << std::endl;
      // The code below inserted temporarily for test.
      // === Begin ===
      //amsChain->Add("root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1307005643.01501575.root");
      // === END ===
      return true;
    }
  }
  else if( argc == 2 ) // Process as many events as user requested. (second argument is number of events to be processed)
  {
    KNUOUT << "RUN MODE : Single Test Run (Cat. 2)" << std::endl;

    if(amsChain->Add(skirmishRunPath) != 1)
    {
      KNUERR << "ERROR    : File open error, [" << skirmishRunPath << "] can not be found!" << std::endl;
      return false;
    }
    else
    {
      KNUOUT << "File [" << skirmishRunPath << "] opened successfully." << std::endl;
      return true;
    }
  }
  else if( argc == 3 )
  {
    KNUOUT << "RUN MODE : Batch-job Mode" << std::endl;

    char listFileName[256];       // The name of the list file.
    char inputFileName[256];

    FILE* fp;                     // File pointer to read list file.

    strcpy(listFileName, argv[1]);     // First argument is list file

    if( ( fp = fopen( listFileName, "r") ) == NULL )
    {
      KNUERR << "ERROR   : Failed to open file [" << listFileName << "]!" << std::endl;
      return false;
    }

    char* line_p;                 // Character pointer to filter out CRLF at the end of each line.
    while( fgets( inputFileName, 256, fp ) != NULL )
    {
      if( inputFileName[31] == 'M' ) // If the 32-th character of the file path is 'M' then it is MC run.
      {
        isMC = true;
        KNUOUT << "==== THE PROGRAM IS RUNNING AS MONTE CARLO MODE ====" << std::endl;
      }
      else
      {
        isMC = false;
        KNUOUT << "==== THE PROGRAM IS RUNNING AS ISS DATA PROCESSING MODE ====" << std::endl;
      }
      if( ( line_p = strchr(inputFileName, '\n') ) != NULL ) // Search \n
        *line_p = 0;  // Filter out \n at the end of lines.

      if(amsChain->Add( inputFileName ) != 1)
      {
        KNUERR << "ERROR     : Failed to open file [" << inputFileName << "]!" << std::endl;
        return false;
      }
      else
      {
        KNUOUT << "The file [" << inputFileName << "] is added to the chain." << std::endl;
        KNUOUT << "Currently loaded events : " << amsChain->GetEntries() << std::endl;
        nTreesInChain++;
      }
    }

    fclose(fp);
    return true;
  }
  else if( argc >= 4 )
  {
    KNUOUT << "RUN MODE : Single Test Run (Cat. 3)" << std::endl;

    strcpy(inputFileName, argv[1]);

    if(amsChain->Add(inputFileName) != 1)
    {
      KNUERR << "ERROR   : File open error, [" << inputFileName << "] can not found!" << std::endl;
      return false;
    }
    else
    {
      KNUOUT << "The file [" << inputFileName << "] is added to the chain." << std::endl;
      return true;
    }
  }

  return false;
}/*}}}*/

/**
 * @brief  This function initializes all variables relavant to this analysis.
 * @return
 */
void KNUTree::Init()
{/*{{{*/
  KNUOUT << "Initialize...";

  nEvents             = pChain->GetEntries();
  nCut                = 0;
  nExaminedParticles  = 0;
  nSelected           = 0;
  nSelectedParticles  = 0;
  nRun                = 0;
  nEvent              = 0;
  nLevel1             = 0;
  nParticle           = 0;
  nCharge             = 0;
  nTrTrack            = 0;
  nTrdTrack           = 0;
  nAntiCluster        = 0;
  nTofClustersInTime  = 0;
  nRichRing           = 0;
  nRichRingB          = 0;
  nBeta               = 0;
  nBetaB              = 0;
  nBetaH              = 0;
  nEcalShower         = 0;
  nVertex             = 0;
  survived_tof        = 0;
  nAntiMCCluster      = 0;
  nTrMCCluster        = 0;
  nTofMCCluster       = 0;
  nTofMCPmtHit        = 0;
  nEcalMCHit          = 0;
  nTrdMCCluster       = 0;
  nRichMCCluster      = 0;
  nMCTrack            = 0;
  nMCEventg           = 0;
  particleID          = 0;
  trkID               = 0;
  parentID            = 0;
  genPosition[0]      = 0;
  genPosition[1]      = 0;
  genPosition[2]      = 0;
  genDirCosine[0]     = 0;
  genDirCosine[1]     = 0;
  genDirCosine[2]     = 0;
  genMass             = 0;
  genCharge           = 0;
  genMomentum         = 0;
  conversionId        = 0;
  conversionIdContainer.clear();
  conversionIdSet.clear();
  particleType        = 0;
  liveTime            = 0;
  utcTime             = 0;
  orbitAltitude       = 0;
  orbitLatitude       = 0;
  orbitLongitude      = 0;
  orbitLatitudeM      = 0;
  orbitLongitudeM     = 0;
  velR                = 0;
  velTheta            = 0;
  velPhi              = 0;
  yaw                 = 0;
  pitch               = 0;
  roll                = 0;
  gLongitude          = 0;
  gLatitude           = 0;
  gCoordCalcResult    = 0;
  gCoordCalcResultBit = 0;
  sunPosAzimuth       = 0;
  sunPosElevation     = 0;
  sunPosCalcResult    = 0;
  unixTime            = 0;
  acEventTime         = 0;
  solarArrayCoord[0]  = 0;
  solarArrayCoord[1]  = 0;
  solarArrayCoord[2]  = 0;
  isInShadow          = 0;
  zenithAngle         = 0;
  isInSAA             = 0;
  //ptlCharge = 0;
  /*
  ptlMomentum = 0;
  ptlTheta = 0;
  ptlPhi = 0;
  ptlCoo[0]           = 0;
  ptlCoo[1]           = 0;
  ptlCoo[2]           = 0;
  ptlTofCoo[0][0]     = 0;
  ptlTofCoo[0][1]     = 0;
  ptlTofCoo[0][2]     = 0;
  ptlTofCoo[1][0]     = 0;
  ptlTofCoo[1][1]     = 0;
  ptlTofCoo[1][2]     = 0;
  ptlTofCoo[2][0]     = 0;
  ptlTofCoo[2][1]     = 0;
  ptlTofCoo[2][2]     = 0;
  ptlTofCoo[3][0]     = 0;
  ptlTofCoo[3][1]     = 0;
  ptlTofCoo[3][2]     = 0;
  ptlTofTrLength[0]   = 0;
  ptlTofTrLength[1]   = 0;
  ptlTofTrLength[2]   = 0;
  ptlTofTrLength[3]   = 0;
  ptlAntiCoo[0][0]    = 0;
  ptlAntiCoo[0][1]    = 0;
  ptlAntiCoo[0][2]    = 0;
  ptlAntiCoo[0][3]    = 0;
  ptlAntiCoo[0][4]    = 0;
  ptlAntiCoo[1][0]    = 0;
  ptlAntiCoo[1][1]    = 0;
  ptlAntiCoo[1][2]    = 0;
  ptlAntiCoo[1][3]    = 0;
  ptlAntiCoo[1][4]    = 0;
  ptlEcalCoo[0][0]    = 0;
  ptlEcalCoo[0][1]    = 0;
  ptlEcalCoo[0][2]    = 0;
  ptlEcalCoo[1][0]    = 0;
  ptlEcalCoo[1][1]    = 0;
  ptlEcalCoo[1][2]    = 0;
  ptlEcalCoo[2][0]    = 0;
  ptlEcalCoo[2][1]    = 0;
  ptlEcalCoo[2][2]    = 0;
  ptlTrCoo[0][0]      = 0;
  ptlTrCoo[0][1]      = 0;
  ptlTrCoo[0][2]      = 0;
  ptlTrCoo[1][0]      = 0;
  ptlTrCoo[1][1]      = 0;
  ptlTrCoo[1][2]      = 0;
  ptlTrCoo[2][0]      = 0;
  ptlTrCoo[2][1]      = 0;
  ptlTrCoo[2][2]      = 0;
  ptlTrCoo[3][0]      = 0;
  ptlTrCoo[3][1]      = 0;
  ptlTrCoo[3][2]      = 0;
  ptlTrCoo[4][0]      = 0;
  ptlTrCoo[4][1]      = 0;
  ptlTrCoo[4][2]      = 0;
  ptlTrCoo[5][0]      = 0;
  ptlTrCoo[5][1]      = 0;
  ptlTrCoo[5][2]      = 0;
  ptlTrCoo[6][0]      = 0;
  ptlTrCoo[6][1]      = 0;
  ptlTrCoo[6][2]      = 0;
  ptlTrCoo[7][0]      = 0;
  ptlTrCoo[7][1]      = 0;
  ptlTrCoo[7][2]      = 0;
  ptlTrCoo[8][0]      = 0;
  ptlTrCoo[8][1]      = 0;
  ptlTrCoo[8][2]      = 0;
  ptlTrdCoo[0][0]     = 0;
  ptlTrdCoo[0][1]     = 0;
  ptlTrdCoo[0][2]     = 0;
  ptlTrdCoo[1][0]     = 0;
  ptlTrdCoo[1][1]     = 0;
  ptlTrdCoo[1][2]     = 0;
  ptlRichCoo[0][0]    = 0;
  ptlRichCoo[0][1]    = 0;
  ptlRichCoo[0][2]    = 0;
  ptlRichCoo[1][0]    = 0;
  ptlRichCoo[1][1]    = 0;
  ptlRichCoo[1][2]    = 0;
  ptlRichPath[0]      = 0;
  ptlRichPath[1]      = 0;
  ptlRichPathBeta[0]  = 0;
  ptlRichPathBeta[1]  = 0; 
  ptlRichLength       = 0;
  ptlRichParticles    = 0;
  ptlCutOffStoermer   = 0;
  ptlCutOffDipole     = 0;
  */
  //ptlCutOffMax[0] = 0;
  //ptlCutOffMax[1] = 0;
  isEcalAvailable       = 0;
  showerEnergyD         = 0;
  showerEnergyDL[0]     = 0;
  showerEnergyDL[1]     = 0;
  showerEnergyDL[2]     = 0;
  showerEnergyDL[3]     = 0;
  showerEnergyDL[4]     = 0;
  showerEnergyDL[5]     = 0;
  showerEnergyDL[6]     = 0;
  showerEnergyDL[7]     = 0;
  showerEnergyDL[8]     = 0;
  showerEnergyDL[9]     = 0;
  showerEnergyDL[10]    = 0;
  showerEnergyDL[11]    = 0;
  showerEnergyDL[12]    = 0;
  showerEnergyDL[13]    = 0;
  showerEnergyDL[14]    = 0;
  showerEnergyDL[15]    = 0;
  showerEnergyDL[16]    = 0;
  showerEnergyDL[17]    = 0;
  showerEnergyE         = 0;
  showerEnergyCorrected = 0;
  showerBDT = 0;
  showerCofG[0] = 0;
  showerCofG[1] = 0;
  showerCofG[2] = 0;
  showerCofGDist = 0;
  showerCofGdX = 0;
  showerCofGdY = 0;
  tofNCluster = 0;
  tofNClusterH = 0;
  tofNUsedHits = 0;
  tofNUnusedHits = 0;
  tofNUsedLayersForQ = 0;
  tofBeta = 0;
  tofInvBetaErr = 0;
  tofMass = 0;
  tofMassError = 0;
  isGoodBeta = 0;
  isTkTofMatch = 0;
  tofReducedChisqC = 0;
  tofReducedChisqT = 0;
  tofDepositedEnergyOnLayer[0] = 0;
  tofDepositedEnergyOnLayer[1] = 0;
  tofDepositedEnergyOnLayer[2] = 0;
  tofDepositedEnergyOnLayer[3] = 0;
  tofEstimatedChargeOnLayer[0] = 0;
  tofEstimatedChargeOnLayer[1] = 0;
  tofEstimatedChargeOnLayer[2] = 0;
  tofEstimatedChargeOnLayer[3] = 0;
  tofCharge = -1;
  tofChargeRMS = 0;
  tofUpperCharge = -1;
  tofLowerCharge = -1;
  tofChargeOnLayer[0] = -1;
  tofChargeOnLayer[1] = -1;
  tofChargeOnLayer[2] = -1;
  tofChargeOnLayer[3] = -1;
  tofZ = -1;
  trkFitCodeL1Inner = -1;
  trkRigidityL1Inner = -9999;
  trkRigidityInverseErrorL1Inner = -9999;
  trkReducedChisquareL1InnerX = -9999;
  trkReducedChisquareL1InnerY = -9999;
  trkFitCodeFS = -1;
  trkRigidityFS = -9999;
  trkRigidityInverseErrorFS = -9999;
  trkReducedChisquareFSX = -9999;
  trkReducedChisquareFSY = -9999;
  trkFitCodeMS = -1;
  trkRigidityMS = -9999;
  trkRigidityInverseErrorMS = -9999;
  trkReducedChisquareMSY = -9999;
  trkReducedChisquareMSX = -9999;
  trkFitCodeInner = -1;
  trkRigidityInverseErrorInner = -9999;
  trkReducedChisquareInnerX = -9999;
  trkReducedChisquareInnerY = -9999;
  trkFitCodeUpperInner = -1;
  trkRigidityInverseErrorUpperInner = -9999;
  trkReducedChisquareUpperInnerX = -9999;
  trkReducedChisquareUpperInnerY = -9999;
  trkFitCodeLowerInner = -1;
  trkRigidityInverseErrorLowerInner = -9999;
  trkReducedChisquareLowerInnerX = -9999;
  trkReducedChisquareLowerInnerY = -9999;
  trkLayerJQ[0] = -1;
  trkLayerJQ[1] = -1;
  trkLayerJQ[2] = -1;
  trkLayerJQ[3] = -1;
  trkLayerJQ[4] = -1;
  trkLayerJQ[5] = -1;
  trkLayerJQ[6] = -1;
  trkLayerJQ[7] = -1;
  trkLayerJQ[8] = -1;
  trkEdepLayerJXSideOK[0] = 0;
  trkEdepLayerJXSideOK[1] = 0;
  trkEdepLayerJXSideOK[2] = 0;
  trkEdepLayerJXSideOK[3] = 0;
  trkEdepLayerJXSideOK[4] = 0;
  trkEdepLayerJXSideOK[5] = 0;
  trkEdepLayerJXSideOK[6] = 0;
  trkEdepLayerJXSideOK[7] = 0;
  trkEdepLayerJXSideOK[8] = 0;
  trkEdepLayerJYSideOK[0] = 0;
  trkEdepLayerJYSideOK[1] = 0;
  trkEdepLayerJYSideOK[2] = 0;
  trkEdepLayerJYSideOK[3] = 0;
  trkEdepLayerJYSideOK[4] = 0;
  trkEdepLayerJYSideOK[5] = 0;
  trkEdepLayerJYSideOK[6] = 0;
  trkEdepLayerJYSideOK[7] = 0;
  trkEdepLayerJYSideOK[8] = 0;
  trkEdepLayerJ[0] = 0;
  trkEdepLayerJ[1] = 0;
  trkEdepLayerJ[2] = 0;
  trkEdepLayerJ[3] = 0;
  trkEdepLayerJ[4] = 0;
  trkEdepLayerJ[5] = 0;
  trkEdepLayerJ[6] = 0;
  trkEdepLayerJ[7] = 0;
  trkEdepLayerJ[8] = 0;
  for( int i = 0; i < 9; i++ )
  {
    trkPosXLJ[i] = 0;
    trkPosYLJ[i] = 0;
    trkPosZLJ[i] = 0;
    trkDirThetaLJ[i] = 0;
    trkDirPhiLJ[i] = 0;
  }
  trkCharge = -1;
  trkInnerCharge = -1;
  trkZ = -1;
  trkInnerZ = -1;
  for( int i = 0; i < 9; i++ ) trkLayerJZ[i] = -1;
  trkHasExtLayers = 0;
  isRichAvailable = 0;
  richRebuild = 0;
  richIsGood = 0;
  richIsClean = 0;
  richIsNaF = 0;
  richTileIndex = -1;
  richDistanceTileBorder = -1;
  richChargeCorrections = 0;
  richPMTCorrectionsFailed = -10;
  richWidth = -1;
  richBetaConsistency = -999;
  richReflectedHits = -1;
  richPMTChargeConsistency = -999;
  richUsedHits = 0;
  richRingWidth = 0;
  richNHits = 0;
  richNPMTsOnRing = 0;
  richBeta = 0;
  richBetaError = 0;
  richChargeSquared = 0;
  richKolmogorovProbability = 0;
  richPhotoelectrons = 0;
  richPhotoElectrons = 0;
  richExpectedPhotoelectrons = 0;
  richExpectedPhotoElectrons = 0;
  richNUsedHits = 0;
  richTheta = 0;
  richPhi = 0;
  richTrackerTrackOnRadiatorX = -999;
  richTrackerTrackOnRadiatorY = -999;
  trdNClusters = 0;
  trdNUnusedHits = 0;
  trdNTracks = 0;
  trdNVertex = 0;
  trdTrackPhi = 0;
  trdTrackTheta = 0;
  trdTrackChi2 = 0;
  trdTrackPattern = 0;
  trdTrackCharge = 0;
  trdTrackEdepL[0] = 0;
  trdTrackEdepL[1] = 0;
  trdTrackEdepL[2] = 0;
  trdTrackEdepL[3] = 0;
  trdTrackEdepL[4] = 0;
  trdTrackEdepL[5] = 0;
  trdTrackEdepL[6] = 0;
  trdTrackEdepL[7] = 0;
  trdTrackEdepL[8] = 0;
  trdTrackEdepL[9] = 0;
  trdTrackEdepL[10] = 0;
  trdTrackEdepL[11] = 0;
  trdTrackEdepL[12] = 0;
  trdTrackEdepL[13] = 0;
  trdTrackEdepL[14] = 0;
  trdTrackEdepL[15] = 0;
  trdTrackEdepL[16] = 0;
  trdTrackEdepL[17] = 0;
  trdTrackEdepL[18] = 0;
  trdTrackEdepL[19] = 0;
  trdTrackMeanDepositedEnergy = 0;
  trdTrackTotalDepositedEnergy = 0;
  trdTrackDeviationXWithInnerTrk = 0;
  trdTrackDeviationYWithInnerTrk = 0;
  trdTrackDistanceBetweenInnerTrk = 0;
  ResidualXBetweenInnerTrackAndSplineTrack = -9999;
  ResidualYBetweenInnerTrackAndSplineTrack = -9999;
  trdQtNActiveLayer = 0;
  trdQtIsCalibrationGood = 0;
  trdQtIsInsideTrdGeometricalAcceptance = 0;
  trdQtIsValid = 0;
  trdQtNActiveStraws = 0;
  trdQtNActiveLayers = 0;
  trdQtNTRDVertex = 0;
  trdQtElectronToProtonLogLikelihoodRatio = 0;
  trdQtHeliumToElectronLogLikelihoodRatio = 0;
  trdQtHeliumToProtonLogLikelihoodRatio = 0;
  trdKNRawHits = 0;
  trdKIsReadAlignmentOK = 0;
  trdKIsReadCalibOK = 0;
  trdKNHits = 0;
  trdKIsValid = 0;
  trdKElectronToProtonLogLikelihoodRatio = 0;
  trdKHeliumToElectronLogLikelihoodRatio = 0;
  trdKHeliumToProtonLogLikelihoodRatio = 0;
  trdKCharge = -1;
  trdKChargeError = 0;
  trdKNUsedHitsForCharge = 0;
  trdKAmpLayer[0] = 0;
  trdKAmpLayer[1] = 0;
  trdKAmpLayer[2] = 0;
  trdKAmpLayer[3] = 0;
  trdKAmpLayer[4] = 0;
  trdKAmpLayer[5] = 0;
  trdKAmpLayer[6] = 0;
  trdKAmpLayer[7] = 0;
  trdKAmpLayer[8] = 0;
  trdKAmpLayer[9] = 0;
  trdKAmpLayer[10] = 0;
  trdKAmpLayer[11] = 0;
  trdKAmpLayer[12] = 0;
  trdKAmpLayer[13] = 0;
  trdKAmpLayer[14] = 0;
  trdKAmpLayer[15] = 0;
  trdKAmpLayer[16] = 0;
  trdKAmpLayer[17] = 0;
  trdKAmpLayer[18] = 0;
  trdKAmpLayer[19] = 0;
  trdKTotalPathLength = 0;
  trdKTotalAmp = 0;
  trdKElectronLikelihood = 0;
  trdKHeliumLikelihood = 0;
  trdPElectronToProtonLogLikelihoodRatio = 0;
  trdPHeliumToElectronLogLikelihoodRatio = 0;
  trdPHeliumToProtonLogLikelihoodRatio = 0;
  trdPElectronLikelihood = 0;
  trdPProtonLikelihood = 0;
  trdPHeliumLikelihood = 0;
  trdPNHybridHits = 0;
  trdPNActiveLayers = 0;
  trdPNLowAdcHits = 0;
  trdPNLowDxHits = 0;
  acTrackerTrdTrackDistanceAtTrdCenterInSigmas = 0;
  acTrackerTrdTrackAngularDistanceInSigmas = 0;
  accNHits = 0;
  accNRecoClPG = 0;
  accSector.clear();
  accTime.clear();
  accHitPosZ.clear();
  accChi2.clear();
  accUnfoldedHitPosZ.clear();
  accRawCharge.clear();
  accHitPosZfromADC.clear();
  accUnfoldedHitPosZfromADC.clear();
  accCrossingPhiAngle = 0;
  accEdep.clear();
  accTimePG.clear();

  outTree->Branch("nRun",                                    &nRun,                                    "nRun/i");
  outTree->Branch("nEvent",                                  &nEvent,                                  "nEvent/i");
  outTree->Branch("nLevel1",                                 &nLevel1,                                 "nLevel1/i");
  outTree->Branch("nParticle",                               &nParticle,                               "nParticle/i");
  outTree->Branch("nCharge",                                 &nCharge,                                 "nCharge/i");
  outTree->Branch("nTrTrack",                                &nTrTrack,                                "nTrTrack/i");
  outTree->Branch("nTrdTrack",                               &nTrdTrack,                               "nTrdTrack/i");
  outTree->Branch("nAntiCluster",                            &nAntiCluster,                            "nAntiCluster/i");
  outTree->Branch("nTofClustersInTime",                      &nTofClustersInTime,                      "nTofClustersInTime/I");
  outTree->Branch("nRichRing",                               &nRichRing,                               "nRichRing/i");
  outTree->Branch("nRichRingB",                              &nRichRingB,                              "nRichRingB/i");
  outTree->Branch("nBeta",                                   &nBeta,                                   "nBeta/i");
  outTree->Branch("nBetaB",                                  &nBetaB,                                  "nBetaB/i");
  outTree->Branch("nBetaH",                                  &nBetaH,                                  "nBetaH/i");
  outTree->Branch("nShower",                                 &nEcalShower,                             "nShower/i");
  outTree->Branch("nVertex",                                 &nVertex,                                 "nVertex/i");
  if( isMC == true )
  {
    outTree->Branch("survived_tof",                          &survived_tof,                            "survived_tof/I");
    outTree->Branch("nAntiMCCluster",                        &nAntiMCCluster,                          "nAntiMCCluster/i");
    outTree->Branch("nTrMCCluster",                          &nTrMCCluster,                            "nTrMCCluster/i");
    outTree->Branch("nTofMCCluster",                         &nTofMCCluster,                           "nTofMCCluster/i");
    outTree->Branch("nTofMCPmtHit",                          &nTofMCPmtHit,                            "nTofMCPmtHit/i");
    outTree->Branch("nEcalMCHit",                            &nEcalMCHit,                              "nEcalMCHit/i");
    outTree->Branch("nTrdMCCluster",                         &nTrdMCCluster,                           "nTrdMCCluster/i");
    outTree->Branch("nRichMCCluster",                        &nRichMCCluster,                          "nRichMCCluster/i");
    outTree->Branch("nMCTrack",                              &nMCTrack,                                "nMCTrack/i");
    outTree->Branch("nMCEventg",                             &nMCEventg,                               "nMCEventg/i");
    outTree->Branch("particleID",                            &particleID,                              "particleID/I");
    outTree->Branch("trkID",                                 &trkID,                                   "trkID/I");
    outTree->Branch("parentID",                              &parentID,                                "parentID/I");
    outTree->Branch("genMomentum",                           &genMomentum,                             "genMomentum/F");
    outTree->Branch("genPosition",                           &genPosition,                             "genPosition[3]/F");
    outTree->Branch("genDirCosine",                          &genDirCosine,                            "genDirCosine[3]/F");
    outTree->Branch("genMass",                               &genMass,                                 "genMass/F");
    outTree->Branch("genCharge",                             &genCharge,                               "genCharge/F");
    outTree->Branch("conversionId",                          &conversionId,                            "conversionId/I");
    outTree->Branch("conversionIdContainer",                 &conversionIdContainer);
    outTree->Branch("conversionIdSet",                       &conversionIdSet);
  }
  outTree->Branch("particleType",                            &particleType,                            "particleType/I");
  if( isMC == false )
  {
    outTree->Branch("liveTime",                                &liveTime,                                "liveTime/F");
    outTree->Branch("utcTime",                                 &utcTime,                                 "utcTime/F");
    outTree->Branch("orbitAltitude",                           &orbitAltitude,                           "orbitAltitude/F");
    outTree->Branch("orbitLatitude",                           &orbitLatitude,                           "orbitLatitude/F");
    outTree->Branch("orbitLongitude",                          &orbitLongitude,                          "orbitLongitude/F");
    outTree->Branch("orbitLatitudeM",                          &orbitLatitudeM,                          "orbitLatitudeM/F");
    outTree->Branch("orbitLongitudeM",                         &orbitLongitudeM,                         "orbitLongitudeM/F");
    outTree->Branch("velR",                                    &velR,                                    "velR/F");
    outTree->Branch("velTheta",                                &velTheta,                                "velTheta/F");
    outTree->Branch("velPhi",                                  &velPhi,                                  "velPhi/F");
    outTree->Branch("yaw",                                     &yaw,                                     "yaw/F");
    outTree->Branch("pitch",                                   &pitch,                                   "pitch/F");
    outTree->Branch("roll",                                    &roll,                                    "roll/F");
    outTree->Branch("gLongitude",                              &gLongitude,                              "gLongitude/F");
    outTree->Branch("gLatitude",                               &gLatitude,                               "gLatitude/F");
    outTree->Branch("gCoordCalcResult",                        &gCoordCalcResult,                        "gCoordCalcResult/I");
    outTree->Branch("gCoordCalcResultBit",                     &gCoordCalcResultBit,                     "gCoordCalcResultBit/I");
    outTree->Branch("sunPosAzimuth",                           &sunPosAzimuth,                           "sunPosAzimuth/F");
    outTree->Branch("sunPosElevation",                         &sunPosElevation,                         "sunPosElevation/F");
    outTree->Branch("sunPosCalcResult",                        &sunPosCalcResult,                        "sunPosCalcResult/I");
    outTree->Branch("unixTime",                                &unixTime,                                "unixTime/i");
    outTree->Branch("acEventTime",                             &acEventTime,                             "acEventTime/F");
    outTree->Branch("solarArrayCoord",                         &solarArrayCoord,                         "solarArrayCoord[3]/F");
    outTree->Branch("isInShadow",                              &isInShadow,                              "isInShadow/I");
    outTree->Branch("zenithAngle",                             &zenithAngle,                             "zenithAngle/F");
    outTree->Branch("isInSAA",                                 &isInSAA,                                 "isInSAA/I");
  }
  //outTree->Branch("ptlCharge",                               &ptlCharge,                               "ptlCharge/i");
  /*
  outTree->Branch("ptlMomentum",                             &ptlMomentum,                             "ptlMomentum/F");
  outTree->Branch("ptlTheta",                                &ptlTheta,                                "ptlTheta/F");
  outTree->Branch("ptlPhi",                                  &ptlPhi,                                  "ptlPhi/F");
  outTree->Branch("ptlCoo",                                  &ptlCoo,                                  "ptlCoo[3]/F");
  outTree->Branch("ptlTofCoo",                               &ptlTofCoo,                               "ptlTofCoo[4][3]/F");
  outTree->Branch("ptlTofTrLength",                          &ptlTofTrLength,                          "ptlTofTrLength[4]/F");
  outTree->Branch("ptlAntiCoo",                              &ptlAntiCoo,                              "ptlAntiCoo[2][5]/F");
  outTree->Branch("ptlEcalCoo",                              &ptlEcalCoo,                              "ptlEcalCoo[3][3]/F");
  outTree->Branch("ptlTrCoo",                                &ptlTrCoo,                                "ptlTrCoo[9][3]/F");
  outTree->Branch("ptlTrdCoo",                               &ptlTrdCoo,                               "ptlTrdCoo[2][3]/F");
  outTree->Branch("ptlRichCoo",                              &ptlRichCoo,                              "ptlRichCoo[2][3]/F");
  outTree->Branch("ptlRichPath",                             &ptlRichPath,                             "ptlRichPath[2]/F");
  outTree->Branch("ptlRichPathBeta",                         &ptlRichPathBeta,                         "ptlRichPathBeta[2]/F");
  outTree->Branch("ptlRichLength",                           &ptlRichLength,                           "ptlRichLength/F");
  outTree->Branch("ptlRichParticles",                        &ptlRichParticles,                        "ptlRichParticles/I");
  outTree->Branch("ptlCutOffStoermer",                       &ptlCutOffStoermer,                       "ptlCutOffStoermer/F");
  outTree->Branch("ptlCutOffDipole",                         &ptlCutOffDipole,                         "ptlCutOffDipole/F");
  */
  //outTree->Branch("ptlCutOffMax",                            &ptlCutOffMax,                            "ptlCutOffMax[2]/F");
  outTree->Branch("isEcalAvailable",                         &isEcalAvailable,                         "isEcalAvailable/I");
  outTree->Branch("showerEnergyD",                           &showerEnergyD,                           "showerEnergyD/F");
  outTree->Branch("showerEnergyDL",                          &showerEnergyDL,                          "showerEnergyDL[18]/F");
  outTree->Branch("showerEnergyE",                           &showerEnergyE,                           "showerEnergyE/F");
  outTree->Branch("showerEnergyCorrected",                   &showerEnergyCorrected,                   "showerEnergyCorrected/F");
  outTree->Branch("showerBDT",                               &showerBDT,                               "showerBDT/F");
  outTree->Branch("showerCofG",                              &showerCofG,                              "showerCofG[3]/F");
  outTree->Branch("showerCofGDist",                          &showerCofGDist,                          "showerCofGDist/F");
  outTree->Branch("showerCofGdX",                            &showerCofGdX,                            "showerCofGdX/F");
  outTree->Branch("showerCofGdY",                            &showerCofGdY,                            "showerCofGdY/F");
  outTree->Branch("tofNCluster",                             &tofNCluster,                             "tofNCluster/I");
  outTree->Branch("tofNClusterH",                            &tofNClusterH,                            "tofNClusterH/I");
  outTree->Branch("tofNUsedHits",                            &tofNUsedHits,                            "tofNUsedHits/I");
  outTree->Branch("tofNUnusedHits",                          &tofNUnusedHits,                          "tofNUnusedHits/I");
  outTree->Branch("tofNUsedLayersForQ",                      &tofNUsedLayersForQ,                      "tofNUsedLayersForQ/I");
  outTree->Branch("tofBeta",                                 &tofBeta,                                 "tofBeta/F");
  outTree->Branch("tofInvBetaErr",                           &tofInvBetaErr,                           "tofInvBetaErr/F");
  outTree->Branch("tofNormEBetaV",                           &tofNormEBetaV,                           "tofNormEBetaV/F");
  outTree->Branch("tofBetaS",                                &tofBetaS,                                "tofBetaS/F");
  outTree->Branch("tofBetaC",                                &tofBetaC,                                "tofBetaC/F");
  outTree->Branch("tofEBetaCV",                              &tofEBetaCV,                              "tofEBetaCV/F");
  outTree->Branch("tofMass",                                 &tofMass,                                 "tofMass/F");
  outTree->Branch("tofMassError",                            &tofMassError,                            "tofMassError/F");
  outTree->Branch("isGoodBeta",                              &isGoodBeta,                              "isGoodBeta/I");
  outTree->Branch("isTkTofMatch",                            &isTkTofMatch,                            "isTkTofMatch/I");
  outTree->Branch("tofChisqT",                               &tofChisqT,                               "tofChisqT/F");
  outTree->Branch("tofChisqC",                               &tofChisqC,                               "tofChisqC/F");
  outTree->Branch("tofReducedChisqT",                        &tofReducedChisqT,                        "tofReducedChisqT/F");
  outTree->Branch("tofReducedChisqC",                        &tofReducedChisqC,                        "tofReducedChisqC/F");
  outTree->Branch("tofSumHit",                               &tofSumHit,                               "tofSumHit/I");
  outTree->Branch("tofNUsedClusterH",                        &tofNUsedClusterH,                        "tofNUsedClusterH/I");
  outTree->Branch("tofPatternOnLayer",                       &tofPatternOnLayer,                       "tofPatternOnLayer[4]/I");
  outTree->Branch("tofTimeOnLayer",                          &tofTimeOnLayer,                          "tofTimeOnLayer[4]/F");
  outTree->Branch("tofETimeOnLayer",                         &tofETimeOnLayer,                         "tofETimeOnLayer[4]/F");
  outTree->Branch("tofETCooOnLayer",                         &tofETCooOnLayer,                         "tofETCooOnLayer[4]/F");
  outTree->Branch("tofTResidualOnLayer",                     &tofTResidualOnLayer,                     "tofTResidualOnLayer[4]/F");
  outTree->Branch("tofT0",                                   &tofT0,                                   "tofT0/F");
  outTree->Branch("tofTkTFLenOnLayer",                       &tofTkTFLenOnLayer,                       "tofTkTFLenOnLayer[4]/F");
  outTree->Branch("tofCResidualXOnLayer",                    &tofCResidualXOnLayer,                    "tofCResidualXOnLayer[4]/F");
  outTree->Branch("tofCResidualYOnLayer",                    &tofCResidualYOnLayer,                    "tofCResidualYOnLayer[4]/F");
  outTree->Branch("tofDepositedEnergyOnLayer",               &tofDepositedEnergyOnLayer,               "tofDepositedEnergyOnLayer[4]/F");
  outTree->Branch("tofEstimatedChargeOnLayer",               &tofEstimatedChargeOnLayer,               "tofEstimatedChargeOnLayer[4]/F");
  outTree->Branch("tofCharge",                               &tofCharge,                               "tofCharge/F");
  outTree->Branch("tofChargeRMS",                            &tofChargeRMS,                            "tofChargeRMS/F");
  outTree->Branch("tofAcCharge",                             &tofAcCharge,                             "tofAcCharge/F");
  outTree->Branch("tofAcUpperCharge",                        &tofAcUpperCharge,                        "tofAcUpperCharge/F");
  outTree->Branch("tofAcLowerCharge",                        &tofAcLowerCharge,                        "tofAcLowerCharge/F");
  outTree->Branch("tofChargeOnLayer",                        &tofChargeOnLayer,                        "tofChargeOnLayer[4]/F");
  outTree->Branch("tofZ",                                    &tofZ,                                    "tofZ/I");
  outTree->Branch("nTofLUsedForZ",                           &nTofLUsedForZ,                           "nTofLUsedForZ/I");
  outTree->Branch("probTOFZ",                                &probTOFZ,                                "probTOFZ/F");

  outTree->Branch("trkFitCodeFS",                            &trkFitCodeFS,                            "trkFitCodeFS/I");
  outTree->Branch("trkRigidityFS",                           &trkRigidityFS,                           "trkRigidityFS/F");
  outTree->Branch("trkRigidityInverseErrorFS",               &trkRigidityInverseErrorFS,               "trkRigidityInverseErrorFS/F");
  outTree->Branch("trkReducedChisquareFSX",                  &trkReducedChisquareFSX,                  "trkReducedChisquareFSX/F");
  outTree->Branch("trkReducedChisquareFSY",                  &trkReducedChisquareFSY,                  "trkReducedChisquareFSY/F");

  outTree->Branch("trkFitCodeL1Inner",                       &trkFitCodeL1Inner,                       "trkFitCodeL1Inner/I");
  outTree->Branch("trkRigidityL1Inner",                      &trkRigidityL1Inner,                      "trkRigidityL1Inner/F");
  outTree->Branch("trkRigidityInverseErrorL1Inner",          &trkRigidityInverseErrorL1Inner,          "trkRigidityInverseErrorL1Inner/F");
  outTree->Branch("trkReducedChisquareL1InnerX",             &trkReducedChisquareL1InnerX,             "trkReducedChisquareL1InnerX/F");
  outTree->Branch("trkReducedChisquareL1InnerY",             &trkReducedChisquareL1InnerY,             "trkReducedChisquareL1InnerY/F");

  outTree->Branch("trkFitCodeL9Inner",                       &trkFitCodeL9Inner,                       "trkFitCodeL9Inner/I");
  outTree->Branch("trkRigidityL9Inner",                      &trkRigidityL9Inner,                      "trkRigidityL9Inner/F");
  outTree->Branch("trkRigidityInverseErrorL9Inner",          &trkRigidityInverseErrorL9Inner,          "trkRigidityInverseErrorL9Inner/F");
  outTree->Branch("trkReducedChisquareL9InnerX",             &trkReducedChisquareL9InnerX,             "trkReducedChisquareL9InnerX/F");
  outTree->Branch("trkReducedChisquareL9InnerY",             &trkReducedChisquareL9InnerY,             "trkReducedChisquareL9InnerY/F");

  outTree->Branch("trkFitCodeUpperInner",                    &trkFitCodeUpperInner,                    "trkFitCodeUpperInner/I");
  outTree->Branch("trkRigidityUpperInner",                   &trkRigidityUpperInner,                   "trkRigidityUpperInner/F");
  outTree->Branch("trkRigidityInverseErrorUpperInner",       &trkRigidityInverseErrorUpperInner,       "trkRigidityInverseErrorUpperInner/F");
  outTree->Branch("trkReducedChisquareUpperInnerX",          &trkReducedChisquareUpperInnerX,          "trkReducedChisquareUpperInnerX/F");
  outTree->Branch("trkReducedChisquareUpperInnerY",          &trkReducedChisquareUpperInnerY,          "trkReducedChisquareUpperInnerY/F");

  outTree->Branch("trkFitCodeLowerInner",                    &trkFitCodeLowerInner,                    "trkFitCodeLowerInner/I");
  outTree->Branch("trkRigidityLowerInner",                   &trkRigidityLowerInner,                   "trkRigidityLowerInner/F");
  outTree->Branch("trkRigidityInverseErrorLowerInner",       &trkRigidityInverseErrorLowerInner,       "trkRigidityInverseErrorLowerInner/F");
  outTree->Branch("trkReducedChisquareLowerInnerX",          &trkReducedChisquareLowerInnerX,          "trkReducedChisquareLowerInnerX/F");
  outTree->Branch("trkReducedChisquareLowerInnerY",          &trkReducedChisquareLowerInnerY,          "trkReducedChisquareLowerInnerY/F");

  outTree->Branch("trkFitCodeMS",                            &trkFitCodeMS,                            "trkFitCodeMS/I");
  outTree->Branch("trkRigidityMS",                           &trkRigidityMS,                           "trkRigidityMS/F");
  outTree->Branch("trkRigidityInverseErrorMS",               &trkRigidityInverseErrorMS,               "trkRigidityInverseErrorMS/F");
  outTree->Branch("trkReducedChisquareMSX",                  &trkReducedChisquareMSX,                  "trkReducedChisquareMSX/F");
  outTree->Branch("trkReducedChisquareMSY",                  &trkReducedChisquareMSY,                  "trkReducedChisquareMSY/F");

  outTree->Branch("trkFitCodeInner",                         &trkFitCodeInner,                         "trkFitCodeInner/I");
  outTree->Branch("trkRigidityInner",                        &trkRigidityInner,                        "trkRigidityInner/F");
  outTree->Branch("trkRigidityInverseErrorInner",            &trkRigidityInverseErrorInner,            "trkRigidityInverseErrorInner/F");
  outTree->Branch("trkReducedChisquareInnerX",               &trkReducedChisquareInnerX,               "trkReducedChisquareInnerX/F");
  outTree->Branch("trkReducedChisquareInnerY",               &trkReducedChisquareInnerY,               "trkReducedChisquareInnerY/F");

  outTree->Branch("trkLayerJQ",                              &trkLayerJQ,                              "trkLayerJQ[9]/F");
  outTree->Branch("trkEdepLayerJXSideOK",                    &trkEdepLayerJXSideOK,                    "trkEdepLayerJXSideOK[9]/I");
  outTree->Branch("trkEdepLayerJYSideOK",                    &trkEdepLayerJYSideOK,                    "trkEdepLayerJYSideOK[9]/I");
  outTree->Branch("trkEdepLayerJ",                           &trkEdepLayerJ,                           "trkEdepLayerJ[9]/F");
  outTree->Branch("trkPosXLJ",                               &trkPosXLJ,                               "trkPosXLJ[9]/F");
  outTree->Branch("trkPosYLJ",                               &trkPosYLJ,                               "trkPosYLJ[9]/F");
  outTree->Branch("trkPosZLJ",                               &trkPosZLJ,                               "trkPosZLJ[9]/F");
  outTree->Branch("trkDirThetaLJ",                           &trkDirThetaLJ,                           "trkDirThetaLJ[9]/F");
  outTree->Branch("trkDirPhiLJ",                             &trkDirPhiLJ,                             "trkDirPhiLJ[9]/F");
  outTree->Branch("trkCharge",                               &trkCharge,                               "trkCharge/F");
  outTree->Branch("trkInnerCharge",                          &trkInnerCharge,                          "trkInnerCharge/F");
  outTree->Branch("trkZ",                                    &trkZ,                                    "trkZ/I");
  outTree->Branch("trkInnerZ",                               &trkInnerZ,                               "trkInnerZ/I");
  outTree->Branch("trkLayerJZ",                              &trkLayerJZ,                              "trkLayerJZ[9]/I");
  outTree->Branch("trkHasExtLayers",                         &trkHasExtLayers,                         "trkHasExtLayers/F");

  outTree->Branch("isRichAvailable",                         &isRichAvailable,                         "isRichAvailable/I");
  outTree->Branch("richRebuild",                             &richRebuild,                             "richRebuild/I");
  outTree->Branch("richIsGood",                              &richIsGood,                              "richIsGood/I");
  outTree->Branch("richIsClean",                             &richIsClean,                             "richIsClean/I");
  outTree->Branch("richIsNaF",                               &richIsNaF,                               "richIsNaF/I");
  outTree->Branch("richTileIndex",                           &richTileIndex,                           "richTileIndex/I");
  outTree->Branch("richDistanceTileBorder",                  &richDistanceTileBorder,                  "richDistanceTileBorder/F");
  outTree->Branch("richChargeCorrections",                   &richChargeCorrections,                   "richChargeCorrections/I");
  outTree->Branch("richPMTCorrectionsFailed",                &richPMTCorrectionsFailed,                "richPMTCorrectionsFailed/I");
  outTree->Branch("richTrackEmissionPoint",                  &richTrackEmissionPoint,                  "richTrackEmissionPoint[5]/F");
  outTree->Branch("richUsedHits",                            &richUsedHits,                            "richUsedHits/I");
  outTree->Branch("richRingWidth",                           &richRingWidth,                           "richRingWidth/F");
  outTree->Branch("richWidth",                               &richWidth,                               "richRingWidth/F");
  outTree->Branch("richNHits",                               &richNHits,                               "richNHits/I");
  outTree->Branch("richNPMTsOnRing",                         &richNPMTsOnRing,                         "richNPMTsOnRing/I");
  outTree->Branch("richBeta",                                &richBeta,                                "richBeta/F");
  outTree->Branch("richBetaError",                           &richBetaError,                           "richBetaError/F");
  outTree->Branch("richChargeSquared",                       &richChargeSquared,                       "richChargeSquared/F");
  outTree->Branch("richKolmogorovProbability",               &richKolmogorovProbability,               "richKolmogorovProbability/F");
  outTree->Branch("richPhotoelectrons",                      &richPhotoelectrons,                      "richPhotoelectrons/F");
  outTree->Branch("richExpectedPhotoelectrons",              &richExpectedPhotoelectrons,              "richExpectedPhotoelectrons/F");
  outTree->Branch("richTheta",                               &richTheta,                               "richTheta/F");
  outTree->Branch("richPhi",                                 &richPhi,                                 "richPhi/F");
  outTree->Branch("richBetaConsistency",                     &richBetaConsistency,                     "richBetaConsistency/F");
  outTree->Branch("richReflectedHits",                       &richReflectedHits,                       "richReflectedHits/I");
  outTree->Branch("richPMTChargeConsistency",                &richPMTChargeConsistency,                "richPMTChargeConsistency/F");
  outTree->Branch("richTrackerTrackOnRadiatorX",             &richTrackerTrackOnRadiatorX,             "richTrackerTrackOnRadiatorX/F");
  outTree->Branch("richTrackerTrackOnRadiatorY",             &richTrackerTrackOnRadiatorY,             "richTrackerTrackOnRadiatorY/F");
  outTree->Branch("trdNClusters",                            &trdNClusters,                            "trdNClusters/I");
  outTree->Branch("trdNUnusedHits",                          &trdNUnusedHits,                          "trdNUnusedHits/I");
  outTree->Branch("trdNTracks",                              &trdNTracks,                              "trdNTracks/I");
  outTree->Branch("trdNVertex",                              &trdNVertex,                              "trdNVertex/I");
  outTree->Branch("trdTrackTheta",                           &trdTrackTheta,                           "trdTrackTheta/F");
  outTree->Branch("trdTrackPhi",                             &trdTrackPhi,                             "trdTrackPhi/F");
  outTree->Branch("trdTrackChi2",                            &trdTrackChi2,                            "trdTrackChi2/F");
  outTree->Branch("trdTrackPattern",                         &trdTrackPattern,                         "trdTrackPattern/I");
  outTree->Branch("trdTrackCharge",                          &trdTrackCharge,                          "trdTrackCharge/F");
  outTree->Branch("trdTrackEdepL",                           &trdTrackEdepL,                           "trdTrackEdepL[20]/F");
  outTree->Branch("trdTrackMeanDepositedEnergy",             &trdTrackMeanDepositedEnergy,             "trdTrackMeanDepositedEnergy/F");
  outTree->Branch("trdTrackTotalDepositedEnergy",            &trdTrackTotalDepositedEnergy,            "trdTrackTotalDepositedEnergy/F");
  outTree->Branch("trdTrackDeviationXWithInnerTrk",          &trdTrackDeviationXWithInnerTrk,          "trdTrackDeviationXWithInnerTrk/F");
  outTree->Branch("trdTrackDeviationYWithInnerTrk",          &trdTrackDeviationYWithInnerTrk,          "trdTrackDeviationYWithInnerTrk/F");
  outTree->Branch("trdTrackDistanceBetweenInnerTrk",         &trdTrackDistanceBetweenInnerTrk,         "trdTrackDistanceBetweenInnerTrk/F");
  outTree->Branch("ResidualXBetweenInnerTrackAndSplineTrack",&ResidualXBetweenInnerTrackAndSplineTrack,"ResidualXBetweenInnerTrackAndSplineTrack/F");
  outTree->Branch("ResidualYBetweenInnerTrackAndSplineTrack",&ResidualYBetweenInnerTrackAndSplineTrack,"ResidualYBetweenInnerTrackAndSplineTrack/F");
  /*
  outTree->Branch("trdQtNActiveLayer",                       &trdQtNActiveLayer,                       "trdQtNActiveLayer/I");
  outTree->Branch("trdQtIsCalibrationGood",                  &trdQtIsCalibrationGood,                  "trdQtIsCalibrationGood/I");
  outTree->Branch("trdQtIsSlowControlDataGood",              &trdQtIsSlowControlDataGood,              "trdQtIsSlowControlDataGood/I");
  outTree->Branch("trdQtIsInsideTrdGeometricalAcceptance",   &trdQtIsInsideTrdGeometricalAcceptance,   "trdQtIsInsideTrdGeometricalAcceptance/I");
  outTree->Branch("trdQtIsValid",                            &trdQtIsValid,                            "trdQtIsValid/I");
  outTree->Branch("trdQtNActiveStraws",                      &trdQtNActiveStraws,                      "trdQtNActiveStraws/I");
  outTree->Branch("trdQtNActiveLayers",                      &trdQtNActiveLayers,                      "trdQtNActiveLayers/I");
  outTree->Branch("trdQtNTRDVertex",                         &trdQtNTRDVertex,                         "trdQtNTRDVertex/I");
  outTree->Branch("trdQtElectronToProtonLogLikelihoodRatio", &trdQtElectronToProtonLogLikelihoodRatio, "trdQtElectronToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdQtHeliumToElectronLogLikelihoodRatio", &trdQtHeliumToElectronLogLikelihoodRatio, "trdQtHeliumToElectronLogLikelihoodRatio/F");
  outTree->Branch("trdQtHeliumToProtonLogLikelihoodRatio",   &trdQtHeliumToProtonLogLikelihoodRatio,   "trdQtHeliumToProtonLogLikelihoodRatio/F");
  */
  outTree->Branch("trdKNRawHits",                            &trdKNRawHits,                            "trdKNRawHits/I");
  outTree->Branch("trdKIsReadAlignmentOK",                   &trdKIsReadAlignmentOK,                   "trdKIsReadAlignmentOK/I");
  outTree->Branch("trdKIsReadCalibOK",                       &trdKIsReadCalibOK,                       "trdKIsReadCalibOK/I");
  outTree->Branch("trdKNHits",                               &trdKNHits,                               "trdKNHits/I");
  outTree->Branch("trdKIsValid",                             &trdKIsValid,                             "trdKIsValid/I");
  outTree->Branch("trdKElectronToProtonLogLikelihoodRatio",  &trdKElectronToProtonLogLikelihoodRatio,  "trdKElectronToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdKHeliumToElectronLogLikelihoodRatio",  &trdKHeliumToElectronLogLikelihoodRatio,  "trdKHeliumToElectronLogLikelihoodRatio/F");
  outTree->Branch("trdKHeliumToProtonLogLikelihoodRatio",    &trdKHeliumToProtonLogLikelihoodRatio,    "trdKHeliumToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdKCharge",                              &trdKCharge,                              "trdCharge/F");
  outTree->Branch("trdKChargeError",                         &trdKChargeError,                         "trdKChargeError/F");
  outTree->Branch("trdKNUsedHitsForCharge",                  &trdKNUsedHitsForCharge,                  "trdKNUsedHitsForCharge/I");
  outTree->Branch("trdKAmpLayer",                            &trdKAmpLayer,                            "trdKAmpLayer[20]/F");
  outTree->Branch("trdKTotalPathLength",                     &trdKTotalPathLength,                     "trdKTotalPathLength/F");
  outTree->Branch("trdKTotalAmp",                            &trdKTotalAmp,                            "trdKTotalAmp/F");
  outTree->Branch("trdKElectronLikelihood",                  &trdKElectronLikelihood,                  "trdKElectronLikelihood/F");
  outTree->Branch("trdKProtonLikelihood",                    &trdKProtonLikelihood,                    "trdKProtonLikelihood/F");
  outTree->Branch("trdKHeliumLikelihood",                    &trdKHeliumLikelihood,                    "trdKHeliumLikelihood/F");
  outTree->Branch("trdPElectronToProtonLogLikelihoodRatio",  &trdPElectronToProtonLogLikelihoodRatio,  "trdPElectronToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdPHeliumToElectronLogLikelihoodRatio",  &trdPHeliumToElectronLogLikelihoodRatio,  "trdPHeliumToElectronLogLikelihoodRatio/F");
  outTree->Branch("trdPHeliumToProtonLogLikelihoodRatio",    &trdPHeliumToProtonLogLikelihoodRatio,    "trdPHeliumToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdPElectronLikelihood",                  &trdPElectronLikelihood,                  "trdPElectronLikelihood/F");
  outTree->Branch("trdPProtonLikelihood",                    &trdPProtonLikelihood,                    "trdPProtonLikelihood/F");
  outTree->Branch("trdPHeliumLikelihood",                    &trdPHeliumLikelihood,                    "trdPHeliumLikelihood/F");
  outTree->Branch("trdPNHybridHits",                         &trdPNHybridHits,                         "trdPNHybridHits/I");
  outTree->Branch("trdPNActiveLayers",                       &trdPNActiveLayers,                       "trdPNActiveLayers/I");
  outTree->Branch("trdPNLowAdcHits",                         &trdPNLowAdcHits,                         "trdPNLowAdcHits/I");
  outTree->Branch("trdPNLowDxHits",                          &trdPNLowDxHits,                          "trdPNLowDxHits/I");
  outTree->Branch("acTrackerTrdTrackDistanceAtTrdCenterInSigmas", &acTrackerTrdTrackDistanceAtTrdCenterInSigmas, "acTrackerTrdTrackDistanceAtTrdCenterInSigmas/F");
  outTree->Branch("acTrackerTrdTrackAngularDistanceInSigmas", &acTrackerTrdTrackAngularDistanceInSigmas, "acTrackerTrdTrackAngularDistanceInSigmas/F");
//  outTree->Branch("accNHits",                                &accNHits,                                "accNHits/I");
//  outTree->Branch("accNRecoClPG", &accNRecoClPG, "accNRecoClPG/I");
//  outTree->Branch("accSector", &accSector);
//  outTree->Branch("accTime", &accTime);
//  outTree->Branch("accHitPosZ", &accHitPosZ);
//  outTree->Branch("accChi2", &accChi2);
//  outTree->Branch("accNPairs", &accNPairs);
//  outTree->Branch("accUnfoldedHitPosZ", &accUnfoldedHitPosZ);
//  outTree->Branch("accRawCharge", &accRawCharge);
//  outTree->Branch("accHitPosZfromADC", &accHitPosZfromADC);
//  outTree->Branch("accUnfoldedHitPosZfromADC", &accUnfoldedHitPosZfromADC);
//  outTree->Branch("accCrossingPhiAngle", &accCrossingPhiAngle, "accCrossingPhiAngle/F");
//  outTree->Branch("accEdep", &accEdep);
//  outTree->Branch("accTimePG", &accTimePG);
  cout << "done." << endl;

}/*}}}*/

/**
 * @brief  This function do some preparation of the beginning process. 1. Bad run check., 2. Initialize TFile instance for output., 3. Initialize event counting histogram.
 * @return
 */
void KNUTree::Begin()
{/*{{{*/
  KNUOUT << "Checking bad run: " << endl;
  if( this->IsBadRun(pChain) ) this->End(1);
  if( isMC == false )
  {
    KNUOUT << "Checking science run: " << endl;
    if( !this->IsScienceRun(pChain) ) this->End(2);
  }
  KNUOUT << "Output file name: " << outputFileName << endl;
  pOutputTFile = new TFile(outputFileName, "RECREATE");
  outTree = new TTree("KNUTree", "KNUTree");

  if( isMC == false )
  {
    nCut = 7;
    hEvtCounter = new TH1F("hEvtCounter", "Number of survived events after cut", nCut, 0., (float)nCut);
    hEvtCounter->GetXaxis()->SetBinLabel(1, "No AMSEventR");
    hEvtCounter->GetXaxis()->SetBinLabel(2, "Bad run check");
    hEvtCounter->GetXaxis()->SetBinLabel(3, "Science-run check");
    hEvtCounter->GetXaxis()->SetBinLabel(4, "DAQ H/W status check");
    hEvtCounter->GetXaxis()->SetBinLabel(5, "Unbiased physics trigger check");
    hEvtCounter->GetXaxis()->SetBinLabel(6, "Tracker alignment test");
    hEvtCounter->GetXaxis()->SetBinLabel(7, "Good particle test");
    //hEvtCounter->GetXaxis()->SetBinLabel(8, "Survived");
    //hEvtCounter->GetXaxis()->SetBinLabel(9, "SAA rejection");
  }
  else
  {
    nCut = 3;
    hEvtCounter = new TH1F("hEvtCounter", "Number of survived events after cut", nCut, 0., (float)nCut);
    hEvtCounter->GetXaxis()->SetBinLabel(1, "No AMSEventR");
    hEvtCounter->GetXaxis()->SetBinLabel(2, "Bad run check");
    //hEvtCounter->GetXaxis()->SetBinLabel(3, "Unbiased physics trigger check");
    //hEvtCounter->GetXaxis()->SetBinLabel(4, "Multiple particle event rejection");
    //hEvtCounter->GetXaxis()->SetBinLabel(5, "Tracker alignment test");
    hEvtCounter->GetXaxis()->SetBinLabel(3, "Good particle test");
    //hEvtCounter->GetXaxis()->SetBinLabel(7, "Survived");
  }
  hGoodParticleCounter = new TH1F("hGoodParticleCounter", "Number of particle cuts", 7, 0, 7);
  TAxis* axisGoodParticleCounterX = (TAxis*)hGoodParticleCounter->GetXaxis();
  axisGoodParticleCounterX->SetBinLabel(1, "No pParticle");
  axisGoodParticleCounterX->SetBinLabel(2, "No pBetaH");
  axisGoodParticleCounterX->SetBinLabel(3, "No pTrTrack");
  axisGoodParticleCounterX->SetBinLabel(4, "No pTrdTrack");
  axisGoodParticleCounterX->SetBinLabel(5, "False GoodBeta");
  axisGoodParticleCounterX->SetBinLabel(6, "No GoodTrTrack");
  axisGoodParticleCounterX->SetBinLabel(7, "Multiple Good Particle");
}/*}}}*/

/**
 * @brief 
 * @return
 */
void KNUTree::End()
{/*{{{*/
  KNUOUT << "Finalizing..." << std::endl << std::endl;
  pOutputTFile->cd("/");
  outTree->Write();
  pOutputTFile->Write();
  pOutputTFile->Close();
  KNUOUT << "The program exited successfully" << std::endl << std::endl;
  KNUOUT << nEvents << " events were analyzed and " << nSelected << " events are selected" << std::endl;
}/*}}}*/

/**
 * @brief 
 * @return
 */
void KNUTree::End(int killarg)
{/*{{{*/
  KNUOUT << "The program exited successfully but it seems there was problem." << std::endl << std::endl;
  KNUOUT << "Finalizing due to kill sign... killarg: " << killarg << std::endl << std::endl;
  if( killarg == 1 || killarg == 2)
    exit(killarg);
  else
    pOutputTFile->Close();
  exit(killarg);
}/*}}}*/

/**
 * @brief 
 * @return
 */
void KNUTree::InitLoop()
{/*{{{*/
  pAMSEvent = NULL;
  pHeader = NULL;
  pParticle = NULL;
  pBetaH = NULL;
  pTrTrack = NULL;
  pTrRecHit = NULL;
  pTrdTrack = NULL;
  pTrdSegment = NULL;
  pEcalShower = NULL;
  pRichRing = NULL;
  pTrdKCluster = NULL;
  pTrdKHit = NULL;
  pAntiCluster = NULL;

  nLevel1 = 0;
  nParticle = 0;
  nCharge = 0;
  nTrTrack = 0;
  nTrdTrack = 0;
  nAntiCluster = 0;
  nTofClustersInTime = 0;
  nRichRing = 0;
  nRichRingB = 0;
  nBeta = 0;
  nBetaB = 0;
  nBetaH = 0;
  nEcalShower = 0;
  nVertex = 0;
  particleID = 0;
  trkID = 0;
  parentID = 0;
  for( int i = 0; i < 3; i++ )
  {
    genPosition[i] = 0;
    genDirCosine[i] = 0;
  }
  genMass = 0;
  genCharge = 0;
  genMomentum = 0;
  conversionId = 0;
  conversionIdContainer.clear();
  conversionIdSet.clear();
  particleType = 0;
  liveTime = 0;
  utcTime = 0;
  orbitAltitude = 0;
  orbitLatitude = 0;
  orbitLongitude = 0;
  orbitLatitudeM = 0;
  orbitLongitudeM = 0;
  velR = 0;
  velTheta = 0;
  velPhi = 0;
  yaw = 0;
  pitch = 0;
  roll = 0;
  gLongitude = 0;
  gLatitude = 0;
  gCoordCalcResult = 0;
  gCoordCalcResultBit = 0;
  sunPosAzimuth = 0;
  sunPosElevation = 0;
  sunPosCalcResult = 0;
  unixTime = 0;
  acEventTime = 0.;
  solarArrayCoord[0] = 0;
  solarArrayCoord[1] = 0;
  solarArrayCoord[2] = 0;
  isInShadow = 0;
  zenithAngle = 0;
  isInSAA = 0;
  //ptlCharge = 0;
  /*
  ptlMomentum = 0;
  ptlTheta = 0;
  ptlPhi = 0;
  ptlCoo[0]           = 0;
  ptlCoo[1]           = 0;
  ptlCoo[2]           = 0;
  ptlTofCoo[0][0]     = 0;
  ptlTofCoo[0][1]     = 0;
  ptlTofCoo[0][2]     = 0;
  ptlTofCoo[1][0]     = 0;
  ptlTofCoo[1][1]     = 0;
  ptlTofCoo[1][2]     = 0;
  ptlTofCoo[2][0]     = 0;
  ptlTofCoo[2][1]     = 0;
  ptlTofCoo[2][2]     = 0;
  ptlTofCoo[3][0]     = 0;
  ptlTofCoo[3][1]     = 0;
  ptlTofCoo[3][2]     = 0;
  ptlTofTrLength[0]   = 0;
  ptlTofTrLength[1]   = 0;
  ptlTofTrLength[2]   = 0;
  ptlTofTrLength[3]   = 0;
  ptlAntiCoo[0][0]    = 0;
  ptlAntiCoo[0][1]    = 0;
  ptlAntiCoo[0][2]    = 0;
  ptlAntiCoo[0][3]    = 0;
  ptlAntiCoo[0][4]    = 0;
  ptlAntiCoo[1][0]    = 0;
  ptlAntiCoo[1][1]    = 0;
  ptlAntiCoo[1][2]    = 0;
  ptlAntiCoo[1][3]    = 0;
  ptlAntiCoo[1][4]    = 0;
  ptlEcalCoo[0][0]    = 0;
  ptlEcalCoo[0][1]    = 0;
  ptlEcalCoo[0][2]    = 0;
  ptlEcalCoo[1][0]    = 0;
  ptlEcalCoo[1][1]    = 0;
  ptlEcalCoo[1][2]    = 0;
  ptlEcalCoo[2][0]    = 0;
  ptlEcalCoo[2][1]    = 0;
  ptlEcalCoo[2][2]    = 0;
  ptlTrCoo[0][0]      = 0;
  ptlTrCoo[0][1]      = 0;
  ptlTrCoo[0][2]      = 0;
  ptlTrCoo[1][0]      = 0;
  ptlTrCoo[1][1]      = 0;
  ptlTrCoo[1][2]      = 0;
  ptlTrCoo[2][0]      = 0;
  ptlTrCoo[2][1]      = 0;
  ptlTrCoo[2][2]      = 0;
  ptlTrCoo[3][0]      = 0;
  ptlTrCoo[3][1]      = 0;
  ptlTrCoo[3][2]      = 0;
  ptlTrCoo[4][0]      = 0;
  ptlTrCoo[4][1]      = 0;
  ptlTrCoo[4][2]      = 0;
  ptlTrCoo[5][0]      = 0;
  ptlTrCoo[5][1]      = 0;
  ptlTrCoo[5][2]      = 0;
  ptlTrCoo[6][0]      = 0;
  ptlTrCoo[6][1]      = 0;
  ptlTrCoo[6][2]      = 0;
  ptlTrCoo[7][0]      = 0;
  ptlTrCoo[7][1]      = 0;
  ptlTrCoo[7][2]      = 0;
  ptlTrCoo[8][0]      = 0;
  ptlTrCoo[8][1]      = 0;
  ptlTrCoo[8][2]      = 0;
  ptlTrdCoo[0][0]     = 0;
  ptlTrdCoo[0][1]     = 0;
  ptlTrdCoo[0][2]     = 0;
  ptlTrdCoo[1][0]     = 0;
  ptlTrdCoo[1][1]     = 0;
  ptlTrdCoo[1][2]     = 0;
  ptlRichCoo[0][0]    = 0;
  ptlRichCoo[0][1]    = 0;
  ptlRichCoo[0][2]    = 0;
  ptlRichCoo[1][0]    = 0;
  ptlRichCoo[1][1]    = 0;
  ptlRichCoo[1][2]    = 0;
  ptlRichPath[0]      = 0;
  ptlRichPath[1]      = 0;
  ptlRichPathBeta[0]  = 0;
  ptlRichPathBeta[1]  = 0; 
  ptlRichLength = 0;
  ptlRichParticles = 0;
  ptlCutOffStoermer = 0;
  ptlCutOffDipole = 0;
  */
  //ptlCutOffMax[0] = 0;
  //ptlCutOffMax[1] = 0;
  isEcalAvailable = 0;
  showerEnergyD = -1.;
  showerEnergyDL[0] = 0;
  showerEnergyDL[1] = 0;
  showerEnergyDL[2] = 0;
  showerEnergyDL[3] = 0;
  showerEnergyDL[4] = 0;
  showerEnergyDL[5] = 0;
  showerEnergyDL[6] = 0;
  showerEnergyDL[7] = 0;
  showerEnergyDL[8] = 0;
  showerEnergyDL[9] = 0;
  showerEnergyDL[10] = 0;
  showerEnergyDL[11] = 0;
  showerEnergyDL[12] = 0;
  showerEnergyDL[13] = 0;
  showerEnergyDL[14] = 0;
  showerEnergyDL[15] = 0;
  showerEnergyDL[16] = 0;
  showerEnergyDL[17] = 0;
  showerEnergyE = -1.;
  showerEnergyCorrected = -1.;
  showerBDT = -99.;
  showerCofG[0] = -99;
  showerCofG[1] = -99;
  showerCofG[2] = -99;
  showerCofGDist = -99;
  showerCofGdX = -99;
  showerCofGdY = -99;
  tofNCluster = 0;
  tofNClusterH = 0;
  tofNUsedHits = 0;
  tofNUnusedHits = 0;
  tofNUsedLayersForQ = 0;
  tofBeta = 0;
  tofInvBetaErr = 0;
  tofMass = 0;
  tofMassError = 0;
  isGoodBeta = 0;
  isTkTofMatch = 0;
  tofReducedChisqC = 0;
  tofReducedChisqT = 0;
  tofDepositedEnergyOnLayer[0] = 0;
  tofDepositedEnergyOnLayer[1] = 0;
  tofDepositedEnergyOnLayer[2] = 0;
  tofDepositedEnergyOnLayer[3] = 0;
  tofEstimatedChargeOnLayer[0] = 0;
  tofEstimatedChargeOnLayer[1] = 0;
  tofEstimatedChargeOnLayer[2] = 0;
  tofEstimatedChargeOnLayer[3] = 0;
  tofCharge = 0;
  tofChargeRMS = 0;
  tofUpperCharge = 0;
  tofLowerCharge = 0;
  tofChargeOnLayer[0] = 0;
  tofChargeOnLayer[1] = 0;
  tofChargeOnLayer[2] = 0;
  tofChargeOnLayer[3] = 0;
  trkFitCodeFS = -1;
  trkRigidityFS = -9999;
  trkRigidityInverseErrorFS = -9999;
  trkReducedChisquareFSX = -9999;
  trkReducedChisquareFSY = -9999;
  trkFitCodeMS = -1;
  trkRigidityMS = -9999;
  trkRigidityInverseErrorMS = -9999;
  trkReducedChisquareMSY = -9999;
  trkReducedChisquareMSX = -9999;
  trkFitCodeInner = -1;
  trkRigidityInverseErrorInner = -9999;
  trkReducedChisquareInnerX = -9999;
  trkReducedChisquareInnerY = -9999;
  trkLayerJQ[0] = -1;
  trkLayerJQ[1] = -1;
  trkLayerJQ[2] = -1;
  trkLayerJQ[3] = -1;
  trkLayerJQ[4] = -1;
  trkLayerJQ[5] = -1;
  trkLayerJQ[6] = -1;
  trkLayerJQ[7] = -1;
  trkLayerJQ[8] = -1;
  trkEdepLayerJXSideOK[0] = -1;
  trkEdepLayerJXSideOK[1] = -1;
  trkEdepLayerJXSideOK[2] = -1;
  trkEdepLayerJXSideOK[3] = -1;
  trkEdepLayerJXSideOK[4] = -1;
  trkEdepLayerJXSideOK[5] = -1;
  trkEdepLayerJXSideOK[6] = -1;
  trkEdepLayerJXSideOK[7] = -1;
  trkEdepLayerJXSideOK[8] = -1;
  trkEdepLayerJYSideOK[0] = -1;
  trkEdepLayerJYSideOK[1] = -1;
  trkEdepLayerJYSideOK[2] = -1;
  trkEdepLayerJYSideOK[3] = -1;
  trkEdepLayerJYSideOK[4] = -1;
  trkEdepLayerJYSideOK[5] = -1;
  trkEdepLayerJYSideOK[6] = -1;
  trkEdepLayerJYSideOK[7] = -1;
  trkEdepLayerJYSideOK[8] = -1;
  trkEdepLayerJ[0] = -1;
  trkEdepLayerJ[1] = -1;
  trkEdepLayerJ[2] = -1;
  trkEdepLayerJ[3] = -1;
  trkEdepLayerJ[4] = -1;
  trkEdepLayerJ[5] = -1;
  trkEdepLayerJ[6] = -1;
  trkEdepLayerJ[7] = -1;
  trkEdepLayerJ[8] = -1;
  trkCharge = -1;
  trkInnerCharge = -1;
  trkHasExtLayers = -1;
  isRichAvailable = 0;
  richRebuild = -1;
  richIsGood = -1;
  richIsClean = -1;
  richIsNaF = -1;
  richUsedHits = -1;
  richRingWidth = -1;
  richNHits = -1;
  richNPMTsOnRing = -1;
  richBeta = -1;
  richBetaError = -1;
  richChargeSquared = -1;
  richKolmogorovProbability = -1;
  richPhotoelectrons = -1;
  richPhotoElectrons = -1;
  richExpectedPhotoelectrons = -1;
  richExpectedPhotoElectrons = -1;
  richNUsedHits = -1;
  richTheta = -9;
  richPhi = -9;
  trdNClusters = 0;
  trdNUnusedHits = 0;
  trdNTracks = 0;
  trdNVertex = 0;
  trdTrackPhi = -9.;
  trdTrackTheta = -9.;
  trdTrackChi2 = -9.;
  trdTrackPattern = -9.;
  trdTrackCharge = -9.;
  trdTrackEdepL[0] = -1;
  trdTrackEdepL[1] = -1;
  trdTrackEdepL[2] = -1;
  trdTrackEdepL[3] = -1;
  trdTrackEdepL[4] = -1;
  trdTrackEdepL[5] = -1;
  trdTrackEdepL[6] = -1;
  trdTrackEdepL[7] = -1;
  trdTrackEdepL[8] = -1;
  trdTrackEdepL[9] = -1;
  trdTrackEdepL[10] = -1;
  trdTrackEdepL[11] = -1;
  trdTrackEdepL[12] = -1;
  trdTrackEdepL[13] = -1;
  trdTrackEdepL[14] = -1;
  trdTrackEdepL[15] = -1;
  trdTrackEdepL[16] = -1;
  trdTrackEdepL[17] = -1;
  trdTrackEdepL[18] = -1;
  trdTrackEdepL[19] = -1;
  trdTrackMeanDepositedEnergy  = -1;
  trdTrackTotalDepositedEnergy = -1;
  trdQtNActiveLayer = 0;
  trdQtIsCalibrationGood = 0;
  trdQtIsInsideTrdGeometricalAcceptance = 0;
  trdQtIsValid = 0;
  trdQtNActiveStraws = 0;
  trdQtNActiveLayers = 0;
  trdQtNTRDVertex = 0;
  trdQtElectronToProtonLogLikelihoodRatio = 0;
  trdQtHeliumToElectronLogLikelihoodRatio = 0;
  trdQtHeliumToProtonLogLikelihoodRatio = 0;
  trdKNRawHits                                = -1;
  trdKIsReadAlignmentOK                       = -1;
  trdKIsReadCalibOK                           = -1;
  trdKNHits                                   = -1;
  trdKIsValid                                 = -1;
  trdKElectronToProtonLogLikelihoodRatio      = -1;
  trdKHeliumToElectronLogLikelihoodRatio      = -1;
  trdKHeliumToProtonLogLikelihoodRatio        = -1;
  trdKCharge                                  = -1;
  trdKChargeError                             = -1;
  trdKNUsedHitsForCharge                      = -1;
  for(int k = 0; k < 20; k++) trdKAmpLayer[k] = -1;
  trdKTotalPathLength                         = -1;
  trdKElectronLikelihood                      = -1;
  trdKProtonLikelihood                        = -1;
  trdKHeliumLikelihood                        = -1;
  trdPElectronToProtonLogLikelihoodRatio = 0;
  trdPHeliumToElectronLogLikelihoodRatio = 0;
  trdPHeliumToProtonLogLikelihoodRatio = 0;
  trdPElectronLikelihood = 0;
  trdPProtonLikelihood = 0;
  trdPHeliumLikelihood = 0;
  trdPNHybridHits = 0;
  accNHits = 0;
  accNRecoClPG = 0;
  accSector.clear();
  accTime.clear();
  accHitPosZ.clear();
  accChi2.clear();
  accUnfoldedHitPosZ.clear();
  accRawCharge.clear();
  accHitPosZfromADC.clear();
  accUnfoldedHitPosZfromADC.clear();
  accCrossingPhiAngle = 0;
  accEdep.clear();
  accTimePG.clear();
}/*}}}*/

/**
 * @brief Initialize class members to prepare new particle object.
 * @return no return
 */
void KNUTree::InitParticle()
{/*{{{*/
  /*
  particleID = 0;
  trkID = 0;
  parentID = 0;
  for( int i = 0; i < 3; i++ )
  {
    genPosition[i] = 0;
    genDirCosine[i] = 0;
  }
  genMass = 0;
  genCharge = 0;
  genMomentum = 0;
  particleType = 0;
  liveTime = 0;
  utcTime = 0;
  orbitAltitude = 0;
  orbitLatitude = 0;
  orbitLongitude = 0;
  orbitLatitudeM = 0;
  orbitLongitudeM = 0;
  velR = 0;
  velTheta = 0;
  velPhi = 0;
  yaw = 0;
  pitch = 0;
  roll = 0;
  gLongitude = 0;
  gLatitude = 0;
  gCoordCalcResult = 0;
  gCoordCalcResultBit = 0;
  sunPosAzimuth = 0;
  sunPosElevation = 0;
  sunPosCalcResult = 0;
  unixTime = 0;
  solarArrayCoord[0] = 0;
  solarArrayCoord[1] = 0;
  solarArrayCoord[2] = 0;
  isInShadow = 0;
  zenithAngle = 0;
  isInSAA = 0;
  */
  //ptlCharge = 0;
  /*
  ptlMomentum = 0;
  ptlTheta = 0;
  ptlPhi = 0;
  ptlCoo[0]           = 0;
  ptlCoo[1]           = 0;
  ptlCoo[2]           = 0;
  ptlTofCoo[0][0]     = 0;
  ptlTofCoo[0][1]     = 0;
  ptlTofCoo[0][2]     = 0;
  ptlTofCoo[1][0]     = 0;
  ptlTofCoo[1][1]     = 0;
  ptlTofCoo[1][2]     = 0;
  ptlTofCoo[2][0]     = 0;
  ptlTofCoo[2][1]     = 0;
  ptlTofCoo[2][2]     = 0;
  ptlTofCoo[3][0]     = 0;
  ptlTofCoo[3][1]     = 0;
  ptlTofCoo[3][2]     = 0;
  ptlTofTrLength[0]   = 0;
  ptlTofTrLength[1]   = 0;
  ptlTofTrLength[2]   = 0;
  ptlTofTrLength[3]   = 0;
  ptlAntiCoo[0][0]    = 0;
  ptlAntiCoo[0][1]    = 0;
  ptlAntiCoo[0][2]    = 0;
  ptlAntiCoo[0][3]    = 0;
  ptlAntiCoo[0][4]    = 0;
  ptlAntiCoo[1][0]    = 0;
  ptlAntiCoo[1][1]    = 0;
  ptlAntiCoo[1][2]    = 0;
  ptlAntiCoo[1][3]    = 0;
  ptlAntiCoo[1][4]    = 0;
  ptlEcalCoo[0][0]    = 0;
  ptlEcalCoo[0][1]    = 0;
  ptlEcalCoo[0][2]    = 0;
  ptlEcalCoo[1][0]    = 0;
  ptlEcalCoo[1][1]    = 0;
  ptlEcalCoo[1][2]    = 0;
  ptlEcalCoo[2][0]    = 0;
  ptlEcalCoo[2][1]    = 0;
  ptlEcalCoo[2][2]    = 0;
  ptlTrCoo[0][0]      = 0;
  ptlTrCoo[0][1]      = 0;
  ptlTrCoo[0][2]      = 0;
  ptlTrCoo[1][0]      = 0;
  ptlTrCoo[1][1]      = 0;
  ptlTrCoo[1][2]      = 0;
  ptlTrCoo[2][0]      = 0;
  ptlTrCoo[2][1]      = 0;
  ptlTrCoo[2][2]      = 0;
  ptlTrCoo[3][0]      = 0;
  ptlTrCoo[3][1]      = 0;
  ptlTrCoo[3][2]      = 0;
  ptlTrCoo[4][0]      = 0;
  ptlTrCoo[4][1]      = 0;
  ptlTrCoo[4][2]      = 0;
  ptlTrCoo[5][0]      = 0;
  ptlTrCoo[5][1]      = 0;
  ptlTrCoo[5][2]      = 0;
  ptlTrCoo[6][0]      = 0;
  ptlTrCoo[6][1]      = 0;
  ptlTrCoo[6][2]      = 0;
  ptlTrCoo[7][0]      = 0;
  ptlTrCoo[7][1]      = 0;
  ptlTrCoo[7][2]      = 0;
  ptlTrCoo[8][0]      = 0;
  ptlTrCoo[8][1]      = 0;
  ptlTrCoo[8][2]      = 0;
  ptlTrdCoo[0][0]     = 0;
  ptlTrdCoo[0][1]     = 0;
  ptlTrdCoo[0][2]     = 0;
  ptlTrdCoo[1][0]     = 0;
  ptlTrdCoo[1][1]     = 0;
  ptlTrdCoo[1][2]     = 0;
  ptlRichCoo[0][0]    = 0;
  ptlRichCoo[0][1]    = 0;
  ptlRichCoo[0][2]    = 0;
  ptlRichCoo[1][0]    = 0;
  ptlRichCoo[1][1]    = 0;
  ptlRichCoo[1][2]    = 0;
  ptlRichPath[0]      = 0;
  ptlRichPath[1]      = 0;
  ptlRichPathBeta[0]  = 0;
  ptlRichPathBeta[1]  = 0; 
  ptlRichLength = 0;
  ptlRichParticles = 0;
  ptlCutOffStoermer = 0;
  ptlCutOffDipole = 0;
  */
  //ptlCutOffMax[0] = 0;
  //ptlCutOffMax[1] = 0;
  isEcalAvailable = 0;
  showerEnergyD = -1.;
  showerEnergyDL[0] = 0;
  showerEnergyDL[1] = 0;
  showerEnergyDL[2] = 0;
  showerEnergyDL[3] = 0;
  showerEnergyDL[4] = 0;
  showerEnergyDL[5] = 0;
  showerEnergyDL[6] = 0;
  showerEnergyDL[7] = 0;
  showerEnergyDL[8] = 0;
  showerEnergyDL[9] = 0;
  showerEnergyDL[10] = 0;
  showerEnergyDL[11] = 0;
  showerEnergyDL[12] = 0;
  showerEnergyDL[13] = 0;
  showerEnergyDL[14] = 0;
  showerEnergyDL[15] = 0;
  showerEnergyDL[16] = 0;
  showerEnergyDL[17] = 0;
  showerEnergyE = -1.;
  showerEnergyCorrected = -1.;
  showerBDT = -99.;
  showerCofG[0] = -99;
  showerCofG[1] = -99;
  showerCofG[2] = -99;
  showerCofGDist = -99;
  showerCofGdX = -99;
  showerCofGdY = -99;
  tofNCluster = 0;
  tofNClusterH = 0;
  tofNUsedHits = 0;
  tofNUnusedHits = 0;
  tofNUsedLayersForQ = 0;
  tofBeta = 0;
  tofInvBetaErr = 0;
  tofMass = 0;
  tofMassError = 0;
  isGoodBeta = 0;
  isTkTofMatch = 0;
  tofReducedChisqC = 0;
  tofReducedChisqT = 0;
  tofDepositedEnergyOnLayer[0] = 0;
  tofDepositedEnergyOnLayer[1] = 0;
  tofDepositedEnergyOnLayer[2] = 0;
  tofDepositedEnergyOnLayer[3] = 0;
  tofEstimatedChargeOnLayer[0] = 0;
  tofEstimatedChargeOnLayer[1] = 0;
  tofEstimatedChargeOnLayer[2] = 0;
  tofEstimatedChargeOnLayer[3] = 0;
  tofCharge = 0;
  tofChargeRMS = 0;
  tofUpperCharge = 0;
  tofLowerCharge = 0;
  tofChargeOnLayer[0] = 0;
  tofChargeOnLayer[1] = 0;
  tofChargeOnLayer[2] = 0;
  tofChargeOnLayer[3] = 0;
  trkFitCodeFS = -1;
  trkRigidityFS = 0;
  trkRigidityInverseErrorFS = 0;
  trkReducedChisquareFSX = 0;
  trkReducedChisquareFSY = 0;
  trkFitCodeMS = -1;
  trkRigidityMS = 0;
  trkRigidityInverseErrorMS = 0;
  trkReducedChisquareMSY = 0;
  trkReducedChisquareMSX = 0;
  trkFitCodeInner = -1;
  trkRigidityInverseErrorInner = 0;
  trkReducedChisquareInnerX = 0;
  trkReducedChisquareInnerY = 0;
  trkLayerJQ[0] = 0;
  trkLayerJQ[1] = 0;
  trkLayerJQ[2] = 0;
  trkLayerJQ[3] = 0;
  trkLayerJQ[4] = 0;
  trkLayerJQ[5] = 0;
  trkLayerJQ[6] = 0;
  trkLayerJQ[7] = 0;
  trkLayerJQ[8] = 0;
  trkEdepLayerJXSideOK[0] = 0;
  trkEdepLayerJXSideOK[1] = 0;
  trkEdepLayerJXSideOK[2] = 0;
  trkEdepLayerJXSideOK[3] = 0;
  trkEdepLayerJXSideOK[4] = 0;
  trkEdepLayerJXSideOK[5] = 0;
  trkEdepLayerJXSideOK[6] = 0;
  trkEdepLayerJXSideOK[7] = 0;
  trkEdepLayerJXSideOK[8] = 0;
  trkEdepLayerJYSideOK[0] = 0;
  trkEdepLayerJYSideOK[1] = 0;
  trkEdepLayerJYSideOK[2] = 0;
  trkEdepLayerJYSideOK[3] = 0;
  trkEdepLayerJYSideOK[4] = 0;
  trkEdepLayerJYSideOK[5] = 0;
  trkEdepLayerJYSideOK[6] = 0;
  trkEdepLayerJYSideOK[7] = 0;
  trkEdepLayerJYSideOK[8] = 0;
  trkEdepLayerJ[0] = 0;
  trkEdepLayerJ[1] = 0;
  trkEdepLayerJ[2] = 0;
  trkEdepLayerJ[3] = 0;
  trkEdepLayerJ[4] = 0;
  trkEdepLayerJ[5] = 0;
  trkEdepLayerJ[6] = 0;
  trkEdepLayerJ[7] = 0;
  trkEdepLayerJ[8] = 0;
  trkCharge = 0;
  trkInnerCharge = 0;
  trkHasExtLayers = 0;
  isRichAvailable = 0;
  richRebuild = -1;
  richIsGood = -1;
  richIsClean = -1;
  richIsNaF = -1;
  richUsedHits = -1;
  richRingWidth = -1;
  richNHits = -1;
  richNPMTsOnRing = -1;
  richBeta = -1;
  richBetaError = -1;
  richChargeSquared = -1;
  richKolmogorovProbability = -1;
  richPhotoelectrons = -1;
  richPhotoElectrons = -1;
  richExpectedPhotoelectrons = -1;
  richExpectedPhotoElectrons = -1;
  richNUsedHits = -1;
  richTheta = -9;
  richPhi = -9;
  trdNClusters = 0;
  trdNUnusedHits = 0;
  trdNTracks = 0;
  trdNVertex = 0;
  trdTrackPhi = -9.;
  trdTrackTheta = -9.;
  trdTrackChi2 = -9.;
  trdTrackPattern = -9.;
  trdTrackCharge = -9.;
  trdTrackEdepL[0] = -1;
  trdTrackEdepL[1] = -1;
  trdTrackEdepL[2] = -1;
  trdTrackEdepL[3] = -1;
  trdTrackEdepL[4] = -1;
  trdTrackEdepL[5] = -1;
  trdTrackEdepL[6] = -1;
  trdTrackEdepL[7] = -1;
  trdTrackEdepL[8] = -1;
  trdTrackEdepL[9] = -1;
  trdTrackEdepL[10] = -1;
  trdTrackEdepL[11] = -1;
  trdTrackEdepL[12] = -1;
  trdTrackEdepL[13] = -1;
  trdTrackEdepL[14] = -1;
  trdTrackEdepL[15] = -1;
  trdTrackEdepL[16] = -1;
  trdTrackEdepL[17] = -1;
  trdTrackEdepL[18] = -1;
  trdTrackEdepL[19] = -1;
  trdTrackMeanDepositedEnergy  = -1;
  trdTrackTotalDepositedEnergy = -1;
  trdQtNActiveLayer = 0;
  trdQtIsCalibrationGood = 0;
  trdQtIsInsideTrdGeometricalAcceptance = 0;
  trdQtIsValid = 0;
  trdQtNActiveStraws = 0;
  trdQtNActiveLayers = 0;
  trdQtNTRDVertex = 0;
  trdQtElectronToProtonLogLikelihoodRatio = 0;
  trdQtHeliumToElectronLogLikelihoodRatio = 0;
  trdQtHeliumToProtonLogLikelihoodRatio = 0;
  trdKNRawHits                                = -1;
  trdKIsReadAlignmentOK                       = -1;
  trdKIsReadCalibOK                           = -1;
  trdKNHits                                   = -1;
  trdKIsValid                                 = -1;
  trdKElectronToProtonLogLikelihoodRatio      = -1;
  trdKHeliumToElectronLogLikelihoodRatio      = -1;
  trdKHeliumToProtonLogLikelihoodRatio        = -1;
  trdKCharge                                  = -1;
  trdKChargeError                             = -1;
  trdKNUsedHitsForCharge                      = -1;
  for(int k = 0; k < 20; k++) trdKAmpLayer[k] = -1;
  trdKTotalPathLength                         = -1;
  trdKElectronLikelihood                      = -1;
  trdKProtonLikelihood                        = -1;
  trdKHeliumLikelihood                        = -1;
  trdPElectronToProtonLogLikelihoodRatio = 0;
  trdPHeliumToElectronLogLikelihoodRatio = 0;
  trdPHeliumToProtonLogLikelihoodRatio = 0;
  trdPElectronLikelihood = 0;
  trdPProtonLikelihood = 0;
  trdPHeliumLikelihood = 0;
  trdPNHybridHits = 0;
  trdPNActiveLayers = 0;
  trdPNLowAdcHits = 0;
  trdPNLowDxHits = 0;
  acTrackerTrdTrackDistanceAtTrdCenterInSigmas = 0;
  acTrackerTrdTrackAngularDistanceInSigmas = 0;
  accNHits = 0;
  accNRecoClPG = 0;
  accSector.clear();
  accTime.clear();
  accHitPosZ.clear();
  accChi2.clear();
  accUnfoldedHitPosZ.clear();
  accRawCharge.clear();
  accHitPosZfromADC.clear();
  accUnfoldedHitPosZfromADC.clear();
  accCrossingPhiAngle = 0;
  accEdep.clear();
  accTimePG.clear();
}/*}}}*/
