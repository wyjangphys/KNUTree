#ifndef _KNUTREE_H_
#define _KNUTREE_H_

#include "amschain.h"
#include "TrdKCluster.h"
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

#define KNUOUT std::cout << "[KNUTree::Loop] "
#define KNUERR std::cerr << "[KNUTree::Loop] Error: "

class KNUTree
{
  public:
    KNUTree(char* _name, int argc, char* argv[]);
    ~KNUTree();
    char* GetName();
    bool SetName(char* _name);
    bool RegisterFileToChain(int argc, char* argv[], AMSChain* amsChain);
    void Init();
    void Begin();
    void Loop();
    void End();
    void End(int killarg);

    bool IsBadRun(AMSEventR* thisEvent);
    bool IsScienceRun(AMSEventR* thisEvent);
    bool IsHardwareStatusGood(AMSEventR* thisEvent);
    bool IsUnbiasedPhysicsTriggerEvent(AMSEventR* thisEvent);
    bool IsACCPatternGood(AMSEventR* thisEvent);
    bool IsGoodBeta(BetaR* thisBeta);
    bool IsGoodLiveTime(AMSEventR* thisEvent);
    bool IsInSouthAtlanticAnomaly(AMSEventR* thisEvent);
    bool IsInSolarArrays(AMSEventR* thisEvent);
    bool IsGoodTrTrack(TrTrackR* thisTrack);
    bool IsShowerTrackMatched(EcalShowerR* thisShower, TrTrackR* thisTrack);
    bool IsTrackInsideEcalFiducialVolume(TrTrackR* thisTrack);
    bool IsTrkAlignmentGood(AMSEventR* thisEvent);
    int GetGoodParticleIndex(AMSEventR* thisEvent);

  private:
    char* name;
    char* outputFileName;
    AMSChain* pChain;
    TTree* outTree; /// TTree for output
    TFile* pOutputTFile;
    TH1F* hEvtCounter;
    unsigned int nTreesInChain;   /// Number of TTree's in current AMSChain
    unsigned long int nevt;

  public:
    unsigned int  nCut;                       // Number of cuts
    unsigned int  nProcessed;                 // Number of processed events
    unsigned int  nRun;                       // Run number
    unsigned int  nEvent;                     // Event number
    unsigned int  nProcessedNumber;           // Number of processed events
    unsigned int  nLevel1;                    // Number of Level1 triggers
    unsigned int  nParticle;                  // Number of particles
    unsigned int  nCharge;                    // Number of charges
    unsigned int  nTrTrack;                   // Number of the Tracker tracks which are successfully reconstructed
    unsigned int  nTrdTrack;                  // Number of the TRD tracks which are successfully reconstructed
    unsigned int  nAntiCluster;               // Number of clusters on the ACC
    unsigned int  nTofClustersInTime;         // Number of in-time clusters on the TOF
    unsigned int  nRichRing;                  // Number of successfully reconstructed the RICH rings
    unsigned int  nRichRingB;                 // Number of successfully reconstructed the RICH rings with algorithm B
    unsigned int  nBeta;                      // Number of successfully estimated beta(v/c) values
    unsigned int  nBetaB;                     // Number of successfully estimated beta(v/c) values with algorithm B
    unsigned int  nBetaH;                     // Number of successfully estimated beta(v/c) values with algorithm H
    unsigned int  nShower;                    // Number of the ECAL shower objects
    unsigned int  nVertex;                    // Number of vertices in this event
    unsigned int  particleType;               // Type of ParticleR
    float         liveTime;                   // Livetime fraction
    float         utcTime;                    // UTC time
    float         orbitAltitude;              // (cm) in GTOD coordinates system.
    float         orbitLatitude;              // (rad) in GTOD coordinates system.
    float         orbitLongitude;             // (rad) in GTOD coordinates system.
    float         orbitLatitudeM;             // (rad) in eccentric dipole coordinate system.
    float         orbitLongitudeM;            // (rad) in eccentric dipole coordinate system.
    float         velR;                       // Speed of the ISS in radial direction
    float         velTheta;                   // Angular speed of the ISS in polar angle direction
    float         velPhi;                     // Angular speed of the ISS in azimuthal angle direction
    float         yaw;                        // A parameter describing ISS attitude (tilted angle with respect to the x-axis)
    float         pitch;                      // A parameter describing ISS attitude (tilted angle with respect to the y-axis)
    float         roll;                       // A parameter describing ISS attitude (tilted angle with respect to the z-axis)
    float         gLongitude;                 // Galactic longitude of the incoming particle.
    float         gLatitude;                  // Galactic latitude of the incoming particle.
    int           gCoordCalcResult;           // Return value for galactic coordinate calculation.
    int           gCoordCalcResultBit;        // Result of GetGalCoo written in bit form
    float         sunPosAzimuth;              // Azimuthal angle of the position of the Sun.
    float         sunPosElevation;            // Elevation angle of the position of the Sun.
    int           sunPosCalcResult;           // Return value for the Sun's position calculation.
    unsigned int  unixTime;                   // UNIX time
    float         solarArrayCoord[3];
    int           isInShadow;                 // Value for check whether the AMS is in ISS solar panel shadow or not.
    float         zenithAngle;
    int           isInSAA;
    unsigned int  ptlCharge;                  // ParticleR::Charge value
    float         ptlMomentum;                // ParticleR::Momentum value
    float         ptlTheta;                   // Direction of the incoming particle (polar angle)
    float         ptlPhi;                     // Direction of the incoming particle (azimuthal angle)
    float         ptlCoo[3];                  // Coo (1st[last tracker plane) cm
    float         ptlTofCoo[4][3];            // Track extrapolated point [TOF Layer][x/y/z]
    float         ptlTofTrLength[4];          // Track length till TOF planes crossing
    float         ptlAntiCoo[2][5];           // Track extrapolated in ACC [Direction][x/y/z/theta/phi]
    float         ptlEcalCoo[3][3];           // Track extrapolated in ECAL [Entering position/ Center of gravity/ Exiting position][x/y/z]
    float         ptlTrCoo[9][3];             // Track extrapolated in Tracker [Layer][x/y/z]
    float         ptlTrdCoo[2][3];            // Track extrapolated in TRD [at center / at top][x/y/z]
    float         ptlRichCoo[2][3];           // Track extrapolated in RICH [at radiator / at PMT][x/y/z]
    float         ptlRichPath[2];
    float         ptlRichPathBeta[2];
    float         ptlRichLength;              // Estimated pathlength of particle within RICH radiator (cm).
    int           ptlRichParticles;
    float         ptlCutOffStoermer;
    float         ptlCutOffDipole;
    float         ptlCutOffMax[2];
    float         showerEnergyD;
    float         showerEnergyDL[18];
    float         showerEnergyE;
    float         showerEnergyCorrected;
    float         showerBDT;
    float         showerCofG[3];
    float         showerCofGDist;
    float         showerCofGdX;
    float         showerCofGdY;
    int           tofNCluster;
    int           tofNClusterH;
    int           tofNUsedHits;
    int           tofNUnusedHits;
    int           tofNUsedLayersForQ;
    float         tofBeta;
    float         tofInvBetaErr;
    float         tofMass;
    float         tofMassError;
    int           isGoodBeta;
    int           isTkTofMatch;
    float         tofReducedChisqT;
    float         tofReducedChisqC;
    float         tofDepositedEnergyOnLayer[4];
    float         tofEstimatedChargeOnLayer[4];
    float         tofCharge;
    float         tofUpperCharge;
    float         tofLowerCharge;
    float         tofChargeOnLayer[4];
    int           trkFitCodeFS;
    float         trkRigidityFS;
    float         trkRigidityInverseErrorFS;
    float         trkReducedChisquareFSX;
    float         trkReducedChisquareFSY;
    int           trkFitCodeMS;
    float         trkRigidityMS;
    float         trkRigidityInverseErrorMS;
    float         trkReducedChisquareMSX;
    float         trkReducedChisquareMSY;
    int           trkFitCodeInner;
    float         trkRigidityInner;
    float         trkRigidityInverseErrorInner;
    float         trkReducedChisquareInnerX;
    float         trkReducedChisquareInnerY;
    int           trkEdepLayerJXSideOK[9];
    int           trkEdepLayerJYSideOK[9];
    float         trkEdepLayerJ[9];
    float         trkCharge;
    float         trkInnerCharge;
    int           trkHasExtLayers;
    int           richRebuild;
    int           richIsGood;
    int           richIsClean;
    int           richIsNaF;
    int           richUsedHits;
    float         richRingWidth;
    int           richNHits;
    int           richNPMTsOnRing;
    float         richBeta;
    float         richBetaError;
    float         richChargeSquared;
    float         richKolmogorovProbability;
    float         richPhotoelectrons;                         // Return value from RichRingR::getPhotoelectrons(bool corr=true)
    float         richPhotoElectrons;                         // Return value from RichRingR::getPhotoElectrons(int pmt, bool corr)
    float         richExpectedPhotoelectrons;                 // Return value from RichRingR::getExpectedPhotoElectrons(bool corr=true)
    float         richExpectedPhotoElectrons;                 // Return value from RichRingR::getExpectedPhotoElectrons(int pmt, bool corr)
    int           richNUsedHits;
    float         richTheta;
    float         richPhi;
    int           trdNClusters;
    int           trdNUnusedHits;
    int           trdNTracks;
    float         trdTrackPhi;
    float         trdTrackTheta;
    float         trdTrackChi2;
    int           trdTrackPattern;
    float         trdTrackCharge;
    float         trdTrackEdepL[20];
    float         trdDepositedEnergyOnLayer[20];
    float         trdTrackMeanDepositedEnergy;
    float         trdTrackTotalDepositedEnergy;
    int           trdQtNActiveLayer;
    int           trdQtIsCalibrationGood;
    int           trdQtIsSlowControlDataGood;
    int           trdQtIsInsideTrdGeometricalAcceptance;
    int           trdQtIsValid;
    int           trdQtNActiveStraws;
    int           trdQtNActiveLayers;
    int           trdQtNTRDVertex;
    float         trdQtElectronToProtonLogLikelihoodRatio;
    float         trdQtHeliumToElectronLogLikelihoodRatio;
    float         trdQtHeliumToProtonLogLikelihoodRatio;
    int           trdKNRawHits;
    int           trdKIsReadAlignmentOK;
    int           trdKIsReadCalibOK;
    int           trdKNHits;
    int           trdKIsValid;
    float         trdKElectronToProtonLogLikelihoodRatio;
    float         trdKHeliumToElectronLogLikelihoodRatio;
    float         trdKHeliumToProtonLogLikelihoodRatio;
    float         trdKCharge;
    float         trdKChargeError;
    int           trdKNUsedHitsForCharge;
    float         trdKAmpLayer[20];
    float         trdKTotalPathLength;
    float         trdKTotalAmp;
    float         trdKElectronLikelihood;
    float         trdKProtonLikelihood;
    float         trdKHeliumLikelihood;
    float         trdPElectronToProtonLogLikelihoodRatio;
    float         trdPHeliumToElectronLogLikelihoodRatio;
    float         trdPHeliumToProtonLogLikelihoodRatio;
    int           accNHits;
    int           accNRecoClPG;
    vector<int>   accSector;
    vector<float> accTime;
    vector<float> accHitPosZ;
    vector<float> accChi2;
    vector<int>   accNPairs;
    vector<float> accUnfoldedHitPosZ;
    vector<float> accRawCharge;
    vector<float> accHitPosZfromADC;
    vector<float> accUnfoldedHitPosZfromADC;
    float         accCrossingPhiAngle;
    vector<float> accEdep;
    vector<float> accTimePG;
  protected:
};
#endif
