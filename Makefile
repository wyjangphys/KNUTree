CXX=g++
CXXFLAGS = -W -Wall -Wno-write-strings -Wno-extra -Wno-unused-parameter -Wno-unknown-pragmas -DAMS_ACQT_INTERFACE -D_PGTRACK_ -Wno-unused-variable -Wno-extra -Wno-unused-function

# AMS Global Environment
CVMFS_AMS_OFFLINE = /cvmfs/ams.cern.ch/Offline

############################# CERN libraries ##################################
CERNLIBS = -lgeant321 -lpacklib -lmathlib -lkernlib
CERN_LEVEL = 2005.gcc64

ifndef CERNDIR
CERNDIR = $(CVMFS_AMS_OFFLINE)/CERN/$(CERN_LEVEL)
endif

CERNSRCDIR = $(CERNDIR)

ifndef AMSLIB
AMSLIB = /$(CVMFS_AMS_OFFLINE)/lib/linux/gcc64
endif

ifndef NAGDIR
NAGDIR = $(CVMFS_AMS_OFFLINE)/CERN/NagLib
endif
######################### End of CERN library settings ########################

# AMS Offline Software Related Includes
INCLUDES = -I${ROOTSYS}/include -I${AMSWD}/include -I./include
NTUPLE_PG = $(AMSWD)/lib/linuxx8664gcc5.34/ntuple_slc6_PG.so
############ End of AMS Offline Software related includes

# ROOT Related Settings
ROOTLIBS = $(shell root-config --libs) -lASImage -lRIO -lNet -lNetx -lMinuit -lTMVA -lMLP -lXMLIO -lTreePlayer
############ End of ROOT related settings

# ACSOFT Flags
#ACSOFTFLAGS = `acsoft-config --definitions` `acsoft-config --cflags` `acsoft-config --auxcflags` `acsoft-config --libs` `acsoft-config --auxlibs`
ACSOFTFLAGS = `acsoft-config --definitions` `acsoft-config --cflags` `acsoft-config --cflags-include` `acsoft-config --auxcflags-include` `acsoft-config --libs` `acsoft-config --auxlibs`

TARGET = main

all : $(TARGET)

$(TARGET) : main.o KNUTree.o KNUTreeLoop.o KNUTreeGetGoodParticleIndex.o KNUTreeIsACCPatternGood.o KNUTreeIsBadRun.o KNUTreeIsGoodBeta.o KNUTreeIsGoodTrTrack.o KNUTreeIsHardwareStatusGood.o KNUTreeIsScienceRun.o KNUTreeIsShowerTrackMatched.o KNUTreeIsTrackInsideEcalFiducialVolume.o KNUTreeIsTrkAlignmentGood.o KNUTreeIsUnbiasedPhysicsTriggerEvent.o
	$(CXX) $(CXXFLAGS) $(ACSOFTFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -o $@ $^ $(NTUPLE_PG)

main.o : main.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTree.o : KNUTree.cxx
	$(CXX) $(CXXFLAGS) $(ACSOFTFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeLoop.o : KNUTreeLoop.cxx
	$(CXX) $(CXXFLAGS) $(ACSOFTFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsBadRun.o : KNUTreeIsBadRun.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeGetGoodParticleIndex.o : KNUTreeGetGoodParticleIndex.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsACCPatternGood.o : KNUTreeIsACCPatternGood.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsGoodBeta.o : KNUTreeIsGoodBeta.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsGoodTrTrack.o : KNUTreeIsGoodTrTrack.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsHardwareStatusGood.o : KNUTreeIsHardwareStatusGood.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsScienceRun.o : KNUTreeIsScienceRun.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsShowerTrackMatched.o : KNUTreeIsShowerTrackMatched.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsTrackInsideEcalFiducialVolume.o : KNUTreeIsTrackInsideEcalFiducialVolume.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsTrkAlignmentGood.o : KNUTreeIsTrkAlignmentGood.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

KNUTreeIsUnbiasedPhysicsTriggerEvent.o : KNUTreeIsUnbiasedPhysicsTriggerEvent.cxx
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -c $^ -o $@

clean :
	rm -rf main
	rm -rf main.o
	rm -rf KNUTree.o
	rm -rf KNUTreeLoop.o
	rm -rf KNUTreeIsBadRun.o
	rm -rf KNUTreeGetGoodParticleIndex.o
	rm -rf KNUTreeIsACCPatternGood.o
	rm -rf KNUTreeIsGoodBeta.o
	rm -rf KNUTreeIsGoodTrTrack.o
	rm -rf KNUTreeIsHardwareStatusGood.o
	rm -rf KNUTreeIsScienceRun.o
	rm -rf KNUTreeIsShowerTrackMatched.o
	rm -rf KNUTreeIsTrackInsideEcalFiducialVolume.o
	rm -rf KNUTreeIsTrkAlignmentGood.o
	rm -rf KNUTreeIsUnbiasedPhysicsTriggerEvent.o
# DO NOT DELETE
#KNUTree.o: KNUTree.cxx KNUTree.h
#KNUTreeLoop.o: KNUTreeLoop.cxx KNUTree.h
#KNUTreeIsBadRun.o: KNUTreeIsBadRun.cxx KNUTree.h
#KNUTreeGetGoodParticleIndex.o: KNUTreeGetGoodParticleIndex.cxx KNUTree.h
#KNUTreeIsACCPatternGood.o: KNUTreeIsACCPatternGood.cxx KNUTree.h
#KNUTreeIsGoodBeta.o: KNUTreeIsGoodBeta.cxx KNUTree.h
#KNUTreeIsGoodTrTrack.o: KNUTreeIsGoodTrTrack.cxx KNUTree.h
#KNUTreeIsHardwareStatusGood.o: KNUTreeIsHardwareStatusGood.cxx KNUTree.h
#KNUTreeIsScienceRun.o: KNUTreeIsScienceRun.cxx KNUTree.h
#KNUTreeIsShowerTrackMatched.o: KNUTreeIsShowerTrackMatched.cxx KNUTree.h
#KNUTreeIsTrackInsideEcalFiducialVolume.o: KNUTreeIsTrackInsideEcalFiducialVolume.cxx KNUTree.h
#KNUTreeIsTrkAlignmentGood.o: KNUTreeIsTrkAlignmentGood.cxx KNUTree.h
#KNUTreeIsUnbiasedPhysicsTriggerEvent.o: KNUTreeIsUnbiasedPhysicsTriggerEvent.cxx KNUTree.h
