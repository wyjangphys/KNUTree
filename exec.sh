#!/bin/bash
source /afs/cern.ch/user/w/wyjang/amsvar_gcc.sh
source /afs/cern.ch/user/w/wyjang/AMS-ACsoft/scripts/thisacsoft.sh

# example command
# /afs/cern.ch/user/w/wyjang/work/KNUTree/main FileList/set_1105/1305853512 Output/set_1105/1305853512.root
/afs/cern.ch/user/w/wyjang/work/KNUTree/main $1 $2
