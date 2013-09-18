#!/bin/bash
cd /afs/cern.ch/work/l/lpernie/ECALpro/CMSSW_4_2_4/src/Analysis/Modules/test
export SCRAM_ARCH=slc5_amd64_gcc434
eval `scramv1 runtime -sh`
echo 'python /afs/cern.ch/work/l/lpernie/ECALpro/CMSSW_4_2_4/src/Analysis/Modules/test/RunOnBatch.py'
/afs/cern.ch/work/l/lpernie/ECALpro/CMSSW_4_2_4/src/Analysis/Modules/test/RunOnBatch.py
