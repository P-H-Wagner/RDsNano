#!/bin/bash
cd /work/pahwagne/CMSSW_10_6_37/src/PhysicsTools/RDsNano/test
scramv1 runtime -sh
mkdir -p /scratch/pahwagne/nanoAOD
ls /scratch/pahwagne/
cmsRun run_nano_rds_cfg.py
xrdcp /scratch/pahwagne/nanoAOD/test.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/test.root
rm /scratch/pahwagne/nanoAOD/test.root
