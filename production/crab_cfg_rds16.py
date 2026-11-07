#!/usr/bin/env python3
"""CRAB submission for the RDs16 production (CMSSW_16_X, test/run_crab.py). One request per sample.

  python3 crab_cfg_rds16.py --pilot  [--dry] [--samples data,sig,hb]   # pilot (few files per sample)
  python3 crab_cfg_rds16.py --full   [--dry] [--samples data,sig,hb]   # full production (only after Ric's OK)

Global tags: PAT-step tag from edmProvDump of the first file of each sample (see LOG.md).
Job size from the stage-4 timing (CPU/event data 0.029 s, MC 0.041 s; RSS <= 1.6 GB):
~3 h nominal per job, limit 8 h.
"""
import argparse, datetime, os
from CRABClient.UserUtilities import config as crab_config
from CRABAPI.RawCommand import crabCommand

HERE = os.path.dirname(os.path.abspath(__file__))
PSET = os.path.join(HERE, '..', 'test', 'run_crab.py')
GOLDEN = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

SAMPLES = {
  #  name : (dataset, isMC, globalTag, files/job, pilot files)
  'data': ('/ParkingBPH1/Run2018D-UL2018_MiniAODv2-v1/MINIAOD', False, '106X_dataRun2_v35', 6, 24),
  'sig' : ('/BsToDsL_DsToPhiKKPi_TauToMuNuNu_TuneCP5_13TeV_pythia8-evtgen/RunIISummer20UL18MiniAODv2-Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM',
           True, '106X_upgrade2018_realistic_v16_L1v1', 8, 8),
  'hb'  : ('/HbToDsMu_DsToPhiKKPi_TuneCP5_13TeV_pythia8-evtgen/RunIISummer20UL18MiniAODv2-Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',
           True, '106X_upgrade2018_realistic_v16_L1v1', 10, 10),
}

def make_config(name, pilot, tag):
    dataset, is_mc, gt, files_per_job, pilot_files = SAMPLES[name]
    c = crab_config()
    c.General.requestName     = 'rds16_%s_%s%s' % (name, tag, '_pilot' if pilot else '')
    c.General.workArea        = '/work/manzoni/rds/rds16/crab_projects'
    c.General.transferOutputs = True
    c.General.transferLogs    = True

    c.JobType.pluginName       = 'Analysis'
    c.JobType.psetName         = PSET
    c.JobType.pyCfgParams      = ['isMC=%s' % is_mc, 'globalTag=%s' % gt]
    c.JobType.outputFiles      = ['rds_nano.root']      # = outputFile default of run_crab.py
    c.JobType.maxJobRuntimeMin = 480
    c.JobType.maxMemoryMB      = 2000
    c.JobType.numCores         = 1
    c.JobType.allowUndistributedCMSSW = True

    c.Data.inputDataset  = dataset
    c.Data.inputDBS      = 'global'
    c.Data.splitting     = 'FileBased'
    c.Data.unitsPerJob   = files_per_job
    c.Data.totalUnits    = pilot_files if pilot else -1
    c.Data.publication   = False
    c.Data.outLFNDirBase = '/store/user/manzoni/RDs/%s_rds16' % tag
    if not is_mc:
        c.Data.lumiMask = GOLDEN

    c.Site.storageSite = 'T2_CH_CSCS'
    return c

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--pilot', action='store_true')
    g.add_argument('--full', action='store_true')
    p.add_argument('--samples', default='data,sig,hb')
    p.add_argument('--tag', default=datetime.date.today().strftime('%Y%m%d'), help='YYYYMMDD for outLFNDirBase / requestName')
    p.add_argument('--dry', action='store_true', help='crab submit --dryrun (no jobs)')
    a = p.parse_args()
    for s in a.samples.split(','):
        cfg = make_config(s, a.pilot, a.tag)
        print('submitting', cfg.General.requestName, cfg.Data.inputDataset, 'unitsPerJob', cfg.Data.unitsPerJob, 'totalUnits', cfg.Data.totalUnits)
        print(crabCommand('submit', config=cfg, dryrun=a.dry))
