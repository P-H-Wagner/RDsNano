#!/usr/bin/env python3
"""CRAB submission for the RDs16 production (CMSSW_16_X, test/run_crab.py). One request per sample.

  python3 crab_cfg_rds16.py --pilot  [--dry] [--samples data,sig,hb]   # pilot (few files per sample)
  --samples qcd = all QCD MuEnrichedPt5 bins (qcd15To20, ..., qcd1000)
  python3 crab_cfg_rds16.py --full   [--dry] [--samples data,sig,hb]   # full production (only after Ric's OK)

Global tags: PAT-step tag from edmProvDump of the first file of each sample (see LOG.md).
Job size from the V8 timing with the final code (val/v8, 5000 ev; CPU/event data 0.030 s, signal/Hb 0.043 s,
QCD 0.024-0.029 s; local RSS 1.5-1.6 GB, QCD 1.8 GB) and the pilot (grid ~2.5x slower than local, grid RSS ~1.3x local):
~6 h median per job on the grid, limit 20 h, 2500 MB.
"""
import argparse, datetime, os
from CRABClient.UserUtilities import config as crab_config
from CRABAPI.RawCommand import crabCommand

HERE = os.path.dirname(os.path.abspath(__file__))
PSET = os.path.join(HERE, '..', 'test', 'run_crab.py')
GOLDEN = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

SAMPLES = {
  #  name : (dataset, isMC, globalTag, files/job, pilot files)
  'data': ('/ParkingBPH1/Run2018D-UL2018_MiniAODv2-v1/MINIAOD', False, '106X_dataRun2_v35', 6, 24),   # 53k ev/file
  'sig' : ('/BsToDsL_DsToPhiKKPi_TauToMuNuNu_TuneCP5_13TeV_pythia8-evtgen/RunIISummer20UL18MiniAODv2-Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM',
           True, '106X_upgrade2018_realistic_v16_L1v1', 6, 6),     # 33k ev/file
  'hb'  : ('/HbToDsMu_DsToPhiKKPi_TuneCP5_13TeV_pythia8-evtgen/RunIISummer20UL18MiniAODv2-Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',
           True, '106X_upgrade2018_realistic_v16_L1v1', 8, 8),     # 26k ev/file
}
# QCD MuEnrichedPt5 (no Custom_RDStarPU_BParking version exists in DAS): standard UL18 MiniAODv2, same PAT global tag
# as signal/Hb (106X_upgrade2018_realistic_v16_L1v1), BParking HLT paths in the menu, standard 2018 pile-up profile.
QCD_BINS = ['15To20', '20To30', '30To50', '50To80', '80To120', '120To170', '170To300',
            '300To470', '470To600', '600To800', '800To1000', '1000']
# ~290k events per job (0.075 s/ev on the grid); DAS events/file: 15To20 94k, 20To120 77-81k, 120To170 74k, >=170 53-55k
QCD_FILES_PER_JOB = {'15To20': 3, '20To30': 4, '30To50': 4, '50To80': 4, '80To120': 4, '120To170': 4,
                     '170To300': 5, '300To470': 5, '470To600': 5, '600To800': 5, '800To1000': 5, '1000': 5}
for b in QCD_BINS:
    SAMPLES['qcd' + b] = ('/QCD_Pt-%s_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' % b,
                          True, '106X_upgrade2018_realistic_v16_L1v1', QCD_FILES_PER_JOB.get(b, 2), QCD_FILES_PER_JOB.get(b, 2))

def make_config(name, pilot, tag):
    dataset, is_mc, gt, files_per_job, pilot_files = SAMPLES[name]
    c = crab_config()
    c.General.requestName     = 'rds16_%s_%s%s' % (name, tag, '_pilot2' if pilot else '')
    c.General.workArea        = '/work/manzoni/rds/rds16/crab_projects'
    c.General.transferOutputs = True
    c.General.transferLogs    = True

    c.JobType.pluginName       = 'Analysis'
    c.JobType.psetName         = PSET
    c.JobType.pyCfgParams      = ['isMC=%s' % is_mc, 'globalTag=%s' % gt]
    c.JobType.outputFiles      = ['rds_nano.root']      # = outputFile default of run_crab.py
    c.JobType.maxJobRuntimeMin = 1200
    c.JobType.maxMemoryMB      = 2500
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

def submit(cfg, dry):
    print('submitting', cfg.General.requestName, cfg.Data.inputDataset, 'unitsPerJob', cfg.Data.unitsPerJob, 'totalUnits', cfg.Data.totalUnits, flush=True)
    print(crabCommand('submit', config=cfg, dryrun=dry), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--pilot', action='store_true')
    g.add_argument('--full', action='store_true')
    p.add_argument('--samples', default='data,sig,hb')
    p.add_argument('--tag', default=datetime.date.today().strftime('%Y%m%d'), help='YYYYMMDD for outLFNDirBase / requestName')
    p.add_argument('--dry', action='store_true', help='crab submit --dryrun (no jobs)')
    a = p.parse_args()
    # one process per request: the CRAB API caches the CMSSW config inside a process
    from multiprocessing import Process
    names = []
    for s in a.samples.split(','):
        names += ['qcd' + b for b in QCD_BINS] if s == 'qcd' else [s]
    for s in names:
        proc = Process(target=submit, args=(make_config(s, a.pilot, a.tag), a.dry))
        proc.start(); proc.join()
