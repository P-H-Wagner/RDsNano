from CRABClient.UserUtilities import config, ClientException
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser


test = False 

date_time = datetime.date.today().strftime('%Y%b%d')

# create instace
config = config()

# set general parameters
config.section_('General')
config.General.transferOutputs = True        # transfer to T2
config.General.transferLogs    = True        # save logs
config.General.workArea        = date_time   # save logs at

# set in- and output data parameters
if not test:
  config.section_('Data')
  config.Data.inputDataset       = '/ParkingBPH1/Run2018D-UL2018_MiniAODv2-v1/MINIAOD' # dataset
  #config.Data.inputDataset      = '/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIISummer20UL18MiniAODv2-Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' # signal
  #config.Data.inputBlocks       = ['/ParkingBPH1/Run2018D-UL2018_MiniAODv2-v1/MINIAOD#007c47d2-b6e8-433a-9b32-5effbf767314'] #process only these files (format: ['dataset#id'])
  config.Data.publication        = False       # dont save on the DAS
  config.Data.outLFNDirBase      = '/store/user/pahwagne/%s' % date_time # destination on T2
  config.Data.inputDBS           = 'global'
  config.Data.allowNonValidInputDataset = True # allow to process data which is not in CMS valid state yet
  config.Data.splitting          = 'FileBased' # split via files
  config.Data.unitsPerJob        = 1   # 5 file per job
  config.Data.totalUnits         = 10000 

else:
  config.Data.inputDataset = '/ParkingBPH1/Run2018D-UL2018_MiniAODv2-v1/MINIAOD'
  config.Data.splitting    = 'FileBased'
  config.Data.unitsPerJob  = 1      # 1 file per job
  config.Data.totalUnits   = 1      # only 1 file in total -> 1 job
  config.Data.publication  = False
  config.Data.outLFNDirBase = '/store/user/pahwagne/%s' % date_time
  config.Data.inputDBS     = 'global'
  config.Data.allowNonValidInputDataset = True


# job settings
config.section_('JobType')
config.JobType.pluginName       = 'Analysis'       #naming
config.JobType.psetName         =  '../test/run_crab.py' #cmssw cfg file to run
config.JobType.inputFiles       = ['../test/run_crab.py' ] #'file1.py','file2.txt','file3.yml'] #specify here additional input files
#config.JobType.psetName        = 'PSet.py'
#config.JobType.scriptExe       = 'crab_script.sh' # only needed for more complex scrips (more input files and diff. environment)

#config.JobType.maxJobRuntimeMin = 2000 #set job time to 1.5 * default for data processing ( ~ 30h, i.e. 3h per data file, since we submit in batches of 10)
config.JobType.maxJobRuntimeMin = 2000 #set job time to 1.5 * default for data processing ( ~ 30h, i.e. 3h per data file, since we submit in batches of 10)

config.JobType.allowUndistributedCMSSW = True #allow to run on single node
config.JobType.outputFiles      = ["test_crab.root"]

config.section_('User')

# specify output site
config.section_('Site')
config.Site.storageSite = 'T2_CH_CSCS'



