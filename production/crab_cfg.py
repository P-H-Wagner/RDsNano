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
  config.Data.inputDataset       = '/ParkingBPH5/Run2018D-UL2018_MiniAODv2-v1/MINIAOD' # dataset
  #config.Data.inputDataset      = '/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIISummer20UL18MiniAODv2-Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' # signal
  #config.Data.inputBlocks       = ['/ParkingBPH1/Run2018D-UL2018_MiniAODv2-v1/MINIAOD#007c47d2-b6e8-433a-9b32-5effbf767314'] #process only these files (format: ['dataset#id'])
  config.Data.publication        = False       # dont save on the DAS
  config.Data.outLFNDirBase      = '/store/user/pahwagne/%s' % date_time # destination on T2
  config.Data.inputDBS           = 'global'
  config.Data.allowNonValidInputDataset = True # allow to process data which is not in CMS valid state yet
  #config.Data.totalUnits        = 200 #200 events when splitting automatic
  config.Data.splitting          = 'FileBased' # split via files
  config.Data.unitsPerJob        = 10           # 1 file per job
else:
  # ONLY WHEN TESTING SINGLE FILE
  config.Data.inputDataset = None
  config.Data.userInputFiles = ['file:root://eoscms.cern.ch//eos/cms/store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2530009/33EBB973-79DA-114A-9C72-1CC48E0ED7C6.root']  # Replace with your input file
  config.Data.splitting = 'FileBased'
  config.Data.unitsPerJob = 1
  config.Data.totalUnits = 1
  config.Data.publication        = False       # dont save on the DAS
  config.Data.outLFNDirBase      = '/store/user/pahwagne/%s' % date_time # destination on T2
  config.Data.inputDBS           = 'global'
  config.Data.allowNonValidInputDataset = True # allow to process data which is not in CMS valid state yet

# job settings
config.section_('JobType')
config.JobType.pluginName       = 'Analysis'       #naming
config.JobType.psetName         =  '../test/run_crab.py' #cmssw cfg file to run
config.JobType.inputFiles       = ['../test/run_crab.py' ] #'file1.py','file2.txt','file3.yml'] #specify here additional input files
#config.JobType.psetName        = 'PSet.py'
#config.JobType.scriptExe       = 'crab_script.sh' # only needed for more complex scrips (more input files and diff. environment)

config.JobType.maxJobRuntimeMin = 1972  #set job time to 1.5 * default for data processing ( ~ 30h, i.e. 3h per data file, since we submit in batches of 10)
config.JobType.allowUndistributedCMSSW = True #allow to run on single node
config.JobType.outputFiles      = ["test_crab.root"]

config.section_('User')

# specify output site
config.section_('Site')
config.Site.storageSite = 'T2_CH_CSCS'



