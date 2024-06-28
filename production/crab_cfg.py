from CRABClient.UserUtilities import config, ClientException
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser

date_time = datetime.date.today().strftime('%Y%b%d')

# create instace
config = config()

# set general parameters
config.section_('General')
config.General.transferOutputs = True        # transfer to T2
config.General.transferLogs    = True        # save logs
config.General.workArea        = date_time   # save logs at

# set in- and output data parameters
config.section_('Data')
config.Data.inputDataset       = '/store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/5D552231-6643-F042-9DD8-4CFC0CFD1B26.root' 
config.Data.publication        = False       # dont save on the DAS
config.Data.outLFNDirBase      = '/store/user/pahwagne/%s' % date_time # destination on T2
config.Data.inputDBS           = 'global'
config.Data.allowNonValidInputDataset = True # allow to process data which is not in CMS valid state yet
#config.Data.splitting          = 'FileBased' # split via files
#config.Data.unitsPerJob        = 1           # 1 file per job

# job settings
config.section_('JobType')
config.JobType.pluginName      = 'Analysis'       #naming
config.JobType.psetName        = '../test/run_crab.py' #cmssw cfg file to run
#config.JobType.inputFiles      = ['../test/run_crab.py' ] #'file1.py','file2.txt','file3.yml'] #specify here additional input files
#config.JobType.psetName       = 'PSet.py'
config.JobType.scriptExe      = 'crab_script.sh' # only needed for more complex scrips (more input files and diff. environment)
#config.JobType.maxJobRuntimeMin = 3000
config.JobType.allowUndistributedCMSSW = True #allow to run on single node
config.JobType.outputFiles     = ["test_crab.root"]

config.section_('User')

# specify output site
config.section_('Site')
config.Site.storageSite = 'T2_CH_CSCS'


"""
if __name__ == '__main__':

  from CRABAPI.RawCommand import crabCommand
  from CRABClient.ClientExceptions import ClientException
  from httplib import HTTPException
  from multiprocessing import Process

  def submit(config):
      try:
          crabCommand('submit', config = config)
      except HTTPException as hte:
          print "Failed submitting task: %s" % (hte.headers)
      except ClientException as cle:
          print "Failed submitting task: %s" % (cle)


  parser = ArgumentParser()
  parser.add_argument('-y', '--yaml', default = 'samples_mc_rjpsi.yml', help = 'File with dataset descriptions')
  parser.add_argument('-f', '--filter', default='*', help = 'filter samples, POSIX regular expressions allowed')
  args = parser.parse_args()

  with open(args.yaml) as f:
    doc = yaml.load(f) # Parse YAML file
    common = doc['common'] if 'common' in doc else {'data' : {}, 'mc' : {}}
    
    # loop over samples
    for sample, info in doc['samples'].iteritems():
      # Given we have repeated datasets check for different parts
      parts = info['parts'] if 'parts' in info else [None]
      for part in parts:
        name = sample % part if part is not None else sample
        
        # filter names according to what we need
        if not fnmatch(name, args.filter): continue
        print 'submitting', name

        isMC = info['isMC']
        config.Data.inputDataset = info['dataset'] % part \
                                   if part is not None else \
                                      info['dataset']

        #config.General.requestName = name + '2'
        config.General.requestName = name 
        common_branch = 'mc' if isMC else 'data'
        config.Data.splitting = 'FileBased' if isMC else 'LumiBased'
        if not isMC:
            config.Data.lumiMask = info.get(
                'lumimask', 
                common[common_branch].get('lumimask', None)
            )
        else:
            config.Data.lumiMask = ''

        config.Data.unitsPerJob = info.get(
            'splitting',
            common[common_branch].get('splitting', None)
        )
        globaltag = info.get(
            'globaltag',
            common[common_branch].get('globaltag', None)
        )
        
        config.JobType.pyCfgParams = [
            'isMC=%s' % isMC, 'reportEvery=1000',
            'tag=%s' % production_tag,
            'globalTag=%s' % globaltag,
        ]
        
        config.JobType.outputFiles = ['_'.join(['RJPsi', 'mc' if isMC else 'data','pu', production_tag])+'.root']
        #config.JobType.outputFiles = ['_'.join(['RJPsi', 'mc' if isMC else 'data', production_tag])+'.root']
        
        print config
        submit(config)
"""
