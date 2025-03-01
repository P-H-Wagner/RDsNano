import os
from glob import glob
from datetime import datetime
import sys 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('channel') # sig or hb or data
parser.add_argument('nFiles')  # nr of miniAOD files to process
parser.add_argument('-f', '--failedCrab') #date_time of crab, if given, the ntuplizer will run on minis which failed during crab
args = parser.parse_args()


######################################

#800 jobs per user on short queue
nMaxJobs = 1000

#default
filesPerJob = 3

if (int(args.nFiles) < nMaxJobs) and (int(args.nFiles) != -1):
  filesPerJob = 1 #below 500 jobs we can take 1 file per job and thus short
  queue = 'standard' 
  time = 3*60
else:
  queue = 'standard'
  time = 8*60
print("========> processing ", filesPerJob, " files per job")

######################################

#queue = 'standard'
#time = 60
nevents = -1
nFiles = int(args.nFiles)
now = datetime.now()
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")

######################################

#binned qcd samples
#from https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-23-016
qcd_samples = ["50to100", "100to200", "200to300", "300to500", "500to700", "700to1000", "1000to1500", "1500to2000", "2000toInf"]

#####################################


def filesFromFolder(direc):
  filenames = os.listdir(direc)
  return ['file:' + direc + filename for filename in filenames ]

def filesFromTxt(txtFile):
  with open(txtFile) as dataFiles: 
    filenames = [line[0:-1] for line in dataFiles] #-2 to remove the newline character \n
  return ['file:root://cms-xrd-global.cern.ch//' + filename for filename in filenames ]

#######################################

# Input source

if args.channel == 'sig':
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/' #old signal from MA thesis
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/signals/all_signals_request_21_11_23.txt' # new signals!!
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1_PV_REFITTED/'
  #inputfiles = filesFromFolder(directory)
  inputfiles = filesFromTxt(directory)
  naming = 'all_signals'

if args.channel == 'hb':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1/' #hb 
  inputfiles = filesFromFolder(directory)
  naming = 'hb_inclusive'

if args.channel == 'b0':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/b0/b0.txt' #b0 
  inputfiles = filesFromTxt(directory)

  naming = 'b0'

if args.channel == 'bplus':

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bplus/bplus_updated.txt' #b0 
  inputfiles = filesFromTxt(directory)
  naming = 'bplus'

if args.channel == 'bs':

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bs/bs_updated.txt' #bs 
  inputfiles = filesFromTxt(directory)
  naming = 'bs'

if args.channel == 'lambdab':

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/lambdab/lambdab_updated.txt' #lambda 
  inputfiles = filesFromTxt(directory)
  naming = 'lambdab'

if args.channel == 'bc':

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bc/bc.txt' #bc 
  inputfiles = filesFromTxt(directory)
  naming = 'bc'

if args.channel == 'data':
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_D/BPark_2018D.txt' #data bParking 2018 part D
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_D/mu9.txt' #data bParking 2018 part D
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/runTest/files_run_325117_lumi_172.txt' # only one run
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_D/ten_parts/BPark_2018D_part5.txt' #data bParking 2018 part D
  #txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/dataTest/single.txt' # test
  inputfiles = filesFromTxt(directory)
  #inputfiles = [inputfiles[0]] #debugging only first file
  naming = 'data'

if args.channel in qcd_samples:

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/qcd/binned/binned_' + args.channel + '.txt'
  inputfiles = filesFromTxt(directory)
  naming = 'binned_' + args.channel

######################################

if args.failedCrab is not None:
  print("Re-batching failed crab jobs!")

  dt_string  = args.failedCrab # take the same naming as in crab
  directory  = "/work/pahwagne/CMSSW_10_6_37/src/PhysicsTools/RDsNano/production/failedJobMinis/"+dt_string+"/all_minis_job.txt" 
  inputfiles = filesFromTxt(directory)
  naming     = 'data_toRetry'


#####################################

if nFiles != -1:
  #process not the whole dataset but only nFiles
  inputfiles = inputfiles[0:nFiles] #50 files give ca 200k events

os.makedirs("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/"+dt_string)
os.makedirs(dt_string+"/logs")
os.makedirs(dt_string+"/errs")

for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
  #print(fin)
  #template
  temp = open("temp_cfg.py", "rt")
  #file to write to
  cfg = open(dt_string+"/cfg_chunk_{0}.py".format(i),"wt")
  #file to save things
  fout = "/scratch/pahwagne/{0}/{1}_chunk_{2}.root".format(dt_string,naming,i)  

  for line in temp:
    if   "HOOK_CHANNEL" in line: cfg.write(line.replace("HOOK_CHANNEL", args.channel))
    elif "HOOK_INPUT" in line: cfg.write(line.replace("HOOK_INPUT", directory))
    elif "HOOK_N_EVENTS" in line: cfg.write(line.replace("HOOK_N_EVENTS", str(nevents)))
    elif "HOOK_FILE_IN" in line: cfg.write(line.replace("HOOK_FILE_IN", str(fin)))
    elif "HOOK_FILE_OUT" in line: cfg.write(line.replace("HOOK_FILE_OUT", fout))
    else: cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/CMSSW_10_6_37/src/PhysicsTools/RDsNano/test/'+dt_string ,
         'scramv1 runtime -sh',
         'mkdir -p /scratch/pahwagne/'+dt_string,
         'ls /scratch/pahwagne/',
         'cmsRun cfg_chunk_{1}.py'.format(dt_string,i),
         'xrdcp /scratch/pahwagne/{0}/{1}_chunk_{2}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/{0}/{1}_chunk_{2}.root'.format(dt_string,naming,i),
         'rm /scratch/pahwagne/{0}/{1}_chunk_{2}.root'.format(dt_string,naming,i),
         '',
     ])

  with open("{0}/submitter_chunk_{1}.sh".format(dt_string,i), "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        '-p '+queue,
        '--account=t3',
        '-o {0}/logs/chunk_{1}.log'.format(dt_string,i),
        '-e {0}/errs/chunk_{1}.err'.format(dt_string,i),
        #'--mem=1200M',
        '--job-name=MINI_{0}_{1}'.format(i,args.channel),
        '--time={0}'.format(time),
        '{0}/submitter_chunk_{1}.sh'.format(dt_string,i),
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







