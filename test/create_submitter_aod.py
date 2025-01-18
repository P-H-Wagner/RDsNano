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
nMaxJobs = 800

#default
filesPerJob = 3

if (int(args.nFiles) < nMaxJobs) and (int(args.nFiles) != -1):
  filesPerJob = 1 #below 500 jobs we can take 1 file per job and thus short
  queue = 'short' 
  time = 60
else:
  queue = 'standard'
  time = 3*60
print("========> processing ", filesPerJob, " files per job")

######################################

#queue = 'standard'
#time = 60
nevents = -1
nFiles = int(args.nFiles)
now = datetime.now()
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")

######################################


def filesFromFolder(direc):
  filenames = os.listdir(direc)
  return ['file:' + direc + filename for filename in filenames ], filenames

def filesFromTxt(txtFile):
  with open(txtFile) as dataFiles: 
    direcs = [line[0:-1] for line in dataFiles] #-2 to remove the newline character \n

    filenames = []
    for direc in direcs: 
      path_parts = direc.split("/")
      # Extract the last two parts
      last_two_parts = path_parts[-2:]
      # Merge them with an underscore
      new_name = "_".join(last_two_parts)
      filenames.append(new_name)

  return ['file:root://cms-xrd-global.cern.ch//' + direc for direc in direcs], filenames

#######################################

# Input source

if args.channel == 'sig':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/signals/all_signals_request_21_11_23.txt' # new signals!!
  inputfiles, filenames = filesFromTxt(directory)
  naming = 'all_signals'
  folder = "mc/signals"

if args.channel == 'hb':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1/' #hb 
  inputfiles, filenames = filesFromFolder(directory)
  naming = 'hb_inclusive'
  folder = "mc/hb/inclusive"

if args.channel == 'b0':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/b0/b0.txt' #b0 
  inputfiles, filenames = filesFromTxt(directory)
  naming = 'b0'
  folder = "mc/hb/b0"

if args.channel == 'bplus':

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bplus/bplus.txt' #b0 
  inputfiles, filenames = filesFromTxt(directory)
  naming = 'bplus'
  folder = "mc/hb/bplus"

if args.channel == 'bs':

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bs/bs.txt' #bs 
  inputfiles, filenames = filesFromTxt(directory)
  naming = 'bs'
  folder = "mc/hb/bs"

if args.channel == 'lambdab':

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/lambdab/lambdab.txt' #bs 
  inputfiles, filenames = filesFromTxt(directory)
  naming = 'lambdab'
  folder = "mc/hb/lambdab"

if args.channel == 'data':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_D/BPark_2018D.txt' #data bParking 2018 part D
  inputfiles , filenames= filesFromTxt(directory)
  naming = 'data'
  folder = "data"


#####################################

if nFiles != -1:
  #process not the whole dataset but only nFiles
  inputfiles = inputfiles[0:nFiles] #50 files give ca 200k events

save_path = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/" + folder + "/refit/"

if os.path.exists(save_path):
  print("Attention! possible files are already in path " + save_path)

else:
  os.makedirs(save_path)

log_path = "aod_adder_" + dt_string
os.makedirs( log_path + "/logs")
os.makedirs( log_path + "/errs")

for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
  #print(fin)
  #template
  temp = open("addAODTracks_cfg_temp.py", "rt")
  #file to write to
  cfg = open(log_path + "/cfg_chunk_{0}.py".format(i),"wt")
  #file to save things
  fout = "/scratch/pahwagne/{0}/{1}".format(log_path,filenames[i])  

  for line in temp:
    if   "HOOK_INPUT" in line: cfg.write(line.replace("HOOK_INPUT", directory))
    elif "HOOK_N_EVENTS" in line: cfg.write(line.replace("HOOK_N_EVENTS", str(nevents)))
    elif "HOOK_FILE_IN" in line: cfg.write(line.replace("HOOK_FILE_IN", str(fin)))
    elif "HOOK_FILE_OUT" in line: cfg.write(line.replace("HOOK_FILE_OUT", fout))
    else: cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/CMSSW_10_6_37/src/PhysicsTools/RDsNano/test/' + log_path ,
         'scramv1 runtime -sh',
         'mkdir -p /scratch/pahwagne/' + log_path,
         'ls /scratch/pahwagne/',
         'cmsRun cfg_chunk_{0}.py'.format(i),
         'xrdcp /scratch/pahwagne/{0}/{1} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/{2}/refit/{1}'.format(log_path,filenames[i],folder),
         'rm /scratch/pahwagne/{0}/{1}'.format(log_path,filenames[i]),
         '',
     ])

  with open("{0}/submitter_chunk_{1}.sh".format(log_path,i), "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        '-p '+queue,
        '--account=t3',
        '-o {0}/logs/chunk_{1}.log'.format(log_path,i),
        '-e {0}/errs/chunk_{1}.err'.format(log_path,i),
        #'--mem=1200M',
        '--job-name=AOD_{0}_{1}'.format(i,args.channel),
        '--time={0}'.format(time),
        '{0}/submitter_chunk_{1}.sh'.format(log_path,i),
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







