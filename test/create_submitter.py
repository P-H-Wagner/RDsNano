import os
from glob import glob
from datetime import datetime
import sys 


queue = 'short' 
time = 60
nevents = -1
now = datetime.now()
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")


filenames = os.listdir('/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/')
inputfiles = ['file:/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/' + filename for filename in filenames ]
inputfiles = inputfiles[0:2]

os.makedirs("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/"+dt_string)
os.makedirs(dt_string+"/logs")
os.makedirs(dt_string+"/errs")

for i, fin in enumerate(inputfiles):

  #template
  temp = open("temp_cfg.py", "rt")
  #file to write to
  cfg = open(dt_string+"/cfg_chunk_{0}.py".format(i),"wt")
  #file to save things
  fout = "/scratch/pahwagne/{0}/all_signals_chunk_{1}.root".format(dt_string,i)  

  for line in temp:
    if "HOOK_N_EVENTS" in line: cfg.write(line.replace("HOOK_N_EVENTS", str(nevents)))
    elif "HOOK_FILE_IN" in line: cfg.write(line.replace("HOOK_FILE_IN", fin))
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
         'xrdcp /scratch/pahwagne/{0}/all_signals_chunk_{1}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/{0}/all_signals_chunk_{1}.root'.format(dt_string,i),
         'rm /scratch/pahwagne/{0}/all_signals_chunk_{1}.root'.format(dt_string,i),
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

        '--job-name=NANO_{0}'.format(i),
        '--time={0}'.format(time),
        '{0}/submitter_chunk_{1}.sh'.format(dt_string,i),
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)






