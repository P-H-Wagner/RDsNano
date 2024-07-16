import os
from glob import glob
from datetime import datetime
import sys
import argparse
import re #for pattern matchin
import csv
import pandas as pd
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument('date_time') # crab submit to inspect
args = parser.parse_args()

######################################
## Collects the failed jobs for a   ##
## crab submit and collects the     ##
## corresponding data files         ##
######################################


## Preparation: 
## 1. Go on Graphana and filter for failed jobs (select the Status 'failed' in the panel)
## 2. Open 'Jobs Table - last retry only' and click on the three dots in the upper right corner, select 'Inspect'
## 3. Download the table as CSV and save it in the 'failedJobLogs' folder


# Convert Date-Time into crab naming:

months = {"06": "Jun", "07": "Jul"} # add more if needed

year  = args.date_time[0:4]
month = months[args.date_time[4:6]]
day   = args.date_time[6:8]

crab_date = year+month+day

print(" ====> Investigating File of:", crab_date, "with id:", args.date_time)

# input csv and project directory
path        = "/work/pahwagne/CMSSW_10_6_37/src/PhysicsTools/RDsNano/production/failedJobCsv/"
project_dir = "/work/pahwagne/CMSSW_10_6_37/src/PhysicsTools/RDsNano/production/"+crab_date+"/crab_"+args.date_time

#output
out         = "/work/pahwagne/CMSSW_10_6_37/src/PhysicsTools/RDsNano/production/failedJobLogs/"+args.date_time+"/"
outMini     = "/work/pahwagne/CMSSW_10_6_37/src/PhysicsTools/RDsNano/production/failedJobMinis/"+args.date_time+"/"
os.makedirs(out)
os.makedirs(outMini)

# Search csv file
logs = os.listdir(path)

log = ''

for l in logs: 
  if args.date_time in l: 
    log = l

if log == '': print("No log csv file found!"); sys.exit();

# Dump Ids of failed jobs into list
jobDetails = pd.read_csv(path + log)
jobIDs     = list(jobDetails["data.CRAB_Id"])
jobIDs     = jobIDs[0:3] #debug
jobIDs     = [str(jobID) for jobID in jobIDs] #convert elements to str

jobIDsStr = ','.join(jobIDs) # comma separate 

# Now get the log files

command = "crab getlog " + " -d " + project_dir + " --jobids=" + jobIDsStr + " --short" + " --outputpath=" + out

print(command)
os.system(command)

#result = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# Search the miniAOD files within the log files:

  
allminis = open(outMini + "all_minis_job.txt","wt")

for name in jobIDs:

  temp  = open(out+"job_out."+name+".0.txt","rt")
  minis = open(outMini + "minis_job_"+name+".txt","wt")

  for line in temp:
    if "/store/data/Run2018D/ParkingBPH" in line and "eos" not in line and "cms" not in line:
      minis.write(line.strip().strip(",").strip("'") + "\n")    #remove whitespaces, etc.
      allminis.write(line.strip().strip(",").strip("'") + "\n") #remove whitespaces, etc.
 
  temp.close()
  minis.close()

allminis.close()



