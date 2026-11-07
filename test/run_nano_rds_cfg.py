from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

channel = 'sig'
refit = True #False 


import os

#globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
globaltag = '106X_upgrade2018_realistic_v11_L1v1'

#what's the purpose of this 
#annotation = '%s nevts:%d' % ('file:/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanotest.root', 100)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('RDsNANO',eras.Run2_2018)

# import of standard configurations (do we need all of them?)
process.load('Configuration.StandardSequences.Services_cff') #
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')      
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')#
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.RDsNano.nanoRDs_cff') #
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') #

#create the collection
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

#load vtx 
process.load("RecoVertex.Configuration.RecoVertex_cff")



#prints the time report
process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True),
                             useJobReport = cms.untracked.bool(True)
)

#load all the chosen options
process.MessageLogger.cerr.FwkReport.reportEvery =1
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(300)
)

def filesFromFolder(directory):
  filenames = os.listdir(directory)
  filenames = [file for file in filenames if "root" in file]
  return ['file:' + directory + filename for filename in filenames ]

def filesFromTxt(txtFile):
  with open(txtFile) as dataFiles: 
    filenames = [line for line in dataFiles]
  return ['file:' + 'root://cms-xrd-global.cern.ch//' + filename for filename in filenames ]

# Input source

if channel == 'sig':
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/' #old signals MA Thesis
  #if refit: directory = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/signals/refit/"
  #if refit: directory = "/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1_PV_REFITTED/"
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/signals/all_signals_request_21_11_23.txt' # new signals!!
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/signals/all_signals_request_18_07_25.txt' # new signals!!
  #inputfiles = filesFromFolder(directory)
  inputfiles = filesFromTxt(directory)
  #inputfiles = inputfiles[0]

if channel == 'hb':

  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1/' #hb 
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/hbInclusiveToDsPhiKKPi_20250829/hbInclusiveToDsPhiKKPi/hbInclusiveToDsPhiKKPi_MINI/250829_074208/0000/' #new hb
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/hb_request_18_07_25.txt' # new hb!!
  inputfiles = filesFromTxt(directory)

if channel == 'bplus':
  txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bplus/bplus.txt' #data
  inputfiles = filesFromTxt(txtFile)

if channel == 'data':
  #txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_D/BPark_2018D.txt' #data
  #txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/runTest/run_325117.txt' # only one run
  #txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_B/BPH1_B.txt' #data bParking 2018 part D
  txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPH1/BPH1_partA.txt' #data bParking 2018 part D
  inputfiles = filesFromTxt(txtFile)

if channel == 'qcd':
  txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/qcd/pt_50_to_80.txt' 
  inputfiles = filesFromTxt(txtFile)


process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring('file:/pnfs/psi.ch/cms/trivcat/store/user/manzoni/inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1/inclusive_HbToDsPhiKKPiMuNu_0.root'), # run=321436 lumi=592 to check ds vtx cos
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/50003/199766A2-70D0-674D-A9CE-D5EB255BD87A.root'), # data file which is always empty, investigate this 
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch////store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2530009/665CB12D-9B84-144F-9811-CB4953DA52CF.root'), # contains mu9 only
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/5D552231-6643-F042-9DD8-4CFC0CFD1B26.root'), #data file to compare with riccs MA code
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/91AA7B5A-A0CA-B845-825B-7428EC79F833.root'), #data file to compare with riccs MA code
    #fileNames = cms.untracked.vstring(inputfiles),# all_signals_HbToDsPhiKKPiMuNu_MT_0.root'), #automized case
    #fileNames = cms.untracked.vstring('file:root://eoscms.cern.ch//eos/cms/store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2530009/33EBB973-79DA-114A-9C72-1CC48E0ED7C6.root'), #fails with crab
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/HbToDsMu_DsToPhiKKPi_TuneCP5_13TeV_pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/90000/9EA8B1BE-5BF7-9D4D-BCDF-D432326847DB.root'), #fails with crab
    fileNames = cms.untracked.vstring(inputfiles), #fails with crab
    #eventsToProcess = cms.untracked.VEventRange('1:8708', "1:17743", "1:16520", "1:30437", "1:26990", "1:58035"),
    #eventsToProcess = cms.untracked.VEventRange('1:17883749', '1:17909804', '1:17946276', '1:17958917','1:17964025','1:17986950','1:18009262','1:18039706', '1:18060883','1:28672802'),#, "1:17743", "1:16520", "1:30437", "1:26990", "1:58035"),
    skipEvents=cms.untracked.uint32(25000) # skip first n events   

)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    TryToContinue = cms.untracked.vstring('ProductNotFound')
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    #annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    #fileName = cms.untracked.string('file:/scratch/pahwagne/nanoAOD/test.root' ),
    fileName = cms.untracked.string('file:nanotest_sig.root'), #used for local tests
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')


# add all sequences as addtributes
from PhysicsTools.RDsNano.nanoRDs_cff import *
process = nanoAOD_customizeStart(process)
process = nanoAOD_customizeMuonTriggerBPark(process)  
process = nanoAOD_customizeBsToDsPhiKKPiMu(process) #comment this out to run only Trigger.cc for debugging


## put the refit process somewhere in the beginning before all the EDAnalyzer
process.primaryVertexRefit = process.unsortedOfflinePrimaryVertices.clone()
process.primaryVertexRefit.TrackLabel = cms.InputTag("unpackedTracksAndVertices")
# Only the beam-spot-constrained collection is read (BsToDsPhiKKPiMu.vtxWithBsCand = primaryVertexRefit:WithBS).
# Do not fit the unconstrained "" collection as well: it is computed for every event and never used.
process.primaryVertexRefit.vertexCollections = cms.VPSet(
    [pset for pset in process.primaryVertexRefit.vertexCollections if pset.label.value() == 'WithBS']
)
#process.primaryVertexRefitSequence = cms.Sequence(process.primaryVertexRefit)

if channel != 'data':
  #can only gen match on mc
  process = nanoAOD_customizeGenMatching(process) 
  # Path and EndPath definitions
  print(" ----------- Processing MC files -------------- ")
  process.nanoAOD_Bs_step= cms.Path( process.nanoSequence + process.unpackedTracksAndVertices +  process.primaryVertexRefit + process.triggerSequence  + process.nanoBsToDsPhiKKPiMuSequence + process.nanoGenMatchingSequence)

else:
  print(" ----------- Processing data files -------------- ")
  process.nanoAOD_Bs_step= cms.Path( process.nanoSequence + process.unpackedTracksAndVertices +  process.primaryVertexRefit + process.triggerSequence  + process.nanoBsToDsPhiKKPiMuSequence )
 


#process.nanoAOD_Bs_step= cms.Path(process.triggerSequence) ## to run only Trigger.cc for debugging
process.endjob_step = cms.EndPath(process.endOfProcess)

process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.nanoAOD_Bs_step,
    process.endjob_step, 
    process.NANOAODoutput_step # commented out ? not saving !!
                               )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('nanoAOD_Bs_step')
)

## multicore, I hope this does not screw everything up
## must be consistent with the cpu number of the batch submission!
#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)


# ?? 
### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
