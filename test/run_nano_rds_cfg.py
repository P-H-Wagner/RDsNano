from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

# TODO: put different samples into parser (flag from command line)
# channel = 'sig'
channel = 'sig'

import os

#globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
globaltag = '106X_upgrade2018_realistic_v11_L1v1'

#what's the purpose of this 
annotation = '%s nevts:%d' % ('file:/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanotest.root', 100)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('RDsNANO',eras.Run2_2018)

# import of standard configurations (do we need all of them?)
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.RDsNano.nanoRDs_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


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
  return ['file:' + directory + filename for filename in filenames ]

def filesFromTxt(txtFile):
  with open(txtFile) as dataFiles: 
    filenames = [line for line in dataFiles]
  return ['file:' + 'root://cms-xrd-global.cern.ch//' + filename for filename in filenames ]

# Input source

if channel == 'sig':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/' #old signals MA Thesis
  #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/signals/all_signals_request_21_11_23.txt' # new signals!!
  inputfiles = filesFromFolder(directory)[0]
  #inputfiles = filesFromTxt(directory)
  #inputfiles = inputfiles[0]

if channel == 'hb':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1/' #hb 
  inputfiles = filesFromFolder(directory)

if channel == 'bplus':
  txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bplus/bplus.txt' #data
  inputfiles = filesFromTxt(txtFile)

if channel == 'data':
  txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_D/BPark_2018D.txt' #data
  txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/runTest/run_325117.txt' # only one run
  inputfiles = filesFromTxt(txtFile)

process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring('file:/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/all_signals_HbToDsPhiKKPiMuNu_MT_97.root'),
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/5D552231-6643-F042-9DD8-4CFC0CFD1B26.root'),
    #fileNames = cms.untracked.vstring("file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH3/MINIAOD/UL2018_MiniAODv2-v1/50000/0B26935C-81C7-1D4F-8994-23CBB40DA54C.root"),
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/50003/199766A2-70D0-674D-A9CE-D5EB255BD87A.root'), # data file which is always empty, investigate this 
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/40000/56D888ED-EB2C-B24F-A5A0-8D162DAFFA25.root'), #data file to compare with riccs MA code
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/50002/D0AE1369-0D7B-554C-BBB9-7B324AACCABD.root'), #data file to compare with riccs MA code
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/50002/36CD4F31-A249-DF49-A3FF-32DCA7223D09.root'), #10
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/40000/56D888ED-EB2C-B24F-A5A0-8D162DAFFA25.root'), #7
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2530000/521F0FA1-D91E-ED4D-8920-5BB1E7CDDE38.root'), #7
    fileNames = cms.untracked.vstring(inputfiles),# all_signals_HbToDsPhiKKPiMuNu_MT_0.root'), #automized case
    #fileNames = cms.untracked.vstring(inputfiles),
    #fileNames = cms.untracked.vstring('file:root://eoscms.cern.ch//eos/cms/store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/2530009/33EBB973-79DA-114A-9C72-1CC48E0ED7C6.root'), #fails with crab
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/CB5855BC-D585-7243-B46E-40130B3B16C9.root'), #lxy is nan here
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/50303941-50B8-C84A-9D02-D5000919ED8A.root'), #lxy is nan here
    #fileNames = cms.untracked.vstring(''), #lxy is nan here
    #fileNames = cms.untracked.vstring(['file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/CB5855BC-D585-7243-B46E-40130B3B16C9.root', 'file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/18444F62-3117-C244-8DCB-A065FB62C65D.root']),
    #fileNames = cms.untracked.vstring(['file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/CB5855BC-D585-7243-B46E-40130B3B16C9.root', 'file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDs_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/50303941-50B8-C84A-9D02-D5000919ED8A.root', 'file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDsMuNu_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v2/40000/18444F62-3117-C244-8DCB-A065FB62C65D.root']),


    secondaryFileNames = cms.untracked.vstring(),
    #eventsToProcess = cms.untracked.VEventRange('325117:316225000', '325117:317459950', '325117:316794666', '325117:316199657'),
    #eventsToProcess = cms.untracked.VEventRange('325117:316306200','325117:316362075','325117:317317084','325117:316229017','325117:317682842','325117:316762974','325117:317803805','325117:316457342','325117:317726844'), #constrained vs ma
    #eventsToProcess = cms.untracked.VEventRange('325117:316306200','325117:316362075','325117:316762974','325117:316164638','325117:317419163'), #unconstrained vs ma
    #eventsToProcess = cms.untracked.VEventRange('325117:316306200','325117:316362075','325117:316762974','325117:316164638','325117:317309145', '325117:317587931','325117:317002588'), #unconstrained vs new ma (same events :))
    #eventsToProcess = cms.untracked.VEventRange('1:3316932894'),
    #eventsToProcess = cms.untracked.VEventRange('1:3314083441'),

    #eventsToProcess  = cms.untracked.VEventRange('325117:316065411','325117:316065411', '325117:316066604', '325117:316067368', '325117:316070783', '325117:339068458', '325117:339074399', '325117:339076186', '325117:339078422', '325117:339078591'),
    skipEvents=cms.untracked.uint32(0) # skip first n events   

)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

#purpose?
process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
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
    fileName = cms.untracked.string('file:/work/pahwagne/test/nanotest.root'), #used for local tests
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
process = nanoAOD_customizeMuonTriggerBPark(process)  
process = nanoAOD_customizeBsToDsPhiKKPiMu(process) #comment this out to run only Trigger.cc for debugging

if channel != 'data':
  #can only gen match on mc
  process = nanoAOD_customizeGenMatching(process) 
  # Path and EndPath definitions
  process.nanoAOD_Bs_step= cms.Path(process.triggerSequence  + process.nanoBsToDsPhiKKPiMuSequence + process.nanoGenMatchingSequence)

else:
  process.nanoAOD_Bs_step= cms.Path(process.triggerSequence  + process.nanoBsToDsPhiKKPiMuSequence )
 


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
