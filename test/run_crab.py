# RDsNano config, used both locally and on CRAB (rds16, CMSSW_16_X).
#   cmsRun run_crab.py isMC=False globalTag=106X_dataRun2_v35 inputFiles=file:a.root maxEvents=500 outputFile=rds_nano.root
# On CRAB: pyCfgParams = ['isMC=...', 'globalTag=...'] (CRAB sets the input files itself).
from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing()   # no preset: plain outputFile (no _numEventN tag), our own option set
options.register('isMC',       False,           VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'MC (adds the gen matching) or data')
options.register('globalTag',  '',              VarParsing.multiplicity.singleton, VarParsing.varType.string, 'global tag (required; take it from the provenance of the sample)')
options.register('maxEvents',  -1,              VarParsing.multiplicity.singleton, VarParsing.varType.int,    'events to process (-1 = all)')
options.register('inputFiles', '',              VarParsing.multiplicity.list,      VarParsing.varType.string, 'input files')
options.register('outputFile', 'rds_nano.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'output NanoAOD file')
options.register('reportEvery', 5000,           VarParsing.multiplicity.singleton, VarParsing.varType.int,    'FwkReport frequency')
options.parseArguments()

if not options.globalTag:
  raise RuntimeError('run_crab.py: globalTag=... is required (from edmProvDump of the input sample, PAT step)')

channel = 'mc' if options.isMC else 'data'
globaltag = options.globalTag

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
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(0)    
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    TryToContinue = cms.untracked.vstring('ProductNotFound')
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
    fileName = cms.untracked.string(options.outputFile),
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


if channel != 'data':
  #can only gen match on mc
  process = nanoAOD_customizeGenMatching(process) 
  # Path and EndPath definitions

  process.nanoAOD_Bs_step= cms.Path(  process.nanoSequence
                                    + process.unpackedTracksAndVertices
                                    + process.primaryVertexRefit
                                    + process.triggerSequence
                                    + process.nanoBsToDsPhiKKPiMuSequence
                                    + process.nanoGenMatchingSequence)

else:

  process.nanoAOD_Bs_step= cms.Path(  process.nanoSequence
                                    + process.unpackedTracksAndVertices
                                    + process.primaryVertexRefit
                                    + process.triggerSequence
                                    + process.nanoBsToDsPhiKKPiMuSequence) 


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
