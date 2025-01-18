import FWCore.ParameterSet.Config as cms
print('========> Processing the miniAOD from: HOOK_INPUT')


process = cms.Process("VertexRefit")

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.MessageLogger.cerr.FwkReport.reportEvery =1
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "auto:run2_mc"
process.GlobalTag.globaltag = "106X_upgrade2018_realistic_v11_L1v1"

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      HOOK_FILE_IN
    )
)

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(HOOK_N_EVENTS))

# create a collection of tracks 
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')



# Load vertex reconstruction
process.load("RecoVertex.Configuration.RecoVertex_cff")

process.primaryVertexRefit = process.unsortedOfflinePrimaryVertices.clone()
process.primaryVertexRefit.TrackLabel = cms.InputTag("unpackedTracksAndVertices")
# process.primaryVertexRefit.TrackLabel = cms.InputTag("generalTracks")
# process.primaryVertexRefit.TrackLabel = cms.InputTag("packedPFCandidates")





process.path = cms.Path(process.unpackedTracksAndVertices + process.primaryVertexRefit)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:HOOK_FILE_OUT"),
    outputCommands = cms.untracked.vstring(
        "keep *",
    )
)

process.end = cms.EndPath(process.output)








# 
# 
# cms.EDProducer("RecoChargedRefCandidatePrimaryVertexSorter",
#     assignment = cms.PSet(
#         maxDistanceToJetAxis = cms.double(0.07),
#         maxDtSigForPrimaryAssignment = cms.double(4.0),
#         maxDxyForJetAxisAssigment = cms.double(0.1),
#         maxDxyForNotReconstructedPrimary = cms.double(0.01),
#         maxDxySigForNotReconstructedPrimary = cms.double(2),
#         maxDzErrorForPrimaryAssignment = cms.double(0.05),
#         maxDzForJetAxisAssigment = cms.double(0.1),
#         maxDzForPrimaryAssignment = cms.double(0.1),
#         maxDzSigForPrimaryAssignment = cms.double(5.0),
#         maxJetDeltaR = cms.double(0.5),
#         minJetPt = cms.double(25),
#         preferHighRanked = cms.bool(False),
#         useTiming = cms.bool(False)
#     ),
#     jets = cms.InputTag("ak4CaloJetsForTrk"),
#     particles = cms.InputTag("trackRefsForJetsBeforeSorting"),
#     produceAssociationToOriginalVertices = cms.bool(False),
#     produceNoPileUpCollection = cms.bool(False),
#     producePileUpCollection = cms.bool(False),
#     produceSortedVertices = cms.bool(True),
#     qualityForPrimary = cms.int32(3),
#     sorting = cms.PSet(
# 
#     ),
#     trackTimeResoTag = cms.InputTag(""),
#     trackTimeTag = cms.InputTag(""),
#     usePVMET = cms.bool(True),
#     vertices = cms.InputTag("unsortedOfflinePrimaryVertices")
# )
