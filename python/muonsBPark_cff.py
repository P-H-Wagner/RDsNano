import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

muonTrgSelector = cms.EDProducer("Trigger",
                                 muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                 trgResultsCollection = cms.InputTag("TriggerResults", "", "HLT"),
                                 trgObjectsCollection = cms.InputTag("slimmedPatTrigger"),
                                 trgPrescaleCollection = cms.InputTag("patTrigger"),
                                 vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 #cuts and selections
                                 maxdR_matching = cms.double(0.02),                             
)

#countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
#                            minNumber = cms.uint32(0),
#                            maxNumber = cms.uint32(999999),
#                            src = cms.InputTag("muonTrgSelector", "trgMuons")
#                         )


muonBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("muonTrgSelector:trgMuons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
        ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        xErr = Var("bestTrack().covariance(0,0)", float, doc = "xError of the muon track", precision=10),
        yErr = Var("bestTrack().covariance(1,1)", float, doc = "xError of the muon track", precision=10),
        zErr = Var("bestTrack().covariance(2,2)", float, doc = "xError of the muon track", precision=10),
        xyErr = Var("bestTrack().covariance(0,1)", float, doc = "xError of the muon track", precision=10),
        yzErr = Var("bestTrack().covariance(0,2)", float, doc = "xError of the muon track", precision=10),
        xzErr = Var("bestTrack().covariance(1,2)", float, doc = "xError of the muon track", precision=10),
        isPFcand = Var("isPFMuon",bool,doc="muon is PF candidate"),
        isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon"),
        isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon"),
    ),
)

muonBParkSequence = cms.Sequence(muonTrgSelector)
muonBParkTables = cms.Sequence(muonBParkTable)
