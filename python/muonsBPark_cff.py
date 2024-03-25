import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

muonTrgSelector = cms.EDProducer("Trigger",
                                 muonCollection        = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                 trgResultsCollection  = cms.InputTag("TriggerResults", "", "HLT"),
                                 trgObjectsCollection  = cms.InputTag("slimmedPatTrigger"),
                                 trgPrescaleCollection = cms.InputTag("patTrigger"),
                                 vtxCollection         = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 #cuts and selections
                                 maxdR_matching        = cms.double(0.02),                             
)

#countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
#                            minNumber = cms.uint32(0),
#                            maxNumber = cms.uint32(999999),
#                            src = cms.InputTag("muonTrgSelector", "trgMuons")
#                         )


muonTrgTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("muonTrgSelector:trgMuons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("muon"), #must be different than the name of BsToDs... Builder!!
    doc  = cms.string("slimmedMuons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
        isPFcand  = Var("isPFMuon",      bool, doc= "muon is PF candidate" ),
        isGlobal  = Var("isGlobalMuon",  bool, doc= "muon is global muon"  ),
        isTracker = Var("isTrackerMuon", bool, doc= "muon is tracker muon" ),
    ),
)

muonTrgSequence = cms.Sequence(muonTrgSelector + muonTrgTable)
