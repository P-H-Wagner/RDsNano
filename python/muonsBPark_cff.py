import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

muonTrgSelector = cms.EDProducer("Trigger",

muonCollection        = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
trgResultsCollection  = cms.InputTag("TriggerResults", "", "HLT"),
trgObjectsCollection  = cms.InputTag("slimmedPatTrigger"),
trgPrescaleCollection = cms.InputTag("patTrigger"),
vtxCollection         = cms.InputTag("offlineSlimmedPrimaryVertices"),
#cuts and selections
trgFilterLabel = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q"),
muSelection    = cms.string(' &&  '.join([
'pt > 7.0', 
'eta > -1.5',
'eta < 1.5',
'isPFMuon',
'isGlobalMuon'])), # pre-selection of hadrons (k1,k2 and pion)
hlt_7_4_p0   = cms.string("HLT_Mu7_IP4_part0_v2"), # trigger menue
hlt_7_4_p1   = cms.string("HLT_Mu7_IP4_part1_v2"), # "
hlt_7_4_p2   = cms.string("HLT_Mu7_IP4_part2_v2"), # "
hlt_7_4_p3   = cms.string("HLT_Mu7_IP4_part3_v2"), # "
hlt_7_4_p4   = cms.string("HLT_Mu7_IP4_part4_v2"), # "
maxdR_matching = cms.double(0.05), #muon trg object matching                             
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

print( " ========> Parameters used:")
print(muonTrgSelector.dumpPython)

muonTrgSequence = cms.Sequence(muonTrgSelector)
