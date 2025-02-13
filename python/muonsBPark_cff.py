import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.RDsNano.variables_cff import triggerVariables

muonTrgSelector = cms.EDProducer("Trigger",

  muonCollection        = cms.InputTag("slimmedMuons"), 
  trgResultsCollection  = cms.InputTag("TriggerResults", "", "HLT"),
  trgObjectsCollection  = cms.InputTag("slimmedPatTrigger"),
  trgPrescaleCollection = cms.InputTag("patTrigger"),
  vtxCollection         = cms.InputTag("offlineSlimmedPrimaryVertices"),

  #cuts and selections
  trgFilterLabelMu7 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q"),
  trgFilterLabelMu9 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q"),
  muSelection    = cms.string(' &&  '.join([
  'pt > 7.0', 
  'eta > -1.5',
  'eta < 1.5',
  'isPFMuon',
  #'isGlobalMuon'
   ])), 
  
  hlt_7_4_p0   = cms.string("HLT_Mu7_IP4_part0_v2"), # trigger menue 
  hlt_7_4_p1   = cms.string("HLT_Mu7_IP4_part1_v2"), # "
  hlt_7_4_p2   = cms.string("HLT_Mu7_IP4_part2_v2"), # "
  hlt_7_4_p3   = cms.string("HLT_Mu7_IP4_part3_v2"), # "
  hlt_7_4_p4   = cms.string("HLT_Mu7_IP4_part4_v2"), # "
  
  hlt_9_6_p0   = cms.string("HLT_Mu9_IP6_part0_v3"), # trigger menue
  hlt_9_6_p1   = cms.string("HLT_Mu9_IP6_part1_v3"), # "
  hlt_9_6_p2   = cms.string("HLT_Mu9_IP6_part2_v3"), # "
  hlt_9_6_p3   = cms.string("HLT_Mu9_IP6_part3_v3"), # "
  hlt_9_6_p4   = cms.string("HLT_Mu9_IP6_part4_v3"), # "
  
  maxdR_matching = cms.double(0.05), #muon trg object matching                             
)

muonTrgTable = cms.EDProducer("SimpleCandidateFlatTableProducer",

    src  = cms.InputTag("muonTrgSelector:trgMuons"),
    cut  = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("trg"), #must be different than the name of BsToDs... Builder!!
    doc  = cms.string("muonTrgSelector"),
    singleton = cms.bool(False), 
    extension = cms.bool(False), 
    variables = triggerVariables

)

print( " ========> Parameters used:")
print(muonTrgSelector.dumpPython)

muonTrgSequence = cms.Sequence(muonTrgSelector + muonTrgTable)
