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
  trgFilterLabelMu7_4 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q"),
  trgFilterLabelMu8_3 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q"),
  trgFilterLabelMu8_5 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP5Q"),
  trgFilterLabelMu8_6 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP6Q"),
  trgFilterLabelMu9_4 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP4Q"),
  trgFilterLabelMu9_5 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q"),
  trgFilterLabelMu9_6 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q"),
  trgFilterLabelMu12_6= cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q"),
  muSelection    = cms.string(' &&  '.join([
  'pt > 7.0', 
  'eta > -1.5',
  'eta < 1.5',
  'isPFMuon',
  #'isGlobalMuon'
   ])), 

  #versions check with https://raw.githubusercontent.com/cms-sw/cmssw/6cb4ba7dea7cedb955469404e84e73063f69aa8f/HLTrigger/Configuration/python/HLT_GRun_cff.py
  
  hlt_7_4_p0   = cms.string("HLT_Mu7_IP4_part0_v2"), # trigger menue 
  hlt_7_4_p1   = cms.string("HLT_Mu7_IP4_part1_v2"), # "
  hlt_7_4_p2   = cms.string("HLT_Mu7_IP4_part2_v2"), # "
  hlt_7_4_p3   = cms.string("HLT_Mu7_IP4_part3_v2"), # "
  hlt_7_4_p4   = cms.string("HLT_Mu7_IP4_part4_v2"), # "

  #take v1 for periods A,B,C
  hlt_8_3_p0   = cms.string("HLT_Mu8_IP3_part0_v1"), # trigger menue 
  hlt_8_3_p1   = cms.string("HLT_Mu8_IP3_part1_v1"), # "
  hlt_8_3_p2   = cms.string("HLT_Mu8_IP3_part2_v1"), # "
  hlt_8_3_p3   = cms.string("HLT_Mu8_IP3_part3_v1"), # "
  hlt_8_3_p4   = cms.string("HLT_Mu8_IP3_part4_v1"), # "

  #hlt_8_3_p0   = cms.string("HLT_Mu8_IP3_part0_v3"), # trigger menue 
  #hlt_8_3_p1   = cms.string("HLT_Mu8_IP3_part1_v3"), # "
  #hlt_8_3_p2   = cms.string("HLT_Mu8_IP3_part2_v3"), # "
  #hlt_8_3_p3   = cms.string("HLT_Mu8_IP3_part3_v3"), # "
  #hlt_8_3_p4   = cms.string("HLT_Mu8_IP3_part4_v3"), # "

  hlt_8_5_p0   = cms.string("HLT_Mu8_IP5_part0_v2"), # trigger menue 
  hlt_8_5_p1   = cms.string("HLT_Mu8_IP5_part1_v2"), # "
  hlt_8_5_p2   = cms.string("HLT_Mu8_IP5_part2_v2"), # "
  hlt_8_5_p3   = cms.string("HLT_Mu8_IP5_part3_v2"), # "
  hlt_8_5_p4   = cms.string("HLT_Mu8_IP5_part4_v2"), # "

  hlt_8_6_p0   = cms.string("HLT_Mu8_IP6_part0_v2"), # trigger menue 
  hlt_8_6_p1   = cms.string("HLT_Mu8_IP6_part1_v2"), # "
  hlt_8_6_p2   = cms.string("HLT_Mu8_IP6_part2_v2"), # "
  hlt_8_6_p3   = cms.string("HLT_Mu8_IP6_part3_v2"), # "
  hlt_8_6_p4   = cms.string("HLT_Mu8_IP6_part4_v2"), # "

  hlt_9_4_p0   = cms.string("HLT_Mu9_IP4_part0_v2"), # trigger menue
  hlt_9_4_p1   = cms.string("HLT_Mu9_IP4_part1_v2"), # "
  hlt_9_4_p2   = cms.string("HLT_Mu9_IP4_part2_v2"), # "
  hlt_9_4_p3   = cms.string("HLT_Mu9_IP4_part3_v2"), # "
  hlt_9_4_p4   = cms.string("HLT_Mu9_IP4_part4_v2"), # "

  hlt_9_5_p0   = cms.string("HLT_Mu9_IP5_part0_v2"), # trigger menue
  hlt_9_5_p1   = cms.string("HLT_Mu9_IP5_part1_v2"), # "
  hlt_9_5_p2   = cms.string("HLT_Mu9_IP5_part2_v2"), # "
  hlt_9_5_p3   = cms.string("HLT_Mu9_IP5_part3_v2"), # "
  hlt_9_5_p4   = cms.string("HLT_Mu9_IP5_part4_v2"), # "
 
  #take v1 for periods A,B,C
  hlt_9_6_p0   = cms.string("HLT_Mu9_IP6_part0_v1"), # trigger menue
  hlt_9_6_p1   = cms.string("HLT_Mu9_IP6_part1_v1"), # "
  hlt_9_6_p2   = cms.string("HLT_Mu9_IP6_part2_v1"), # "
  hlt_9_6_p3   = cms.string("HLT_Mu9_IP6_part3_v1"), # "
  hlt_9_6_p4   = cms.string("HLT_Mu9_IP6_part4_v1"), # "

  #hlt_9_6_p0   = cms.string("HLT_Mu9_IP6_part0_v3"), # trigger menue
  #hlt_9_6_p1   = cms.string("HLT_Mu9_IP6_part1_v3"), # "
  #hlt_9_6_p2   = cms.string("HLT_Mu9_IP6_part2_v3"), # "
  #hlt_9_6_p3   = cms.string("HLT_Mu9_IP6_part3_v3"), # "
  #hlt_9_6_p4   = cms.string("HLT_Mu9_IP6_part4_v3"), # "

  hlt_12_6_p0   = cms.string("HLT_Mu12_IP6_part0_v2"), # trigger menue
  hlt_12_6_p1   = cms.string("HLT_Mu12_IP6_part1_v2"), # "
  hlt_12_6_p2   = cms.string("HLT_Mu12_IP6_part2_v2"), # "
  hlt_12_6_p3   = cms.string("HLT_Mu12_IP6_part3_v2"), # "
  hlt_12_6_p4   = cms.string("HLT_Mu12_IP6_part4_v2"), # "

  
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

muonTrgSequence = cms.Sequence(muonTrgSelector)
