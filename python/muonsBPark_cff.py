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
  trgFilterLabelMu8p5_3p5 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8p5Q"),
  trgFilterLabelMu8_5 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP5Q"),
  trgFilterLabelMu8_6 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP6Q"),
  trgFilterLabelMu9_4 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP4Q"),
  trgFilterLabelMu9_5 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q"),
  trgFilterLabelMu9_6 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q"),
  trgFilterLabelMu10p5_3p5 = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered10p5Q"),
  trgFilterLabelMu12_6= cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q"),
  muSelection    = cms.string(' &&  '.join([
  'pt > 7.0', 
  'eta > -1.5',
  'eta < 1.5',
  'isPFMuon',
  #'isGlobalMuon'
   ])), 

  #versions check with https://raw.githubusercontent.com/cms-sw/cmssw/6cb4ba7dea7cedb955469404e84e73063f69aa8f/HLTrigger/Configuration/python/HLT_GRun_cff.py
  
  hlt_7_4   = cms.string("HLT_Mu7_IP4"), # trigger menue 
  hlt_8_3   = cms.string("HLT_Mu8_IP3"), # trigger menue 
  hlt_8p5_3p5   = cms.string("HLT_Mu8p5_IP3p5"), # trigger menue 
  hlt_8_5   = cms.string("HLT_Mu8_IP5"), # trigger menue 
  hlt_8_6   = cms.string("HLT_Mu8_IP6"), # trigger menue 
  hlt_9_4   = cms.string("HLT_Mu9_IP4"), # trigger menue
  hlt_9_5   = cms.string("HLT_Mu9_IP5"), # trigger menue
  hlt_9_6   = cms.string("HLT_Mu9_IP6"), # trigger menue
  hlt_10p5_3p5   = cms.string("HLT_Mu10p5_IP3p5"), # trigger menue 
  hlt_12_6  = cms.string("HLT_Mu12_IP6"), # trigger menue
  
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
