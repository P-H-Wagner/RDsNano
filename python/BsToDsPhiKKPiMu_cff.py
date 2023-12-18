import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RDsNano.common_cff import RDsCandVars, ufloat, uint, ubool
from PhysicsTools.RDsNano.rds_common_cff import TableDefaultVariables, TableDefault
from PhysicsTools.RDsNano.variables_cff import BsToDsPhiKKPiMuVariables
#from PhysicsTools.RDsNano.primaryVertices_cff import *

#BsToDsPhiKKPiMuCfg = BuilderDefaultCfg.clone()
#BTo3MuCfg.dileptons             = cms.InputTag('JpsiMuonPairs')
#BTo3MuCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc

BsToDsPhiKKPiMu = cms.EDProducer(
    'BsToDsPhiKKPiMuBuilder',
    pfCand = cms.InputTag('packedPFCandidates'),
    muCand = cms.InputTag('muonTrgSelector', 'trgMuons'),
    pvCand = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

#BsToDsPhiKKPiMuTableVariables = TableDefaultVariables.clone()

"""
BsToDsPhiKKPiMuTable = TableDefault.clone()
BsToDsPhiKKPiMuTable.src       = cms.InputTag("KKPair")
BsToDsPhiKKPiMuTable.cut       = cms.string("")
BsToDsPhiKKPiMuTable.name      = cms.string("KKPair")
BsToDsPhiKKPiMuTable.doc       = cms.string("KKPair Variables")
BsToDsPhiKKPiMuTable.singleton = cms.bool(False)
BsToDsPhiKKPiMuTable.extension = cms.bool(False)
BsToDsPhiKKPiMuTable.variables = BsToDsPhiKKPiMuTableVariables
"""


BsToDsPhiKKPiMuTable = cms.EDProducer('SimpleCompositeCandidateFlatTableProducer',

    src = cms.InputTag("BsToDsPhiKKPiMu:bs"),
    cut = cms.string(""),
    name = cms.string(""),
    doc  = cms.string("BsToDsPhiKKPiMu"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(BsToDsPhiKKPiMuVariables)   #BsToDsPhiKKPiMuTableVariables),

)


BsToDsPhiKKPiMuSequence = cms.Sequence(BsToDsPhiKKPiMu + BsToDsPhiKKPiMuTable)

#?? why needed
CountBsToDsPhiKKPiMu = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BsToDsPhiKKPiMu")
)
