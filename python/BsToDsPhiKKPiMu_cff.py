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
    prunedCand = cms.InputTag("prunedGenParticles"),
    packedCand = cms.InputTag("packedGenParticles"),
    hadSelection = cms.string(' &&  '.join([
'pdgId != 11',
'pdgId != 13',
'charge != 0',
'pt > 1.0', 
'eta > -2.4',
'eta < 2.4',
'hasTrackDetails'])), # pre-selection of hadrons (k1,k2 and pion)
    hadSelectionGen = cms.string(' &&  '.join([
'charge != 0',
'pt > 1.0', 
'eta > -2.4',
'eta < 2.4'])), # pre-selection of hadrons (k1,k2 and pion)
maxdRHadMuon = cms.double(1.2),       # max dR between hadron and muon
mindRHadMuon = cms.double(0.005),     # min dR "
maxdzDiffHadMuon = cms.double(0.5),   # difference in dz between muon/pv and had/pv
phiMassAllowance = cms.double(0.015), # allow 15 MeV when collecting candidates for phi 
dsMassAllowance = cms.double(0.05),   # allow 50 MeV when collecting candidates for ds
drMatchGen = cms.double(0.5),        # allow 0.05 in dR when gen matching

piMass = cms.double(0.13957039),      # pi mass
kMass = cms.double(0.493677),         # kaon mass
phiMass = cms.double(1.019461),       # phi mass
dsMass = cms.double(1.96834),         # ds mass
muMass = cms.double(0.105658),        # mu mass
bsMass = cms.double(5.36688)          # bs mass

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
