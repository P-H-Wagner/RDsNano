import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RDsNano.common_cff import RDsCandVars, ufloat, uint, ubool
from PhysicsTools.RDsNano.rds_common_cff import TableDefaultVariables, TableDefault
from PhysicsTools.RDsNano.variables_cff import * 
#from PhysicsTools.RDsNano.primaryVertices_cff import *

#BsToDsPhiKKPiMuCfg = BuilderDefaultCfg.clone()
#BTo3MuCfg.dileptons             = cms.InputTag('JpsiMuonPairs')
#BTo3MuCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc

#cuts for ma thesis at:https://github.com/rmanzoni/Bmmm/blob/main/Bmmm/Analysis/test/rds/cuts.py

BsToDsPhiKKPiMu = cms.EDProducer(
    'BsToDsPhiKKPiMuBuilder',
    offBeamSpot = cms.InputTag('offlineBeamSpot'),
    pfCand = cms.InputTag('packedPFCandidates'),
    muCand = cms.InputTag('muonTrgSelector', 'trgMuons'),
    pvCand = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tracks = cms.InputTag('packedPFCandidates'), # for isolation
    lostTracks = cms.InputTag("lostTracks"),     # consider also lost pf cands!
    prunedCand = cms.InputTag("prunedGenParticles"),
    packedCand = cms.InputTag("packedGenParticles"),
    photonCand = cms.InputTag("reducedEgamma"),
    unpackedTracksCand = cms.InputTag("unpackedTracksAndVertices"),
    vtxWithBsCand = cms.InputTag("primaryVertexRefit","WithBS", "RDsNANO"),
    hadSelection = cms.string(' &&  '.join([
'abs(pdgId) != 11',
'abs(pdgId) != 13',
'charge != 0',
'pt > 1.0', 
'eta > -2.4',
'eta < 2.4',
'hasTrackDetails'])), # pre-selection of hadrons (k1,k2 and pion)

    gSelection = cms.string(' &&  '.join([
'(abs(pdgId) == 22)'])), #pre-selection of photons

    hadSelectionGen = cms.string(' &&  '.join([
'charge != 0',
'pt > 0.7', 
'eta > -2.1',
'eta < 2.1'])), # pre-selection of Gen hadrons (k1,k2 and pion), allow some tolerance w.r.t. hadSelection
maxdRHadMuon     = cms.double( 1.2),       # max dR between hadron and muon
mindRHadMuon     = cms.double( 0.005),     # min dR "
maxdRPhotonDs    = cms.double( 1.0),       # max dR between photon and Ds
maxdzDiffHadMuon = cms.double( 0.5),   # difference in dz between muon/pv and had/pv
maxdxyHadPv      = cms.double( 0.6),
phiMassAllowance = cms.double( 0.015),  # allow 15 MeV when collecting candidates for phi 
dsMassAllowance  = cms.double( 0.06),   # allow 60 MeV when collecting candidates for ds
dsStarMassAllowance  = cms.double( 0.1), # allow 100 MeV when collecting candidates for dsstar
drMatchGen       = cms.double( 0.1),         # allow 0.1, (0.05 would also be reasonable) in dR when gen matching
maxBsMass        = cms.double( 8.0 ),  
piMass           = cms.double( 0.13957039),      # pi mass
piMassSigma      = cms.double( 0.00000018),      # pi mass
kMass            = cms.double( 0.493677),         # kaon mass
kMassSigma       = cms.double( 0.000016),         # kaon mass
phiMass          = cms.double( 1.019461),       # phi mass
constrainPhiMass = cms.bool(True),    # constrain phi mass in the vtx fit?
minPhiVtxProb    = cms.double(0.01),
dsMass           = cms.double( 1.96834),         # ds mass
constrainDsMass  = cms.bool(True),     # constrain Ds mass in the vtx fit?
minDsVtxProb     = cms.double(0.01),
dsStarMass       = cms.double( 2.112204),    # ds star mass
muMass           = cms.double( 0.105658),        # mu ma
muMassSigma      = cms.double( 0.0000000023),        # mu mass
bsMass           = cms.double( 5.36688),         # bs mass
isoCone          = cms.double( 0.5)             # cut on dR for the mu isolation cone
)

print( " ========> Parameters used:")
print(BsToDsPhiKKPiMu.dumpPython)

#BsToDsPhiKKPiMuVariables.extend(vertexVariables)
combined_variables = cms.PSet(
  triggerVariables,
  prefitBasicVariables,
  #genVariables,
  vertexVariables,
  postfitBasicVariables,
  helicityVariables,
  bsMomentumVariables,
  photonVariables

)
  
BsToDsPhiKKPiMuTable = cms.EDProducer('SimpleCompositeCandidateFlatTableProducer',

    src = cms.InputTag("BsToDsPhiKKPiMu:bs"),
    cut = cms.string(""),
    name = cms.string(""),
    doc  = cms.string("BsToDsPhiKKPiMu"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = combined_variables

)



BsToDsPhiKKPiMuSequence = cms.Sequence(BsToDsPhiKKPiMu + BsToDsPhiKKPiMuTable)

