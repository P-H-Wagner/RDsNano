import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RDsNano.common_cff import RDsCandVars, ufloat, uint, ubool
from PhysicsTools.RDsNano.rds_common_cff import TableDefaultVariables, TableDefault
from PhysicsTools.RDsNano.variables_cff import genVariables

#from PhysicsTools.RDsNano.primaryVertices_cff import *

#BsToDsPhiKKPiMuCfg = BuilderDefaultCfg.clone()
#BTo3MuCfg.dileptons             = cms.InputTag('JpsiMuonPairs')
#BTo3MuCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc

genMatcher = cms.EDProducer(
    'genMatching',
    bs           = cms.InputTag('BsToDsPhiKKPiMu', 'bs'),
    prunedCand   = cms.InputTag("prunedGenParticles"),
    packedCand   = cms.InputTag("packedGenParticles"),
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
'eta < 2.4'])), # pre-selection of Gen hadrons (k1,k2 and pion), allow some tolerance w.r.t. hadSelection
maxdRHadMuon = cms.double(1.2),       # max dR between hadron and muon
mindRHadMuon = cms.double(0.005),     # min dR "
maxdzDiffHadMuon = cms.double(0.5),   # difference in dz between muon/pv and had/pv
phiMassAllowance = cms.double(0.015), # allow 15 MeV when collecting candidates for phi 
dsMassAllowance = cms.double(0.05),   # allow 50 MeV when collecting candidates for ds
drMatchGen = cms.double(0.1),         # allow 0.1, (0.05 would also be reasonable) in dR when gen matching

piMass = cms.double(0.13957039),      # pi mass
kMass = cms.double(0.493677),         # kaon mass
phiMass = cms.double(1.019461),       # phi mass
dsMass = cms.double(1.96834),         # ds mass
dsStarMass = cms.double(2.112204),    # ds star mass
muMass = cms.double(0.105658),        # mu mass
bsMass = cms.double(5.36688),         # bs mass
isoCone = cms.double(0.5)             # cut on dR for the mu isolation cone
)

#BsToDsPhiKKPiMuVariables.extend(vertexVariables)
combined_variables = cms.PSet(
  genVariables,
)
  
genMatcherTable= cms.EDProducer('SimpleCompositeCandidateFlatTableProducer',

    src = cms.InputTag("genMatcher:gen"),
    cut = cms.string(""),
    name = cms.string("gen"),
    doc  = cms.string("genMatcher"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), 
    variables = combined_variables

)

#tables = cms.Sequence(BsToDsPhiKKPiMuTable) # * vertexTable)
genMatchingSequence = cms.Sequence(genMatcher + genMatcherTable)

