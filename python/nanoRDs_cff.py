import FWCore.ParameterSet.Config as cms

# port (CMSSW_16_X): import only what is used from PhysicsTools.NanoAOD instead of nano_cff wholesale.
# In 16_X globals_cff has no 'globalTables' Sequence any more (globalTablesTask = beamSpotTable + rhoTable),
# and its rhoTable is named "Rho" (branches Rho_*) with an extra variable (fixedGridRhoAll).
# Keep the 10_6 content: globalTables == rhoTable with the 10_6 variables and unprefixed branch names.
from PhysicsTools.NanoAOD.common_cff import ExtVar
from PhysicsTools.NanoAOD.nano_cff import nanoMetadata
from PhysicsTools.NanoAOD.globals_cff import rhoTable as _rhoTable16
# MC: event-level pileup (Pileup_nTrueInt, Pileup_nPU, ...) for every event, including events without candidates
from PhysicsTools.NanoAOD.globals_cff import puTable

rhoTable = _rhoTable16.clone(
    name = "",
    variables = cms.PSet(
        fixedGridRhoFastjetAll = ExtVar( cms.InputTag("fixedGridRhoFastjetAll"), "double", doc = "rho from all PF Candidates, used e.g. for JECs" ),
        fixedGridRhoFastjetCentralNeutral = ExtVar( cms.InputTag("fixedGridRhoFastjetCentralNeutral"), "double", doc = "rho from neutral PF Candidates with |eta| < 2.5, used e.g. for rho corrections of some lepton isolations" ),
        fixedGridRhoFastjetCentralCalo = ExtVar( cms.InputTag("fixedGridRhoFastjetCentralCalo"), "double", doc = "rho from calo towers with |eta| < 2.5, used e.g. egamma PFCluster isolation" ),
        fixedGridRhoFastjetCentral = ExtVar( cms.InputTag("fixedGridRhoFastjetCentral"), "double", doc = "rho from all PF Candidates for central region, used e.g. for JECs" ),
        fixedGridRhoFastjetCentralChargedPileUp = ExtVar( cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"), "double", doc = "rho from charged PF Candidates for central region, used e.g. for JECs" ),
    )
)
globalTables = cms.Sequence(rhoTable)

# for trigger muon
from PhysicsTools.RDsNano.muonsBPark_cff import * 

# Bs chain collection
from PhysicsTools.RDsNano.BsToDsPhiKKPiMu_cff import *

# L1 seeds of the parking paths: decision + brilcalc prescale (plugins/L1SeedsTableProducer.cc, as Bmmm)
l1SeedsTable = cms.EDProducer("L1SeedsTableProducer",
    src     = cms.InputTag("gtStage2Digis"),
    name    = cms.string("L1"),
    payload = cms.string("PhysicsTools/RDsNano/data/l1_seeds_bparking2018.json"),
)

# gen matching 
from PhysicsTools.RDsNano.genMatching_cff import *


def nanoAOD_customizeStart(process):
    # from PhysiscsTools.NanoAOD
    process.nanoSequence = cms.Sequence(nanoMetadata + globalTables + l1SeedsTable)
    return process

def nanoAOD_customizeMuonTriggerBPark(process):  #delete trigger inside
    process.triggerSequence = cms.Sequence( muonTrgSequence)
    return process

def nanoAOD_customizeBsToDsPhiKKPiMu(process):
    process.nanoBsToDsPhiKKPiMuSequence = cms.Sequence( BsToDsPhiKKPiMuSequence ) #+ CountBsToDsPhiKKPiMu)#+HighMassLowMassFlagsTables )
    return process

def nanoAOD_customizeGenMatching(process):
    process.nanoGenMatchingSequence = cms.Sequence( genMatchingSequence + puTable )

    return process
