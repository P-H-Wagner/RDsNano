from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *   #why do we need all these nanoaod functions?
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *

# for trigger muon
from PhysicsTools.RDsNano.muonsBPark_cff import * 

# Bs chain collection
from PhysicsTools.RDsNano.BsToDsPhiKKPiMu_cff import *

# gen matching 
from PhysicsTools.RDsNano.genMatching_cff import *

#G: nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)  #purpose?


def nanoAOD_customizeStart(process):
    # from PhysiscsTools.NanoAOD
    process.nanoSequence = cms.Sequence(nanoMetadata + globalTables)
    return process

def nanoAOD_customizeMuonTriggerBPark(process):  #delete trigger inside
    process.triggerSequence = cms.Sequence( muonTrgSequence)
    return process

def nanoAOD_customizeBsToDsPhiKKPiMu(process):
    process.nanoBsToDsPhiKKPiMuSequence = cms.Sequence( BsToDsPhiKKPiMuSequence ) #+ CountBsToDsPhiKKPiMu)#+HighMassLowMassFlagsTables )
    return process

def nanoAOD_customizeGenMatching(process):
    process.nanoGenMatchingSequence = cms.Sequence( genMatchingSequence )

    return process
