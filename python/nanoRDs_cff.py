from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *   #why do we need all these nanoaod functions?
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *

##for gen and trigger muon
from PhysicsTools.RDsNano.muonsBPark_cff import * 

## B collections
from PhysicsTools.RDsNano.BsToDsPhiKKPiMu_cff import *


#G: nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)  #purpose?

# from PhysiscsTools.NanoAOD
nanoSequence = cms.Sequence(nanoMetadata + globalTables)


def nanoAOD_customizeMuonTriggerBPark(process):  #delete trigger inside
    process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonBParkTables)
    return process

def nanoAOD_customizeBsToDsPhiKKPiMu(process):
    process.nanoBsToDsPhiKKPiMuSequence = cms.Sequence( BsToDsPhiKKPiMuSequence ) #+ CountBsToDsPhiKKPiMu)#+HighMassLowMassFlagsTables )
    return process

