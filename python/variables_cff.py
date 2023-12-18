import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RDsNano.common_cff import RDsCandVars, ufloat, uint, ubool
 
BsToDsPhiKKPiMuVariables = cms.PSet(

        #k1_pt = Var("userCand('k1').pt()"),
        #k1_eta = Var("userCand('k1').eta()"),
        #k1_phi = Var("userCand('k1').phi()"),
        #k1_mass = Var("userCand('k1').mass()"),
        #k1_charge = Var("userCand('k1').charge()"),
        #k1_pdgId = Var("userCand('k1').pdgId()"),

        #k2_pt = Var("userCand('k2').pt()"),
        #k2_eta = Var("userCand('k2').eta()"),
        #k2_phi = Var("userCand('k2').phi()"),
        #k2_mass = Var("userCand('k2').mass()"),
        #k2_charge = Var("userCand('k2').charge()"),
        #k2_pdgId = Var("userCand('k2').pdgId()"),

        bs_pt = Var("pt()",float),
        bs_eta = Var("eta()",float),
        bs_phi = Var("phi()",float),
        bs_mass = Var("mass()",float),
        bs_charge = Var("charge()",float),
        bs_pdgId = Var("pdgId()",int),

        k1_pt = Var("userCand('k1').pt()",float),
        k1_eta = Var("userCand('k1').eta()",float),
        k1_phi = Var("userCand('k1').phi()",float),
        k1_mass = ufloat('k1_mass'),
        k1_charge = Var("userCand('k1').charge()",float),
        k1_pdgId = Var("userCand('k1').pdgId()",int),

        k2_pt = Var("userCand('k2').pt()",float),
        k2_eta = Var("userCand('k2').eta()",float),
        k2_phi = Var("userCand('k2').phi()",float),
        k2_mass = ufloat('k2_mass'),
        k2_charge = Var("userCand('k2').charge()",float),
        k2_pdgId = Var("userCand('k2').pdgId()",int),

        pi_pt = Var("userCand('pi').pt()",float),
        pi_eta = Var("userCand('pi').eta()",float),
        pi_phi = Var("userCand('pi').phi()",float),
        pi_mass = Var("userCand('pi').mass()",float),
        pi_charge = Var("userCand('pi').charge()",float),
        pi_pdgId = Var("userCand('pi').pdgId()",int),

        mu_pt = Var("userCand('mu').pt()",float),
        mu_eta = Var("userCand('mu').eta()",float),
        mu_phi = Var("userCand('mu').phi()",float),
        mu_mass = Var("userCand('mu').mass()",float),
        mu_charge = Var("userCand('mu').charge()",float),
        mu_pdgId = Var("userCand('mu').pdgId()",int),

        phi_pt = Var("userFloat('phi_pt')",float),
        phi_eta = Var("userFloat('phi_eta')",float),
        phi_phi = Var("userFloat('phi_phi')",float),
        phi_mass = Var("userFloat('phi_mass')",float),
        phi_charge = Var("userFloat('phi_charge')",float),

        ds_pt = Var("userFloat('ds_pt')",float),
        ds_eta = Var("userFloat('ds_eta')",float),
        ds_phi = Var("userFloat('ds_phi')",float),
        ds_mass = Var("userFloat('ds_mass')",float),
        ds_charge = Var("userFloat('ds_charge')",float),

        dsMu_pt = Var("userFloat('dsMu_pt')",float),
        dsMu_eta = Var("userFloat('dsMu_eta')",float),
        dsMu_phi = Var("userFloat('dsMu_phi')",float),
        dsMu_mass = Var("userFloat('dsMu_mass')",float),
        dsMu_charge = Var("userFloat('dsMu_charge')",float),


        piMu_charge = uint("piMuCharge"),
        kk_charge = uint("kkCharge"),

)
