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

       
        # bs particle, can directly access its member functions 
        bs_pt = Var("pt()",float),
        bs_eta = Var("eta()",float),
        bs_phi = Var("phi()",float),
        bs_mass = Var("mass()",float),
        bs_charge = Var("charge()",float),
        bs_pdgId = Var("pdgId()",int),
      
        #save the indices of the final states in the pruned collection
        #mu_idx = uint("mu_idx"), #always 0 :)
        k1_idx = uint("k1_idx"),
        k2_idx = uint("k2_idx"),
        pi_idx = uint("pi_idx"),

        # final state particle info
        mu_pt = Var("userCand('mu').pt()",float),
        mu_eta = Var("userCand('mu').eta()",float),
        mu_phi = Var("userCand('mu').phi()",float),
        mu_mass = Var("userCand('mu').mass()",float),
        mu_charge = Var("userCand('mu').charge()",float),
        mu_pdgId = Var("userCand('mu').pdgId()",int),

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

        # gen particle information
        gen_mu_pt = Var("userCand('gen_mu').pt()",float),
        gen_mu_eta = Var("userCand('gen_mu').eta()",float),
        gen_mu_phi = Var("userCand('gen_mu').phi()",float),
        gen_mu_mass = Var("userCand('gen_mu').mass()",float),
        gen_mu_charge = Var("userCand('gen_mu').charge()",float),
        gen_mu_pdgId = Var("userCand('gen_mu').pdgId()",int),

        gen_k1_pt = Var("userCand('gen_k1').pt()",float),
        gen_k1_eta = Var("userCand('gen_k1').eta()",float),
        gen_k1_phi = Var("userCand('gen_k1').phi()",float),
        gen_k1_mass = Var("userCand('gen_k1').mass()",float),
        gen_k1_charge = Var("userCand('gen_k1').charge()",float),
        gen_k1_pdgId = Var("userCand('gen_k1').pdgId()",int),

        gen_k2_pt = Var("userCand('gen_k2').pt()",float),
        gen_k2_eta = Var("userCand('gen_k2').eta()",float),
        gen_k2_phi = Var("userCand('gen_k2').phi()",float),
        gen_k2_mass = Var("userCand('gen_k2').mass()",float),
        gen_k2_charge = Var("userCand('gen_k2').charge()",float),
        gen_k2_pdgId = Var("userCand('gen_k2').pdgId()",int),

        gen_pi_pt = Var("userCand('gen_pi').pt()",float),
        gen_pi_eta = Var("userCand('gen_pi').eta()",float),
        gen_pi_phi = Var("userCand('gen_pi').phi()",float),
        gen_pi_mass = Var("userCand('gen_pi').mass()",float),
        gen_pi_charge = Var("userCand('gen_pi').charge()",float),
        gen_pi_pdgId = Var("userCand('gen_pi').pdgId()",int),

        gen_phi_pt  = Var("userFloat('gen_phi_pt')",float),
        gen_phi_eta = Var("userFloat('gen_phi_eta')",float),
        gen_phi_phi = Var("userFloat('gen_phi_phi')",float),
        gen_phi_vx = Var("userFloat('gen_phi_vx')",float),
        gen_phi_vy = Var("userFloat('gen_phi_vy')",float),
        gen_phi_vz = Var("userFloat('gen_phi_vz')",float),

        gen_ds_pt  = Var("userFloat('gen_ds_pt')",float),
        gen_ds_eta = Var("userFloat('gen_ds_eta')",float),
        gen_ds_phi = Var("userFloat('gen_ds_phi')",float),
        gen_ds_vx = Var("userFloat('gen_ds_vx')",float),
        gen_ds_vy = Var("userFloat('gen_ds_vy')",float),
        gen_ds_vz = Var("userFloat('gen_ds_vz')",float),

        gen_bs_pt  = Var("userFloat('gen_bs_pt')",float),
        gen_bs_eta = Var("userFloat('gen_bs_eta')",float),
        gen_bs_phi = Var("userFloat('gen_bs_phi')",float),
        gen_bs_vx = Var("userFloat('gen_bs_vx')",float),
        gen_bs_vy = Var("userFloat('gen_bs_vy')",float),
        gen_bs_vz = Var("userFloat('gen_bs_vz')",float),

        # signal id
        sig = Var("userFloat('sig')",int),
  
        # from global fit
        sv_x = Var("userFloat('sv_x')",float),
        sv_y = Var("userFloat('sv_y')",float),
        sv_z = Var("userFloat('sv_z')",float),
        bs_fitted_px = Var("userFloat('bs_fitted_px')",float),
        bs_fitted_py = Var("userFloat('bs_fitted_py')",float),
        bs_fitted_pz = Var("userFloat('bs_fitted_pz')",float),
        bs_fitted_m =  Var("userFloat('bs_fitted_m')",float),

        phi_pt = Var("userFloat('phi_pt')",float),
        phi_eta = Var("userFloat('phi_eta')",float),
        phi_phi = Var("userFloat('phi_phi')",float),
        phi_mass = Var("userFloat('phi_mass')",float),
        phi_charge = Var("userFloat('phi_charge')",float),

        fv_x = Var("userFloat('fv_x')",float),
        fv_y = Var("userFloat('fv_y')",float),
        fv_z = Var("userFloat('fv_z')",float),

        #easy_bs_vtx_x = Var("userFloat('easy_bs_vtx_x')",float),
        #easy_bs_vtx_y = Var("userFloat('easy_bs_vtx_y')",float),
        #easy_bs_vtx_z = Var("userFloat('easy_bs_vtx_z')",float),


        phi_fitted_px = Var("userFloat('phi_fitted_px')",float),
        phi_fitted_py = Var("userFloat('phi_fitted_py')",float),
        phi_fitted_pz = Var("userFloat('phi_fitted_pz')",float),
        phi_fitted_pt = Var("userFloat('phi_fitted_pt')",float),
        phi_fitted_m =  Var("userFloat('phi_fitted_m')",float),

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

        tv_x = Var("userFloat('tv_x')",float),
        tv_y = Var("userFloat('tv_y')",float),
        tv_z = Var("userFloat('tv_z')",float),
        ds_fitted_px = Var("userFloat('ds_fitted_px')",float),
        ds_fitted_py = Var("userFloat('ds_fitted_py')",float),
        ds_fitted_pz = Var("userFloat('ds_fitted_pz')",float),
        ds_fitted_pt = Var("userFloat('ds_fitted_pt')",float),
        ds_fitted_m =  Var("userFloat('ds_fitted_m')",float),

        k1_refitted_vx = Var("userFloat('k1_refitted_vx')",float),
        k1_refitted_vy = Var("userFloat('k1_refitted_vy')",float),
        k1_refitted_vz = Var("userFloat('k1_refitted_vz')",float),
        k1_refitted_px = Var("userFloat('k1_refitted_px')",float),
        k1_refitted_py = Var("userFloat('k1_refitted_py')",float),
        k1_refitted_pz = Var("userFloat('k1_refitted_pz')",float),
        k1_refitted_m =  Var("userFloat('k1_refitted_m')",float),

        k2_refitted_vx = Var("userFloat('k2_refitted_vx')",float),
        k2_refitted_vy = Var("userFloat('k2_refitted_vy')",float),
        k2_refitted_vz = Var("userFloat('k2_refitted_vz')",float),
        k2_refitted_px = Var("userFloat('k2_refitted_px')",float),
        k2_refitted_py = Var("userFloat('k2_refitted_py')",float),
        k2_refitted_pz = Var("userFloat('k2_refitted_pz')",float),
        k2_refitted_m =  Var("userFloat('k2_refitted_m')",float),

        phi_refitted_vx = Var("userFloat('phi_refitted_vx')",float),
        phi_refitted_vy = Var("userFloat('phi_refitted_vy')",float),
        phi_refitted_vz = Var("userFloat('phi_refitted_vz')",float),
        phi_refitted_px = Var("userFloat('phi_refitted_px')",float),
        phi_refitted_py = Var("userFloat('phi_refitted_py')",float),
        phi_refitted_pz = Var("userFloat('phi_refitted_pz')",float),
        phi_refitted_m =  Var("userFloat('phi_refitted_m')",float),

        pi_refitted_vx = Var("userFloat('pi_refitted_vx')",float),
        pi_refitted_vy = Var("userFloat('pi_refitted_vy')",float),
        pi_refitted_vz = Var("userFloat('pi_refitted_vz')",float),
        pi_refitted_px = Var("userFloat('pi_refitted_px')",float),
        pi_refitted_py = Var("userFloat('pi_refitted_py')",float),
        pi_refitted_pz = Var("userFloat('pi_refitted_pz')",float),
        pi_refitted_m =  Var("userFloat('pi_refitted_m')",float),

        mu_refitted_vx = Var("userFloat('mu_refitted_vx')",float),
        mu_refitted_vy = Var("userFloat('mu_refitted_vy')",float),
        mu_refitted_vz = Var("userFloat('mu_refitted_vz')",float),
        mu_refitted_px = Var("userFloat('mu_refitted_px')",float),
        mu_refitted_py = Var("userFloat('mu_refitted_py')",float),
        mu_refitted_pz = Var("userFloat('mu_refitted_pz')",float),
        mu_refitted_m =  Var("userFloat('mu_refitted_m')",float),

        ds_refitted_vx = Var("userFloat('ds_refitted_vx')",float),
        ds_refitted_vy = Var("userFloat('ds_refitted_vy')",float),
        ds_refitted_vz = Var("userFloat('ds_refitted_vz')",float),
        ds_refitted_px = Var("userFloat('ds_refitted_px')",float),
        ds_refitted_py = Var("userFloat('ds_refitted_py')",float),
        ds_refitted_pz = Var("userFloat('ds_refitted_pz')",float),
        ds_refitted_m =  Var("userFloat('ds_refitted_m')",float),

        cosPiK1 = Var("userFloat('cosPiK1')",float),

        cosPiK2 = Var("userFloat('cosPiK2')",float),

        cosMuWColl = Var("userFloat('cosMuWColl')",float),
        cosMuWLhcb = Var("userFloat('cosMuWLhcb')",float),
        cosMuWLhcbAlt = Var("userFloat('cosMuWLhcbAlt')",float),
        cosMuWReco1 = Var("userFloat('cosMuWReco1')",float),
        cosMuWReco2 = Var("userFloat('cosMuWReco1')",float),
        cosMuWGen = Var("userFloat('cosMuWGen')",float),

        cosPhiDs = Var("userFloat('cosPhiDs')",float),
        cosPiDs = Var("userFloat('cosPiDs')",float),
        angPiK1 = Var("userFloat('angPiK1')",float),
        angPiK2 = Var("userFloat('angPiK2')",float),

        angMuWColl = Var("userFloat('angMuWColl')",float),

        angPhiDs = Var("userFloat('angPhiDs')",float),
        angPiDs = Var("userFloat('angPiDs')",float),

        pv_x = Var("userFloat('pv_x')",float),        
        pv_y = Var("userFloat('pv_y')",float),        
        pv_z = Var("userFloat('pv_z')",float),        

        lxy_bs = Var("userFloat('lxy_bs')",float),        
        lxyz_bs = Var("userFloat('lxyz_bs')",float),        

        lxy_ds = Var("userFloat('lxy_ds')",float),        
        lxyz_ds = Var("userFloat('lxyz_ds')",float),        

        lxy_phi = Var("userFloat('lxy_phi')",float),        
        lxyz_phi = Var("userFloat('lxyz_phi')",float),        

        dxy_mu = Var("userFloat('dxy_mu')",float),
        dz_mu = Var("userFloat('dz_mu')",float),
        dxy_mu_err = Var("userFloat('dxy_mu_err')",float),
        dz_mu_err = Var("userFloat('dz_mu_err')",float),
        dxy_mu_sig = Var("userFloat('dxy_mu_sig')",float),
        dz_mu_sig = Var("userFloat('dz_mu_sig')",float),

        dxy_pi = Var("userFloat('dxy_pi')",float),
        dz_pi = Var("userFloat('dz_pi')",float),
        dxy_pi_err = Var("userFloat('dxy_pi_err')",float),
        dz_pi_err = Var("userFloat('dz_pi_err')",float),
        dxy_pi_sig = Var("userFloat('dxy_pi_sig')",float),
        dz_pi_sig = Var("userFloat('dz_pi_sig')",float),

        dxy_k1 = Var("userFloat('dxy_k1')",float),
        dz_k1 = Var("userFloat('dz_k1')",float),
        dxy_k1_err = Var("userFloat('dxy_k1_err')",float),
        dz_k1_err = Var("userFloat('dz_k1_err')",float),
        dxy_k1_sig = Var("userFloat('dxy_k1_sig')",float),
        dz_k1_sig = Var("userFloat('dz_k1_sig')",float),

        dxy_k2 = Var("userFloat('dxy_k2')",float),
        dz_k2 = Var("userFloat('dz_k2')",float),
        dxy_k2_err = Var("userFloat('dxy_k2_err')",float),
        dz_k2_err = Var("userFloat('dz_k2_err')",float),
        dxy_k2_sig = Var("userFloat('dxy_k2_sig')",float),
        dz_k2_sig = Var("userFloat('dz_k2_sig')",float),

        bsMassCorr = Var("userFloat('bs_mass_corr')",float), 

        pi_mu_charge = Var("userInt('pi_mu_charge')",int),
        kk_charge = Var("userInt('kk_charge')",int),

        bs_pt_coll = Var("userFloat('bs_pt_coll')",float),
        bs_eta_coll = Var("userFloat('bs_eta_coll')",float),
        bs_phi_coll = Var("userFloat('bs_phi_coll')",float),

        bs_pt_lhcb = Var("userFloat('bs_pt_lhcb')",float),
        bs_eta_lhcb = Var("userFloat('bs_eta_lhcb')",float),
        bs_phi_lhcb = Var("userFloat('bs_phi_lhcb')",float),

        bs_pt_lhcb_alt = Var("userFloat('bs_pt_lhcb_alt')",float),
        bs_eta_lhcb_alt = Var("userFloat('bs_eta_lhcb_alt')",float),
        bs_phi_lhcb_alt = Var("userFloat('bs_phi_lhcb_alt')",float),

        bs_pt_reco_1 = Var("userFloat('bs_pt_reco_1')",float),
        bs_pt_reco_2 = Var("userFloat('bs_pt_reco_2')",float),
        bs_eta_reco_1 = Var("userFloat('bs_eta_reco_1')",float),
        bs_eta_reco_2 = Var("userFloat('bs_eta_reco_2')",float),
        bs_phi_reco_1 = Var("userFloat('bs_phi_reco_1')",float),
        bs_phi_reco_2 = Var("userFloat('bs_phi_reco_2')",float),

        lhcb_pz = Var("userFloat('lhcb_pz')",float),
        theta = Var("userFloat('theta')",float), 

        m2_miss_gen = Var("userFloat('m2_miss_gen')",float),
        m2_miss_coll = Var("userFloat('m2_miss_coll')",float),
        m2_miss_lhcb = Var("userFloat('m2_miss_lhcb')",float),
        m2_miss_lhcb_alt = Var("userFloat('m2_miss_lhcb_alt')",float),
        m2_miss_reco_1 = Var("userFloat('m2_miss_reco_1')",float),
        m2_miss_reco_2 = Var("userFloat('m2_miss_reco_2')",float),

        q2_gen = Var("userFloat('q2_gen')",float),
        q2_coll = Var("userFloat('q2_coll')",float),
        q2_lhcb = Var("userFloat('q2_lhcb')",float),
        q2_lhcb_alt = Var("userFloat('q2_lhcb_alt')",float),
        q2_reco_1 = Var("userFloat('q2_reco_1')",float),
        q2_reco_2 = Var("userFloat('q2_reco_2')",float),


        gen_match_success = Var("userInt('gen_match_success')",int),
        arrived = Var("userInt('arrived')",int),
)

empty = cms.PSet(arrived = Var("userInt('arrived')",int))

