import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RDsNano.common_cff import RDsCandVars, ufloat, uint, ubool

###################################################
## Here we define all variables we want to store ##
## in the output file. Remark that in python, we ##
## we can not give more than 255 arguments to a  ##
## function (here PSet (..)). Therefore, we need ##
## several PSets and combine them later in the   ##
## BsToDsPhiKKPiMu_cff.py                        ##
###################################################

print("i am in the table!!")
## TODO check variable types (int/float)

##################################################

triggerVariables = cms.PSet(

        muon_idx             = Var("userInt('muIdx')", int),
        trg_obj_idx          = Var("userInt('trgObjIdx')", int),
        
        mu7_ip4              = Var("userInt('mu7_ip4')", int),
        prescale_mu7_ip4     = Var("userFloat('prescale_mu7_ip4')"   , float),

        mu8_ip3              = Var("userInt('mu8_ip3')", int),
        prescale_mu8_ip3     = Var("userFloat('prescale_mu8_ip3')"   , float),

        mu8p5_3p5            = Var("userInt('mu8p5_ip3p5')", int),
        prescale_mu8p5_3p5   = Var("userFloat('prescale_mu8p5_ip3p5')"   , float),

        mu8_ip5              = Var("userInt('mu8_ip5')", int),
        prescale_mu8_ip5     = Var("userFloat('prescale_mu8_ip5')"   , float),

        mu8_ip6              = Var("userInt('mu8_ip6')", int),
        prescale_mu8_ip6     = Var("userFloat('prescale_mu8_ip6')"   , float),

        mu9_ip4              = Var("userInt('mu9_ip4')", int),
        prescale_mu9_ip4     = Var("userFloat('prescale_mu9_ip4')"   , float),

        mu9_ip5              = Var("userInt('mu9_ip5')", int),
        prescale_mu9_ip5     = Var("userFloat('prescale_mu9_ip5')"   , float),

        mu9_ip6              = Var("userInt('mu9_ip6')", int),
        prescale_mu9_ip6     = Var("userFloat('prescale_mu9_ip6')"   , float),

        mu10p5_3p5           = Var("userInt('mu10p5_ip3p5')", int),
        prescale_mu10p5_3p5  = Var("userFloat('prescale_mu10p5_ip3p5')"   , float),

        mu12_ip6             = Var("userInt('mu12_ip6')", int),
        prescale_mu12_ip6    = Var("userFloat('prescale_mu12_ip6')"   , float),

        filter_mu7_4         = Var("userInt('filter_mu7_4')"       , int),
        filter_mu8_3         = Var("userInt('filter_mu8_3')"       , int),
        filter_mu8p5_3p5     = Var("userInt('filter_mu8p5_3p5')"   , int),
        filter_mu8_5         = Var("userInt('filter_mu8_5')"       , int),
        filter_mu8_6         = Var("userInt('filter_mu8_6')"       , int),
        filter_mu9_4         = Var("userInt('filter_mu9_4')"       , int),
        filter_mu9_5         = Var("userInt('filter_mu9_5')"       , int),
        filter_mu9_6         = Var("userInt('filter_mu9_6')"       , int),
        filter_mu10p5_3p5    = Var("userInt('filter_mu10p5_3p5')"   , int),
        filter_mu12_6        = Var("userInt('filter_mu12_6')"      , int),

)
 
prefitBasicVariables = cms.PSet(

        #save the indices of the final states in the pruned collection
        #mu_idx = uint("mu_idx"), #always 0 :) 
        k1_idx = uint("k1_idx"),
        k2_idx = uint("k2_idx"),
        pi_idx = uint("pi_idx"),

        # final state particle info, access via userCand()
        mu_px_old         = Var("userCand('mu').px()",float),
        mu_py_old         = Var("userCand('mu').py()",float),
        mu_pz_old         = Var("userCand('mu').pz()",float),
        mu_pt_old         = Var("userCand('mu').pt()",float),
        mu_eta_old        = Var("userCand('mu').eta()",float),
        mu_phi_old        = Var("userCand('mu').phi()",float),

        k1_px_old         = Var("userCand('k1').px()",float),
        k1_py_old         = Var("userCand('k1').py()",float),
        k1_pz_old         = Var("userCand('k1').pz()",float),
        k1_pt_old         = Var("userCand('k1').pt()",float),
        k1_eta_old        = Var("userCand('k1').eta()",float),
        k1_phi_old        = Var("userCand('k1').phi()",float),

        k2_px_old         = Var("userCand('k2').px()",float),
        k2_py_old         = Var("userCand('k2').py()",float),
        k2_pz_old         = Var("userCand('k2').pz()",float),
        k2_pt_old         = Var("userCand('k2').pt()",float),
        k2_eta_old        = Var("userCand('k2').eta()",float),
        k2_phi_old        = Var("userCand('k2').phi()",float),

        pi_px_old         = Var("userCand('pi').px()",float),
        pi_py_old         = Var("userCand('pi').py()",float),
        pi_pz_old         = Var("userCand('pi').pz()",float),
        pi_pt_old         = Var("userCand('pi').pt()",float),
        pi_eta_old        = Var("userCand('pi').eta()",float),
        pi_phi_old        = Var("userCand('pi').phi()",float),
 

        # CORRECTED MOMENTA (NEW) 
        mu_px         = Var("userFloat('mu_px')",float),
        mu_py         = Var("userFloat('mu_py')",float),
        mu_pz         = Var("userFloat('mu_pz')",float),
        mu_pt         = Var("userFloat('mu_pt')",float),
        mu_eta        = Var("userFloat('mu_eta')",float),
        mu_phi        = Var("userFloat('mu_phi')",float),
        mu_m          = Var("userFloat('mu_m')",float),
        mu_charge     = Var("userFloat('mu_charge')",float),
        mu_pdgId      = Var("userInt('mu_pdgId')",int),
        
        k1_px         = Var("userFloat('k1_px')",float),
        k1_py         = Var("userFloat('k1_py')",float),
        k1_pz         = Var("userFloat('k1_pz')",float),
        k1_pt         = Var("userFloat('k1_pt')",float),
        k1_eta        = Var("userFloat('k1_eta')",float),
        k1_phi        = Var("userFloat('k1_phi')",float),
        k1_m          = Var("userFloat('k1_m')",float),
        k1_charge     = Var("userFloat('k1_charge')",float),
        k1_pdgId      = Var("userInt('k1_pdgId')",int),
        
        k2_px         = Var("userFloat('k2_px')",float),
        k2_py         = Var("userFloat('k2_py')",float),
        k2_pz         = Var("userFloat('k2_pz')",float),
        k2_pt         = Var("userFloat('k2_pt')",float),
        k2_eta        = Var("userFloat('k2_eta')",float),
        k2_phi        = Var("userFloat('k2_phi')",float),
        k2_m          = Var("userFloat('k2_m')",float),
        k2_charge     = Var("userFloat('k2_charge')",float),
        k2_pdgId      = Var("userInt('k2_pdgId')",int),
        
        pi_px         = Var("userFloat('pi_px')",float),
        pi_py         = Var("userFloat('pi_py')",float),
        pi_pz         = Var("userFloat('pi_pz')",float),
        pi_pt         = Var("userFloat('pi_pt')",float),
        pi_eta        = Var("userFloat('pi_eta')",float),
        pi_phi        = Var("userFloat('pi_phi')",float),
        pi_m          = Var("userFloat('pi_m')",float),
        pi_charge     = Var("userFloat('pi_charge')",float),
        pi_pdgId      = Var("userInt('pi_pdgId')",int),

        mu_cov_00     = Var("userFloat('mu_cov_00')",float),
        mu_cov_01     = Var("userFloat('mu_cov_01')",float),
        mu_cov_02     = Var("userFloat('mu_cov_02')",float),
        mu_cov_03     = Var("userFloat('mu_cov_03')",float),
        mu_cov_04     = Var("userFloat('mu_cov_04')",float),
        mu_cov_11     = Var("userFloat('mu_cov_11')",float),
        mu_cov_12     = Var("userFloat('mu_cov_12')",float),
        mu_cov_13     = Var("userFloat('mu_cov_13')",float),
        mu_cov_14     = Var("userFloat('mu_cov_14')",float),
        mu_cov_22     = Var("userFloat('mu_cov_22')",float),
        mu_cov_23     = Var("userFloat('mu_cov_23')",float),
        mu_cov_24     = Var("userFloat('mu_cov_24')",float),
        mu_cov_33     = Var("userFloat('mu_cov_33')",float),
        mu_cov_34     = Var("userFloat('mu_cov_34')",float),
        mu_cov_44     = Var("userFloat('mu_cov_44')",float),

        k1_cov_00     = Var("userFloat('k1_cov_00')",float),
        k1_cov_01     = Var("userFloat('k1_cov_01')",float),
        k1_cov_02     = Var("userFloat('k1_cov_02')",float),
        k1_cov_03     = Var("userFloat('k1_cov_03')",float),
        k1_cov_04     = Var("userFloat('k1_cov_04')",float),
        k1_cov_11     = Var("userFloat('k1_cov_11')",float),
        k1_cov_12     = Var("userFloat('k1_cov_12')",float),
        k1_cov_13     = Var("userFloat('k1_cov_13')",float),
        k1_cov_14     = Var("userFloat('k1_cov_14')",float),
        k1_cov_22     = Var("userFloat('k1_cov_22')",float),
        k1_cov_23     = Var("userFloat('k1_cov_23')",float),
        k1_cov_24     = Var("userFloat('k1_cov_24')",float),
        k1_cov_33     = Var("userFloat('k1_cov_33')",float),
        k1_cov_34     = Var("userFloat('k1_cov_34')",float),
        k1_cov_44     = Var("userFloat('k1_cov_44')",float),

        k2_cov_00     = Var("userFloat('k2_cov_00')",float),
        k2_cov_01     = Var("userFloat('k2_cov_01')",float),
        k2_cov_02     = Var("userFloat('k2_cov_02')",float),
        k2_cov_03     = Var("userFloat('k2_cov_03')",float),
        k2_cov_04     = Var("userFloat('k2_cov_04')",float),
        k2_cov_11     = Var("userFloat('k2_cov_11')",float),
        k2_cov_12     = Var("userFloat('k2_cov_12')",float),
        k2_cov_13     = Var("userFloat('k2_cov_13')",float),
        k2_cov_14     = Var("userFloat('k2_cov_14')",float),
        k2_cov_22     = Var("userFloat('k2_cov_22')",float),
        k2_cov_23     = Var("userFloat('k2_cov_23')",float),
        k2_cov_24     = Var("userFloat('k2_cov_24')",float),
        k2_cov_33     = Var("userFloat('k2_cov_33')",float),
        k2_cov_34     = Var("userFloat('k2_cov_34')",float),
        k2_cov_44     = Var("userFloat('k2_cov_44')",float),

        pi_cov_00     = Var("userFloat('pi_cov_00')",float),
        pi_cov_01     = Var("userFloat('pi_cov_01')",float),
        pi_cov_02     = Var("userFloat('pi_cov_02')",float),
        pi_cov_03     = Var("userFloat('pi_cov_03')",float),
        pi_cov_04     = Var("userFloat('pi_cov_04')",float),
        pi_cov_11     = Var("userFloat('pi_cov_11')",float),
        pi_cov_12     = Var("userFloat('pi_cov_12')",float),
        pi_cov_13     = Var("userFloat('pi_cov_13')",float),
        pi_cov_14     = Var("userFloat('pi_cov_14')",float),
        pi_cov_22     = Var("userFloat('pi_cov_22')",float),
        pi_cov_23     = Var("userFloat('pi_cov_23')",float),
        pi_cov_24     = Var("userFloat('pi_cov_24')",float),
        pi_cov_33     = Var("userFloat('pi_cov_33')",float),
        pi_cov_34     = Var("userFloat('pi_cov_34')",float),
        pi_cov_44     = Var("userFloat('pi_cov_44')",float),


        mu_is_tracker    = Var("userInt('mu_is_tracker')",int),
        mu_is_pf         = Var("userInt('mu_is_pf')",int),
        mu_is_global     = Var("userInt('mu_is_global')",int),
        mu_is_standalone = Var("userInt('mu_is_standalone')",int),

        mu_id_loose      = Var("userInt('mu_id_loose')",int),
        mu_id_medium     = Var("userInt('mu_id_medium')",int),
        mu_id_tight      = Var("userInt('mu_id_tight')",int),

        mu_id_medium_promt        = Var("userInt('mu_id_medium_prompt')",       int),
        mu_id_global_high_pt      = Var("userInt('mu_id_global_high_pt')",      int),
        mu_id_trk_high_pt         = Var("userInt('mu_id_trk_high_pt')",         int),
        mu_pf_iso_very_loose      = Var("userInt('mu_pf_iso_very_loose')",      int),
        mu_pf_iso_loose           = Var("userInt('mu_pf_iso_loose')",           int),
        mu_pf_iso_medium          = Var("userInt('mu_pf_iso_medium')",          int),
        mu_pf_iso_tight           = Var("userInt('mu_pf_iso_tight')",           int),
        mu_pf_iso_very_tight      = Var("userInt('mu_pf_iso_very_tight')",      int),

        mu_tk_iso_loose           = Var("userInt('mu_tk_iso_loose')",           int),
        mu_tk_iso_tight           = Var("userInt('mu_tk_iso_tight')",           int),
        mu_id_soft                = Var("userInt('mu_id_soft')",                int),
        mu_id_soft_mva            = Var("userInt('mu_id_soft_mva')",            int),
        mu_mva_loose              = Var("userInt('mu_mva_loose')",              int),
        mu_mva_medium             = Var("userInt('mu_mva_medium')",             int),
        mu_mvs_tight              = Var("userInt('mu_mva_tight')",              int),
        mu_mini_iso_loose         = Var("userInt('mu_mini_iso_loose')",         int),
        mu_mini_iso_medium        = Var("userInt('mu_mini_iso_medium')",        int),
        mu_mini_iso_tight         = Var("userInt('mu_mini_iso_tight')",         int),

        mu_mini_iso_very_tight    = Var("userInt('mu_mini_iso_very_tight')",    int),
        mu_trigger_id_loose       = Var("userInt('mu_trigger_id_loose')",       int),
        mu_in_time_muon           = Var("userInt('mu_in_time_muon')",           int),
        mu_pf_iso_very_very_tight = Var("userInt('mu_pf_iso_very_very_tight')", int),
        mu_multi_iso_loose        = Var("userInt('mu_multi_iso_loose')",        int),
        mu_multi_iso_medium       = Var("userInt('mu_multi_iso_medium')",       int),
        mu_puppi_iso_loose        = Var("userInt('mu_puppi_iso_loose')",        int),
        mu_puppi_iso_medium       = Var("userInt('mu_puppi_iso_medium')",       int),
        mu_puppi_iso_tight        = Var("userInt('mu_puppi_iso_tight')",        int),
        mu_mva_v_tight            = Var("userInt('mu_mva_v_tight')",            int),

        mu_mva_vv_tight           = Var("userInt('mu_mva_vv_tight')",           int),
        mu_low_pt_mva_loose       = Var("userInt('mu_low_pt_mva_loose')",       int),
        mu_low_pt_mva_medium      = Var("userInt('mu_low_pt_mva_medium')",      int),
        mu_mv_id_wp_medium        = Var("userInt('mu_mv_id_wp_medium')",        int),
        mu_mv_id_wp_tight         = Var("userInt('mu_mv_id_wp_tight')",         int),

        ## add prefit resonances (f.e. Ds is simply phi + pi, without any fit)

        kk_pt         = Var("userFloat('kk_pt')",float),
        kk_eta        = Var("userFloat('kk_eta')",float),
        kk_phi        = Var("userFloat('kk_phi')",float),
        kk_m       = Var("userFloat('kk_m')",float),
        kk_charge     = Var("userFloat('kk_charge')",float),
        kk_deltaR     = Var("userFloat('kk_deltaR')",float),

        phiPi_pt      = Var("userFloat('phiPi_pt')",float),
        phiPi_eta     = Var("userFloat('phiPi_eta')",float),
        phiPi_phi     = Var("userFloat('phiPi_phi')",float),
        phiPi_m    = Var("userFloat('phiPi_m')",float),
        phiPi_charge  = Var("userFloat('phiPi_charge')",float),
        phiPi_deltaR  = Var("userFloat('phiPi_deltaR')",float),
       
        dsMu_pt       = Var("userFloat('dsMu_pt')",float),
        dsMu_eta      = Var("userFloat('dsMu_eta')",float),
        dsMu_phi      = Var("userFloat('dsMu_phi')",float),
        dsMu_m     = Var("userFloat('dsMu_m')",float),
        dsMu_charge   = Var("userFloat('dsMu_charge')",float),
        dsMu_deltaR   = Var("userFloat('dsMu_deltaR')",float),

        # bs particle, can directly access its member functions 
        # this bs is obtained from the dsMu + coll approximation -> simple approach!
        bs_px         = Var("px()",float),
        bs_py         = Var("py()",float),
        bs_pz         = Var("pz()",float),
        bs_pt         = Var("pt()",float),
        bs_eta        = Var("eta()",float),
        bs_phi        = Var("phi()",float),
        bs_m       = Var("mass()",float),
        bs_charge     = Var("charge()",float),
        bs_pdgId      = Var("pdgId()",int),

        # rel charge, important for flattening
        pi_mu_charge  = Var("userFloat('pi_mu_charge')",float),

        matchedMu     = Var("userInt('matchedMu')",int),
        matchedK1     = Var("userInt('matchedK1')",int),
        matchedK2     = Var("userInt('matchedK2')",int),
        matchedPi     = Var("userInt('matchedPi')",int),

)

##################################################

genVariables = cms.PSet(
        ## gen particle information
        mu_px     = Var("userFloat('mu_gen_px')",float),
        mu_py     = Var("userFloat('mu_gen_py')",float),
        mu_pz     = Var("userFloat('mu_gen_pz')",float),
        mu_pt     = Var("userFloat('mu_gen_pt')",float),
        mu_eta    = Var("userFloat('mu_gen_eta')",float),
        mu_phi    = Var("userFloat('mu_gen_phi')",float),
        mu_m      = Var("userFloat('mu_gen_m')",float),
        mu_charge = Var("userFloat('mu_gen_charge')",float),
        mu_pdgid  = Var("userInt('mu_gen_pdgid')", int),

        k1_px     = Var("userFloat('k1_gen_px')",float),
        k1_py     = Var("userFloat('k1_gen_py')",float),
        k1_pz     = Var("userFloat('k1_gen_pz')",float),
        k1_pt     = Var("userFloat('k1_gen_pt')",float),
        k1_eta    = Var("userFloat('k1_gen_eta')",float),
        k1_phi    = Var("userFloat('k1_gen_phi')",float),
        k1_m   = Var("userFloat('k1_gen_m')",float),
        k1_charge = Var("userFloat('k1_gen_charge')",float),
        k1_pdgid  = Var("userInt('k1_gen_pdgid')",int),

        k2_px     = Var("userFloat('k2_gen_px')",float),
        k2_py     = Var("userFloat('k2_gen_py')",float),
        k2_pz     = Var("userFloat('k2_gen_pz')",float),
        k2_pt     = Var("userFloat('k2_gen_pt')",float),
        k2_eta    = Var("userFloat('k2_gen_eta')",float),
        k2_phi    = Var("userFloat('k2_gen_phi')",float),
        k2_m   = Var("userFloat('k2_gen_m')",float),
        k2_charge = Var("userFloat('k2_gen_charge')",float),
        k2_pdgid  = Var("userInt('k2_gen_pdgid')",int),

        pi_px     = Var("userFloat('pi_gen_px')",float),
        pi_py     = Var("userFloat('pi_gen_py')",float),
        pi_pz     = Var("userFloat('pi_gen_pz')",float),
        pi_pt     = Var("userFloat('pi_gen_pt')",float),
        pi_eta    = Var("userFloat('pi_gen_eta')",float),
        pi_phi    = Var("userFloat('pi_gen_phi')",float),
        pi_m   = Var("userFloat('pi_gen_m')",float),
        pi_charge = Var("userFloat('pi_gen_charge')",float),
        pi_pdgid  = Var("userInt('pi_gen_pdgid')",int),

        phi_px    = Var("userFloat('phi_gen_px')",float),
        phi_py    = Var("userFloat('phi_gen_py')",float),
        phi_pz    = Var("userFloat('phi_gen_pz')",float),
        phi_pt    = Var("userFloat('phi_gen_pt')",float),
        phi_eta   = Var("userFloat('phi_gen_eta')",float),
        phi_phi   = Var("userFloat('phi_gen_phi')",float),
        tv_x      = Var("userFloat('tv_x_gen')",float),
        tv_y      = Var("userFloat('tv_y_gen')",float),
        tv_z      = Var("userFloat('tv_z_gen')",float),
        phi_charge= Var("userFloat('phi_gen_charge')",float),
        phi_pdgid = Var("userInt('phi_gen_pdgid')",int),

        ds_px     = Var("userFloat('ds_gen_px')",float),
        ds_py     = Var("userFloat('ds_gen_py')",float),
        ds_pz     = Var("userFloat('ds_gen_pz')",float),
        ds_pt     = Var("userFloat('ds_gen_pt')",float),
        ds_eta    = Var("userFloat('ds_gen_eta')",float),
        ds_phi    = Var("userFloat('ds_gen_phi')",float),
        ds_boost  = Var("userFloat('ds_gen_boost')",float),
        ds_m      = Var("userFloat('ds_gen_m')",float),

        sv_x      = Var("userFloat('sv_x_gen')",float),
        sv_y      = Var("userFloat('sv_y_gen')",float),
        sv_z      = Var("userFloat('sv_z_gen')",float),
        ds_charge = Var("userFloat('ds_gen_charge')",float),
        ds_pdgid  = Var("userInt('ds_gen_pdgid')",int),

        bs_px     = Var("userFloat('bs_gen_px')",float),
        bs_py     = Var("userFloat('bs_gen_py')",float),
        bs_pz     = Var("userFloat('bs_gen_pz')",float),
        bs_pt     = Var("userFloat('bs_gen_pt')",float),
        bs_eta    = Var("userFloat('bs_gen_eta')",float),
        bs_phi    = Var("userFloat('bs_gen_phi')",float),
        bs_m      = Var("userFloat('bs_gen_m')",float),

        pv_x      = Var("userFloat('pv_x_gen')",float),
        pv_y      = Var("userFloat('pv_y_gen')",float),
        pv_z      = Var("userFloat('pv_z_gen')",float),

        scnd_pv_x      = Var("userFloat('scnd_pv_x_gen')",float),
        scnd_pv_y      = Var("userFloat('scnd_pv_y_gen')",float),
        scnd_pv_z      = Var("userFloat('scnd_pv_z_gen')",float),
        scnd_pv_idx    = Var("userInt('scnd_pv_idx_gen')",int),

        bs_charge = Var("userFloat('bs_gen_charge')",float),
        bs_pdgid  = Var("userInt('bs_gen_pdgid')",int),
        bs_boost  = Var("userFloat('b_boost_gen')",float),
        bs_boost_pt   = Var("userFloat('b_boost_gen_pt')",float),
        bs_boost_eta  = Var("userFloat('b_boost_gen_eta')",float),
        bs_boost_phi  = Var("userFloat('b_boost_gen_phi')",float),

        
        disc_is_negative = Var("userInt('disc_is_negative_gen')",int),
        disc_negativity  = Var("userFloat('disc_negativity_gen')",float),


        fv_x          = Var("userFloat('fv_x_gen')",float),
        fv_y          = Var("userFloat('fv_y_gen')",float),
        fv_z          = Var("userFloat('fv_z_gen')",float),

        numTrueInts   = Var("userFloat('numTrueInts')", float),
        # ---- Bmmm candidate truth (plugins/genBmmm.h; own track->gen matching, dR < 0.03)
        cand            = Var("userFloat('cand')", float, doc="Bmmm gen_cand: 0 Ds mu, 1 Ds* mu, 2 Ds tau, 3 Ds* tau, 10 Bs feed-down, 11 Bs double charm, 12-15 Ds from B0/B+/Lb/other b, 20 Ds and mu from different b, 21 mu not from b, 30 KKpi not from one Ds, -1 <4 tracks matched"),
        ds_true         = Var("userFloat('ds_true')", float, doc="K K pi from one Ds"),
        charm_pdgid     = Var("userFloat('charm_pdgid')", float, doc="common charm ancestor of K K pi"),
        ds_mother_pdgid = Var("userFloat('ds_mother_pdgid')", float),
        ds_b_pdgid      = Var("userFloat('ds_b_pdgid')", float, doc="b hadron of the Ds (first copy)"),
        mu_b_pdgid      = Var("userFloat('mu_b_pdgid')", float, doc="b hadron of the muon (first copy)"),
        mu_mother_pdgid = Var("userFloat('mu_mother_pdgid')", float),
        mu_from_tau     = Var("userFloat('mu_from_tau')", float),
        same_b          = Var("userFloat('same_b')", float, doc="Ds and muon from the same b hadron"),
        n_matched       = Var("userFloat('n_matched')", float, doc="tracks matched to gen (0-4)"),
        ds_b_ct         = Var("userFloat('ds_b_ct')", float, doc="c*t [cm] of the b hadron of the Ds (lifetime reweighting)"),
        mu_b_ct         = Var("userFloat('mu_b_ct')", float, doc="c*t [cm] of the b hadron of the muon"),
        numPUInts     = Var("userInt('numPUInts')"  , int),

        m2_miss       = ufloat('m2_miss_gen'),
        pt_miss       = ufloat('pt_miss_gen'),
        q2            = ufloat('q2_gen'),
        e_star        = ufloat('e_star_gen'),

        cosMuW        = Var("userFloat('cosMuWGen')",float),
        cosMuWReco1   = Var("userFloat('cosMuWGenReco1')",float),
        cosMuWReco2   = Var("userFloat('cosMuWGenReco2')",float),

        cosPiK1       = Var("userFloat('cosPiK1Gen')",float),

        cosPiK2       = Var("userFloat('cosPiK2Gen')",float),

        cosPiDs       = Var("userFloat('cosPiDsGen')",float),

        cosPhiDs      = Var("userFloat('cosPhiDsGen')",float),

        cosPlaneBs    = Var("userFloat('cosPlaneBsGen')",float),
 
        cosPlaneDs    = Var("userFloat('cosPlaneDsGen')",float),

        match_success = uint('gen_match_success'),

        #mu_iso_03     = Var("userFloat('mu_iso_03_gen')", float),
        #mu_iso_04     = Var("userFloat('mu_iso_04_gen')", float),
        #mu_rel_iso_03 = Var("userFloat('mu_rel_iso_03_gen')", float),
        #mu_rel_iso_04 = Var("userFloat('mu_rel_iso_04_gen')", float),

        e_gamma       = Var("userFloat('e_gamma_gen')",float),

        ## signal id
        sig               = Var("userInt('sig')",int),
        b_mother_id       = Var("userInt('b_mother_id')",int),

        tau_pt      = Var("userFloat('tau_gen_pt')",float),
        tau_eta     = Var("userFloat('tau_gen_eta')",float),
        tau_phi     = Var("userFloat('tau_gen_phi')",float),
        tau_m       = Var("userFloat('tau_gen_m')",float),
        tau_pdgid   = Var("userInt('tau_gen_pdgid')",int),

        dsStar_pt      = Var("userFloat('dsStar_gen_pt')",float),
        dsStar_eta     = Var("userFloat('dsStar_gen_eta')",float),
        dsStar_phi     = Var("userFloat('dsStar_gen_phi')",float),
        dsStar_m       = Var("userFloat('dsStar_gen_m')",float),
        dsStar_pdgid   = Var("userInt('dsStar_gen_pdgid')",int),

        photon_pt      = Var("userFloat('photon_gen_pt')",float),
        photon_eta     = Var("userFloat('photon_gen_eta')",float),
        photon_phi     = Var("userFloat('photon_gen_phi')",float),
        photon_pdgid   = Var("userInt('photon_gen_pdgid')",int),

        dr_photon_ds   = Var("userFloat('dr_gen_photon_ds')",float),

        same_mother    = Var("userInt('same_mother')",int),
)

##################################################

vertexVariables = cms.PSet(

        bs_x0        = Var("userFloat('bs_x0')",float),        
        bs_y0        = Var("userFloat('bs_y0')",float),        
        bs_z0        = Var("userFloat('bs_z0')",float),        

        pv_general_x = Var("userFloat('pv_general_x')",float),        
        pv_general_y = Var("userFloat('pv_general_y')",float),        
        pv_general_z = Var("userFloat('pv_general_z')",float),        

        pv_x         = Var("userFloat('pv_x')",float),        
        pv_y         = Var("userFloat('pv_y')",float),        
        pv_x_w_bs    = Var("userFloat('pv_x_w_bs')",float),        
        pv_y_w_bs    = Var("userFloat('pv_y_w_bs')",float),        

        pv_z         = Var("userFloat('pv_z')",float),        
        pv_x_wout_bs = Var("userFloat('pv_x_wout_bs')", float),
        pv_y_wout_bs = Var("userFloat('pv_y_wout_bs')", float),

        pv_dz_x      = Var("userFloat('pv_dz_x')",float),        
        pv_dz_y      = Var("userFloat('pv_dz_y')",float),        
        pv_dz_z      = Var("userFloat('pv_dz_z')",float),        

        pv_chi2      = Var("userFloat('pv_chi2')",float),        
        pv_ndof      = Var("userFloat('pv_ndof')",float),        
        pv_redchi2   = Var("userFloat('pv_redchi2')",float),        
        pv_prob      = Var("userFloat('pv_prob')",float),        
        pv_idx       = Var("userInt('pv_idx')",int),
        pv_cov_xx    = Var("userFloat('pv_cov_xx')",float),
        pv_cov_xy    = Var("userFloat('pv_cov_xy')",float),
        pv_cov_xz    = Var("userFloat('pv_cov_xz')",float),
        pv_cov_yy    = Var("userFloat('pv_cov_yy')",float),
        pv_cov_yz    = Var("userFloat('pv_cov_yz')",float),
        pv_cov_zz    = Var("userFloat('pv_cov_zz')",float),
        pv_ntrk      = Var("userInt('pv_ntrk')",int),        
        # ---- ported from Bmmm (same definitions; plugins/bmmmExtras.h)
        mu_cov_scale_dxy          = Var("userFloat('mu_cov_scale_dxy')", float),
        dxy_mu_err_pv_raw         = Var("userFloat('dxy_mu_err_pv_raw')", float),
        dz_mu_err_pv_raw          = Var("userFloat('dz_mu_err_pv_raw')", float),
        dxy_mu_err_beamSpot_raw   = Var("userFloat('dxy_mu_err_beamSpot_raw')", float),
        dxy_k1_pv                 = Var("userFloat('dxy_k1_pv')", float),
        dxy_k1_err_pv             = Var("userFloat('dxy_k1_err_pv')", float),
        dxy_k1_sig_pv             = Var("userFloat('dxy_k1_sig_pv')", float),
        dz_k1_pv                  = Var("userFloat('dz_k1_pv')", float),
        dz_k1_err_pv              = Var("userFloat('dz_k1_err_pv')", float),
        dz_k1_sig_pv              = Var("userFloat('dz_k1_sig_pv')", float),
        dxy_k2_pv                 = Var("userFloat('dxy_k2_pv')", float),
        dxy_k2_err_pv             = Var("userFloat('dxy_k2_err_pv')", float),
        dxy_k2_sig_pv             = Var("userFloat('dxy_k2_sig_pv')", float),
        dz_k2_pv                  = Var("userFloat('dz_k2_pv')", float),
        dz_k2_err_pv              = Var("userFloat('dz_k2_err_pv')", float),
        dz_k2_sig_pv              = Var("userFloat('dz_k2_sig_pv')", float),
        dxy_k1_beamSpot           = Var("userFloat('dxy_k1_beamSpot')", float),
        dxy_k1_err_beamSpot       = Var("userFloat('dxy_k1_err_beamSpot')", float),
        dxy_k1_sig_beamSpot       = Var("userFloat('dxy_k1_sig_beamSpot')", float),
        dxy_k2_beamSpot           = Var("userFloat('dxy_k2_beamSpot')", float),
        dxy_k2_err_beamSpot       = Var("userFloat('dxy_k2_err_beamSpot')", float),
        dxy_k2_sig_beamSpot       = Var("userFloat('dxy_k2_sig_beamSpot')", float),
        dxy_pi_beamSpot           = Var("userFloat('dxy_pi_beamSpot')", float),
        dxy_pi_err_beamSpot       = Var("userFloat('dxy_pi_err_beamSpot')", float),
        dxy_pi_sig_beamSpot       = Var("userFloat('dxy_pi_sig_beamSpot')", float),
        fv_cov_xx                 = Var("userFloat('fv_cov_xx')", float),
        fv_cov_xy                 = Var("userFloat('fv_cov_xy')", float),
        fv_cov_xz                 = Var("userFloat('fv_cov_xz')", float),
        fv_cov_yy                 = Var("userFloat('fv_cov_yy')", float),
        fv_cov_yz                 = Var("userFloat('fv_cov_yz')", float),
        fv_cov_zz                 = Var("userFloat('fv_cov_zz')", float),
        tv_cov_xx                 = Var("userFloat('tv_cov_xx')", float),
        tv_cov_xy                 = Var("userFloat('tv_cov_xy')", float),
        tv_cov_xz                 = Var("userFloat('tv_cov_xz')", float),
        tv_cov_yy                 = Var("userFloat('tv_cov_yy')", float),
        tv_cov_yz                 = Var("userFloat('tv_cov_yz')", float),
        tv_cov_zz                 = Var("userFloat('tv_cov_zz')", float),
        sv_cov_xx                 = Var("userFloat('sv_cov_xx')", float),
        sv_cov_xy                 = Var("userFloat('sv_cov_xy')", float),
        sv_cov_xz                 = Var("userFloat('sv_cov_xz')", float),
        sv_cov_yy                 = Var("userFloat('sv_cov_yy')", float),
        sv_cov_yz                 = Var("userFloat('sv_cov_yz')", float),
        sv_cov_zz                 = Var("userFloat('sv_cov_zz')", float),
        ds_vtx_cosine_xyz_sv      = Var("userFloat('ds_vtx_cosine_xyz_sv')", float),
        ds_vtx_cosine_xy_sv       = Var("userFloat('ds_vtx_cosine_xy_sv')", float),
        kk_m_err                  = Var("userFloat('kk_m_err')", float),
        phiPi_m_err               = Var("userFloat('phiPi_m_err')", float),
        phiPi_m_err_cal           = Var("userFloat('phiPi_m_err_cal')", float),
        mu_trig_dr                = Var("userFloat('mu_trig_dr')", float),
        mu_iso_ch_03              = Var("userFloat('mu_iso_ch_03')", float),
        mu_iso_ch_clean_03        = Var("userFloat('mu_iso_ch_clean_03')", float),
        mu_iso_pu_03              = Var("userFloat('mu_iso_pu_03')", float),
        mu_iso_nh_03              = Var("userFloat('mu_iso_nh_03')", float),
        mu_iso_ph_03              = Var("userFloat('mu_iso_ph_03')", float),
        mu_iso_03                 = Var("userFloat('mu_iso_03')", float),
        mu_iso_clean_03           = Var("userFloat('mu_iso_clean_03')", float),
        mu_reliso_03              = Var("userFloat('mu_reliso_03')", float),
        mu_reliso_clean_03        = Var("userFloat('mu_reliso_clean_03')", float),
        mu_iso_ch_04              = Var("userFloat('mu_iso_ch_04')", float),
        mu_iso_ch_clean_04        = Var("userFloat('mu_iso_ch_clean_04')", float),
        mu_iso_pu_04              = Var("userFloat('mu_iso_pu_04')", float),
        mu_iso_nh_04              = Var("userFloat('mu_iso_nh_04')", float),
        mu_iso_ph_04              = Var("userFloat('mu_iso_ph_04')", float),
        mu_iso_04                 = Var("userFloat('mu_iso_04')", float),
        mu_iso_clean_04           = Var("userFloat('mu_iso_clean_04')", float),
        mu_reliso_04              = Var("userFloat('mu_reliso_04')", float),
        mu_reliso_clean_04        = Var("userFloat('mu_reliso_clean_04')", float),
        ds_iso_ch_03              = Var("userFloat('ds_iso_ch_03')", float),
        ds_iso_ch_clean_03        = Var("userFloat('ds_iso_ch_clean_03')", float),
        ds_iso_pu_03              = Var("userFloat('ds_iso_pu_03')", float),
        ds_iso_nh_03              = Var("userFloat('ds_iso_nh_03')", float),
        ds_iso_ph_03              = Var("userFloat('ds_iso_ph_03')", float),
        ds_iso_03                 = Var("userFloat('ds_iso_03')", float),
        ds_iso_clean_03           = Var("userFloat('ds_iso_clean_03')", float),
        ds_reliso_03              = Var("userFloat('ds_reliso_03')", float),
        ds_reliso_clean_03        = Var("userFloat('ds_reliso_clean_03')", float),
        ds_iso_ch_04              = Var("userFloat('ds_iso_ch_04')", float),
        ds_iso_ch_clean_04        = Var("userFloat('ds_iso_ch_clean_04')", float),
        ds_iso_pu_04              = Var("userFloat('ds_iso_pu_04')", float),
        ds_iso_nh_04              = Var("userFloat('ds_iso_nh_04')", float),
        ds_iso_ph_04              = Var("userFloat('ds_iso_ph_04')", float),
        ds_iso_04                 = Var("userFloat('ds_iso_04')", float),
        ds_iso_clean_04           = Var("userFloat('ds_iso_clean_04')", float),
        ds_reliso_04              = Var("userFloat('ds_reliso_04')", float),
        ds_reliso_clean_04        = Var("userFloat('ds_reliso_clean_04')", float),
        b_iso_ch_03               = Var("userFloat('b_iso_ch_03')", float),
        b_iso_ch_clean_03         = Var("userFloat('b_iso_ch_clean_03')", float),
        b_iso_pu_03               = Var("userFloat('b_iso_pu_03')", float),
        b_iso_nh_03               = Var("userFloat('b_iso_nh_03')", float),
        b_iso_ph_03               = Var("userFloat('b_iso_ph_03')", float),
        b_iso_03                  = Var("userFloat('b_iso_03')", float),
        b_iso_clean_03            = Var("userFloat('b_iso_clean_03')", float),
        b_reliso_03               = Var("userFloat('b_reliso_03')", float),
        b_reliso_clean_03         = Var("userFloat('b_reliso_clean_03')", float),
        b_iso_ch_04               = Var("userFloat('b_iso_ch_04')", float),
        b_iso_ch_clean_04         = Var("userFloat('b_iso_ch_clean_04')", float),
        b_iso_pu_04               = Var("userFloat('b_iso_pu_04')", float),
        b_iso_nh_04               = Var("userFloat('b_iso_nh_04')", float),
        b_iso_ph_04               = Var("userFloat('b_iso_ph_04')", float),
        b_iso_04                  = Var("userFloat('b_iso_04')", float),
        b_iso_clean_04            = Var("userFloat('b_iso_clean_04')", float),
        b_reliso_04               = Var("userFloat('b_reliso_04')", float),
        b_reliso_clean_04         = Var("userFloat('b_reliso_clean_04')", float),
        npv                       = Var("userInt('npv')", int),
        mu_n_pix_hit              = Var("userInt('mu_n_pix_hit')", int),
        mu_n_pix_b_hit            = Var("userInt('mu_n_pix_b_hit')", int),
        mu_n_pix_e_hit            = Var("userInt('mu_n_pix_e_hit')", int),
        mu_n_pix_layer            = Var("userInt('mu_n_pix_layer')", int),
        mu_n_trk_layer            = Var("userInt('mu_n_trk_layer')", int),
        mu_pix_first_b_layer      = Var("userInt('mu_pix_first_b_layer')", int),
        mu_pix_first_e_disk       = Var("userInt('mu_pix_first_e_disk')", int),
        mu_pix_first_layer        = Var("userInt('mu_pix_first_layer')", int),
        mu_n_pix_miss_inner       = Var("userInt('mu_n_pix_miss_inner')", int),
        mu_n_pix_inact_inner      = Var("userInt('mu_n_pix_inact_inner')", int),
        mu_pix_valid_mask         = Var("userInt('mu_pix_valid_mask')", int),
        mu_pix_miss_mask          = Var("userInt('mu_pix_miss_mask')", int),
        mu_pix_inact_mask         = Var("userInt('mu_pix_inact_mask')", int),
        mu_pix_hit_count          = Var("userInt('mu_pix_hit_count')", int),
        k1_n_pix_hit              = Var("userInt('k1_n_pix_hit')", int),
        k1_n_pix_b_hit            = Var("userInt('k1_n_pix_b_hit')", int),
        k1_n_pix_e_hit            = Var("userInt('k1_n_pix_e_hit')", int),
        k1_n_pix_layer            = Var("userInt('k1_n_pix_layer')", int),
        k1_n_trk_layer            = Var("userInt('k1_n_trk_layer')", int),
        k1_pix_first_b_layer      = Var("userInt('k1_pix_first_b_layer')", int),
        k1_pix_first_e_disk       = Var("userInt('k1_pix_first_e_disk')", int),
        k1_pix_first_layer        = Var("userInt('k1_pix_first_layer')", int),
        k1_n_pix_miss_inner       = Var("userInt('k1_n_pix_miss_inner')", int),
        k1_n_pix_inact_inner      = Var("userInt('k1_n_pix_inact_inner')", int),
        k1_pix_valid_mask         = Var("userInt('k1_pix_valid_mask')", int),
        k1_pix_miss_mask          = Var("userInt('k1_pix_miss_mask')", int),
        k1_pix_inact_mask         = Var("userInt('k1_pix_inact_mask')", int),
        k1_pix_hit_count          = Var("userInt('k1_pix_hit_count')", int),
        k2_n_pix_hit              = Var("userInt('k2_n_pix_hit')", int),
        k2_n_pix_b_hit            = Var("userInt('k2_n_pix_b_hit')", int),
        k2_n_pix_e_hit            = Var("userInt('k2_n_pix_e_hit')", int),
        k2_n_pix_layer            = Var("userInt('k2_n_pix_layer')", int),
        k2_n_trk_layer            = Var("userInt('k2_n_trk_layer')", int),
        k2_pix_first_b_layer      = Var("userInt('k2_pix_first_b_layer')", int),
        k2_pix_first_e_disk       = Var("userInt('k2_pix_first_e_disk')", int),
        k2_pix_first_layer        = Var("userInt('k2_pix_first_layer')", int),
        k2_n_pix_miss_inner       = Var("userInt('k2_n_pix_miss_inner')", int),
        k2_n_pix_inact_inner      = Var("userInt('k2_n_pix_inact_inner')", int),
        k2_pix_valid_mask         = Var("userInt('k2_pix_valid_mask')", int),
        k2_pix_miss_mask          = Var("userInt('k2_pix_miss_mask')", int),
        k2_pix_inact_mask         = Var("userInt('k2_pix_inact_mask')", int),
        k2_pix_hit_count          = Var("userInt('k2_pix_hit_count')", int),
        pi_n_pix_hit              = Var("userInt('pi_n_pix_hit')", int),
        pi_n_pix_b_hit            = Var("userInt('pi_n_pix_b_hit')", int),
        pi_n_pix_e_hit            = Var("userInt('pi_n_pix_e_hit')", int),
        pi_n_pix_layer            = Var("userInt('pi_n_pix_layer')", int),
        pi_n_trk_layer            = Var("userInt('pi_n_trk_layer')", int),
        pi_pix_first_b_layer      = Var("userInt('pi_pix_first_b_layer')", int),
        pi_pix_first_e_disk       = Var("userInt('pi_pix_first_e_disk')", int),
        pi_pix_first_layer        = Var("userInt('pi_pix_first_layer')", int),
        pi_n_pix_miss_inner       = Var("userInt('pi_n_pix_miss_inner')", int),
        pi_n_pix_inact_inner      = Var("userInt('pi_n_pix_inact_inner')", int),
        pi_pix_valid_mask         = Var("userInt('pi_pix_valid_mask')", int),
        pi_pix_miss_mask          = Var("userInt('pi_pix_miss_mask')", int),
        pi_pix_inact_mask         = Var("userInt('pi_pix_inact_mask')", int),
        pi_pix_hit_count          = Var("userInt('pi_pix_hit_count')", int),
        mu_best_trk_type          = Var("userInt('mu_best_trk_type')", int),
        k1_from_pv                = Var("userInt('k1_from_pv')", int),
        k1_pv_assoc_q             = Var("userInt('k1_pv_assoc_q')", int),
        k1_high_purity            = Var("userInt('k1_high_purity')", int),
        k1_lost_inner_hits        = Var("userInt('k1_lost_inner_hits')", int),
        k1_is_lost                = Var("userInt('k1_is_lost')", int),
        k2_from_pv                = Var("userInt('k2_from_pv')", int),
        k2_pv_assoc_q             = Var("userInt('k2_pv_assoc_q')", int),
        k2_high_purity            = Var("userInt('k2_high_purity')", int),
        k2_lost_inner_hits        = Var("userInt('k2_lost_inner_hits')", int),
        k2_is_lost                = Var("userInt('k2_is_lost')", int),
        pi_from_pv                = Var("userInt('pi_from_pv')", int),
        pi_pv_assoc_q             = Var("userInt('pi_pv_assoc_q')", int),
        pi_high_purity            = Var("userInt('pi_high_purity')", int),
        pi_lost_inner_hits        = Var("userInt('pi_lost_inner_hits')", int),
        pi_is_lost                = Var("userInt('pi_is_lost')", int),

        sv_x         = Var("userFloat('sv_x')",float),
        sv_y         = Var("userFloat('sv_y')",float),
        sv_z         = Var("userFloat('sv_z')",float),
        sv_chi2      = Var("userFloat('sv_chi2')",float),        
        sv_ndof      = Var("userFloat('sv_ndof')",float),        
        sv_redchi2   = Var("userFloat('sv_redchi2')",float),        
        sv_prob      = Var("userFloat('sv_prob')",float),        

        tv_x         = Var("userFloat('tv_x')",float),
        tv_y         = Var("userFloat('tv_y')",float),
        tv_z         = Var("userFloat('tv_z')",float),
        tv_chi2      = Var("userFloat('tv_chi2')",float),        
        tv_ndof      = Var("userFloat('tv_ndof')",float),        
        tv_redchi2   = Var("userFloat('tv_redchi2')",float),        
        tv_prob      = Var("userFloat('tv_prob')",float),        

        fv_x         = Var("userFloat('fv_x')",float),
        fv_y         = Var("userFloat('fv_y')",float),
        fv_z         = Var("userFloat('fv_z')",float),
        fv_chi2      = Var("userFloat('fv_chi2')",float),        
        fv_ndof      = Var("userFloat('fv_ndof')",float),        
        fv_redchi2   = Var("userFloat('fv_redchi2')",float),        
        fv_prob      = Var("userFloat('fv_prob')",float),        


        lxy_bs             = Var("userFloat('lxy_bs')",float),        
        lxy_bs_err         = Var("userFloat('lxy_bs_err')",float),        
        lxy_bs_sig         = Var("userFloat('lxy_bs_sig')",float),        

        lxyz_bs            = Var("userFloat('lxyz_bs')",float),        
        lxyz_bs_err        = Var("userFloat('lxyz_bs_err')",float),        
        lxyz_bs_sig        = Var("userFloat('lxyz_bs_sig')",float),        

        lxy_ds             = Var("userFloat('lxy_ds')",float),        
        lxy_ds_err         = Var("userFloat('lxy_ds_err')",float),        
        lxy_ds_sig         = Var("userFloat('lxy_ds_sig')",float),        

        lxyz_ds            = Var("userFloat('lxyz_ds')",float),        
        lxyz_ds_err        = Var("userFloat('lxyz_ds_err')",float),        
        lxyz_ds_sig        = Var("userFloat('lxyz_ds_sig')",float),        

        lxy_phi            = Var("userFloat('lxy_phi')",float),        
        lxy_phi_err        = Var("userFloat('lxy_phi_err')",float),        
        lxy_phi_sig        = Var("userFloat('lxy_phi_sig')",float),        

        lxyz_phi           = Var("userFloat('lxyz_phi')",float),        
        lxyz_phi_err       = Var("userFloat('lxyz_phi_err')",float),        
        lxyz_phi_sig       = Var("userFloat('lxyz_phi_sig')",float),        

        signed_decay_ip3d_mu_ds_sv      = Var("userFloat('signed_decay_ip3d_mu_ds_sv')",float),
        signed_decay_ip3d_mu_ds_err_sv  = Var("userFloat('signed_decay_ip3d_mu_ds_err_sv')",float),
        signed_decay_ip3d_mu_ds_sig_sv  = Var("userFloat('signed_decay_ip3d_mu_ds_sig_sv')",float),

        signed_decay_ip3d_mu_bs_sv      = Var("userFloat('signed_decay_ip3d_mu_bs_sv')",float),
        signed_decay_ip3d_mu_bs_err_sv  = Var("userFloat('signed_decay_ip3d_mu_bs_err_sv')",float),
        signed_decay_ip3d_mu_bs_sig_sv  = Var("userFloat('signed_decay_ip3d_mu_bs_sig_sv')",float),

        signed_impact_ip3d_mu_ds_sv      = Var("userFloat('signed_impact_ip3d_mu_ds_sv')",float),
        signed_impact_ip3d_mu_ds_err_sv  = Var("userFloat('signed_impact_ip3d_mu_ds_err_sv')",float),
        signed_impact_ip3d_mu_ds_sig_sv  = Var("userFloat('signed_impact_ip3d_mu_ds_sig_sv')",float),

        signed_impact_ip3d_mu_bs_sv      = Var("userFloat('signed_impact_ip3d_mu_bs_sv')",float),
        signed_impact_ip3d_mu_bs_err_sv  = Var("userFloat('signed_impact_ip3d_mu_bs_err_sv')",float),
        signed_impact_ip3d_mu_bs_sig_sv  = Var("userFloat('signed_impact_ip3d_mu_bs_sig_sv')",float),

        btv_ip3d_mu_ds_sv      = Var("userFloat('btv_ip3d_mu_ds_sv')",float),
        btv_ip3d_mu_ds_err_sv  = Var("userFloat('btv_ip3d_mu_ds_err_sv')",float),
        btv_ip3d_mu_ds_sig_sv  = Var("userFloat('btv_ip3d_mu_ds_sig_sv')",float),

        btv_ip3d_mu_bs_sv      = Var("userFloat('btv_ip3d_mu_bs_sv')",float),
        btv_ip3d_mu_bs_err_sv  = Var("userFloat('btv_ip3d_mu_bs_err_sv')",float),
        btv_ip3d_mu_bs_sig_sv  = Var("userFloat('btv_ip3d_mu_bs_sig_sv')",float),

        dxy_mu_pv      = Var("userFloat('dxy_mu_pv')",float),
        dxy_mu_err_pv  = Var("userFloat('dxy_mu_err_pv')",float),
        dxy_mu_sig_pv  = Var("userFloat('dxy_mu_sig_pv')",float),

        dxy_mu_pv_dz      = Var("userFloat('dxy_mu_pv_dz')",float),
        dxy_mu_err_pv_dz  = Var("userFloat('dxy_mu_err_pv_dz')",float),
        dxy_mu_sig_pv_dz  = Var("userFloat('dxy_mu_sig_pv_dz')",float),

        dxy_mu_pv_dz_w_bs      = Var("userFloat('dxy_mu_pv_dz_w_bs')",float),
        dxy_mu_err_pv_dz_w_bs  = Var("userFloat('dxy_mu_err_pv_dz_w_bs')",float),
        dxy_mu_sig_pv_dz_w_bs  = Var("userFloat('dxy_mu_sig_pv_dz_w_bs')",float),


        dxy_mu_pv_general      = Var("userFloat('dxy_mu_pv_general')",float),
        dxy_mu_err_pv_general  = Var("userFloat('dxy_mu_err_pv_general')",float),
        dxy_mu_sig_pv_general  = Var("userFloat('dxy_mu_sig_pv_general')",float),

        dxy_mu_beamSpot      = Var("userFloat('dxy_mu_beamSpot')",float),
        dxy_mu_err_beamSpot  = Var("userFloat('dxy_mu_err_beamSpot')",float),
        dxy_mu_sig_beamSpot  = Var("userFloat('dxy_mu_sig_beamSpot')",float),


        dxy_mu_sv      = Var("userFloat('dxy_mu_sv')",float),
        dxy_mu_err_sv  = Var("userFloat('dxy_mu_err_sv')",float),
        dxy_mu_sig_sv  = Var("userFloat('dxy_mu_sig_sv')",float),

        dz_mu_pv       = Var("userFloat('dz_mu_pv')",float),
        dz_mu_sig_pv   = Var("userFloat('dz_mu_sig_pv')",float),
        dz_mu_err_pv   = Var("userFloat('dz_mu_err_pv')",float),

        dz_mu_sv       = Var("userFloat('dz_mu_sv')",float),
        dz_mu_sig_sv   = Var("userFloat('dz_mu_sig_sv')",float),
        dz_mu_err_sv   = Var("userFloat('dz_mu_err_sv')",float),

        dxy_pi_pv      = Var("userFloat('dxy_pi_pv')",float),
        dxy_pi_err_pv  = Var("userFloat('dxy_pi_err_pv')",float),
        dxy_pi_sig_pv  = Var("userFloat('dxy_pi_sig_pv')",float),

        dz_pi_pv       = Var("userFloat('dz_pi_pv')",float),
        dz_pi_err_pv   = Var("userFloat('dz_pi_err_pv')",float),
        dz_pi_sig_pv   = Var("userFloat('dz_pi_sig_pv')",float),

        dxy_pi_sv      = Var("userFloat('dxy_pi_sv')",float),
        dxy_pi_err_sv  = Var("userFloat('dxy_pi_err_sv')",float),
        dxy_pi_sig_sv  = Var("userFloat('dxy_pi_sig_sv')",float),

        dz_pi_sv       = Var("userFloat('dz_pi_sv')",float),
        dz_pi_err_sv   = Var("userFloat('dz_pi_err_sv')",float),
        dz_pi_sig_sv   = Var("userFloat('dz_pi_sig_sv')",float),

        ds_vtx_cosine_xy     = Var("userFloat('ds_vtx_cosine_xy')"    ,float),
        ds_vtx_cosine_xy_pv  = Var("userFloat('ds_vtx_cosine_xy_pv')" ,float),
        ds_vtx_cosine_xyz    = Var("userFloat('ds_vtx_cosine_xyz')"   ,float),
        ds_vtx_cosine_xyz_pv = Var("userFloat('ds_vtx_cosine_xyz_pv')",float),

)

##################################################

postfitBasicVariables = cms.PSet(

        # Phi fit
        phi_fitted_px            = Var("userFloat('phi_fitted_px')",float),
        phi_fitted_py            = Var("userFloat('phi_fitted_py')",float),
        phi_fitted_pz            = Var("userFloat('phi_fitted_pz')",float),
        phi_fitted_pt            = Var("userFloat('phi_fitted_pt')",float),
        phi_fitted_eta           = Var("userFloat('phi_fitted_eta')",float),
        phi_fitted_phi           = Var("userFloat('phi_fitted_phi')",float),
        phi_fitted_m             = Var("userFloat('phi_fitted_m')",float),
        phi_fitted_m_err         = Var("userFloat('phi_fitted_m_err')",float),

        # Ds fit
        ds_fitted_px             = Var("userFloat('ds_fitted_px')",float),
        ds_fitted_py             = Var("userFloat('ds_fitted_py')",float),
        ds_fitted_pz             = Var("userFloat('ds_fitted_pz')",float),
        ds_fitted_pt             = Var("userFloat('ds_fitted_pt')",float),
        ds_fitted_eta            = Var("userFloat('ds_fitted_eta')",float),
        ds_fitted_phi            = Var("userFloat('ds_fitted_phi')",float),
        ds_fitted_m              = Var("userFloat('ds_fitted_m')",float),
        ds_fitted_m_err          = Var("userFloat('ds_fitted_m_err')",float),
        ds_custom_m_err          = Var("userFloat('ds_custom_m_err')",float),

        # Bs fit
        bs_fitted_px             = Var("userFloat('bs_fitted_px')",float),
        bs_fitted_py             = Var("userFloat('bs_fitted_py')",float),
        bs_fitted_pz             = Var("userFloat('bs_fitted_pz')",float),
        bs_fitted_pt             = Var("userFloat('bs_fitted_pt')",float),
        bs_fitted_eta            = Var("userFloat('bs_fitted_eta')",float),
        bs_fitted_phi            = Var("userFloat('bs_fitted_phi')",float),
        bs_fitted_m              = Var("userFloat('bs_fitted_m')",float),
        bs_fitted_m_err          = Var("userFloat('bs_fitted_m_err')",float),

        iso_03                   = Var("userFloat('iso_03')", float),
        iso_04                   = Var("userFloat('iso_04')", float),
        rel_iso_03               = Var("userFloat('rel_iso_03')", float),
        rel_iso_04               = Var("userFloat('rel_iso_04')", float),

        iso_03_pv                = Var("userFloat('iso_03_pv')", float),
        iso_04_pv                = Var("userFloat('iso_04_pv')", float),
        rel_iso_03_pv            = Var("userFloat('rel_iso_03_pv')", float),
        rel_iso_04_pv            = Var("userFloat('rel_iso_04_pv')", float),

        iso_03_sv                = Var("userFloat('iso_03_sv')", float),
        iso_04_sv                = Var("userFloat('iso_04_sv')", float),
        rel_iso_03_sv            = Var("userFloat('rel_iso_03_sv')", float),
        rel_iso_04_sv            = Var("userFloat('rel_iso_04_sv')", float),

        iso_03_tv                = Var("userFloat('iso_03_tv')", float),
        iso_04_tv                = Var("userFloat('iso_04_tv')", float),
        rel_iso_03_tv            = Var("userFloat('rel_iso_03_tv')", float),
        rel_iso_04_tv            = Var("userFloat('rel_iso_04_tv')", float),

        iso_03_ds_sv               = Var("userFloat('iso_03_ds_sv')", float),
        iso_04_ds_sv               = Var("userFloat('iso_04_ds_sv')", float),
        rel_iso_03_ds_sv           = Var("userFloat('rel_iso_03_ds_sv')", float),
        rel_iso_04_ds_sv           = Var("userFloat('rel_iso_04_ds_sv')", float),

        e_gamma                  = Var("userFloat('e_gamma')",float),
        kappa                    = Var("userFloat('kappa')",float),
)

##################################################

helicityVariables = cms.PSet(

        cosPiK1             = Var("userFloat('cosPiK1')",float),
        cosPiK2             = Var("userFloat('cosPiK2')",float),

        cosPhiDs_coll       = Var("userFloat('cosPhiDsColl')",float),
        cosPhiDs_lhcb_alt   = Var("userFloat('cosPhiDsLhcbAlt')",float),
        cosPhiDs_reco_1      = Var("userFloat('cosPhiDsReco1')",float),
        cosPhiDs_reco_2      = Var("userFloat('cosPhiDsReco2')",float),

        cosPiDs_coll        = Var("userFloat('cosPiDsColl')",float),
        cosPiDs_lhcb_alt    = Var("userFloat('cosPiDsLhcbAlt')",float),
        cosPiDs_reco_1       = Var("userFloat('cosPiDsReco1')",float),
        cosPiDs_reco_2       = Var("userFloat('cosPiDsReco2')",float),

        cosMuW_coll         = Var("userFloat('cosMuWColl')",float),
        cosMuW_lhcb_alt     = Var("userFloat('cosMuWLhcbAlt')",float),
        cosMuW_reco_1        = Var("userFloat('cosMuWReco1')",float),
        cosMuW_reco_2        = Var("userFloat('cosMuWReco2')",float),

        cosPlaneBs_coll     = Var("userFloat('cosPlaneBsColl')",float),
        cosPlaneBs_lhcb_alt = Var("userFloat('cosPlaneBsLhcbAlt')",float),
        cosPlaneBs_reco_1    = Var("userFloat('cosPlaneBsReco1')",float),
        cosPlaneBs_reco_2    = Var("userFloat('cosPlaneBsReco2')",float),

        cosPlaneDs_coll     = Var("userFloat('cosPlaneDsColl')",float),
        cosPlaneDs_lhcb_alt = Var("userFloat('cosPlaneDsLhcbAlt')",float),
        cosPlaneDs_reco_1    = Var("userFloat('cosPlaneDsReco1')",float),
        cosPlaneDs_reco_2    = Var("userFloat('cosPlaneDsReco2')",float),


)

##################################################

bsMomentumVariables = cms.PSet(

        bs_mass_corr       = Var("userFloat('bs_m_corr')",float), 
        ds_perp            = Var("userFloat('ds_perp')",float), 
        ds_mu_perp         = Var("userFloat('ds_mu_perp')",float), 

        #pi_mu_charge = Var("userInt('pi_mu_charge')",int),
        bs_px_coll         = Var("userFloat('bs_px_coll')",float),
        bs_py_coll         = Var("userFloat('bs_py_coll')",float),
        bs_pz_coll         = Var("userFloat('bs_pz_coll')",float),
        bs_pt_coll         = Var("userFloat('bs_pt_coll')",float),
        bs_eta_coll        = Var("userFloat('bs_eta_coll')",float),
        bs_phi_coll        = Var("userFloat('bs_phi_coll')",float),
        bs_boost_coll      = Var("userFloat('b_boost_coll')",float),
        bs_boost_x_coll   = Var("userFloat('b_boost_coll_x')",float),
        bs_boost_y_coll  = Var("userFloat('b_boost_coll_y')",float),
        bs_boost_z_coll  = Var("userFloat('b_boost_coll_z')",float),


        bs_px_lhcb_alt        = Var("userFloat('bs_px_lhcb_alt')",float),
        bs_py_lhcb_alt        = Var("userFloat('bs_py_lhcb_alt')",float),
        bs_pz_lhcb_alt        = Var("userFloat('bs_pz_lhcb_alt')",float),
        bs_pt_lhcb_alt        = Var("userFloat('bs_pt_lhcb_alt')",float),
        bs_eta_lhcb_alt       = Var("userFloat('bs_eta_lhcb_alt')",float),
        bs_phi_lhcb_alt       = Var("userFloat('bs_phi_lhcb_alt')",float),
        bs_boost_lhcb_alt     = Var("userFloat('b_boost_lhcb_alt')",float),
        bs_boost_x_lhcb_alt  = Var("userFloat('b_boost_lhcb_alt_x')",float),
        bs_boost_y_lhcb_alt = Var("userFloat('b_boost_lhcb_alt_y')",float),
        bs_boost_z_lhcb_alt = Var("userFloat('b_boost_lhcb_alt_z')",float),

        disc_is_negative    = Var("userInt('disc_is_negative')",int),
        disc_negativity     = Var("userFloat('disc_negativity')",float),

        bs_px_reco_1        = Var("userFloat('bs_px_reco_1')",float),
        bs_py_reco_1        = Var("userFloat('bs_py_reco_1')",float),
        bs_pz_reco_1        = Var("userFloat('bs_pz_reco_1')",float),
        bs_pt_reco_1        = Var("userFloat('bs_pt_reco_1')",float),
        bs_eta_reco_1       = Var("userFloat('bs_eta_reco_1')",float),
        bs_phi_reco_1       = Var("userFloat('bs_phi_reco_1')",float),
        bs_boost_reco_1     = Var("userFloat('b_boost_reco_1')",float),
        bs_boost_x_reco_1  = Var("userFloat('b_boost_reco_1_x')",float),
        bs_boost_y_reco_1 = Var("userFloat('b_boost_reco_1_y')",float),
        bs_boost_z_reco_1 = Var("userFloat('b_boost_reco_1_z')",float),

        bs_px_reco_2         = Var("userFloat('bs_px_reco_2')",float),
        bs_py_reco_2         = Var("userFloat('bs_py_reco_2')",float),
        bs_pz_reco_2         = Var("userFloat('bs_pz_reco_2')",float),
        bs_pt_reco_2         = Var("userFloat('bs_pt_reco_2')",float),
        bs_eta_reco_2        = Var("userFloat('bs_eta_reco_2')",float),
        bs_phi_reco_2        = Var("userFloat('bs_phi_reco_2')",float),
        bs_boost_reco_2      = Var("userFloat('b_boost_reco_2')",float),
        bs_boost_x_reco_2   = Var("userFloat('b_boost_reco_2_x')",float),
        bs_boost_y_reco_2  = Var("userFloat('b_boost_reco_2_y')",float),
        bs_boost_z_reco_2  = Var("userFloat('b_boost_reco_2_z')",float),

        #lhcb_pz = Var("userFloat('lhcb_pz')",float),
        #theta = Var("userFloat('theta')",float), 

        m2_miss_coll     = Var("userFloat('m2_miss_coll')",float),
        m2_miss_lhcb_alt = Var("userFloat('m2_miss_lhcb_alt')",float),
        m2_miss_reco_1   = Var("userFloat('m2_miss_reco_1')",float),
        m2_miss_reco_2   = Var("userFloat('m2_miss_reco_2')",float),

        pt_miss_coll     = Var("userFloat('pt_miss_coll')",float),
        pt_miss_lhcb_alt = Var("userFloat('pt_miss_lhcb_alt')",float),
        pt_miss_reco_1   = Var("userFloat('pt_miss_reco_1')",float),
        pt_miss_reco_2   = Var("userFloat('pt_miss_reco_2')",float),

        q2_coll          = Var("userFloat('q2_coll')",float),
        q2_lhcb_alt      = Var("userFloat('q2_lhcb_alt')",float),
        q2_reco_1        = Var("userFloat('q2_reco_1')",float),
        q2_reco_2        = Var("userFloat('q2_reco_2')",float),

        e_star_coll          = Var("userFloat('e_star_coll')",float),
        e_star_lhcb_alt      = Var("userFloat('e_star_lhcb_alt')",float),
        e_star_reco_1        = Var("userFloat('e_star_reco_1')",float),
        e_star_reco_2        = Var("userFloat('e_star_reco_2')",float),

        #arrived          = Var("userInt('arrived')",int),
)

photonVariables = cms.PSet(

        photon_pt                    = Var("userFloat('photon_pt')",float),
        photon_eta                   = Var("userFloat('photon_eta')",float),
        photon_phi                   = Var("userFloat('photon_phi')",float),
        photon_pdgid                 = Var("userInt('photon_pdgid')",int),

        dr_photon_ds                 = Var("userFloat('dr_photon_ds')",float),
        dsPhoton_m                   = Var("userFloat('dsPhoton_m')",float),
        dsPhotonMu_m                 = Var("userFloat('dsPhotonMu_m')",float),
        found_photon                 = Var("userInt('found_photon')",int),
        photon_multiplicity          = Var("userInt('photon_multiplicity')",int),

        bs_mass_corr_photon          = Var("userFloat('bs_m_corr_photon')",float), 
        ds_perp_photon               = Var("userFloat('ds_perp_photon')",float), 
        ds_mu_perp_photon            = Var("userFloat('ds_mu_perp_photon')",float), 

        bs_px_coll_photon            = Var("userFloat('bs_px_coll_photon')",float),
        bs_py_coll_photon            = Var("userFloat('bs_py_coll_photon')",float),
        bs_pz_coll_photon            = Var("userFloat('bs_pz_coll_photon')",float),
        bs_pt_coll_photon            = Var("userFloat('bs_pt_coll_photon')",float),
        bs_eta_coll_photon           = Var("userFloat('bs_eta_coll_photon')",float),
        bs_phi_coll_photon           = Var("userFloat('bs_phi_coll_photon')",float),
        bs_boost_coll_photon         = Var("userFloat('b_boost_coll_photon')",float),
        bs_boost_x_coll_photon      = Var("userFloat('b_boost_coll_x_photon')",float),
        bs_boost_y_coll_photon     = Var("userFloat('b_boost_coll_y_photon')",float),
        bs_boost_z_coll_photon     = Var("userFloat('b_boost_coll_z_photon')",float),

        bs_px_lhcb_alt_photon        = Var("userFloat('bs_px_lhcb_alt_photon')",float),
        bs_py_lhcb_alt_photon        = Var("userFloat('bs_py_lhcb_alt_photon')",float),
        bs_pz_lhcb_alt_photon        = Var("userFloat('bs_pz_lhcb_alt_photon')",float),
        bs_pt_lhcb_alt_photon        = Var("userFloat('bs_pt_lhcb_alt_photon')",float),
        bs_eta_lhcb_alt_photon       = Var("userFloat('bs_eta_lhcb_alt_photon')",float),
        bs_phi_lhcb_alt_photon       = Var("userFloat('bs_phi_lhcb_alt_photon')",float),
        bs_boost_lhcb_alt_photon     = Var("userFloat('b_boost_lhcb_alt_photon')",float),
        bs_boost_x_lhcb_alt_photon  = Var("userFloat('b_boost_lhcb_alt_x_photon')",float),
        bs_boost_y_lhcb_alt_photon = Var("userFloat('b_boost_lhcb_alt_y_photon')",float),
        bs_boost_z_lhcb_alt_photon = Var("userFloat('b_boost_lhcb_alt_z_photon')",float),

        disc_is_negative_photon      = Var("userInt('disc_is_negative_photon')",int),
        disc_negativity_photon       = Var("userFloat('disc_negativity_photon')",float),

        bs_px_reco_1_photon          = Var("userFloat('bs_px_reco_1_photon')",float),
        bs_py_reco_1_photon          = Var("userFloat('bs_py_reco_1_photon')",float),
        bs_pz_reco_1_photon          = Var("userFloat('bs_pz_reco_1_photon')",float),
        bs_pt_reco_1_photon          = Var("userFloat('bs_pt_reco_1_photon')",float),
        bs_eta_reco_1_photon         = Var("userFloat('bs_eta_reco_1_photon')",float),
        bs_phi_reco_1_photon         = Var("userFloat('bs_phi_reco_1_photon')",float),
        bs_boost_reco_1_photon       = Var("userFloat('b_boost_reco_1_photon')",float),
        bs_boost_x_reco_1_photon    = Var("userFloat('b_boost_reco_1_x_photon')",float),
        bs_boost_y_reco_1_photon   = Var("userFloat('b_boost_reco_1_y_photon')",float),
        bs_boost_z_reco_1_photon   = Var("userFloat('b_boost_reco_1_z_photon')",float),

        bs_px_reco_2_photon          = Var("userFloat('bs_px_reco_2_photon')",float),
        bs_py_reco_2_photon          = Var("userFloat('bs_py_reco_2_photon')",float),
        bs_pz_reco_2_photon          = Var("userFloat('bs_pz_reco_2_photon')",float),
        bs_pt_reco_2_photon          = Var("userFloat('bs_pt_reco_2_photon')",float),
        bs_eta_reco_2_photon         = Var("userFloat('bs_eta_reco_2_photon')",float),
        bs_phi_reco_2_photon         = Var("userFloat('bs_phi_reco_2_photon')",float),
        bs_boost_reco_2_photon       = Var("userFloat('b_boost_reco_2_photon')",float),
        bs_boost_x_reco_2_photon    = Var("userFloat('b_boost_reco_2_x_photon')",float),
        bs_boost_y_reco_2_photon   = Var("userFloat('b_boost_reco_2_y_photon')",float),
        bs_boost_z_reco_2_photon   = Var("userFloat('b_boost_reco_2_z_photon')",float),

        m2_miss_coll_photon          = Var("userFloat('m2_miss_coll_photon')",float),
        m2_miss_lhcb_alt_photon      = Var("userFloat('m2_miss_lhcb_alt_photon')",float),
        m2_miss_reco_1_photon        = Var("userFloat('m2_miss_reco_1_photon')",float),
        m2_miss_reco_2_photon        = Var("userFloat('m2_miss_reco_2_photon')",float),

        pt_miss_coll_photon          = Var("userFloat('pt_miss_coll_photon')",float),
        pt_miss_lhcb_alt_photon      = Var("userFloat('pt_miss_lhcb_alt_photon')",float),
        pt_miss_reco_1_photon        = Var("userFloat('pt_miss_reco_1_photon')",float),
        pt_miss_reco_2_photon        = Var("userFloat('pt_miss_reco_2_photon')",float),

        q2_coll_photon               = Var("userFloat('q2_coll_photon')",float),
        q2_lhcb_alt_photon           = Var("userFloat('q2_lhcb_alt_photon')",float),
        q2_reco_1_photon             = Var("userFloat('q2_reco_1_photon')",float),
        q2_reco_2_photon             = Var("userFloat('q2_reco_2_photon')",float),

        e_star_coll_photon           = Var("userFloat('e_star_coll_photon')",float),
        e_star_lhcb_alt_photon       = Var("userFloat('e_star_lhcb_alt_photon')",float),
        e_star_reco_1_photon         = Var("userFloat('e_star_reco_1_photon')",float),
        e_star_reco_2_photon         = Var("userFloat('e_star_reco_2_photon')",float),

        cosPhiDs_coll_photon         = Var("userFloat('cosPhiDsColl_photon')",float),
        cosPhiDs_lhcb_alt_photon     = Var("userFloat('cosPhiDsLhcbAlt_photon')",float),
        cosPhiDs_reco_1_photon       = Var("userFloat('cosPhiDsReco1_photon')",float),
        cosPhiDs_reco_2_photon       = Var("userFloat('cosPhiDsReco2_photon')",float),

        cosPiDs_coll_photon          = Var("userFloat('cosPiDsColl_photon')",float),
        cosPiDs_lhcb_alt_photon      = Var("userFloat('cosPiDsLhcbAlt_photon')",float),
        cosPiDs_reco_1_photon        = Var("userFloat('cosPiDsReco1_photon')",float),
        cosPiDs_reco_2_photon        = Var("userFloat('cosPiDsReco2_photon')",float),

        cosMuW_coll_photon           = Var("userFloat('cosMuWColl_photon')",float),
        cosMuW_lhcb_alt_photon       = Var("userFloat('cosMuWLhcbAlt_photon')",float),
        cosMuW_reco_1_photon         = Var("userFloat('cosMuWReco1_photon')",float),
        cosMuW_reco_2_photon         = Var("userFloat('cosMuWReco2_photon')",float),

        cosPlaneBs_coll_photon       = Var("userFloat('cosPlaneBsColl_photon')",float),
        cosPlaneBs_lhcb_alt_photon   = Var("userFloat('cosPlaneBsLhcbAlt_photon')",float),
        cosPlaneBs_reco_1_photon     = Var("userFloat('cosPlaneBsReco1_photon')",float),
        cosPlaneBs_reco_2_photon     = Var("userFloat('cosPlaneBsReco2_photon')",float),

        cosPlaneDs_coll_photon       = Var("userFloat('cosPlaneDsColl_photon')",float),
        cosPlaneDs_lhcb_alt_photon   = Var("userFloat('cosPlaneDsLhcbAlt_photon')",float),
        cosPlaneDs_reco_1_photon     = Var("userFloat('cosPlaneDsReco1_photon')",float),
        cosPlaneDs_reco_2_photon     = Var("userFloat('cosPlaneDsReco2_photon')",float),

        iso_03_photon                    = Var("userFloat('iso_03_photon')", float),
        iso_04_photon                    = Var("userFloat('iso_04_photon')", float),
        rel_iso_03_photon                = Var("userFloat('rel_iso_03_photon')", float),
        rel_iso_04_photon                = Var("userFloat('rel_iso_04_photon')", float),

        iso_03_pv_photon                 = Var("userFloat('iso_03_pv_photon')", float),
        iso_04_pv_photon                 = Var("userFloat('iso_04_pv_photon')", float),
        rel_iso_03_pv_photon             = Var("userFloat('rel_iso_03_pv_photon')", float),
        rel_iso_04_pv_photon             = Var("userFloat('rel_iso_04_pv_photon')", float),

        iso_03_sv_photon                 = Var("userFloat('iso_03_sv_photon')", float),
        iso_04_sv_photon                 = Var("userFloat('iso_04_sv_photon')", float),
        rel_iso_03_sv_photon             = Var("userFloat('rel_iso_03_sv_photon')", float),
        rel_iso_04_sv_photon             = Var("userFloat('rel_iso_04_sv_photon')", float),

        iso_03_tv_photon                 = Var("userFloat('iso_03_tv_photon')", float),
        iso_04_tv_photon                 = Var("userFloat('iso_04_tv_photon')", float),
        rel_iso_03_tv_photon             = Var("userFloat('rel_iso_03_tv_photon')", float),
        rel_iso_04_tv_photon             = Var("userFloat('rel_iso_04_tv_photon')", float),

        iso_03_ds_sv_photon              = Var("userFloat('iso_03_ds_sv_photon')", float),
        iso_04_ds_sv_photon              = Var("userFloat('iso_04_ds_sv_photon')", float),
        rel_iso_03_ds_sv_photon          = Var("userFloat('rel_iso_03_ds_sv_photon')", float),
        rel_iso_04_ds_sv_photon          = Var("userFloat('rel_iso_04_ds_sv_photon')", float),

)
##################################################

## this is for the ED Filter
#empty = cms.PSet(arrived = Var("userInt('arrived')",int))

