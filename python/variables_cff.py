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
 
prefitBasicVariables = cms.PSet(

        #save the indices of the final states in the pruned collection
        #mu_idx = uint("mu_idx"), #always 0 :) 
        k1_idx = uint("k1_idx"),
        k2_idx = uint("k2_idx"),
        pi_idx = uint("pi_idx"),

        # final state particle info, access via userCand()
        mu_px         = Var("userCand('mu').px()",float),
        mu_py         = Var("userCand('mu').py()",float),
        mu_pz         = Var("userCand('mu').pz()",float),
        mu_pt         = Var("userCand('mu').pt()",float),
        mu_eta        = Var("userCand('mu').eta()",float),
        mu_phi        = Var("userCand('mu').phi()",float),
        mu_m       = ufloat("mu_m"),
        mu_charge     = Var("userCand('mu').charge()",float),
        mu_pdgId      = Var("userCand('mu').pdgId()",int),

        mu_is_tracker    = Var("userInt('mu_is_tracker')",int),
        mu_is_pf         = Var("userInt('mu_is_pf')",int),
        mu_is_global     = Var("userInt('mu_is_global')",int),
        mu_is_standalone = Var("userInt('mu_is_standalone')",int),

        mu_id_loose      = Var("userInt('mu_id_loose')",int),
        mu_id_medium     = Var("userInt('mu_id_medium')",int),
        mu_id_tight      = Var("userInt('mu_id_tight')",int),


        k1_px         = Var("userCand('k1').px()",float),
        k1_py         = Var("userCand('k1').py()",float),
        k1_pz         = Var("userCand('k1').pz()",float),
        k1_pt         = Var("userCand('k1').pt()",float),
        k1_eta        = Var("userCand('k1').eta()",float),
        k1_phi        = Var("userCand('k1').phi()",float),
        k1_m       = ufloat("k1_m"),
        k1_charge     = Var("userCand('k1').charge()",float),
        k1_pdgId      = Var("userCand('k1').pdgId()",int),

        k2_px         = Var("userCand('k2').px()",float),
        k2_py         = Var("userCand('k2').py()",float),
        k2_pz         = Var("userCand('k2').pz()",float),
        k2_pt         = Var("userCand('k2').pt()",float),
        k2_eta        = Var("userCand('k2').eta()",float),
        k2_phi        = Var("userCand('k2').phi()",float),
        k2_m       = ufloat("k2_m"),
        k2_charge     = Var("userCand('k2').charge()",float),
        k2_pdgId      = Var("userCand('k2').pdgId()",int),

        pi_px         = Var("userCand('pi').px()",float),
        pi_py         = Var("userCand('pi').py()",float),
        pi_pz         = Var("userCand('pi').pz()",float),
        pi_pt         = Var("userCand('pi').pt()",float),
        pi_eta        = Var("userCand('pi').eta()",float),
        pi_phi        = Var("userCand('pi').phi()",float),
        pi_m       = ufloat("pi_m"),
        pi_charge     = Var("userCand('pi').charge()",float),
        pi_pdgId      = Var("userCand('pi').pdgId()",int),

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
        k_k_charge    = Var("userFloat('k_k_charge')",float),

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

        bs_lhcb_pt    = Var("userFloat('bs_gen_lhcb_pt')",float),
        bs_lhcb_eta   = Var("userFloat('bs_gen_lhcb_eta')",float),
        bs_lhcb_phi   = Var("userFloat('bs_gen_lhcb_phi')",float),

        fv_x          = Var("userFloat('fv_x_gen')",float),
        fv_y          = Var("userFloat('fv_y_gen')",float),
        fv_z          = Var("userFloat('fv_z_gen')",float),

        m2_miss       = ufloat('m2_miss_gen'),
        pt_miss       = ufloat('pt_miss_gen'),
        q2            = ufloat('q2_gen'),
        e_star            = ufloat('e_star_gen'),

        angMuW        = Var("userFloat('angMuWGen')",float),
        cosMuW        = Var("userFloat('cosMuWGen')",float),
        cosMuWLhcb    = Var("userFloat('cosMuWGenLhcb')",float),
        cosMuWReco1   = Var("userFloat('cosMuWGenReco1')",float),
        cosMuWReco2   = Var("userFloat('cosMuWGenReco2')",float),

        angPiK1       = Var("userFloat('angPiK1Gen')",float),
        cosPiK1       = Var("userFloat('cosPiK1Gen')",float),

        angPiK2       = Var("userFloat('angPiK2Gen')",float),
        cosPiK2       = Var("userFloat('cosPiK2Gen')",float),

        angPiDs       = Var("userFloat('angPiDsGen')",float),
        cosPiDs       = Var("userFloat('cosPiDsGen')",float),
        cosPiDsLhcb   = Var("userFloat('cosPiDsGenLhcb')",float),

        angPhiDs      = Var("userFloat('angPhiDsGen')",float),
        cosPhiDs      = Var("userFloat('cosPhiDsGen')",float),

        angPlaneBs    = Var("userFloat('angPlaneBsGen')",float),
        cosPlaneBs    = Var("userFloat('cosPlaneBsGen')",float),
 
        angPlaneDs    = Var("userFloat('angPlaneDsGen')",float),
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

)

##################################################

vertexVariables = cms.PSet(

        bs_x0 = Var("userFloat('bs_x0')",float),        
        bs_y0 = Var("userFloat('bs_y0')",float),        
        bs_z0 = Var("userFloat('bs_z0')",float),        

        bs_x_wrt_pv = Var("userFloat('bs_x_wrt_pv_z')",float),        
        bs_y_wrt_pv = Var("userFloat('bs_y_wrt_pv_z')",float),        

        pv_x = Var("userFloat('pv_x')",float),        
        pv_y = Var("userFloat('pv_y')",float),        
        pv_z = Var("userFloat('pv_z')",float),        
        pv_chi2    = Var("userFloat('pv_chi2')",float),        
        pv_ndof    = Var("userFloat('pv_ndof')",float),        
        pv_redchi2 = Var("userFloat('pv_redchi2')",float),        
        pv_prob    = Var("userFloat('pv_prob')",float),        
        pv_idx     = Var("userInt('pv_idx')",int),        

        sv_x = Var("userFloat('sv_x')",float),
        sv_y = Var("userFloat('sv_y')",float),
        sv_z = Var("userFloat('sv_z')",float),
        sv_chi2    = Var("userFloat('sv_chi2')",float),        
        sv_ndof    = Var("userFloat('sv_ndof')",float),        
        sv_redchi2 = Var("userFloat('sv_redchi2')",float),        
        sv_prob    = Var("userFloat('sv_prob')",float),        

        tv_x = Var("userFloat('tv_x')",float),
        tv_y = Var("userFloat('tv_y')",float),
        tv_z = Var("userFloat('tv_z')",float),
        tv_chi2    = Var("userFloat('tv_chi2')",float),        
        tv_ndof    = Var("userFloat('tv_ndof')",float),        
        tv_redchi2 = Var("userFloat('tv_redchi2')",float),        
        tv_prob    = Var("userFloat('tv_prob')",float),        

        fv_x = Var("userFloat('fv_x')",float),
        fv_y = Var("userFloat('fv_y')",float),
        fv_z = Var("userFloat('fv_z')",float),
        fv_chi2    = Var("userFloat('fv_chi2')",float),        
        fv_ndof    = Var("userFloat('fv_ndof')",float),        
        fv_redchi2 = Var("userFloat('fv_redchi2')",float),        
        fv_prob    = Var("userFloat('fv_prob')",float),        

        ##easy_bs_vtx_x = Var("userFloat('easy_bs_vtx_x')",float),
        ##easy_bs_vtx_y = Var("userFloat('easy_bs_vtx_y')",float),
        ##easy_bs_vtx_z = Var("userFloat('easy_bs_vtx_z')",float),

        lxy_bs         = Var("userFloat('lxy_bs')",float),        
        lxy_bs_err     = Var("userFloat('lxy_bs_err')",float),        
        lxy_bs_sig     = Var("userFloat('lxy_bs_sig')",float),        

        lxyz_bs        = Var("userFloat('lxyz_bs')",float),        
        lxyz_bs_err    = Var("userFloat('lxyz_bs_err')",float),        
        lxyz_bs_sig    = Var("userFloat('lxyz_bs_sig')",float),        

        lxy_ds         = Var("userFloat('lxy_ds')",float),        
        lxy_ds_err     = Var("userFloat('lxy_ds_err')",float),        
        lxy_ds_sig     = Var("userFloat('lxy_ds_sig')",float),        

        lxyz_ds        = Var("userFloat('lxyz_ds')",float),        
        lxyz_ds_err    = Var("userFloat('lxyz_ds_err')",float),        
        lxyz_ds_sig    = Var("userFloat('lxyz_ds_sig')",float),        

        lxy_phi        = Var("userFloat('lxy_phi')",float),        
        lxy_phi_err    = Var("userFloat('lxy_phi_err')",float),        
        lxy_phi_sig    = Var("userFloat('lxy_phi_sig')",float),        

        lxyz_phi       = Var("userFloat('lxyz_phi')",float),        
        lxyz_phi_err   = Var("userFloat('lxyz_phi_err')",float),        
        lxyz_phi_sig   = Var("userFloat('lxyz_phi_sig')",float),        

        dxy_mu     = Var("userFloat('dxy_mu')",float),
        dz_mu      = Var("userFloat('dz_mu')",float),
        dxy_mu_err = Var("userFloat('dxy_mu_err')",float),
        dz_mu_err  = Var("userFloat('dz_mu_err')",float),
        dxy_mu_sig = Var("userFloat('dxy_mu_sig')",float),
        dz_mu_sig  = Var("userFloat('dz_mu_sig')",float),

        dxy_pi     = Var("userFloat('dxy_pi')",float),
        dz_pi      = Var("userFloat('dz_pi')",float),
        dxy_pi_err = Var("userFloat('dxy_pi_err')",float),
        dz_pi_err  = Var("userFloat('dz_pi_err')",float),
        dxy_pi_sig = Var("userFloat('dxy_pi_sig')",float),
        dz_pi_sig  = Var("userFloat('dz_pi_sig')",float),

        dxy_k1     = Var("userFloat('dxy_k1')",float),
        dz_k1      = Var("userFloat('dz_k1')",float),
        dxy_k1_err = Var("userFloat('dxy_k1_err')",float),
        dz_k1_err  = Var("userFloat('dz_k1_err')",float),
        dxy_k1_sig = Var("userFloat('dxy_k1_sig')",float),
        dz_k1_sig  = Var("userFloat('dz_k1_sig')",float),

        dxy_k2     = Var("userFloat('dxy_k2')",float),
        dz_k2      = Var("userFloat('dz_k2')",float),
        dxy_k2_err = Var("userFloat('dxy_k2_err')",float),
        dz_k2_err  = Var("userFloat('dz_k2_err')",float),
        dxy_k2_sig = Var("userFloat('dxy_k2_sig')",float),
        dz_k2_sig  = Var("userFloat('dz_k2_sig')",float),

        ds_vtx_cosine = Var("userFloat('ds_vtx_cosine')",float),

)

##################################################

postfitBasicVariables = cms.PSet(

        # Phi fit
        k1_refitted_px   = Var("userFloat('k1_refitted_px')",float),
        k1_refitted_py   = Var("userFloat('k1_refitted_py')",float),
        k1_refitted_pz   = Var("userFloat('k1_refitted_pz')",float),
        k1_refitted_pt   = Var("userFloat('k1_refitted_pt')",float),
        k1_refitted_eta  = Var("userFloat('k1_refitted_eta')",float),
        k1_refitted_phi  = Var("userFloat('k1_refitted_phi')",float),
        k1_refitted_m    = Var("userFloat('k1_refitted_m')",float),

        k2_refitted_px   = Var("userFloat('k2_refitted_px')",float),
        k2_refitted_py   = Var("userFloat('k2_refitted_py')",float),
        k2_refitted_pz   = Var("userFloat('k2_refitted_pz')",float),
        k2_refitted_pt   = Var("userFloat('k2_refitted_pt')",float),
        k2_refitted_eta  = Var("userFloat('k2_refitted_eta')",float),
        k2_refitted_phi  = Var("userFloat('k2_refitted_phi')",float),
        k2_refitted_m    = Var("userFloat('k2_refitted_m')",float),

        phi_fitted_px    = Var("userFloat('phi_fitted_px')",float),
        phi_fitted_py    = Var("userFloat('phi_fitted_py')",float),
        phi_fitted_pz    = Var("userFloat('phi_fitted_pz')",float),
        phi_fitted_pt    = Var("userFloat('phi_fitted_pt')",float),
        phi_fitted_eta   = Var("userFloat('phi_fitted_eta')",float),
        phi_fitted_phi   = Var("userFloat('phi_fitted_phi')",float),
        phi_fitted_m     = Var("userFloat('phi_fitted_m')",float),

        # Ds fit
        phi_refitted_px  = Var("userFloat('phi_refitted_px')",float),
        phi_refitted_py  = Var("userFloat('phi_refitted_py')",float),
        phi_refitted_pz  = Var("userFloat('phi_refitted_pz')",float),
        phi_refitted_pt  = Var("userFloat('phi_refitted_pt')",float),
        phi_refitted_eta = Var("userFloat('phi_refitted_eta')",float),
        phi_refitted_phi = Var("userFloat('phi_refitted_phi')",float),
        phi_refitted_m   = Var("userFloat('phi_refitted_m')",float),

        pi_refitted_px   = Var("userFloat('pi_refitted_px')",float),
        pi_refitted_py   = Var("userFloat('pi_refitted_py')",float),
        pi_refitted_pz   = Var("userFloat('pi_refitted_pz')",float),
        pi_refitted_pt   = Var("userFloat('pi_refitted_pt')",float),
        pi_refitted_eta  = Var("userFloat('pi_refitted_eta')",float),
        pi_refitted_phi  = Var("userFloat('pi_refitted_phi')",float),
        pi_refitted_m    = Var("userFloat('pi_refitted_m')",float),

        ds_fitted_px     = Var("userFloat('ds_fitted_px')",float),
        ds_fitted_py     = Var("userFloat('ds_fitted_py')",float),
        ds_fitted_pz     = Var("userFloat('ds_fitted_pz')",float),
        ds_fitted_pt     = Var("userFloat('ds_fitted_pt')",float),
        ds_fitted_eta    = Var("userFloat('ds_fitted_eta')",float),
        ds_fitted_phi    = Var("userFloat('ds_fitted_phi')",float),
        ds_fitted_m      = Var("userFloat('ds_fitted_m')",float),
        ds_fitted_boost  = Var("userFloat('ds_fitted_boost')",float),

        phiPi_refitted_px     = Var("userFloat('phiPi_refitted_px')",float),
        phiPi_refitted_py     = Var("userFloat('phiPi_refitted_py')",float),
        phiPi_refitted_pz     = Var("userFloat('phiPi_refitted_pz')",float),
        phiPi_refitted_pt     = Var("userFloat('phiPi_refitted_pt')",float),
        phiPi_refitted_eta    = Var("userFloat('phiPi_refitted_eta')",float),
        phiPi_refitted_phi    = Var("userFloat('phiPi_refitted_phi')",float),
        phiPi_refitted_m      = Var("userFloat('phiPi_refitted_m')",float),
        phiPi_refitted_boost  = Var("userFloat('phiPi_refitted_boost')",float),

        # Bs fit
        mu_refitted_px   = Var("userFloat('mu_refitted_px')",float),
        mu_refitted_py   = Var("userFloat('mu_refitted_py')",float),
        mu_refitted_pz   = Var("userFloat('mu_refitted_pz')",float),
        mu_refitted_pt   = Var("userFloat('mu_refitted_pt')",float),
        mu_refitted_eta  = Var("userFloat('mu_refitted_eta')",float),
        mu_refitted_phi  = Var("userFloat('mu_refitted_phi')",float),
        mu_refitted_m    = Var("userFloat('mu_refitted_m')",float),

        ds_refitted_px   = Var("userFloat('ds_refitted_px')",float),
        ds_refitted_py   = Var("userFloat('ds_refitted_py')",float),
        ds_refitted_pz   = Var("userFloat('ds_refitted_pz')",float),
        ds_refitted_pt   = Var("userFloat('ds_refitted_pt')",float),
        ds_refitted_eta  = Var("userFloat('ds_refitted_eta')",float),
        ds_refitted_phi  = Var("userFloat('ds_refitted_phi')",float),
        ds_refitted_m    = Var("userFloat('ds_refitted_m')",float),
        ds_refitted_boost  = Var("userFloat('ds_refitted_boost')",float),

        bs_fitted_px     = Var("userFloat('bs_fitted_px')",float),
        bs_fitted_py     = Var("userFloat('bs_fitted_py')",float),
        bs_fitted_pz     = Var("userFloat('bs_fitted_pz')",float),
        bs_fitted_pt     = Var("userFloat('bs_fitted_pt')",float),
        bs_fitted_eta    = Var("userFloat('bs_fitted_eta')",float),
        bs_fitted_phi    = Var("userFloat('bs_fitted_phi')",float),
        bs_fitted_m      = Var("userFloat('bs_fitted_m')",float),

        mu_iso_03     = Var("userFloat('mu_iso_03')", float),
        mu_iso_04     = Var("userFloat('mu_iso_04')", float),
        mu_rel_iso_03 = Var("userFloat('mu_rel_iso_03')", float),
        mu_rel_iso_04 = Var("userFloat('mu_rel_iso_04')", float),
        mu_rel_iso_03_refitted = Var("userFloat('mu_rel_iso_03_refitted')", float),
        mu_rel_iso_04_refitted = Var("userFloat('mu_rel_iso_04_refitted')", float),

        e_gamma          = Var("userFloat('e_gamma')",float),
)

##################################################

helicityVariables = cms.PSet(

        angPiK1             = Var("userFloat('angPiK1')",float),
        angPiK2             = Var("userFloat('angPiK2')",float),
        cosPiK1             = Var("userFloat('cosPiK1')",float),
        cosPiK2             = Var("userFloat('cosPiK2')",float),

        angPhiDs_coll       = Var("userFloat('angPhiDsColl')",float),
        angPhiDs_lhcb       = Var("userFloat('angPhiDsLhcb')",float),
        angPhiDs_lhcb_alt   = Var("userFloat('angPhiDsLhcbAlt')",float),
        angPhiDs_reco_1      = Var("userFloat('angPhiDsReco1')",float),
        angPhiDs_reco_2      = Var("userFloat('angPhiDsReco2')",float),

        angPiDs_coll        = Var("userFloat('angPiDsColl')",float),
        angPiDs_lhcb        = Var("userFloat('angPiDsLhcb')",float),
        angPiDs_lhcb_alt    = Var("userFloat('angPiDsLhcbAlt')",float),
        angPiDs_reco_1       = Var("userFloat('angPiDsReco1')",float),
        angPiDs_reco_2       = Var("userFloat('angPiDsReco2')",float),

        cosPhiDs_coll       = Var("userFloat('cosPhiDsColl')",float),
        cosPhiDs_lhcb       = Var("userFloat('cosPhiDsLhcb')",float),
        cosPhiDs_lhcb_alt   = Var("userFloat('cosPhiDsLhcbAlt')",float),
        cosPhiDs_reco_1      = Var("userFloat('cosPhiDsReco1')",float),
        cosPhiDs_reco_2      = Var("userFloat('cosPhiDsReco2')",float),

        cosPiDs_coll        = Var("userFloat('cosPiDsColl')",float),
        cosPiDs_lhcb        = Var("userFloat('cosPiDsLhcb')",float),
        cosPiDs_lhcb_alt    = Var("userFloat('cosPiDsLhcbAlt')",float),
        cosPiDs_reco_1       = Var("userFloat('cosPiDsReco1')",float),
        cosPiDs_reco_2       = Var("userFloat('cosPiDsReco2')",float),

        angMuW_coll         = Var("userFloat('angMuWColl')",float),
        angMuW_lhcb         = Var("userFloat('angMuWLhcb')",float),
        angMuW_lhcb_alt     = Var("userFloat('angMuWLhcbAlt')",float),
        angMuW_reco_1        = Var("userFloat('angMuWReco1')",float),
        angMuW_reco_2        = Var("userFloat('angMuWReco1')",float),

        cosMuW_coll         = Var("userFloat('cosMuWColl')",float),
        cosMuW_lhcb         = Var("userFloat('cosMuWLhcb')",float),
        cosMuW_lhcb_alt     = Var("userFloat('cosMuWLhcbAlt')",float),
        cosMuW_reco_1        = Var("userFloat('cosMuWReco1')",float),
        cosMuW_reco_2        = Var("userFloat('cosMuWReco2')",float),

        cosPlaneBs_coll     = Var("userFloat('cosPlaneBsColl')",float),
        cosPlaneBs_lhcb     = Var("userFloat('cosPlaneBsLhcb')",float),
        cosPlaneBs_lhcb_alt = Var("userFloat('cosPlaneBsLhcbAlt')",float),
        cosPlaneBs_reco_1    = Var("userFloat('cosPlaneBsReco1')",float),
        cosPlaneBs_reco_2    = Var("userFloat('cosPlaneBsReco2')",float),

        cosPlaneDs_coll     = Var("userFloat('cosPlaneDsColl')",float),
        cosPlaneDs_lhcb     = Var("userFloat('cosPlaneDsLhcb')",float),
        cosPlaneDs_lhcb_alt = Var("userFloat('cosPlaneDsLhcbAlt')",float),
        cosPlaneDs_reco_1    = Var("userFloat('cosPlaneDsReco1')",float),
        cosPlaneDs_reco_2    = Var("userFloat('cosPlaneDsReco2')",float),


)

##################################################

bsMomentumVariables = cms.PSet(

        bsMassCorr       = Var("userFloat('bs_m_corr')",float), 

        #pi_mu_charge = Var("userInt('pi_mu_charge')",int),
        bs_px_coll         = Var("userFloat('bs_px_coll')",float),
        bs_py_coll         = Var("userFloat('bs_py_coll')",float),
        bs_pz_coll         = Var("userFloat('bs_pz_coll')",float),
        bs_pt_coll         = Var("userFloat('bs_pt_coll')",float),
        bs_eta_coll        = Var("userFloat('bs_eta_coll')",float),
        bs_phi_coll        = Var("userFloat('bs_phi_coll')",float),
        bs_boost_coll      = Var("userFloat('b_boost_coll')",float),
        bs_boost_pt_coll   = Var("userFloat('b_boost_coll_pt')",float),
        bs_boost_eta_coll  = Var("userFloat('b_boost_coll_eta')",float),
        bs_boost_phi_coll  = Var("userFloat('b_boost_coll_phi')",float),

        bs_px_lhcb         = Var("userFloat('bs_px_lhcb')",float),
        bs_py_lhcb         = Var("userFloat('bs_py_lhcb')",float),
        bs_pz_lhcb         = Var("userFloat('bs_pz_lhcb')",float),
        bs_pt_lhcb         = Var("userFloat('bs_pt_lhcb')",float),
        bs_eta_lhcb        = Var("userFloat('bs_eta_lhcb')",float),
        bs_phi_lhcb        = Var("userFloat('bs_phi_lhcb')",float),
        bs_boost_lhcb      = Var("userFloat('b_boost_lhcb')",float),
        bs_boost_pt_lhcb   = Var("userFloat('b_boost_lhcb_pt')",float),
        bs_boost_eta_lhcb  = Var("userFloat('b_boost_lhcb_eta')",float),
        bs_boost_phi_lhcb  = Var("userFloat('b_boost_lhcb_phi')",float),

        bs_px_lhcb_alt        = Var("userFloat('bs_px_lhcb_alt')",float),
        bs_py_lhcb_alt        = Var("userFloat('bs_py_lhcb_alt')",float),
        bs_pz_lhcb_alt        = Var("userFloat('bs_pz_lhcb_alt')",float),
        bs_pt_lhcb_alt        = Var("userFloat('bs_pt_lhcb_alt')",float),
        bs_eta_lhcb_alt       = Var("userFloat('bs_eta_lhcb_alt')",float),
        bs_phi_lhcb_alt       = Var("userFloat('bs_phi_lhcb_alt')",float),
        bs_boost_lhcb_alt     = Var("userFloat('b_boost_lhcb_alt')",float),
        bs_boost_pt_lhcb_alt  = Var("userFloat('b_boost_lhcb_alt_pt')",float),
        bs_boost_eta_lhcb_alt = Var("userFloat('b_boost_lhcb_alt_eta')",float),
        bs_boost_phi_lhcb_alt = Var("userFloat('b_boost_lhcb_alt_phi')",float),

        disc_is_negative    = Var("userInt('disc_is_negative')",int),
        disc_negativity     = Var("userFloat('disc_negativity')",float),

        bs_px_reco_1        = Var("userFloat('bs_px_reco_1')",float),
        bs_py_reco_1        = Var("userFloat('bs_py_reco_1')",float),
        bs_pz_reco_1        = Var("userFloat('bs_pz_reco_1')",float),
        bs_pt_reco_1        = Var("userFloat('bs_pt_reco_1')",float),
        bs_eta_reco_1       = Var("userFloat('bs_eta_reco_1')",float),
        bs_phi_reco_1       = Var("userFloat('bs_phi_reco_1')",float),
        bs_boost_reco_1     = Var("userFloat('b_boost_reco_1')",float),
        bs_boost_pt_reco_1  = Var("userFloat('b_boost_reco_1_pt')",float),
        bs_boost_eta_reco_1 = Var("userFloat('b_boost_reco_1_eta')",float),
        bs_boost_phi_reco_1 = Var("userFloat('b_boost_reco_1_phi')",float),

        bs_px_reco_2         = Var("userFloat('bs_px_reco_2')",float),
        bs_py_reco_2         = Var("userFloat('bs_py_reco_2')",float),
        bs_pz_reco_2         = Var("userFloat('bs_pz_reco_2')",float),
        bs_pt_reco_2         = Var("userFloat('bs_pt_reco_2')",float),
        bs_eta_reco_2        = Var("userFloat('bs_eta_reco_2')",float),
        bs_phi_reco_2        = Var("userFloat('bs_phi_reco_2')",float),
        bs_boost_reco_2      = Var("userFloat('b_boost_reco_2')",float),
        bs_boost_pt_reco_2   = Var("userFloat('b_boost_reco_2_pt')",float),
        bs_boost_eta_reco_2  = Var("userFloat('b_boost_reco_2_eta')",float),
        bs_boost_phi_reco_2  = Var("userFloat('b_boost_reco_2_phi')",float),

        #lhcb_pz = Var("userFloat('lhcb_pz')",float),
        #theta = Var("userFloat('theta')",float), 

        m2_miss_coll     = Var("userFloat('m2_miss_coll')",float),
        m2_miss_lhcb     = Var("userFloat('m2_miss_lhcb')",float),
        m2_miss_lhcb_alt = Var("userFloat('m2_miss_lhcb_alt')",float),
        m2_miss_reco_1   = Var("userFloat('m2_miss_reco_1')",float),
        m2_miss_reco_2   = Var("userFloat('m2_miss_reco_2')",float),

        pt_miss_coll     = Var("userFloat('pt_miss_coll')",float),
        pt_miss_lhcb     = Var("userFloat('pt_miss_lhcb')",float),
        pt_miss_lhcb_alt = Var("userFloat('pt_miss_lhcb_alt')",float),
        pt_miss_reco_1   = Var("userFloat('pt_miss_reco_1')",float),
        pt_miss_reco_2   = Var("userFloat('pt_miss_reco_2')",float),

        q2_coll          = Var("userFloat('q2_coll')",float),
        q2_lhcb          = Var("userFloat('q2_lhcb')",float),
        q2_lhcb_alt      = Var("userFloat('q2_lhcb_alt')",float),
        q2_reco_1        = Var("userFloat('q2_reco_1')",float),
        q2_reco_2        = Var("userFloat('q2_reco_2')",float),

        e_star_coll          = Var("userFloat('e_star_coll')",float),
        e_star_lhcb          = Var("userFloat('e_star_lhcb')",float),
        e_star_lhcb_alt      = Var("userFloat('e_star_lhcb_alt')",float),
        e_star_reco_1        = Var("userFloat('e_star_reco_1')",float),
        e_star_reco_2        = Var("userFloat('e_star_reco_2')",float),

        #arrived          = Var("userInt('arrived')",int),
)

##################################################

## this is for the ED Filter
#empty = cms.PSet(arrived = Var("userInt('arrived')",int))

