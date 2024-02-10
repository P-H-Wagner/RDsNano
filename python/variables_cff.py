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
        mu_mass       = ufloat("mu_mass"),
        mu_charge     = Var("userCand('mu').charge()",float),
        mu_pdgId      = Var("userCand('mu').pdgId()",int),

        k1_px         = Var("userCand('k1').px()",float),
        k1_py         = Var("userCand('k1').py()",float),
        k1_pz         = Var("userCand('k1').pz()",float),
        k1_pt         = Var("userCand('k1').pt()",float),
        k1_eta        = Var("userCand('k1').eta()",float),
        k1_phi        = Var("userCand('k1').phi()",float),
        k1_mass       = ufloat("k1_mass"),
        k1_charge     = Var("userCand('k1').charge()",float),
        k1_pdgId      = Var("userCand('k1').pdgId()",int),

        k2_px         = Var("userCand('k2').px()",float),
        k2_py         = Var("userCand('k2').py()",float),
        k2_pz         = Var("userCand('k2').pz()",float),
        k2_pt         = Var("userCand('k2').pt()",float),
        k2_eta        = Var("userCand('k2').eta()",float),
        k2_phi        = Var("userCand('k2').phi()",float),
        k2_mass       = ufloat("k2_mass"),
        k2_charge     = Var("userCand('k2').charge()",float),
        k2_pdgId      = Var("userCand('k2').pdgId()",int),

        pi_px         = Var("userCand('pi').px()",float),
        pi_py         = Var("userCand('pi').py()",float),
        pi_pz         = Var("userCand('pi').pz()",float),
        pi_pt         = Var("userCand('pi').pt()",float),
        pi_eta        = Var("userCand('pi').eta()",float),
        pi_phi        = Var("userCand('pi').phi()",float),
        pi_mass       = ufloat("pi_mass"),
        pi_charge     = Var("userCand('pi').charge()",float),
        pi_pdgId      = Var("userCand('pi').pdgId()",int),

        ## add prefit resonances (f.e. Ds is simply phi + pi, without any fit)

        kk_pt         = Var("userFloat('kk_pt')",float),
        kk_eta        = Var("userFloat('kk_eta')",float),
        kk_phi        = Var("userFloat('kk_phi')",float),
        kk_mass       = Var("userFloat('kk_mass')",float),
        kk_charge     = Var("userFloat('kk_charge')",float),
        kk_deltaR     = Var("userFloat('kk_deltaR')",float),

        phiPi_pt      = Var("userFloat('phiPi_pt')",float),
        phiPi_eta     = Var("userFloat('phiPi_eta')",float),
        phiPi_phi     = Var("userFloat('phiPi_phi')",float),
        phiPi_mass    = Var("userFloat('phiPi_mass')",float),
        phiPi_charge  = Var("userFloat('phiPi_charge')",float),
        phiPi_deltaR  = Var("userFloat('phiPi_deltaR')",float),
       
        dsMu_pt       = Var("userFloat('dsMu_pt')",float),
        dsMu_eta      = Var("userFloat('dsMu_eta')",float),
        dsMu_phi      = Var("userFloat('dsMu_phi')",float),
        dsMu_mass     = Var("userFloat('dsMu_mass')",float),
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
        bs_mass       = Var("mass()",float),
        bs_charge     = Var("charge()",float),
        bs_pdgId      = Var("pdgId()",int),
)

##################################################

genVariables = cms.PSet(
        ## gen particle information
        mu_gen_px     = Var("userCand('mu_gen').px()",float),
        mu_gen_py     = Var("userCand('mu_gen').py()",float),
        mu_gen_pz     = Var("userCand('mu_gen').pz()",float),
        mu_gen_pt     = Var("userCand('mu_gen').pt()",float),
        mu_gen_eta    = Var("userCand('mu_gen').eta()",float),
        mu_gen_phi    = Var("userCand('mu_gen').phi()",float),
        mu_gen_mass   = Var("userCand('mu_gen').mass()",float),
        mu_gen_charge = Var("userCand('mu_gen').charge()",float),
        mu_gen_pdgId  = Var("userCand('mu_gen').pdgId()", int),

        k1_gen_px     = Var("userCand('k1_gen').px()",float),
        k1_gen_py     = Var("userCand('k1_gen').py()",float),
        k1_gen_pz     = Var("userCand('k1_gen').pz()",float),
        k1_gen_pt     = Var("userCand('k1_gen').pt()",float),
        k1_gen_eta    = Var("userCand('k1_gen').eta()",float),
        k1_gen_phi    = Var("userCand('k1_gen').phi()",float),
        k1_gen_mass   = Var("userCand('k1_gen').mass()",float),
        k1_gen_charge = Var("userCand('k1_gen').charge()",float),
        k1_gen_pdgId  = Var("userCand('k1_gen').pdgId()",int),

        k2_gen_px     = Var("userCand('k2_gen').px()",float),
        k2_gen_py     = Var("userCand('k2_gen').py()",float),
        k2_gen_pz     = Var("userCand('k2_gen').pz()",float),
        k2_gen_pt     = Var("userCand('k2_gen').pt()",float),
        k2_gen_eta    = Var("userCand('k2_gen').eta()",float),
        k2_gen_phi    = Var("userCand('k2_gen').phi()",float),
        k2_gen_mass   = Var("userCand('k2_gen').mass()",float),
        k2_gen_charge = Var("userCand('k2_gen').charge()",float),
        k2_gen_pdgId  = Var("userCand('k2_gen').pdgId()",int),

        pi_gen_px     = Var("userCand('pi_gen').px()",float),
        pi_gen_py     = Var("userCand('pi_gen').py()",float),
        pi_gen_pz     = Var("userCand('pi_gen').pz()",float),
        pi_gen_pt     = Var("userCand('pi_gen').pt()",float),
        pi_gen_eta    = Var("userCand('pi_gen').eta()",float),
        pi_gen_phi    = Var("userCand('pi_gen').phi()",float),
        pi_gen_mass   = Var("userCand('pi_gen').mass()",float),
        pi_gen_charge = Var("userCand('pi_gen').charge()",float),
        pi_gen_pdgId  = Var("userCand('pi_gen').pdgId()",int),

        phi_gen_px    = Var("userFloat('phi_gen_px')",float),
        phi_gen_py    = Var("userFloat('phi_gen_py')",float),
        phi_gen_pz    = Var("userFloat('phi_gen_pz')",float),
        phi_gen_pt    = Var("userFloat('phi_gen_pt')",float),
        phi_gen_eta   = Var("userFloat('phi_gen_eta')",float),
        phi_gen_phi   = Var("userFloat('phi_gen_phi')",float),
        tv_x_gen      = Var("userFloat('tv_x_gen')",float),
        tv_y_gen      = Var("userFloat('tv_y_gen')",float),
        tv_z_gen      = Var("userFloat('tv_z_gen')",float),
        phi_gen_charge= Var("userFloat('phi_gen_charge')",float),
        phi_gen_pdgid = Var("userInt('phi_gen_pdgid')",int),

        ds_gen_px     = Var("userFloat('ds_gen_px')",float),
        ds_gen_py     = Var("userFloat('ds_gen_py')",float),
        ds_gen_pz     = Var("userFloat('ds_gen_pz')",float),
        ds_gen_pt     = Var("userFloat('ds_gen_pt')",float),
        ds_gen_eta    = Var("userFloat('ds_gen_eta')",float),
        ds_gen_phi    = Var("userFloat('ds_gen_phi')",float),
        sv_x_gen      = Var("userFloat('sv_x_gen')",float),
        sv_y_gen      = Var("userFloat('sv_y_gen')",float),
        sv_z_gen      = Var("userFloat('sv_z_gen')",float),
        ds_gen_charge = Var("userFloat('ds_gen_charge')",float),
        ds_gen_pdgid  = Var("userInt('ds_gen_pdgid')",int),

        bs_gen_px     = Var("userFloat('bs_gen_px')",float),
        bs_gen_py     = Var("userFloat('bs_gen_py')",float),
        bs_gen_pz     = Var("userFloat('bs_gen_pz')",float),
        bs_gen_pt     = Var("userFloat('bs_gen_pt')",float),
        bs_gen_eta    = Var("userFloat('bs_gen_eta')",float),
        bs_gen_phi    = Var("userFloat('bs_gen_phi')",float),
        pv_x_gen      = Var("userFloat('pv_x_gen')",float),
        pv_y_gen      = Var("userFloat('pv_y_gen')",float),
        pv_z_gen      = Var("userFloat('pv_z_gen')",float),
        bs_gen_charge = Var("userFloat('bs_gen_charge')",float),
        bs_gen_pdgid  = Var("userInt('bs_gen_pdgid')",int),

        fv_x_gen      = Var("userFloat('fv_x_gen')",float),
        fv_y_gen      = Var("userFloat('fv_y_gen')",float),
        fv_z_gen      = Var("userFloat('fv_z_gen')",float),

        m2_miss_gen   = ufloat('m2_miss_gen'),
        q2_gen        = ufloat('q2_gen'),

        angMuWGen         = Var("userFloat('angMuWGen')",float),
        cosMuWGen         = Var("userFloat('cosMuWGen')",float),

        angPiK1Gen        = Var("userFloat('angPiK1Gen')",float),
        cosPiK1Gen        = Var("userFloat('cosPiK1Gen')",float),

        angPiK2Gen        = Var("userFloat('angPiK2Gen')",float),
        cosPiK2Gen        = Var("userFloat('cosPiK2Gen')",float),

        angPiDsGen        = Var("userFloat('angPiDsGen')",float),
        cosPiDsGen        = Var("userFloat('cosPiDsGen')",float),

        angPhiDsGen       = Var("userFloat('angPhiDsGen')",float),
        cosPhiDsGen       = Var("userFloat('cosPhiDsGen')",float),

        gen_match_success = uint('gen_match_success'),


        ## signal id
        sig               = Var("userFloat('sig')",int),
)

##################################################

vertexVariables = cms.PSet(

        pv_x = Var("userFloat('pv_x')",float),        
        pv_y = Var("userFloat('pv_y')",float),        
        pv_z = Var("userFloat('pv_z')",float),        

        sv_x = Var("userFloat('sv_x')",float),
        sv_y = Var("userFloat('sv_y')",float),
        sv_z = Var("userFloat('sv_z')",float),

        tv_x = Var("userFloat('tv_x')",float),
        tv_y = Var("userFloat('tv_y')",float),
        tv_z = Var("userFloat('tv_z')",float),

        fv_x = Var("userFloat('fv_x')",float),
        fv_y = Var("userFloat('fv_y')",float),
        fv_z = Var("userFloat('fv_z')",float),

        ##easy_bs_vtx_x = Var("userFloat('easy_bs_vtx_x')",float),
        ##easy_bs_vtx_y = Var("userFloat('easy_bs_vtx_y')",float),
        ##easy_bs_vtx_z = Var("userFloat('easy_bs_vtx_z')",float),

        lxy_bs     = Var("userFloat('lxy_bs')",float),        
        lxyz_bs    = Var("userFloat('lxyz_bs')",float),        

        lxy_ds     = Var("userFloat('lxy_ds')",float),        
        lxyz_ds    = Var("userFloat('lxyz_ds')",float),        

        lxy_phi    = Var("userFloat('lxy_phi')",float),        
        lxyz_phi   = Var("userFloat('lxyz_phi')",float),        

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

        bs_fitted_px     = Var("userFloat('bs_fitted_px')",float),
        bs_fitted_py     = Var("userFloat('bs_fitted_py')",float),
        bs_fitted_pz     = Var("userFloat('bs_fitted_pz')",float),
        bs_fitted_pt     = Var("userFloat('bs_fitted_pt')",float),
        bs_fitted_eta    = Var("userFloat('bs_fitted_eta')",float),
        bs_fitted_phi    = Var("userFloat('bs_fitted_phi')",float),
        bs_fitted_m      = Var("userFloat('bs_fitted_m')",float),

)

##################################################

helicityVariables = cms.PSet(

        angPiK1          = Var("userFloat('angPiK1')",float),
        angPiK2          = Var("userFloat('angPiK2')",float),
        cosPiK1          = Var("userFloat('cosPiK1')",float),
        cosPiK2          = Var("userFloat('cosPiK2')",float),

        angPhiDsColl     = Var("userFloat('angPhiDsColl')",float),
        angPhiDsLhcb     = Var("userFloat('angPhiDsLhcb')",float),
        angPhiDsLhcbAlt  = Var("userFloat('angPhiDsLhcbAlt')",float),
        angPhiDsReco1    = Var("userFloat('angPhiDsReco1')",float),
        angPhiDsReco2    = Var("userFloat('angPhiDsReco2')",float),

        angPiDsColl     = Var("userFloat('angPiDsColl')",float),
        angPiDsLhcb     = Var("userFloat('angPiDsLhcb')",float),
        angPiDsLhcbAlt  = Var("userFloat('angPiDsLhcbAlt')",float),
        angPiDsReco1    = Var("userFloat('angPiDsReco1')",float),
        angPiDsReco2    = Var("userFloat('angPiDsReco2')",float),

        cosPhiDsColl     = Var("userFloat('cosPhiDsColl')",float),
        cosPhiDsLhcb     = Var("userFloat('cosPhiDsLhcb')",float),
        cosPhiDsLhcbAlt  = Var("userFloat('cosPhiDsLhcbAlt')",float),
        cosPhiDsReco1    = Var("userFloat('cosPhiDsReco1')",float),
        cosPhiDsReco2    = Var("userFloat('cosPhiDsReco2')",float),

        cosPiDsColl     = Var("userFloat('cosPiDsColl')",float),
        cosPiDsLhcb     = Var("userFloat('cosPiDsLhcb')",float),
        cosPiDsLhcbAlt  = Var("userFloat('cosPiDsLhcbAlt')",float),
        cosPiDsReco1    = Var("userFloat('cosPiDsReco1')",float),
        cosPiDsReco2    = Var("userFloat('cosPiDsReco2')",float),

        angMuWColl       = Var("userFloat('angMuWColl')",float),
        angMuWLhcb       = Var("userFloat('angMuWLhcb')",float),
        angMuWLhcbAlt    = Var("userFloat('angMuWLhcbAlt')",float),
        angMuWReco1      = Var("userFloat('angMuWReco1')",float),
        angMuWReco2      = Var("userFloat('angMuWReco1')",float),

        cosMuWColl       = Var("userFloat('cosMuWColl')",float),
        cosMuWLhcb       = Var("userFloat('cosMuWLhcb')",float),
        cosMuWLhcbAlt    = Var("userFloat('cosMuWLhcbAlt')",float),
        cosMuWReco1      = Var("userFloat('cosMuWReco1')",float),
        cosMuWReco2      = Var("userFloat('cosMuWReco2')",float),


)

##################################################

bsMomentumVariables = cms.PSet(

        bsMassCorr       = Var("userFloat('bs_mass_corr')",float), 

        #pi_mu_charge = Var("userInt('pi_mu_charge')",int),
        #kk_charge = Var("userInt('kk_charge')",int),

        bs_px_coll       = Var("userFloat('bs_px_coll')",float),
        bs_py_coll       = Var("userFloat('bs_py_coll')",float),
        bs_pz_coll       = Var("userFloat('bs_pz_coll')",float),
        bs_pt_coll       = Var("userFloat('bs_pt_coll')",float),
        bs_eta_coll      = Var("userFloat('bs_eta_coll')",float),
        bs_phi_coll      = Var("userFloat('bs_phi_coll')",float),

        bs_px_lhcb       = Var("userFloat('bs_px_lhcb')",float),
        bs_py_lhcb       = Var("userFloat('bs_py_lhcb')",float),
        bs_pz_lhcb       = Var("userFloat('bs_pz_lhcb')",float),
        bs_pt_lhcb       = Var("userFloat('bs_pt_lhcb')",float),
        bs_eta_lhcb      = Var("userFloat('bs_eta_lhcb')",float),
        bs_phi_lhcb      = Var("userFloat('bs_phi_lhcb')",float),

        bs_px_lhcb_alt   = Var("userFloat('bs_px_lhcb_alt')",float),
        bs_py_lhcb_alt   = Var("userFloat('bs_py_lhcb_alt')",float),
        bs_pz_lhcb_alt   = Var("userFloat('bs_pz_lhcb_alt')",float),
        bs_pt_lhcb_alt   = Var("userFloat('bs_pt_lhcb_alt')",float),
        bs_eta_lhcb_alt  = Var("userFloat('bs_eta_lhcb_alt')",float),
        bs_phi_lhcb_alt  = Var("userFloat('bs_phi_lhcb_alt')",float),

        bs_px_reco_1     = Var("userFloat('bs_px_reco_1')",float),
        bs_py_reco_1     = Var("userFloat('bs_py_reco_1')",float),
        bs_pz_reco_1     = Var("userFloat('bs_pz_reco_1')",float),
        bs_pt_reco_1     = Var("userFloat('bs_pt_reco_1')",float),
        bs_eta_reco_1    = Var("userFloat('bs_eta_reco_1')",float),
        bs_phi_reco_1    = Var("userFloat('bs_phi_reco_1')",float),

        bs_px_reco_2     = Var("userFloat('bs_px_reco_2')",float),
        bs_py_reco_2     = Var("userFloat('bs_py_reco_2')",float),
        bs_pz_reco_2     = Var("userFloat('bs_pz_reco_2')",float),
        bs_pt_reco_2     = Var("userFloat('bs_pt_reco_2')",float),
        bs_eta_reco_2    = Var("userFloat('bs_eta_reco_2')",float),
        bs_phi_reco_2    = Var("userFloat('bs_phi_reco_2')",float),

        #lhcb_pz = Var("userFloat('lhcb_pz')",float),
        #theta = Var("userFloat('theta')",float), 

        m2_miss_coll     = Var("userFloat('m2_miss_coll')",float),
        m2_miss_lhcb     = Var("userFloat('m2_miss_lhcb')",float),
        m2_miss_lhcb_alt = Var("userFloat('m2_miss_lhcb_alt')",float),
        m2_miss_reco_1   = Var("userFloat('m2_miss_reco_1')",float),
        m2_miss_reco_2   = Var("userFloat('m2_miss_reco_2')",float),

        q2_coll          = Var("userFloat('q2_coll')",float),
        q2_lhcb          = Var("userFloat('q2_lhcb')",float),
        q2_lhcb_alt      = Var("userFloat('q2_lhcb_alt')",float),
        q2_reco_1        = Var("userFloat('q2_reco_1')",float),
        q2_reco_2        = Var("userFloat('q2_reco_2')",float),

        arrived          = Var("userInt('arrived')",int),
)

##################################################

## this is for the ED Filter
empty = cms.PSet(arrived = Var("userInt('arrived')",int))

