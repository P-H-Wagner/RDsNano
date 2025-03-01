#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
// for vertex fitting (both global and sequential)
#include "TrackingTools/Records/interface/TransientTrackRecord.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h" 
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h" 
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h" 
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h" 
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h" 
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h" 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h" 

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // for the tracks!
#include "FWCore/Framework/interface/MakerMacros.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include "KinVtxFitter.h" //--> not needed now 
#include "helper.h" // helper functions
// include "PxPyPzMVector.h" // to new :(
#include "TLorentzVector.h" // use this instead 
#include "TVector3.h" // for boost vector
// for gen matching
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h" 

// B field
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
// for 3DPoint
#include "DataFormats/Math/interface/Point3D.h"
// for the cov matrix correction
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h" 


// counters
int nEventsGen = 0;
int nMuonsGen  = 0;
int nTracksGen = 0;
int nRecoCandidates = 0;
int nGenMatched = 0;

int muSel1CounterGen = 0;
int muSel2CounterGen = 0;
int k1Sel1CounterGen = 0;
int k1Sel2CounterGen = 0;
int k2Sel1CounterGen = 0;
int k2Sel2CounterGen = 0;
int piSel1CounterGen = 0;
int piSel2CounterGen = 0;
int nKKPiMuGen = 0;
int nFoundPhi  = 0;
int nFoundDs   = 0;
int nFoundB    = 0;
int nBMassCut   = 0;

class genMatching : public edm::global::EDProducer<> {

public:

  //define collections which dont exist by default  
  typedef std::vector<reco::GenParticle> GenParticleCollection;
  typedef std::vector<pat::PackedGenParticle> PackedGenParticleCollection;
  //constructor
  explicit genMatching(const edm::ParameterSet&);
  //destructor
  ~genMatching() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  virtual void endJob() override;

  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}


  //must be constant as it takes constant arguments!! otherwise compiler rises an error
  //because it thinks it may modify the input!!
  reco::TransientTrack getTransientTrack(const reco::Track track) const {    
      reco::TransientTrack transientTrack(track, paramField);

      return transientTrack;
    }


private:
  
  //Bfield
  OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");

  //cuts 
  const StringCutObjectSelector<pat::PackedCandidate> hadSelection_; // cut on hadrons
  //const StringCutObjectSelector<pat::PackedGenParticle> hadSelectionGen_; // cut on gen hadrons
  const StringCutObjectSelector<reco::GenParticle> hadSelectionGen_; // cut on gen hadrons for test with pruned only!

  const double minMuPt_;
  const double maxMuEta_;
  const double maxdRHadMuon_;
  const double mindRHadMuon_;
  const double maxdzDiffHadMuon_; 
  const double phiMassAllowance_;
  const double dsMassAllowance_;
  const double drMatchGen_;
  const double maxBsMass_;
  const double piMass_;
  const double kMass_;
  const double phiMass_;
  const double dsMass_;
  const double dsStarMass_;
  const double muMass_;
  const double tauMass_;
  const double bsMass_;
  const double isoCone_;
  //tokens to access data later
  //edm::Input tag can not be directly initialized inside the construcor! Why did it work fro Trigger.cc??
  //anyway ... 

  //for the muons

  // vertices
  const edm::InputTag primaryVtxTag;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVtx_;

  const edm::InputTag bsTag;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> bs_;

  //gen for gen-matching
  const edm::InputTag prunedGenTag; //pruned is a compressed packed format
  const edm::EDGetTokenT<reco::GenParticleCollection> prunedGen_;

  const edm::InputTag packedGenTag; //packed contains much more info->most likely not needed!
  const edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGen_;

};

//define the constructor
genMatching::genMatching(const edm::ParameterSet& iConfig):
    // f.e. hadSelection_ = cfg.getPatameter...

    hadSelection_(iConfig.getParameter<std::string>("hadSelection")),
    hadSelectionGen_(iConfig.getParameter<std::string>("hadSelectionGen")),
    minMuPt_(iConfig.getParameter<double>("minMuPt")),
    maxMuEta_(iConfig.getParameter<double>("maxMuEta")),
    maxdRHadMuon_(iConfig.getParameter<double>("maxdRHadMuon")),
    mindRHadMuon_(iConfig.getParameter<double>("mindRHadMuon")),
    maxdzDiffHadMuon_(iConfig.getParameter<double>("maxdzDiffHadMuon")),
    phiMassAllowance_(iConfig.getParameter<double>("phiMassAllowance")),
    dsMassAllowance_(iConfig.getParameter<double>("dsMassAllowance")),
    drMatchGen_(iConfig.getParameter<double>("drMatchGen")),
    maxBsMass_(iConfig.getParameter<double>("maxBsMass")),

    piMass_(iConfig.getParameter<double>("piMass")),
    kMass_(iConfig.getParameter<double>("kMass")),
    phiMass_(iConfig.getParameter<double>("phiMass")),
    dsMass_(iConfig.getParameter<double>("dsMass")),
    dsStarMass_(iConfig.getParameter<double>("dsStarMass")),
    muMass_(iConfig.getParameter<double>("muMass")),
    tauMass_(iConfig.getParameter<double>("tauMass")),
    bsMass_(iConfig.getParameter<double>("bsMass")),
    isoCone_(iConfig.getParameter<double>("isoCone")),

    primaryVtxTag(iConfig.getParameter<edm::InputTag>("pvCand")),
    primaryVtx_(consumes<reco::VertexCollection>(primaryVtxTag)),

    bsTag(iConfig.getParameter<edm::InputTag>("bs")),
    bs_(consumes<pat::CompositeCandidateCollection>(bsTag)), 

    prunedGenTag(iConfig.getParameter<edm::InputTag>("prunedCand")),
    prunedGen_(consumes<reco::GenParticleCollection>(prunedGenTag)),
    packedGenTag(iConfig.getParameter<edm::InputTag>("packedCand")),
    packedGen_(consumes<pat::PackedGenParticleCollection>(packedGenTag)){
       // output collection
       produces<pat::CompositeCandidateCollection>("gen");
       //produces<pat::CompositeCandidateCollection>("gen");
       //produces<TransientTrackCollection>("kkTransientTracks");
    }

//check const keywords 

// this starts the event loop
void genMatching::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  //input
  edm::Handle<reco::VertexCollection> primaryVtx;
  iEvent.getByToken(primaryVtx_,primaryVtx);

  edm::Handle<pat::CompositeCandidateCollection> bsColl;
  iEvent.getByToken(bs_,bsColl);
 
  edm::Handle<reco::GenParticleCollection> prunedGen;
  iEvent.getByToken(prunedGen_,prunedGen);

  edm::Handle<pat::PackedGenParticleCollection> packedGen;
  iEvent.getByToken(packedGen_,packedGen);

  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  //std::unique_ptr<pat::CompositeCandidateCollection> ret_value_gen(new pat::CompositeCandidateCollection());
  //std::unique_ptr<TransientTrackCollection> kkpi_ttrack(new TransientTrackCollection);

  nEventsGen++;

  //////////////////////////////////////////////////////
  // Match the bsCand (and all its final states) with //
  // a gen particle. We loop over all candidates      //
  // per event!                                       //
  //////////////////////////////////////////////////////

  pat::CompositeCandidate gen; 
  for(size_t bsIdx = 0; bsIdx < bsColl->size(); ++bsIdx){
  

    nRecoCandidates++;
 
    pat::CompositeCandidate gen; 
    edm::Ptr<pat::CompositeCandidate> bsPtr(bsColl, bsIdx);

    //get kinematic info of final states
    auto muBs = bsPtr->userCand("mu");
    auto k1Bs = bsPtr->userCand("k1");
    auto k2Bs = bsPtr->userCand("k2");
    auto piBs = bsPtr->userCand("pi");
 

    // reco cand
    //std::cout << "mu pt: " << muBs->pt() << std::endl;
    //std::cout << "k1 pt: " << k1Bs->pt() << std::endl;


    int sigId = -9999;
    int bId = 0;
    int genMatchSuccess = 0;

    //count the number of gen matches we find, ideally only 1
    int nGenMatches = 0;

    ////////////////////////////////////////////////////
    // find the gen-matched muon                      //
    ////////////////////////////////////////////////////

    for(size_t muIdxGen = 0; muIdxGen < prunedGen->size(); ++muIdxGen){

      nMuonsGen++;

      //define a pointer to the gen muon    
      edm::Ptr<reco::GenParticle> muPtrGen(prunedGen, muIdxGen);

      //select only useful gen muons -> check this selection!
      if((fabs(muPtrGen->pdgId()) != 13) || muPtrGen->pt() < minMuPt_ || fabs(muPtrGen->eta()) > maxMuEta_ || (muBs->charge() * muPtrGen->charge() < 0)) continue; 
      muSel1CounterGen++;

      //now check the dR of the reco muon wrt to the gen Muon 
      float drMuonMatch = reco::deltaR(*muBs,*muPtrGen);

      //TODO:define as variable
      if(drMuonMatch > drMatchGen_) continue;
      muSel2CounterGen++;

      //std::cout << muPtrGen->pt() << std::endl;
      //std::cout << "found a gen matched muon" << std::endl;

      ////////////////////////////////////////////////
      // find gen matched k1                        //
      ////////////////////////////////////////////////

      for(size_t k1IdxGen = 0; k1IdxGen < prunedGen->size(); ++k1IdxGen){
           
        nTracksGen++;
        //define a pointer to the gen kaon    
        edm::Ptr<reco::GenParticle> k1PtrGen(prunedGen, k1IdxGen);

        //select only useful kaons -> check this selection!
        if((fabs(k1PtrGen->pdgId()) != 321) || !hadSelectionGen_(*k1PtrGen) || (k1Bs->charge() * k1PtrGen->charge() < 0)) continue; 
        k1Sel1CounterGen++;

        //std::cout << k1PtrGen->pt() << std::endl;
        //std::cout << "passed selection for k1!" << std::endl;
        //now check the dR of the reco muon wrt to the gen Muon 
        float drK1Match = reco::deltaR(*k1Bs,*k1PtrGen);

        //TODO:define as variable
        if(drK1Match > drMatchGen_) continue;
        k1Sel2CounterGen++;
        //std::cout << "found a gen matched k1!" << std::endl;

        ////////////////////////////////////////////////
        // find gen matched k2                        //
        ////////////////////////////////////////////////
     
        for(size_t k2IdxGen = 0; k2IdxGen < prunedGen->size(); ++k2IdxGen){
     
             //avoid picking the same gen particle as for k1
             if(k2IdxGen == k1IdxGen) continue; 

             //define a pointer to the gen kaon    
             edm::Ptr<reco::GenParticle> k2PtrGen(prunedGen, k2IdxGen);
     
             //select only useful kaons -> check this selection!
             if((fabs(k2PtrGen->pdgId()) != 321) || !hadSelectionGen_(*k2PtrGen) || (k2Bs->charge() * k2PtrGen->charge() < 0 )) continue; 
             //std::cout << "reco  pt " << k2Bs->pt() << "reco eta " << k2Bs->eta() << std::endl;
             //std::cout << "gen pt " << k2PtrGen->pt() << "gen eta " << k2PtrGen->eta() << std::endl;
             k2Sel1CounterGen++;
             //std::cout << "passed selection for k2" << std::endl;
             //now check the dR of the reco muon wrt to the gen Muon 
             float drK2Match = reco::deltaR(*k2Bs,*k2PtrGen);
     
             //TODO:define as variable
             if(drK2Match > drMatchGen_ ) continue;
             k2Sel2CounterGen++;

             //std::cout << "found a gen matched k2!\n" << std::endl;

             ////////////////////////////////////////////////
             // find gen matched pion                      //
             ////////////////////////////////////////////////

             for(size_t piIdxGen = 0; piIdxGen < prunedGen->size(); ++piIdxGen){
      
               //avoid picking the same gen particle as for k1 or k2
               if((piIdxGen == k1IdxGen) || (piIdxGen == k2IdxGen)) continue; 
                       
               //define a pointer to the gen kaon    
               edm::Ptr<reco::GenParticle> piPtrGen(prunedGen, piIdxGen);
     
               //select only useful kaons -> check this selection!
               if((fabs(piPtrGen->pdgId()) != 211) || !hadSelectionGen_(*piPtrGen) || (piBs->charge() * piPtrGen->charge() < 0 )) continue; 
               piSel1CounterGen++;
     
               //std::cout << "gen pt " << piPtrGen->pt() << "gen eta " << piPtrGen->eta() << std::endl;
               //now check the dR of the reco muon wrt to the gen Muon 
               float drPiMatch = reco::deltaR(*piBs,*piPtrGen);
               //std::cout << drPiMatch << std::endl; 
               //TODO:define as variable
               if(drPiMatch > drMatchGen_) continue;
               //std::cout << "found a gen matched pion!" << std::endl;
               piSel2CounterGen++;

               //////////////////////////////////////////////////
               // Find resonances at gen level                 //
               //////////////////////////////////////////////////
               nKKPiMuGen++;     
 
               //Should we pick the best gen match (in terms of dR) only? -> No, like this is better 

               const reco::Candidate* k1Reco = k1PtrGen.get(); 
               const reco::Candidate* k2Reco = k2PtrGen.get(); 
               const reco::Candidate* piReco = piPtrGen.get(); 
               const reco::Candidate* muReco = muPtrGen.get(); 

               // searching for phi resonance 
               auto phiFromK1 = getAncestor(k1Reco,333);
               auto phiFromK2 = getAncestor(k2Reco,333);
               if( (phiFromK1 != phiFromK2) || (phiFromK1 == nullptr) || (phiFromK2 == nullptr)) continue; 
               nFoundPhi++;               
     
               // searching for ds resonance 
               auto dsFromPhi = getAncestor(phiFromK1,431);
               auto dsFromPi  = getAncestor(piReco,431);
               if( (dsFromPhi != dsFromPi) || (dsFromPhi == nullptr) || (dsFromPi == nullptr)) continue; 
               nFoundDs++;               

               // we dont know what b mother we have
               int bMotherId = 0;

               // first search for b baryon (isAncestor also checks neg. Ids)
               for(int bIdx = 5000; bIdx < 6000; bIdx++){
                 if (isAncestor(dsFromPhi, bIdx)){
                   bMotherId = bIdx;
                   break;
                 }
               }

               // Then search for possible B meson coming from B-baryon
               for(int bIdx = 500; bIdx < 600; bIdx++){
                 if (isAncestor(dsFromPhi, bIdx)){
                   bMotherId = bIdx;
                   break;
                 }
               }
              
               if (bMotherId == 0) break; // no b mother found

               // Even if the mu is not required to come brom the b mother directly (would be signal case)
               // if it comes from another D meson (double charm background case), we still want
               // that this D meson is coming from the b mother. So the muon should share
               // the same ancestor as the Ds.
               auto bsFromDs = getAncestor(dsFromPhi,bMotherId);
               auto bsFromMu = getAncestor(muReco,   bMotherId);

               if( (bsFromDs != bsFromMu) || (bsFromDs == nullptr) || (bsFromMu == nullptr)) continue; 
   
               nFoundB++;
               
               if (bsFromMu->mass() > maxBsMass_) continue;
               nBMassCut++;

               //remove oscillations
               auto bsFromMuWOOsc = removeOscillations(bsFromMu);

               nGenMatches++;
               genMatchSuccess = 1;
               nGenMatched++;
               gen.addUserInt("gen_match_success",genMatchSuccess);
               
               //if(nGenMatches > 1) continue; //std::cout <<"there is more than one match!!" << std::endl;

               //get gen 4 momenta
               TLorentzVector genMuTlv; 
               TLorentzVector genK1Tlv; 
               TLorentzVector genK2Tlv; 
               TLorentzVector genPiTlv; 
               TLorentzVector genPhiTlv; 
               TLorentzVector genDsTlv; 
               TLorentzVector genBsTlv; 

               TLorentzVector genMissTlv; //for m2 miss 
               TLorentzVector genQTlv;  // for q2

               genMuTlv.SetXYZM(  muPtrGen->px(),      muPtrGen->py(),      muPtrGen->pz(),      muMass_);
               genK1Tlv.SetXYZM(  k1PtrGen->px(),      k1PtrGen->py(),      k1PtrGen->pz(),      kMass_);
               genK2Tlv.SetXYZM(  k2PtrGen->px(),      k2PtrGen->py(),      k2PtrGen->pz(),      kMass_);
               genPiTlv.SetXYZM(  piPtrGen->px(),      piPtrGen->py(),      piPtrGen->pz(),      piMass_);

               genPhiTlv.SetXYZM( phiFromK1->px(),     phiFromK1->py(),     phiFromK1->pz(),     phiMass_);
               genDsTlv.SetXYZM(  dsFromPi->px(),      dsFromPi->py(),      dsFromPi->pz(),      dsMass_);
               //genBsTlv.SetXYZM(  bsFromMuWOOsc->px(), bsFromMuWOOsc->py(), bsFromMuWOOsc->pz(), bsMass_); //changed
               genBsTlv.SetXYZM(  bsFromMu->px(), bsFromMu->py(), bsFromMu->pz(), bsMass_); //changed

               genMissTlv = genBsTlv - (genDsTlv + genMuTlv); 
               genQTlv    = genBsTlv - (genDsTlv); 

               float m2_miss_gen = genMissTlv.M2();
               float pt_miss_gen = genMissTlv.Pt();
               float q2_gen = genQTlv.M2();
               float e_star_gen   = getEStar(genBsTlv,genMuTlv);
               float e_gamma_gen  = getEGamma(genDsTlv, dsMass_, dsStarMass_);

               gen.addUserFloat("m2_miss_gen",m2_miss_gen);
               gen.addUserFloat("pt_miss_gen",pt_miss_gen);
               gen.addUserFloat("e_star_gen",e_star_gen);
               gen.addUserFloat("e_gamma_gen",e_gamma_gen);
               gen.addUserFloat("q2_gen",q2_gen);

               //vertices
               float pv_x_gen = bsFromMuWOOsc->vx();//This is the bs production vertex!
               float pv_y_gen = bsFromMuWOOsc->vy();
               float pv_z_gen = bsFromMuWOOsc->vz();

               // Do a cross check on gen: Is there another primary vtx in the primary vertex collection which
               // is closer to the gen truth?

               float dummy = 10000;
               int goldenIdxGen = -1;
               reco::Vertex scndPvCandGen;
               for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){
                 edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, vtxIdx);
                 
                 float minDz = std::sqrt(std::pow(vtxPtr->position().x() - pv_x_gen, 2) + std::pow(vtxPtr->position().y() - pv_y_gen, 2) + std::pow(vtxPtr->position().z() - pv_z_gen, 2));
                 if(minDz < dummy){
                   dummy = minDz;
                   goldenIdxGen = vtxIdx;
                 }
               }

               scndPvCandGen = primaryVtx->at(goldenIdxGen);
               float scnd_pv_x_gen = scndPvCandGen.x();//This is the bs production vertex!
               float scnd_pv_y_gen = scndPvCandGen.y();
               float scnd_pv_z_gen = scndPvCandGen.z();
 
               //std::cout << dsFromPi->vtx() << std::endl;

               float sv_x_gen = dsFromPi->vx(); //This is the ds production vertex!
               float sv_y_gen = dsFromPi->vy();
               float sv_z_gen = dsFromPi->vz();

               float tv_x_gen = phiFromK1->vx(); //This is the phi production vertex!
               float tv_y_gen = phiFromK1->vy();
               float tv_z_gen = phiFromK1->vz();

               float fv_x_gen = k1PtrGen->vx(); //This is the k1 production vertex!
               float fv_y_gen = k1PtrGen->vy();
               float fv_z_gen = k1PtrGen->vz();

               //std::cout << "on gen:" << fv_x_gen << std::endl; //This is the k1 production vertex!
               //save the gen info by adding gen candidates of final states

               //gen.addUserCand("mu_gen",muPtrGen);
               //gen.addUserCand("k1_gen",k1PtrGen);
               //gen.addUserCand("k2_gen",k2PtrGen);
               //gen.addUserCand("pi_gen",piPtrGen);


               //printDaughters(bsFromMu); //-> for debugging

               gen.addUserFloat("mu_gen_px"      ,muPtrGen->px());
               gen.addUserFloat("mu_gen_py"      ,muPtrGen->py());
               gen.addUserFloat("mu_gen_pz"      ,muPtrGen->pz());
               gen.addUserFloat("mu_gen_pt"      ,muPtrGen->pt());
               gen.addUserFloat("mu_gen_eta"     ,muPtrGen->eta());
               gen.addUserFloat("mu_gen_phi"     ,muPtrGen->phi());
               gen.addUserFloat("mu_gen_m"    ,muPtrGen->mass());
               gen.addUserFloat("mu_gen_charge"  ,muPtrGen->charge());
               gen.addUserInt(  "mu_gen_pdgid"   ,muPtrGen->pdgId());
     
               gen.addUserFloat("k1_gen_px"      ,k1PtrGen->px());
               gen.addUserFloat("k1_gen_py"      ,k1PtrGen->py());
               gen.addUserFloat("k1_gen_pz"      ,k1PtrGen->pz());
               gen.addUserFloat("k1_gen_pt"      ,k1PtrGen->pt());
               gen.addUserFloat("k1_gen_eta"     ,k1PtrGen->eta());
               gen.addUserFloat("k1_gen_phi"     ,k1PtrGen->phi());
               gen.addUserFloat("k1_gen_m"    ,k1PtrGen->mass());
               gen.addUserFloat("k1_gen_charge"  ,k1PtrGen->charge());
               gen.addUserInt(  "k1_gen_pdgid"   ,k1PtrGen->pdgId());
     
               gen.addUserFloat("k2_gen_px"      ,k2PtrGen->px());
               gen.addUserFloat("k2_gen_py"      ,k2PtrGen->py());
               gen.addUserFloat("k2_gen_pz"      ,k2PtrGen->pz());
               gen.addUserFloat("k2_gen_pt"      ,k2PtrGen->pt());
               gen.addUserFloat("k2_gen_eta"     ,k2PtrGen->eta());
               gen.addUserFloat("k2_gen_phi"     ,k2PtrGen->phi());
               gen.addUserFloat("k2_gen_m"    ,k2PtrGen->mass());
               gen.addUserFloat("k2_gen_charge"  ,k2PtrGen->charge());
               gen.addUserInt(  "k2_gen_pdgid"   ,k2PtrGen->pdgId());
     
               gen.addUserFloat("pi_gen_px"      ,piPtrGen->px());
               gen.addUserFloat("pi_gen_py"      ,piPtrGen->py());
               gen.addUserFloat("pi_gen_pz"      ,piPtrGen->pz());
               gen.addUserFloat("pi_gen_pt"      ,piPtrGen->pt());
               gen.addUserFloat("pi_gen_eta"     ,piPtrGen->eta());
               gen.addUserFloat("pi_gen_phi"     ,piPtrGen->phi());
               gen.addUserFloat("pi_gen_m"    ,piPtrGen->mass());
               gen.addUserFloat("pi_gen_charge"  ,piPtrGen->charge());
               gen.addUserInt(  "pi_gen_pdgid"   ,piPtrGen->pdgId());

               //and gen info from the resonances
               gen.addUserFloat("phi_gen_px"     ,phiFromK1->px());
               gen.addUserFloat("phi_gen_py"     ,phiFromK1->py());
               gen.addUserFloat("phi_gen_pz"     ,phiFromK1->pz());
               gen.addUserFloat("phi_gen_pt"     ,phiFromK1->pt());
               gen.addUserFloat("phi_gen_eta"    ,phiFromK1->eta());
               gen.addUserFloat("phi_gen_phi"    ,phiFromK1->phi());
               gen.addUserFloat("tv_x_gen"       ,tv_x_gen);//This is the phi production vertex!
               gen.addUserFloat("tv_y_gen"       ,tv_y_gen);
               gen.addUserFloat("tv_z_gen"       ,tv_z_gen);
               gen.addUserFloat("phi_gen_charge" ,phiFromK1->charge());
               gen.addUserInt(  "phi_gen_pdgid"  ,phiFromK1->pdgId());

               gen.addUserFloat("ds_gen_px"     ,dsFromPi->px());
               gen.addUserFloat("ds_gen_py"     ,dsFromPi->py());
               gen.addUserFloat("ds_gen_pz"     ,dsFromPi->pz());
               gen.addUserFloat("ds_gen_pt"     ,dsFromPi->pt());
               gen.addUserFloat("ds_gen_eta"    ,dsFromPi->eta());
               gen.addUserFloat("ds_gen_phi"    ,dsFromPi->phi());
               gen.addUserFloat("ds_gen_m"      ,dsFromPi->mass());
               gen.addUserFloat("ds_gen_boost"  ,genDsTlv.BoostVector().Mag());


               gen.addUserFloat("sv_x_gen"      ,sv_x_gen);//This is the ds production vertex!
               gen.addUserFloat("sv_y_gen"      ,sv_y_gen);
               gen.addUserFloat("sv_z_gen"      ,sv_z_gen);
               gen.addUserFloat("ds_gen_charge" ,dsFromPi->charge());
               gen.addUserInt(  "ds_gen_pdgid"  ,dsFromPi->pdgId());

               gen.addUserFloat("bs_gen_px"     ,bsFromMu->px());
               gen.addUserFloat("bs_gen_py"     ,bsFromMu->py());
               gen.addUserFloat("bs_gen_pz"     ,bsFromMu->pz());
               gen.addUserFloat("bs_gen_pt"     ,bsFromMu->pt());
               gen.addUserFloat("bs_gen_eta"    ,bsFromMu->eta());
               gen.addUserFloat("bs_gen_phi"    ,bsFromMu->phi());
               gen.addUserFloat("bs_gen_m"      ,bsFromMu->mass());

               gen.addUserFloat("pv_x_gen"      ,pv_x_gen); //This is the bs production vertex!
               gen.addUserFloat("pv_y_gen"      ,pv_y_gen);
               gen.addUserFloat("pv_z_gen"      ,pv_z_gen);
               gen.addUserFloat("scnd_pv_x_gen"      ,scnd_pv_x_gen); //This is the bs production vertex!
               gen.addUserFloat("scnd_pv_y_gen"      ,scnd_pv_y_gen);
               gen.addUserFloat("scnd_pv_z_gen"      ,scnd_pv_z_gen);
               gen.addUserInt("scnd_pv_idx_gen"      ,goldenIdxGen);

               gen.addUserFloat("bs_gen_charge" ,bsFromMu->charge());
               gen.addUserInt(  "bs_gen_pdgid"  ,bsFromMu->pdgId());
               gen.addUserFloat("b_boost_gen"   ,genBsTlv.BoostVector().Mag());
               gen.addUserFloat("b_boost_gen_pt"   ,genBsTlv.BoostVector().Pt());
               gen.addUserFloat("b_boost_gen_eta"   ,genBsTlv.BoostVector().Eta());
               gen.addUserFloat("b_boost_gen_phi"   ,genBsTlv.BoostVector().Phi());

               //lets also store the fourth vertex ( the k production vertex)
               gen.addUserFloat("fv_x_gen"       ,fv_x_gen);//This is the kaon production vertex!
               gen.addUserFloat("fv_y_gen"       ,fv_y_gen);
               gen.addUserFloat("fv_z_gen"       ,fv_z_gen);



               ///////////////////////////// 
               // now find the channel ID //
               /////////////////////////////
  
               bId = bMotherId;
  
               // Step1: the b Mother fixes the 10s
               switch(abs(bMotherId)){
                 case 521:  sigId = 100;  break;  // B+
                 case 511:  sigId = 200;  break;  // B0
                 case 531:  sigId = 300;  break;  // Bs
                 case 5122: sigId = 400;  break;  // LambdaB
                 default:   sigId = 500;          // anything else
               }
  
               //bool isNotDoubleCharm = false;
               //auto photonPtr = printDirectDaughters(bsFromMu, false);
  
               //std::cout << "candidate nr: " << nRecoCandidates << std::endl;
               int dsID = getDsID(piPtrGen); // get charmed strange ID
               //std::cout << "ds Id is: " << dsID << std::endl;
               int dID  = getSecondCharmID(muPtrGen); // get charmed ID
               //std::cout << "d Id is: " << dID << std::endl;
               bool isTau = isAncestor(muPtrGen, 15);
  
               switch(dsID){
                 case 431:   sigId += 0;  break; // Ds+
                 case 433:   sigId += 10; break; // Ds+*
                 case 10431: sigId += 20; break; // Ds+(2317)*
                 case 20433: sigId += 30; break; // Ds+(2457)*
                 default:    sigId = 500; break; // anything else
               }
  
               int checkSignal = -1;
               int checkKNuMu  = -1;
  
               if (abs(bMotherId) == 531) checkSignal = isSignal(bsFromMu);
               if (abs(bMotherId) == 521) checkKNuMu  = isKNuMu(bsFromMu);
  
               if (checkSignal==-1){
                 std::cout << "this is not tagged as signal" << std::endl;
                 //printDaughters(bsFromMu);
               }
  
               // Signal candidates enter here
               if (checkSignal != -1)    sigId = checkSignal;
               else if(checkKNuMu != -1) sigId = checkKNuMu;
  
               else {
                 switch(dID){
                   case 411:   sigId += 0;  break; // D+
                   case 421:   sigId += 1; break;  // D0
                   case 431:   sigId += 2; break;  // Ds
                   case 413:   sigId += 3; break;  // D+*
                   case 423:   sigId += 4; break;  // D0*
                   case 433:   sigId += 5; break;  // Ds+*
                   case 4122:  sigId += 6; break;  // Lambda c
  
                   default:    sigId = 500; break; // anything else
                 }
               }
  
               // we want to be sure that we dont tag a Hb channel which was not simulated!
               // since the other b meson can decay freely, this can happen! F.e. we can have stuff like
               // Bs -> Double charm + additional pions/kaons/gammas
               std::set<int> bs_channels      = {302, 300, 303, 312, 305, 315, 0, 1, 10 ,11};
               std::set<int> b0_channels      = {200, 203, 210, 213, 212, 205, 215, 220, 223, 230, 233};
               std::set<int> bplus_channels   = { 107, 108, 101, 104, 117, 118, 121, 124, 131, 134, 111, 114};
               std::set<int> lambdab_channels = {406, 416};
  
               if      ((abs(bMotherId) == 531) && (bs_channels.find(sigId)       == bs_channels.end() ))      sigId = 500;
               else if ((abs(bMotherId) == 521) && (bplus_channels.find(sigId)    == bplus_channels.end() ))   sigId = 500;
               else if ((abs(bMotherId) == 511) && (b0_channels.find(sigId)       == b0_channels.end() ))      sigId = 500;
               else if ((abs(bMotherId) == 5122) && (lambdab_channels.find(sigId) == lambdab_channels.end() )) sigId = 500;
  
               //double photonEnergy;
               //if (photonPtr != nullptr) photonEnergy =photonPtr->energy();
               //else photonEnergy = -9999;
  
  
               /////////////////////////////////////////////

               ////////////////////////////////////
               // SPECIAL FOR HAMMER:            //
               // Save also tau and Ds* info     //
               ////////////////////////////////////
  
               float tau_gen_pt;
               float tau_gen_eta;
               float tau_gen_phi;
               float tau_gen_m;
               int   tau_gen_pdgid; 
 
               float dsStar_gen_pt;
               float dsStar_gen_eta;
               float dsStar_gen_phi;
               float dsStar_gen_m;
               int   dsStar_gen_pdgid; 
 
               ////////////////////////////////////
               // SPECIAL FOR STAR SIGNALS:      //
               // Save also photon info          //
               ////////////////////////////////////

               float photon_gen_pt;
               float photon_gen_eta;
               float photon_gen_phi;
               float photon_gen_pdgid;
               float dr_gen_photon_ds;
               float gen_dsPhoton_m;
               float gen_dsPhotonMu_m;
 

               if (sigId == 0){
  
                 //  Bs -> Ds + mu + nu 
                 tau_gen_pt       = -9999;
                 tau_gen_eta      = -9999;
                 tau_gen_phi      = -9999;
                 tau_gen_m        = -9999;
                 tau_gen_pdgid    = -9999;
  
                 dsStar_gen_pt    = -9999;
                 dsStar_gen_eta   = -9999;
                 dsStar_gen_phi   = -9999;
                 dsStar_gen_m     = -9999;
                 dsStar_gen_pdgid = -9999;
 
                 photon_gen_pt     = -9999; 
                 photon_gen_eta    = -9999; 
                 photon_gen_phi    = -9999; 
                 photon_gen_pdgid  = -9999; 

                 dr_gen_photon_ds  = -9999; 
                 gen_dsPhoton_m    = -9999; 
                 gen_dsPhotonMu_m  = -9999; 
 
 
               }
  
  
               else if (sigId == 1){
  
                 //  Bs -> Ds + tau + nu 
  
                 // get the tau (we know its there)
                 auto tauFromMu   = getAncestor(muReco,15);
  
                 tau_gen_pt       = tauFromMu->pt();
                 tau_gen_eta      = tauFromMu->eta();
                 tau_gen_phi      = tauFromMu->phi();
                 tau_gen_m        = tauMass_; 
                 tau_gen_pdgid    = tauFromMu->pdgId();
  
                 dsStar_gen_pt    = -9999;
                 dsStar_gen_eta   = -9999;
                 dsStar_gen_phi   = -9999;
                 dsStar_gen_m     = -9999;
                 dsStar_gen_pdgid = -9999;
 
                 photon_gen_pt     = -9999; 
                 photon_gen_eta    = -9999; 
                 photon_gen_phi    = -9999; 
                 photon_gen_pdgid  = -9999; 
 
                 dr_gen_photon_ds  = -9999; 
                 gen_dsPhoton_m    = -9999; 
                 gen_dsPhotonMu_m  = -9999; 

               }
  
               else if (sigId == 10){
  
                 //  Bs -> Ds* + mu + nu 
                 tau_gen_pt       = -9999;
                 tau_gen_eta      = -9999;
                 tau_gen_phi      = -9999;
                 tau_gen_m        = -9999;
                 tau_gen_pdgid    = -9999;
  
                 // get the Ds* (we know its there)
                 auto dsStarFromDs = getAncestor(dsFromPi,433);
                 auto gFromDs      = getDaughter(dsStarFromDs, 22); 
 
                 dsStar_gen_pt    = dsStarFromDs->pt();
                 dsStar_gen_eta   = dsStarFromDs->eta();
                 dsStar_gen_phi   = dsStarFromDs->phi();
                 dsStar_gen_m     = dsStarMass_; 
                 dsStar_gen_pdgid = dsStarFromDs->pdgId();

                 if (gFromDs != nullptr){ 
                   //can be nullptr if f.e. Ds* -> Ds + pi^0
                   photon_gen_pt     = gFromDs->pt(); 
                   photon_gen_eta    = gFromDs->eta();
                   photon_gen_phi    = gFromDs->phi();
                   photon_gen_pdgid  = gFromDs->pdgId();
 
                   dr_gen_photon_ds  = reco::deltaR(*gFromDs, *dsFromPi); 

                 }

                 else{
                   photon_gen_pt     = -9999; 
                   photon_gen_eta    = -9999; 
                   photon_gen_phi    = -9999; 
                   photon_gen_pdgid  = -9999; 
   
                   dr_gen_photon_ds  = -9999; 
                   gen_dsPhoton_m    = -9999; 
                   gen_dsPhotonMu_m  = -9999; 


                 }


               }
  
               else if (sigId == 11){
  
                 //  Bs -> Ds* + tau + nu 
  
                 // get the tau (we know its there)
                 auto tauFromMu = getAncestor(muReco,15);
  
                 tau_gen_pt       = tauFromMu->pt();
                 tau_gen_eta      = tauFromMu->eta();
                 tau_gen_phi      = tauFromMu->phi();
                 tau_gen_m        = tauMass_; 
                 tau_gen_pdgid    = tauFromMu->pdgId();
  
                 // get the Ds* (we know its there)
                 auto dsStarFromDs = getAncestor(dsFromPi,433);
                 auto gFromDs      = getDaughter(dsStarFromDs, 22);
  
                 dsStar_gen_pt     = dsStarFromDs->pt();
                 dsStar_gen_eta    = dsStarFromDs->eta();
                 dsStar_gen_phi    = dsStarFromDs->phi();
                 dsStar_gen_m      = dsStarMass_; 
                 dsStar_gen_pdgid  = dsStarFromDs->pdgId();
 
                 if (gFromDs != nullptr){ 
                   //can be nullptr if f.e. Ds* -> Ds + pi^0
                   photon_gen_pt     = gFromDs->pt(); 
                   photon_gen_eta    = gFromDs->eta();
                   photon_gen_phi    = gFromDs->phi();
                   photon_gen_pdgid  = gFromDs->pdgId();
 
                   dr_gen_photon_ds  = reco::deltaR(*gFromDs, *dsFromPi); 
                 }
                 else{
                   photon_gen_pt     = -9999; 
                   photon_gen_eta    = -9999; 
                   photon_gen_phi    = -9999; 
                   photon_gen_pdgid  = -9999; 
   
                   dr_gen_photon_ds  = -9999; 
                   gen_dsPhoton_m    = -9999; 
                   gen_dsPhotonMu_m  = -9999; 

                 }
               }
  
               else{
  
                 tau_gen_pt      = -9999;
                 tau_gen_eta     = -9999;
                 tau_gen_phi     = -9999;
                 tau_gen_m       = -9999;
                 tau_gen_pdgid   = -9999;
  
                 dsStar_gen_pt     = -9999;
                 dsStar_gen_eta    = -9999;
                 dsStar_gen_phi    = -9999;
                 dsStar_gen_m      = -9999;
                 dsStar_gen_pdgid  = -9999;

                 photon_gen_pt     = -9999;
                 photon_gen_eta    = -9999;
                 photon_gen_phi    = -9999;
                 photon_gen_pdgid  = -9999;

                 dr_gen_photon_ds  = -9999;
                 gen_dsPhoton_m    = -9999; 
                 gen_dsPhotonMu_m  = -9999; 


  
               }
  
  
               gen.addUserFloat("tau_gen_pt",         tau_gen_pt);
               gen.addUserFloat("tau_gen_eta",        tau_gen_eta);
               gen.addUserFloat("tau_gen_phi",        tau_gen_phi);
               gen.addUserFloat("tau_gen_m",          tau_gen_m);
               gen.addUserInt("tau_gen_pdgid",      tau_gen_pdgid);
  
               gen.addUserFloat("dsStar_gen_pt",      dsStar_gen_pt);
               gen.addUserFloat("dsStar_gen_eta",     dsStar_gen_eta);
               gen.addUserFloat("dsStar_gen_phi",     dsStar_gen_phi);
               gen.addUserFloat("dsStar_gen_m",       dsStar_gen_m);
               gen.addUserInt("dsStar_gen_pdgid",   dsStar_gen_pdgid);

               gen.addUserFloat("photon_gen_pt",      photon_gen_pt);
               gen.addUserFloat("photon_gen_eta",     photon_gen_eta);
               gen.addUserFloat("photon_gen_phi",     photon_gen_phi);
               gen.addUserInt("photon_gen_pdgid",   photon_gen_pdgid);

               gen.addUserFloat("dr_gen_photon_ds",   dr_gen_photon_ds);

               gen.addUserFloat("gen_dsPhoton_m", gen_dsPhoton_m);
               gen.addUserFloat("gen_dsPhotonMu_m", gen_dsPhotonMu_m );


               //define helicity angles

               //test lhcb method on gen
               TLorentzVector lhcbBsTlvGen = lhcbMethod(genDsTlv + genMuTlv, pv_x_gen, pv_y_gen, pv_z_gen, sv_x_gen, sv_y_gen, sv_z_gen, bsMass_);      
               //std::cout << "this is the gen matched case: " << std::endl;
               std::tuple<std::vector<TLorentzVector>,float> recoResultGen = recoMethod(genDsTlv + genMuTlv, pv_x_gen, pv_y_gen, pv_z_gen, sv_x_gen, sv_y_gen, sv_z_gen, bsMass_);      

               std::vector<TLorentzVector> recosGen = std::get<0>(recoResultGen);
               float discNegativityGen                 = std::get<1>(recoResultGen);

               int discIsNegativeGen = 0;
               if (discNegativityGen > 0) discIsNegativeGen = 1;

               TLorentzVector reco1BsTlvGen = recosGen.at(0);
               TLorentzVector reco2BsTlvGen = recosGen.at(1);

               gen.addUserInt("disc_is_negative_gen", discIsNegativeGen); 
               gen.addUserFloat("disc_negativity_gen", discNegativityGen); 

               gen.addUserFloat("bs_gen_lhcb_pt",   lhcbBsTlvGen.Pt());
               gen.addUserFloat("bs_gen_lhcb_eta",  lhcbBsTlvGen.Eta());
               gen.addUserFloat("bs_gen_lhcb_phi",  lhcbBsTlvGen.Phi());

               //angle between Mu and W


               float angMuWGen = angMuW(genDsTlv,genBsTlv,genMuTlv); 

               float angMuWGenLhcb  = angMuW(genDsTlv,lhcbBsTlvGen,genMuTlv); 
               float angMuWGenReco1 = angMuW(genDsTlv,reco1BsTlvGen,genMuTlv); 
               float angMuWGenReco2 = angMuW(genDsTlv,reco2BsTlvGen,genMuTlv); 

               gen.addUserFloat("angMuWGen",angMuWGen);
               gen.addUserFloat("cosMuWGen",cos(angMuWGen));
               gen.addUserFloat("cosMuWGenLhcb",cos(angMuWGenLhcb));
               gen.addUserFloat("cosMuWGenReco1",cos(angMuWGenReco1));
               gen.addUserFloat("cosMuWGenReco2",cos(angMuWGenReco2));


               //angle between k1(k2) and pion in phi rest frame
               float angPiK1Gen  = angDoubleDecay(genPhiTlv, genK1Tlv,  genPiTlv);
               float angPiK2Gen  = angDoubleDecay(genPhiTlv, genK2Tlv,  genPiTlv);
               gen.addUserFloat("angPiK1Gen", angPiK1Gen);
               gen.addUserFloat("angPiK2Gen", angPiK2Gen);
               gen.addUserFloat("cosPiK1Gen", cos(angPiK1Gen));
               gen.addUserFloat("cosPiK2Gen", cos(angPiK2Gen));

               // equivalently, angle between phi(pi) and bs in ds rest frame
               float angPhiDsGen = angDsPi(genDsTlv,  genPhiTlv, genBsTlv);
               float angPiDsGen  = angDsPi(genDsTlv,  genPiTlv,  genBsTlv);
               float angPiDsGenLhcb  = angDsPi(genDsTlv,  genPiTlv,  lhcbBsTlvGen);

               gen.addUserFloat("angPhiDsGen", angPhiDsGen);
               gen.addUserFloat("angPiDsGen",  angPiDsGen);
               gen.addUserFloat("cosPhiDsGen", cos(angPhiDsGen));
               gen.addUserFloat("cosPiDsGen",  cos(angPiDsGen));
               gen.addUserFloat("cosPiDsGenLhcb",  cos(angPiDsGenLhcb));

               //plane angle
               float angPlaneBsGen = angPlane(genDsTlv, genBsTlv, genMuTlv, genPiTlv);
               gen.addUserFloat("angPlaneBsGen", angPlaneBsGen);
               gen.addUserFloat("cosPlaneBsGen", cos(angPlaneBsGen));

               float angPlaneDsGen = angPlane2(genDsTlv, genBsTlv, genK1Tlv, genPiTlv);
               gen.addUserFloat("angPlaneDsGen", angPlaneDsGen);
               gen.addUserFloat("cosPlaneDsGen", cos(angPlaneDsGen));

               //if we reached this point we have found our gen match and we can stop the loop

               goto end;
               //////////////////////////////////////////////////

             }//close gen matching pi loop 
            //break;
            }//close gen matching k2 loop 
          //break;
          } //close gen matching k1 loop
        //break;
        } //close gen matching mu loop

        //if (genMatchSuccess == 0) continue; 
 
        end:
          gen.addUserInt("sig",sigId);
          gen.addUserInt("b_mother_id",bId);
          gen.addUserInt("gen_match_success",genMatchSuccess);
  
          if (genMatchSuccess == 0){
            // no gen match, we store nans
  
            //prepare a dummy (This does not work!! can not add the empty vector as candidate even it compiles.. why??)
            //reco::GenParticle dummy;
            //math::PtEtaPhiMLorentzVector dummyP4(std::nan(""),std::nan("") ,std::nan("") ,std::nan(""));
            //dummy.setP4(dummyP4); 
            //dummy.setCharge(-9999); 
            //dummy.setPdgId( -9999); 
            //edm::Ptr<reco::GenParticle> empty(&dummy, 0);
  
            //gen.addUserCand("mu_gen"          ,empty);
            //gen.addUserCand("k1_gen"          ,empty);
            //gen.addUserCand("k2_gen"          ,empty);
            //gen.addUserCand("pi_gen"          ,empty);
  
            //well, then its a little more tedious
            gen.addUserFloat("m2_miss_gen"    ,std::nan(""));
            gen.addUserFloat("pt_miss_gen"    ,std::nan(""));
            gen.addUserFloat("q2_gen"         ,std::nan(""));
            gen.addUserFloat("e_star_gen"     ,std::nan(""));
            gen.addUserFloat("e_gamma_gen"    ,std::nan(""));
  
            gen.addUserFloat("mu_gen_px"      ,std::nan(""));
            gen.addUserFloat("mu_gen_py"      ,std::nan(""));
            gen.addUserFloat("mu_gen_pz"      ,std::nan(""));
            gen.addUserFloat("mu_gen_pt"      ,std::nan(""));
            gen.addUserFloat("mu_gen_eta"     ,std::nan(""));
            gen.addUserFloat("mu_gen_phi"     ,std::nan(""));
            gen.addUserFloat("mu_gen_m"    ,std::nan(""));
            gen.addUserFloat("mu_gen_charge"  ,std::nan(""));
            gen.addUserInt(  "mu_gen_pdgid"   ,-9999);
  
            gen.addUserFloat("k1_gen_px"      ,std::nan(""));
            gen.addUserFloat("k1_gen_py"      ,std::nan(""));
            gen.addUserFloat("k1_gen_pz"      ,std::nan(""));
            gen.addUserFloat("k1_gen_pt"      ,std::nan(""));
            gen.addUserFloat("k1_gen_eta"     ,std::nan(""));
            gen.addUserFloat("k1_gen_phi"     ,std::nan(""));
            gen.addUserFloat("k1_gen_m"    ,std::nan(""));
            gen.addUserFloat("k1_gen_charge"  ,std::nan(""));
            gen.addUserInt(  "k1_gen_pdgid"   ,-9999);
  
            gen.addUserFloat("k2_gen_px"      ,std::nan(""));
            gen.addUserFloat("k2_gen_py"      ,std::nan(""));
            gen.addUserFloat("k2_gen_pz"      ,std::nan(""));
            gen.addUserFloat("k2_gen_pt"      ,std::nan(""));
            gen.addUserFloat("k2_gen_eta"     ,std::nan(""));
            gen.addUserFloat("k2_gen_phi"     ,std::nan(""));
            gen.addUserFloat("k2_gen_m"    ,std::nan(""));
            gen.addUserFloat("k2_gen_charge"  ,std::nan(""));
            gen.addUserInt(  "k2_gen_pdgid"   ,-9999);
  
            gen.addUserFloat("pi_gen_px"      ,std::nan(""));
            gen.addUserFloat("pi_gen_py"      ,std::nan(""));
            gen.addUserFloat("pi_gen_pz"      ,std::nan(""));
            gen.addUserFloat("pi_gen_pt"      ,std::nan(""));
            gen.addUserFloat("pi_gen_eta"     ,std::nan(""));
            gen.addUserFloat("pi_gen_phi"     ,std::nan(""));
            gen.addUserFloat("pi_gen_m"    ,std::nan(""));
            gen.addUserFloat("pi_gen_charge"  ,std::nan(""));
            gen.addUserInt(  "pi_gen_pdgid"   ,-9999);
  
            gen.addUserFloat("phi_gen_px"     ,std::nan(""));
            gen.addUserFloat("phi_gen_py"     ,std::nan(""));
            gen.addUserFloat("phi_gen_pz"     ,std::nan(""));
            gen.addUserFloat("phi_gen_pt"     ,std::nan(""));
            gen.addUserFloat("phi_gen_eta"    ,std::nan(""));
            gen.addUserFloat("phi_gen_phi"    ,std::nan(""));
            gen.addUserFloat("tv_x_gen"       ,std::nan(""));
            gen.addUserFloat("tv_y_gen"       ,std::nan(""));
            gen.addUserFloat("tv_z_gen"       ,std::nan(""));
            gen.addUserFloat("phi_gen_charge" ,std::nan(""));
            gen.addUserInt(  "phi_gen_pdgid"  ,-9999);
  
            gen.addUserFloat("ds_gen_px"      ,std::nan(""));
            gen.addUserFloat("ds_gen_py"      ,std::nan(""));
            gen.addUserFloat("ds_gen_pz"      ,std::nan(""));
            gen.addUserFloat("ds_gen_pt"      ,std::nan(""));
            gen.addUserFloat("ds_gen_eta"     ,std::nan(""));
            gen.addUserFloat("ds_gen_phi"     ,std::nan(""));
            gen.addUserFloat("ds_gen_m"       ,std::nan(""));
            gen.addUserFloat("sv_x_gen"       ,std::nan(""));
            gen.addUserFloat("sv_y_gen"       ,std::nan(""));
            gen.addUserFloat("sv_z_gen"       ,std::nan(""));
            gen.addUserFloat("ds_gen_charge"  ,std::nan(""));
            gen.addUserInt(  "ds_gen_pdgid"   ,-9999);
            gen.addUserFloat(  "ds_gen_boost"   ,std::nan(""));
  
            gen.addUserFloat("bs_gen_px"      ,std::nan(""));
            gen.addUserFloat("bs_gen_py"      ,std::nan(""));
            gen.addUserFloat("bs_gen_pz"      ,std::nan(""));
            gen.addUserFloat("bs_gen_pt"      ,std::nan(""));
            gen.addUserFloat("bs_gen_eta"     ,std::nan(""));
            gen.addUserFloat("bs_gen_phi"     ,std::nan(""));
            gen.addUserFloat("bs_gen_m"       ,std::nan(""));
  
            gen.addUserFloat("pv_x_gen"       ,std::nan(""));
            gen.addUserFloat("pv_y_gen"       ,std::nan(""));
            gen.addUserFloat("pv_z_gen"       ,std::nan(""));
            gen.addUserFloat("scnd_pv_x_gen"  ,std::nan("")); //This is the bs production vertex!
            gen.addUserFloat("scnd_pv_y_gen"  ,std::nan(""));
            gen.addUserFloat("scnd_pv_z_gen"  ,std::nan(""));
            gen.addUserInt("scnd_pv_idx_gen"  ,-9999);
  
  
            gen.addUserFloat("bs_gen_charge"  ,std::nan(""));
            gen.addUserInt(  "bs_gen_pdgid"   ,-9999);
            gen.addUserFloat("b_boost_gen" ,std::nan(""));
            gen.addUserFloat("b_boost_gen_pt" ,std::nan(""));
            gen.addUserFloat("b_boost_gen_eta" ,std::nan(""));
            gen.addUserFloat("b_boost_gen_phi" ,std::nan(""));
  
            gen.addUserInt("disc_is_negative_gen", -9999); 
            gen.addUserFloat("disc_negativity_gen", std::nan("")); 
  
            gen.addUserFloat("bs_gen_lhcb_pt", std::nan(""));
            gen.addUserFloat("bs_gen_lhcb_eta", std::nan(""));
            gen.addUserFloat("bs_gen_lhcb_phi", std::nan(""));
  
            gen.addUserFloat("fv_x_gen"       ,std::nan(""));
            gen.addUserFloat("fv_y_gen"       ,std::nan(""));
            gen.addUserFloat("fv_z_gen"       ,std::nan(""));
   
            gen.addUserFloat("angMuWGen"      ,std::nan(""));
            gen.addUserFloat("cosMuWGen"      ,std::nan(""));
            gen.addUserFloat("cosMuWGenLhcb"      ,std::nan(""));
            gen.addUserFloat("cosMuWGenReco1"      ,std::nan(""));
            gen.addUserFloat("cosMuWGenReco2"      ,std::nan(""));
  
            gen.addUserFloat("angPiK1Gen"     ,std::nan(""));
            gen.addUserFloat("angPiK2Gen"     ,std::nan(""));
            gen.addUserFloat("cosPiK1Gen"     ,std::nan(""));
            gen.addUserFloat("cosPiK2Gen"     ,std::nan(""));
  
            gen.addUserFloat("angPhiDsGen"    ,std::nan(""));
            gen.addUserFloat("angPiDsGen"     ,std::nan(""));
            gen.addUserFloat("cosPhiDsGen"    ,std::nan(""));
            gen.addUserFloat("cosPiDsGen"     ,std::nan(""));
            gen.addUserFloat("cosPiDsGenLhcb"     ,std::nan(""));
  
            gen.addUserFloat("angPlaneBsGen"  ,std::nan(""));
            gen.addUserFloat("cosPlaneBsGen"  ,std::nan(""));
            gen.addUserFloat("angPlaneDsGen"  ,std::nan(""));
            gen.addUserFloat("cosPlaneDsGen"  ,std::nan(""));
  
            //gen.addUserFloat("mu_iso_03_gen"     ,std::nan(""));
            //gen.addUserFloat("mu_iso_04_gen"     ,std::nan(""));
            //gen.addUserFloat("mu_rel_iso_03_gen" ,std::nan(""));
            //gen.addUserFloat("mu_rel_iso_04_gen" ,std::nan(""));

            gen.addUserFloat("tau_gen_pt",  std::nan(""));
            gen.addUserFloat("tau_gen_eta", std::nan(""));
            gen.addUserFloat("tau_gen_phi", std::nan(""));
            gen.addUserFloat("tau_gen_m",  std::nan(""));
            gen.addUserInt("tau_gen_pdgid",  -9999);

            gen.addUserFloat("dsStar_gen_pt",  std::nan(""));
            gen.addUserFloat("dsStar_gen_eta", std::nan(""));
            gen.addUserFloat("dsStar_gen_phi", std::nan(""));
            gen.addUserFloat("dsStar_gen_m",   std::nan(""));
            gen.addUserInt("dsStar_gen_pdgid",   -9999);
 
            gen.addUserFloat("photon_gen_pt",     std::nan("") );
            gen.addUserFloat("photon_gen_eta",    std::nan("") );
            gen.addUserFloat("photon_gen_phi",    std::nan("") );
            gen.addUserInt("photon_gen_pdgid",   -9999);

            gen.addUserFloat("dr_gen_photon_ds", std::nan("") );
            gen.addUserFloat("gen_dsPhoton_m", std::nan("") );
            gen.addUserFloat("gen_dsPhotonMu_m", std::nan("") );
 
  
          }
  
          /////////////////////// END OF VARIABLE DEFINITION //////////////////////
  
          //append candidate at the end of our return value :)
          //ret_value can be a vector!!
          ret_value->emplace_back(gen);
          //ret_value_gen->emplace_back(gen);
  
    } //closing loop over Bs
    iEvent.put(std::move(ret_value), "gen");
}//closing event loop


void genMatching::endJob(){
// Printouts:
std::cout << "\n--------- GEN MATCHING MODULE ----------\n" << std::endl;
std::cout << "#Events in file                                           : " << nEventsGen  << std::endl;
std::cout << "#Gen Muons in file                                        : " << nMuonsGen   << std::endl;
std::cout << "#Gen Tracks in file                                       : " << nTracksGen  << std::endl;
std::cout << "#Reco candidates found                                    : " << nRecoCandidates << std::endl;
std::cout << "#Gen matched candidates                                   : " << nGenMatched << std::endl;

std::cout << "#Gen Kaon 1 which passed the hadronic selection           : " << k1Sel1CounterGen << std::endl;
std::cout << "#Gen Kaon 1 which passed the dR matching                  : " << k1Sel2CounterGen << std::endl;
std::cout << "#Gen Kaon 2 which passed the hadronic selection           : " << k2Sel1CounterGen << std::endl;
std::cout << "#Gen Kaon 2 which passed the dR matching                  : " << k2Sel2CounterGen << std::endl;
std::cout << "#Gen Pions  which passed the hadronic selection           : " << piSel1CounterGen << std::endl;
std::cout << "#Gen Pions  which passed the dR matching                  : " << piSel2CounterGen << std::endl;

std::cout << "\n#KKPiMu Gen combinations:" << nKKPiMuGen << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a Phi         : " << nFoundPhi  << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a Ds          : " << nFoundDs   << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a B-mom       : " << nFoundB    << std::endl;
std::cout << "#KKPiMu Gen combinations for which the B-mom < B mass cut : " << nBMassCut << std::endl;
}


DEFINE_FWK_MODULE(genMatching);
