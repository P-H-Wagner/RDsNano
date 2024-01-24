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
//for vertex fitting (both global and sequential)
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
//#include "PxPyPzMVector.h" // to new :(
#include "TLorentzVector.h" // use this instead 
#include "TVector3.h" // for boost vector
//for gen matching
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h" 

//B field
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

////////////////////////////////////////////////////////////////////////////////////////////
// TODOS:
//
// - move gen matching in separate module?
// - remove hardcoded numbers
// - helicity plane angles
// - redefine all variables after the fit? save both ?
// - beautify the bs.addUserFloat (...)
// - pos. def. cov matrix 
// - pruned vs packed -> discuss
// - output tree has now empty entries when there is no trigger/signal -> how to avoid this?
// - adapt for Hb background sample
// - how to save kk same sign pair? --> I have an idea, done
// - generally: save only gen matched signals?
// - what if an event has two signals?  
// - divide into submitter chunks
// - save gen information!! 
// - do gen tests, check f.e. refitted p's with gen p's and unfitted p's
///////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////


// function which checks if a genParticle has a certain ancestor 
bool isAncestor(const auto dau, const int id){
  //std::cout << "pdgId = "<< dau->pdgId() << std::endl;
  if (fabs(dau->pdgId()) == id){ 
    return true;
    }
  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    if (isAncestor(dau->mother(momIdx), id)) return true;  
  }
  return false;
}

// function which checks ifdau mom is ancestor of dau
bool hasAncestor(const auto dau, const auto mom){
  //std::cout << "pdgId = "<< dau->pdgId() << std::endl;
  if (dau == mom){ 
    return true;
    }
  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    if (hasAncestor(dau->mother(momIdx), mom)) return true;  
  }
  return false;
}

//function which returns pt eta phi of the ancestor in order to compare
//ancestors.
std::vector<double> infoAncestor(const auto dau, const int id){

  if (fabs(dau->pdgId()) == id){
    return {dau->pt(),dau->eta(),dau->phi(),dau->vx(),dau->vy(),dau->vz()};
  }

  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    if (isAncestor(dau->mother(momIdx), id)) {
      std::vector<double> ptEtaPhiVxVyVz = infoAncestor(dau->mother(momIdx),id);
      return ptEtaPhiVxVyVz;
    }
  }
  std::vector<double> failedVector(6, std::numeric_limits<double>::quiet_NaN());
  return failedVector;
}


// function which returns ancestor, we can match by pointer! :)
auto getAncestor(const auto dau, const int id){

  //the pointer type changes when accessing moms, VERY ANNOYING IN A RECURSIVE FUNCTION
 
  //std::cout << "I am at pdg Id = " << dau->pdgId() << " and vertex vx = " << dau->vx() << std::endl; 
  if ((fabs(dau->pdgId()) == id)){
    //std::cout << "sucess!" << std::endl;
    return dau;
  }

  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    //std::cout << "Now I access mom Nr " << momIdx << std::endl;
    if (isAncestor(dau->mother(momIdx), id)) {
      auto dau2 = getAncestor(dau->mother(momIdx),id);
      return dau2;
    }
  }
  const reco::Candidate* empty = nullptr; 
  return empty;
}

// functino which removes the un-oscillated ancestor of dau
// f.e. dau is -531 and comes from 531 via oscillation, then this function removes oscillation
// and returns a pointer to the 531 particle, which has the correct vertex!
const reco::Candidate* removeOscillations(const auto dau){

  //std::cout << "I am at pdg Id = " << dau->pdgId() << " and vertex vx = " << dau->vx() << std::endl; 

  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){

    //check if dau has a mother with the same pdg Id but opposite sign
    if (dau->mother(momIdx)->pdgId() == (-1 * dau->pdgId())) {
      //std::cout << "oscillation!" << std::endl;

      auto dau2 = removeOscillations(dau->mother(momIdx));
      return dau2;
    }
  }
  return dau;
}


// counters for filters

int nEvents = 0;
int k1Sel1Counter = 0;
int k1Sel2Counter = 0;

int k2Sel1Counter = 0;
int k2Sel2Counter = 0;

int piSel1Counter = 0;
int piSel2Counter = 0;

class BsToDsPhiKKPiMuBuilder : public edm::global::EDProducer<> {

public:

  //define collections which dont exist by default  
  typedef std::vector<reco::GenParticle> GenParticleCollection;
  typedef std::vector<pat::PackedGenParticle> PackedGenParticleCollection;
  //constructor
  explicit BsToDsPhiKKPiMuBuilder(const edm::ParameterSet&);
  //destructor
  ~BsToDsPhiKKPiMuBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}


  //must be constant as it takes constant arguments!! otherwise compiler rises an error
  //because it thinks it may modify the input!!
  reco::TransientTrack getTransientTrack(const reco::Track* track) const {    
      reco::TransientTrack transientTrack(*track, paramField);

      return transientTrack;
    }


private:
  
  //Bfield
  OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");

  //cuts 
  const StringCutObjectSelector<pat::PackedCandidate> hadSelection_; // cut on hadrons
  //const StringCutObjectSelector<pat::PackedGenParticle> hadSelectionGen_; // cut on gen hadrons
  const StringCutObjectSelector<reco::GenParticle> hadSelectionGen_; // cut on gen hadrons for test with pruned only!

  const double maxdRHadMuon_;
  const double mindRHadMuon_;
  const double maxdzDiffHadMuon_; 
  const double phiMassAllowance_;
  const double dsMassAllowance_;
  const double drMatchGen_;
  const double piMass_;
  const double kMass_;
  const double phiMass_;
  const double dsMass_;
  const double muMass_;
  const double bsMass_;

  //tokens to access data later
  //edm::Input tag can not be directly initialized inside the construcor! Why did it work fro Trigger.cc??
  //anyway ... 

  const edm::InputTag srcTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> src_;

  //for the muons

  const edm::InputTag trgMuonTag;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuons_;

  // vertices
  const edm::InputTag primaryVtxTag;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVtx_;
 
  //gen for gen-matching
  const edm::InputTag prunedGenTag; //pruned is a compressed packed format
  const edm::EDGetTokenT<reco::GenParticleCollection> prunedGen_;

  const edm::InputTag packedGenTag; //packed contains much more info->most likely not needed!
  const edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGen_;

};

//define the constructor
BsToDsPhiKKPiMuBuilder::BsToDsPhiKKPiMuBuilder(const edm::ParameterSet& iConfig):
    // f.e. hadSelection_ = cfg.getPatameter...
    hadSelection_(iConfig.getParameter<std::string>("hadSelection")),
    hadSelectionGen_(iConfig.getParameter<std::string>("hadSelectionGen")),
    maxdRHadMuon_(iConfig.getParameter<double>("maxdRHadMuon")),
    mindRHadMuon_(iConfig.getParameter<double>("mindRHadMuon")),
    maxdzDiffHadMuon_(iConfig.getParameter<double>("maxdzDiffHadMuon")),
    phiMassAllowance_(iConfig.getParameter<double>("phiMassAllowance")),
    dsMassAllowance_(iConfig.getParameter<double>("dsMassAllowance")),
    drMatchGen_(iConfig.getParameter<double>("drMatchGen")),
   

    piMass_(iConfig.getParameter<double>("piMass")),
    kMass_(iConfig.getParameter<double>("kMass")),
    phiMass_(iConfig.getParameter<double>("phiMass")),
    dsMass_(iConfig.getParameter<double>("dsMass")),
    muMass_(iConfig.getParameter<double>("muMass")),
    bsMass_(iConfig.getParameter<double>("bsMass")),
 
    srcTag(iConfig.getParameter<edm::InputTag>("pfCand")),
    src_(consumes<pat::PackedCandidateCollection>(srcTag)), 

    trgMuonTag(iConfig.getParameter<edm::InputTag>("muCand")),
    trgMuons_(consumes<pat::MuonCollection>(trgMuonTag)), 

    primaryVtxTag(iConfig.getParameter<edm::InputTag>("pvCand")),
    primaryVtx_(consumes<reco::VertexCollection>(primaryVtxTag)),

    prunedGenTag(iConfig.getParameter<edm::InputTag>("prunedCand")),
    prunedGen_(consumes<reco::GenParticleCollection>(prunedGenTag)),
    packedGenTag(iConfig.getParameter<edm::InputTag>("packedCand")),
    packedGen_(consumes<pat::PackedGenParticleCollection>(packedGenTag)){
       // output collection
       produces<pat::CompositeCandidateCollection>("bs");
       //produces<pat::CompositeCandidateCollection>("gen");
       //produces<TransientTrackCollection>("kkTransientTracks");
    }

//check const keywords 

// this starts the event loop
void BsToDsPhiKKPiMuBuilder::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  //input
  edm::Handle<pat::PackedCandidateCollection> pcand;
  iEvent.getByToken(src_, pcand);
  
  edm::Handle<pat::MuonCollection> trgMuons;
  iEvent.getByToken(trgMuons_,trgMuons);
 
  edm::Handle<reco::VertexCollection> primaryVtx;
  iEvent.getByToken(primaryVtx_,primaryVtx);

  edm::Handle<reco::GenParticleCollection> prunedGen;
  iEvent.getByToken(prunedGen_,prunedGen);

  edm::Handle<pat::PackedGenParticleCollection> packedGen;
  iEvent.getByToken(packedGen_,packedGen);

  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  //std::unique_ptr<pat::CompositeCandidateCollection> ret_value_gen(new pat::CompositeCandidateCollection());
  //std::unique_ptr<TransientTrackCollection> kkpi_ttrack(new TransientTrackCollection);

  std::cout << "---------------- NEW EVENT ---------------" << std::endl;
  nEvents++;
  int arrived = -1;

  //////////////////////////////////////////////////////
  // Match the trigger muon with a muon from the      //
  // packed pat collection. Reminder, we only have one//
  // unique trigger muon candidate from Trigger.cc    //
  //////////////////////////////////////////////////////

  /*
  //define a pointer to the trigger muon called mu_ptr
  edm::Ptr<reco::Muon> trgMuPtr(trgMuons, 0);
  std::cout << trgMuPtr->pt() << std::endl;

  std::vector<float> dzMuPV;
  //Fix the primary vertex to be the one closest to the trg Muon in dz
  // more accurate for mu than for tau signals
    
  for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){

    // when do we use bestTrack() and when do we use TransientTracks()?
    edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, ++vtxIdx);
    dzMuPV.push_back(fabs(trgMuPtr->bestTrack()->dz(vtxPtr->position())));      
  }
    
  // take as the primary vertex the one which has minimal dz with the muon 
  auto dzMuPVMin = std::min_element(std::begin(dzMuPV),std::end(dzMuPV));
  int pvIdx = std::distance(std::begin(dzMuPV), dzMuPVMin); 
  reco::Vertex pv = primaryVtx->at(pvIdx);
  */

  for(size_t muIdx = 0; muIdx < pcand->size(); ++muIdx){

    //define a pointer to the muon called mu_ptr
    edm::Ptr<pat::PackedCandidate> muPtr(pcand, muIdx);
    
    if (!(muPtr->hasTrackDetails()) || (muPtr->pt() < 6) || (fabs(muPtr->eta()) > 2.4) || (fabs(muPtr->pdgId()) != 13) ) continue;

    std::vector<float> dzMuPV;

    //Fix the primary vertex to be the one closest to the trg Muon in dz
    // more accurate for mu than for tau signals
    /*
    for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){

      // when do we use bestTrack() and when do we use TransientTracks()?
      edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, ++vtxIdx);
      dzMuPV.push_back(fabs(muPtr->bestTrack()->dz(vtxPtr->position())));      
    }
    
    // take as the primary vertex the one which has minimal dz with the muon 
    auto dzMuPVMin = std::min_element(std::begin(dzMuPV),std::end(dzMuPV));
    int pvIdx = std::distance(std::begin(dzMuPV), dzMuPVMin); 
    reco::Vertex pv = primaryVtx->at(pvIdx);
    */
    auto pvRef = muPtr->vertexRef();
    GlobalPoint pv(pvRef->x(), pvRef->y(), pvRef->z()); 

    //////////////////////////////////////////////////
    // Loop over k1 and select the good tracks      //
    //////////////////////////////////////////////////

    for(size_t k1Idx = 0; k1Idx < pcand->size(); ++k1Idx) {

      if(k1Idx == muIdx) continue;

      //define a pointer to the kaon at position k1Idx
      edm::Ptr<pat::PackedCandidate> k1Ptr(pcand, k1Idx);

      if (!hadSelection_(*k1Ptr)) continue; 
      k1Sel1Counter++;

      //the PF algorithm assigns a pdgId hypothesis, generall it distinguishes between:
      // photons, electron/muon, charged hadron, neutral hadrons
      // and we trust the algorithm that when it says its an electron (11) or muon (13), that it is not a kaon or pion

      float muonK1dR = reco::deltaR(*k1Ptr,*muPtr);

      bool k1Sel = (( muonK1dR < maxdRHadMuon_ ) && (reco::deltaR(*k1Ptr, *muPtr) > mindRHadMuon_) && (abs(k1Ptr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_));

      if (!k1Sel) continue;
      k1Sel2Counter++;
    //////////////////////////////////////////////////
    // Loop over k2 and select the good tracks      //
    //////////////////////////////////////////////////

      for(size_t k2Idx = k1Idx + 1; k2Idx < pcand->size(); ++k2Idx) {

      //make sure k2 is not k1
      if (k2Idx == k1Idx) continue;

      // pointer to the second kaon candidate
      edm::Ptr<pat::PackedCandidate> k2Ptr(pcand, k2Idx);
     
      // if this kaon does not pass the selection, jump to the next!
      if(!hadSelection_(*k2Ptr)) continue;
      k2Sel1Counter++;

      float muonK2dR = reco::deltaR(*k2Ptr,*muPtr);

      bool k2Sel = (( muonK2dR < maxdRHadMuon_ ) && (reco::deltaR(*k2Ptr, *muPtr) > mindRHadMuon_) && (abs(k2Ptr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_));

      //k1 and k2 must have oppoiste charge -> only for signal tests, later we keep everything
      int kkCharge = k1Ptr->charge() * k2Ptr->charge();

      if (kkCharge > 0) continue; //To be commented out 


      if (!k2Sel) continue;
      k2Sel2Counter++;

      //////////////////////////////////////////////////
      // Loop over pi and select the good tracks      //
      //////////////////////////////////////////////////

      for(size_t piIdx = 0; piIdx < pcand->size(); ++piIdx) {

        //make sure the pion is none of the kaons:
        if((piIdx == k1Idx) || (piIdx == k2Idx)) continue;

        // pointer to the second kaon candidate
        edm::Ptr<pat::PackedCandidate> piPtr(pcand, piIdx);

        // if this pion does not pass the selection, jump to the next!
        if(!hadSelection_(*piPtr)) continue;

        float muonPidR = reco::deltaR(*piPtr,*muPtr);

        bool piSel = ((muonPidR < maxdRHadMuon_) && (reco::deltaR(*piPtr, *muPtr) > mindRHadMuon_) &&(abs(piPtr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_));

        //pi and mu must have opposite charge -> only for signal tests, later we keep everything
        int piMuCharge = piPtr->charge() * muPtr->charge();
        if (piMuCharge > 0) continue; //To be commented out

        if (!piSel) continue;

        //////////////////////////////////////////////////
        // Build Phi resonance                          //
        //////////////////////////////////////////////////

        //define a composite candidate pair 
        pat::CompositeCandidate kkPair;

        //PF canidates are always assigned with the pi mass, we have to force the kaon mass
        math::PtEtaPhiMLorentzVector k1P4(k1Ptr->pt(), k1Ptr->eta(), k1Ptr->phi(), K_MASS);
        math::PtEtaPhiMLorentzVector k2P4(k2Ptr->pt(), k2Ptr->eta(), k2Ptr->phi(), K_MASS);
        math::PtEtaPhiMLorentzVector piP4(piPtr->pt(), piPtr->eta(), piPtr->phi(), PI_MASS); //just to be sure lets also force the pi mass
 
        kkPair.setP4(k1P4 + k2P4);
      
        //only continue when they build a phi resonance, allow 15MeV:
        if (fabs(kkPair.mass() - phiMass_) > phiMassAllowance_) continue;     

        kkPair.setCharge(k1Ptr->charge() + k2Ptr->charge());

        //////////////////////////////////////////////////
        // Build Ds resonance                           //
        // ///////////////////////////////////////////////

        pat::CompositeCandidate ds;
        ds.setP4(kkPair.p4() + piP4); 

        //only continue when they build a ds resonance, allow 50MeV:
        if (fabs(ds.mass() - dsMass_) > dsMassAllowance_) continue;

        ds.setCharge(kkPair.charge() + piPtr->charge());

        //////////////////////////////////////////////////
        // Build Bs resonance                           //
        //////////////////////////////////////////////////

        pat::CompositeCandidate dsMu;
        dsMu.setP4(ds.p4() + muPtr->p4()); 

        dsMu.setCharge(ds.charge() + muPtr->charge()); //sanity check:shoould be 0

        //build bs with collinear approximation
        pat::CompositeCandidate bs;

        bs.setP4(dsMu.p4() * bsMass_ / dsMu.mass()); //the bs_mass will thus be fixed at 536688 (peak in the histo)
        bs.setCharge(dsMu.charge());


        ////////////////////////////////////////////////////
        // Now that we have found resonances at RECO level//
        // check the gen level and assign the signal ID   //
        ////////////////////////////////////////////////////

        // We store the Hb background under the index 4 -> no need for gen match
        // We store the wrong sign pairs for background estimation under the Index 5 -> no need for gen match

        // this value (-1) will change, but it prevents CMSSW from raising unitialized error,
        // but in any case, if we do not find a matched signal, we dont save it!
       
        int sigId = -1; 
        int genMatchSuccess = 0;
        double bsPtGen;
        //pat::CompositeCandidate gen;

        //count the number of gen matches we find, ideally only 1
        int nGenMatches = 0;

        ////////////////////////////////////////////////////
        // find the gen-matched muon                      //
        ////////////////////////////////////////////////////
        int nMuGen = 0;

        for(size_t muIdxGen = 0; muIdxGen < prunedGen->size(); ++muIdxGen){

          //define a pointer to the gen muon    
          edm::Ptr<reco::GenParticle> muPtrGen(prunedGen, muIdxGen);

          //select only useful gen muons -> check this selection!
          if((fabs(muPtrGen->pdgId()) != 13) || muPtrGen->pt() < 6.5 || fabs(muPtrGen->eta()) > 1.5 || (muPtr->charge() * muPtrGen->charge() < 0)) continue; 

          //now check the dR of the reco muon wrt to the gen Muon 
          float drMuonMatch = reco::deltaR(*muPtr,*muPtrGen);

          //TODO:define as variable
          if(drMuonMatch > drMatchGen_) continue;

          std::cout << "found a gen matched muon" << std::endl;
          ++nMuGen;

          ////////////////////////////////////////////////
          // find gen matched k1                        //
          ////////////////////////////////////////////////
          int nK1Gen = 0;

          for(size_t k1IdxGen = 0; k1IdxGen < prunedGen->size(); ++k1IdxGen){
           
            //define a pointer to the gen kaon    
            edm::Ptr<reco::GenParticle> k1PtrGen(prunedGen, k1IdxGen);

            //select only useful kaons -> check this selection!
            if((fabs(k1PtrGen->pdgId()) != 321) || !hadSelectionGen_(*k1PtrGen) || (k1Ptr->charge() * k1PtrGen->charge() < 0)) continue; 

            //now check the dR of the reco muon wrt to the gen Muon 
            float drK1Match = reco::deltaR(*k1Ptr,*k1PtrGen);

            //TODO:define as variable
            if(drK1Match > drMatchGen_) continue;

            std::cout << "found a gen matched k1!" << std::endl;
            ++nK1Gen;

           ////////////////////////////////////////////////
           // find gen matched k2                        //
           ////////////////////////////////////////////////
           int nK2Gen = 0;
     
           for(size_t k2IdxGen = 0; k2IdxGen < prunedGen->size(); ++k2IdxGen){
     
             //avoid picking the same gen particle as for k1
             if(k2IdxGen == k1IdxGen) continue; 

             //define a pointer to the gen kaon    
             edm::Ptr<reco::GenParticle> k2PtrGen(prunedGen, k2IdxGen);
     
             //select only useful kaons -> check this selection!
             if((fabs(k2PtrGen->pdgId()) != 321) || !hadSelectionGen_(*k2PtrGen) || (k2Ptr->charge() * k2PtrGen->charge() < 0 )) continue; 
     
             //now check the dR of the reco muon wrt to the gen Muon 
             float drK2Match = reco::deltaR(*k2Ptr,*k2PtrGen);
     
             //TODO:define as variable
             if(drK2Match > drMatchGen_ ) continue;

             std::cout << "found a gen matched k2!" << std::endl;
             ++nK2Gen;

             ////////////////////////////////////////////////
             // find gen matched pion                      //
             ////////////////////////////////////////////////

             int nPiGen = 0; 
             for(size_t piIdxGen = 0; piIdxGen < prunedGen->size(); ++piIdxGen){
      
               //avoid picking the same gen particle as for k1 or k2
               if((piIdxGen == k1IdxGen) || (piIdxGen == k2IdxGen)) continue; 
                       
               //define a pointer to the gen kaon    
               edm::Ptr<reco::GenParticle> piPtrGen(prunedGen, piIdxGen);
     
               //select only useful kaons -> check this selection!
               if((fabs(piPtrGen->pdgId()) != 211) || !hadSelectionGen_(*piPtrGen) || (piPtr->charge() * piPtrGen->charge() < 0 )) continue; 
     
               //now check the dR of the reco muon wrt to the gen Muon 
               float drPiMatch = reco::deltaR(*piPtr,*piPtrGen);
     
               //TODO:define as variable
               if(drPiMatch > drMatchGen_) continue;
               std::cout << "found a gen matched pion!" << std::endl;
               ++nPiGen;

               //////////////////////////////////////////////////
               // Find resonances at gen level                 //
               //////////////////////////////////////////////////
      
               //Should we pick the best gen match (in terms of dR) only? -> No, like this is better 

               const reco::Candidate* k1Reco = k1PtrGen.get(); 
               const reco::Candidate* k2Reco = k2PtrGen.get(); 
               const reco::Candidate* piReco = piPtrGen.get(); 
               const reco::Candidate* muReco = muPtrGen.get(); 

               //std::cout << "searching for phi resonance .. " << std::endl;
               auto phiFromK1 = getAncestor(k1Reco,333);
               auto phiFromK2 = getAncestor(k2Reco,333);
               if( (phiFromK1 != phiFromK2) || (phiFromK1 == nullptr) || (phiFromK2 == nullptr)) continue; 
                    
               //std::cout << "searching for ds resonance .. " << std::endl;
               auto dsFromPhi = getAncestor(phiFromK1,431);
               auto dsFromPi  = getAncestor(piReco,431);
               if( (dsFromPhi != dsFromPi) || (dsFromPhi == nullptr) || (dsFromPi == nullptr)) continue; 

               //std::cout << "searching for bs resonance .. " << std::endl;
               auto bsFromDs = getAncestor(dsFromPhi,531);
               auto bsFromMu = getAncestor(muReco,531);
               if( (bsFromDs != bsFromMu) || (bsFromDs == nullptr) || (bsFromMu == nullptr)) continue; 

               //remove oscillations
               auto bsFromMuWOOsc = removeOscillations(bsFromMu);
               //bool related = hasAncestor(dsFromK1,bsFromK1);
               
               nGenMatches++;
               genMatchSuccess = 1;

               if(nGenMatches > 1) continue; //std::cout <<"there is more than one match!!" << std::endl;
              
               //save the gen info by adding gen candidates of final states
               bs.addUserCand("gen_mu",muPtrGen);
               bs.addUserCand("gen_k1",k1PtrGen);
               bs.addUserCand("gen_k2",k2PtrGen);
               bs.addUserCand("gen_pi",piPtrGen);

               //and gen info from the resonances
               bs.addUserFloat("gen_phi_pt" ,phiFromK1->pt());
               bs.addUserFloat("gen_phi_eta",phiFromK1->eta());
               bs.addUserFloat("gen_phi_phi",phiFromK1->phi());
               bs.addUserFloat("gen_phi_vx" ,phiFromK1->vx());
               bs.addUserFloat("gen_phi_vy" ,phiFromK1->vy());
               bs.addUserFloat("gen_phi_vz" ,phiFromK1->vz());

               bs.addUserFloat("gen_ds_pt"  ,dsFromPi->pt());
               bs.addUserFloat("gen_ds_eta" ,dsFromPi->eta());
               bs.addUserFloat("gen_ds_phi" ,dsFromPi->phi());
               bs.addUserFloat("gen_ds_vx"  ,dsFromPi->vx());
               bs.addUserFloat("gen_ds_vy"  ,dsFromPi->vy());
               bs.addUserFloat("gen_ds_vz"  ,dsFromPi->vz());

               bs.addUserFloat("gen_bs_pt"  ,bsFromMuWOOsc->pt());
               bs.addUserFloat("gen_bs_eta" ,bsFromMuWOOsc->eta());
               bs.addUserFloat("gen_bs_phi" ,bsFromMuWOOsc->phi());
               bs.addUserFloat("gen_bs_vx"  ,bsFromMuWOOsc->vx());
               bs.addUserFloat("gen_bs_vy"  ,bsFromMuWOOsc->vy());
               bs.addUserFloat("gen_bs_vz"  ,bsFromMuWOOsc->vz());

               // now find the signal ID
               // Ds mu   = 0
               // Ds* mu  = 1
               // Ds tau  = 2
               // Ds* tau = 3
               sigId = 0; //default

               // look for tau in signal 
               if(isAncestor(muPtrGen, 15)) sigId = 2;               

               // now look for possible Ds*
               // TODO: stupid question: this (433) is the only resonance we look at right?
               if(isAncestor(piPtrGen, 433)) sigId += 1; 

               //////////////////////////////////////////////////

             }//close gen matching pi loop 
            }//close gen matching k2 loop 
          } //close gen matching k1 loop
        } //close gen matching mu loop

    
        if (genMatchSuccess == 0) continue; 

        /*
        std::cout << iEvent.id() << std::endl;
        if (genMatchSuccess == 0){
          std::cout << "no gen match, I store nan" << std::endl;

          bs.addUserFloat("gen_phi_pt" ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_phi_eta",std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_phi_phi",std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_phi_vx" ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_phi_vy" ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_phi_vz" ,std::numeric_limits<double>::quiet_NaN());

          bs.addUserFloat("gen_ds_pt"  ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_ds_eta" ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_ds_phi" ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_ds_vx"  ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_ds_vy"  ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_ds_vz"  ,std::numeric_limits<double>::quiet_NaN());

          bs.addUserFloat("gen_bs_pt"  ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_bs_eta" ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_bs_phi" ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_bs_vx"  ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_bs_vy"  ,std::numeric_limits<double>::quiet_NaN());
          bs.addUserFloat("gen_bs_vz"  ,std::numeric_limits<double>::quiet_NaN());
       

        }
        */
        //std::cout << sigId << std::endl; 
        ////std::cout << "1" << std::endl;
        // save signal Id

        //std::cout << "2" << std::endl;

        //easy fit as a sanity check
        KinVtxFitter easyFitter(
        {getTransientTrack(k1Ptr->bestTrack()), getTransientTrack(k2Ptr->bestTrack()),getTransientTrack(piPtr->bestTrack()),getTransientTrack(muPtr->bestTrack())},
        {K_MASS, K_MASS,piMass_,muMass_},
        {0.00005,0.00005,0.00005,0.00005} //some small sigma for the lepton mass
        );
        if(!easyFitter.success()) continue;

        //std::cout << "3" << std::endl;

        ////////////////////////////////////////////////
        // Now we do a proper fit                     //
        ////////////////////////////////////////////////

        //define a factory
        KinematicParticleFactoryFromTransientTrack pFactory;

        //define the vector for the particles to be fitted
        std::vector<RefCountedKinematicParticle> phiToFit;
        std::vector<RefCountedKinematicParticle> dsToFit;
        std::vector<RefCountedKinematicParticle> bsToFit;

        // add masses
        ParticleMass piMass = piMass_;
        ParticleMass kMass = kMass_;
        ParticleMass phiMass = phiMass_;
        ParticleMass dsMass = dsMass_;
        ParticleMass muMass = muMass_;

        float ndf = 0.0;
        float chi = 0.0;
        float sigma = 0.00005;
        //float phiMassSigma = 0.000005; //fine, something small
        //float dsMassSigma = 0.00005;  //discuss
        //float muMassSigma = 0.000005;  //discuss

        phiToFit.push_back(pFactory.particle(getTransientTrack(k1Ptr->bestTrack()),kMass,chi,ndf,sigma));
        phiToFit.push_back(pFactory.particle(getTransientTrack(k2Ptr->bestTrack()),kMass,chi,ndf,sigma));
        dsToFit.push_back(pFactory.particle(getTransientTrack(piPtr->bestTrack()),piMass,chi,ndf,sigma));
        bsToFit.push_back(pFactory.particle(getTransientTrack(muPtr->bestTrack()),muMass,chi,ndf,sigma));

        //std::cout << "4" << std::endl;

        /////////////////////////////// global fitter /////////////////////////////////////////////////////////
        // some remarks: With movePointerToTheFirstChild(), movePointerToTheNextChild() one can
        // access the fitted daughters. In the example of the Ds fit, this is the phi and the pi
        // if one calls finalStateParticles(), one gets *all* the daughters (going down the chain).
        // In the case of the Ds, this is the pi, k1, and k2. Remark that what gets fitted are not
        // these finals states but thr daughters we access directly via : 
        // movePointerToTheFirstChild(), movePointerToTheNextChild()
        // If one wants to go down the chain, this is what is done in finalStateParticles(), one can call
        // movePointerToTheFirstChild() *after* one has alrady called "movePointerToTheFirstChild().
        // This will go one level down.
        // Example in the Ds fit:
        // Call move pointer to the first child -> get pi
        // Call move pointer to the next child -> get phi
        // if you now call move pointer to the first child, it does not go to pi again but it accesses
        // the children of the phi!! So in this case it would return one of the kaons.
        // To go back to the pi one would have to move the pointer to the top again :)
        ///////////////////////////////////////////////////////////////////////////////////////////////////////

        // define fitter
        KinematicConstrainedVertexFitter phiFitter;
        // define cosntraints
        MultiTrackKinematicConstraint* phiConstr = new TwoTrackMassKinematicConstraint(phiMass);
        MultiTrackKinematicConstraint* dsConstr = new TwoTrackMassKinematicConstraint(dsMass);

        // PHI VERTEX FIT
        RefCountedKinematicTree phiTree = phiFitter.fit(phiToFit, phiConstr);
        if (!phiTree->isValid()) continue; //check if fit result is valid

        //access the fitted resonance and the refitted children
        phiTree->movePointerToTheTop();
        RefCountedKinematicParticle phiParticle = phiTree->currentParticle();
        auto phiVtx = phiTree->currentDecayVertex();
        phiTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle phiDau1 = phiTree->currentParticle();
        phiTree->movePointerToTheNextChild();
        RefCountedKinematicParticle phiDau2 = phiTree->currentParticle();

        // get the vectors full of fit information (vertex and momenta and mass)
        AlgebraicVector7 phiParams = phiParticle->currentState().kinematicParameters().vector();     
        AlgebraicVector7 phiDau1Params = phiDau1->currentState().kinematicParameters().vector();     
        AlgebraicVector7 phiDau2Params = phiDau2->currentState().kinematicParameters().vector();
     
        // add the phi to the list of particles (pi) to fit the ds 
        dsToFit.push_back(phiParticle);

        // DS VERTEX FIT
        KinematicConstrainedVertexFitter dsFitter; //check --> define new fitter?

        RefCountedKinematicTree dsTree = dsFitter.fit(dsToFit, dsConstr);
        if (!dsTree->isValid()) continue; //check if fit result is valid

        // access the fitted resonance and the refitted children
        dsTree->movePointerToTheTop();
        RefCountedKinematicParticle dsParticle = dsTree->currentParticle();
        RefCountedKinematicVertex dsDecayVertex = dsTree->currentDecayVertex(); //compare to the access via AlgebraicVector7
        dsTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle dsDau1 = dsTree->currentParticle();
        dsTree->movePointerToTheNextChild();
        RefCountedKinematicParticle dsDau2 = dsTree->currentParticle();

        // get the vectors full of fit information (vertex and momenta and mass)
        AlgebraicVector7 dsParams = dsParticle->currentState().kinematicParameters().vector();     
        AlgebraicVector7 dsDau1Params = dsDau1->currentState().kinematicParameters().vector();     
        AlgebraicVector7 dsDau2Params = dsDau2->currentState().kinematicParameters().vector();

        // add the ds to the list of particles to fit the bs
        bsToFit.push_back(dsParticle);
        
        // BS VERTEX FIT
        KinematicConstrainedVertexFitter bsFitter;

        RefCountedKinematicTree bsTree = bsFitter.fit(bsToFit); // no constraint for bs because missing momentum
        if (!bsTree->isValid()) continue; //check if fit result is valid 

        // access the fitted resonance and the refitted children
        bsTree->movePointerToTheTop();
        RefCountedKinematicParticle bsParticle = bsTree->currentParticle();
        RefCountedKinematicVertex bsDecayVertex = bsTree->currentDecayVertex();

        bsTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle bsDau1 = bsTree->currentParticle();
        bsTree->movePointerToTheNextChild();
        RefCountedKinematicParticle bsDau2 = bsTree->currentParticle();

        // get the vectors full of fit information (vertex and momenta and mass)
        AlgebraicVector7 bsParams = bsParticle->currentState().kinematicParameters().vector();     
        AlgebraicVector7 bsDau1Params = bsDau1->currentState().kinematicParameters().vector();     
        AlgebraicVector7 bsDau2Params = bsDau2->currentState().kinematicParameters().vector();

        //std::cout << "5" << std::endl;

        //////////////////////////////////// end of global fitter /////////////////////////////////////

        ////store basic variables

        bs.addUserInt("gen_match_success",genMatchSuccess);

        //save the indices of the final states in the pruned Collection
        bs.addUserInt("k1_idx",k1Idx);
        bs.addUserInt("k2_idx",k2Idx);
        bs.addUserInt("pi_idx",piIdx);

        //bs.addUserInt("mu_idx",muIdx); //always 0 :)

        //add final states as Candidates
        bs.addUserCand("k1",k1Ptr);  //be aware that this ptr has the wrong mass, need to assign it in variables_cff.py
        bs.addUserFloat("k1_mass",kMass_);

        bs.addUserCand("k2",k2Ptr);  // "
        bs.addUserFloat("k2_mass",kMass_);

        bs.addUserCand("pi",piPtr);
        bs.addUserFloat("pi_mass",piMass_);

        bs.addUserCand("mu",muPtr); //muon mass is correct -> check

        bs.addUserFloat("sig",sigId);

        //add resonances --> are not part of collection and can thus not be
        //automatically access the pt(), .. as I can for k1,k2,pi,mu in the variables_cff.py

        kkPair.addUserFloat("kk_delta_R", reco::deltaR(*k1Ptr, *k2Ptr));

        bs.addUserFloat("phi_pt", kkPair.pt());
        bs.addUserFloat("phi_eta", kkPair.eta());
        bs.addUserFloat("phi_phi", kkPair.phi());
        bs.addUserFloat("phi_mass", kkPair.mass());
        bs.addUserFloat("phi_charge", kkPair.charge());
        bs.addUserFloat("phi_deltaR", kkPair.userFloat("kk_delta_R"));
        bs.addUserInt("kkCharge",kkCharge); 

        ds.addUserFloat("phi_pi_delta_R", reco::deltaR(kkPair, *piPtr));

        bs.addUserFloat("ds_pt", ds.pt());
        bs.addUserFloat("ds_eta", ds.eta());
        bs.addUserFloat("ds_phi", ds.phi());
        bs.addUserFloat("ds_mass", ds.mass());
        bs.addUserFloat("ds_charge", ds.charge());
        bs.addUserFloat("ds_deltaR", ds.userFloat("phi_pi_delta_R"));

        dsMu.addUserFloat("ds_mu_delta_R", reco::deltaR(ds, *muPtr));

        bs.addUserFloat("dsMu_pt", dsMu.pt());
        bs.addUserFloat("dsMu_eta", dsMu.eta());
        bs.addUserFloat("dsMu_phi", dsMu.phi());
        bs.addUserFloat("dsMu_mass", dsMu.mass()); //we miss the neutrino mass
        bs.addUserFloat("dsMu_charge", dsMu.charge());
        bs.addUserFloat("dsMu_deltaR", dsMu.userFloat("ds_mu_delta_R"));

        //rel charges
        bs.addUserInt("kk_charge",kkCharge); 
        bs.addUserInt("pi_mu_charge",piMuCharge); 

        bs.addUserFloat("easy_bs_vtx_x",easyFitter.fitted_vtx().x());
        bs.addUserFloat("easy_bs_vtx_y",easyFitter.fitted_vtx().y());
        bs.addUserFloat("easy_bs_vtx_z",easyFitter.fitted_vtx().z());


        //AlgebraicMatrix77 phiErr = phiParticle->currentState().kinematicParametersError().matrix();
        //AlgebraicMatrix77 dsErr = dsParticle->currentState().kinematicParametersError().matrix();
        //AlgebraicMatrix77 bsErr = bsParticle->currentState().kinematicParametersError().matrix();
        //check if vertex position is truly at 0,1,2: result: Yes it is!
        //double xdecay = phiDecayVertex->position().x();

        // save fitter info
        // Remark: fitted describes the first fit, refitted describes the refit due to daughters.
        // TODO: In the end its sequential, and not fitting all at the same time.. maybe there is a possibility? -> keep looking around
        // TODO: Double check if the order of daughters is kept according to the **ToFit vectors -> include sanity checks


        ////////////////////////////////////// All vertices ///////////////////////////////////////////////
        // Remark: currentDecayVertex() returns the same as the ALgebraic Vector components 0,1,2--> checked :)
        // Remark: the daughters of a fitted vertex return the same vertex via Algebraic vector 0,1,2 like the
        // fitted mother --> checked 
        

        // primary vertex ( = Bs production vertex)
        float pv_x = pv.x();
        float pv_y = pv.y();
        float pv_z = pv.z();

        bs.addUserFloat("pv_x",  pv_x); 
        bs.addUserFloat("pv_y",  pv_y); 
        bs.addUserFloat("pv_z",  pv_z); 
        
        // secondary vertex ( = Bs decay vertex)  

        float sv_x = bsParams(0);
        float sv_y = bsParams(1);
        float sv_z = bsParams(2);

        bs.addUserFloat("sv_x",  sv_x); 
        bs.addUserFloat("sv_y",  sv_y); 
        bs.addUserFloat("sv_z",  sv_z); 
  
        // tertiary vertex ( = Ds decay vertex) 

        float tv_x = dsParams(0);
        float tv_y = dsParams(1);
        float tv_z = dsParams(2);

        bs.addUserFloat("tv_x",  tv_x); 
        bs.addUserFloat("tv_y",  tv_y); 
        bs.addUserFloat("tv_z",  tv_z); 

        // fourth vertex ( = Phi Decay Vertex)

        float fv_x = phiParams(0);
        float fv_y = phiParams(1);
        float fv_z = phiParams(2);

        bs.addUserFloat("fv_x",  fv_x); 
        bs.addUserFloat("fv_y",  fv_y); 
        bs.addUserFloat("fv_z",  fv_z); 

        //std::cout << "6" << std::endl;

        ////////////////////////////////////// lxy(z), dxy(z) //////////////////////////////////////////

        // lxy(z) is the flight distance in the xy(z) plane(space)

        float lxyBs  = std::sqrt(std::pow((pv_x - sv_x),2) + std::pow((pv_y - sv_y),2) ); 
        float lxyzBs = std::sqrt(std::pow((pv_x - sv_x),2) + std::pow((pv_y - sv_y),2) + std::pow((pv_z - sv_z),2) ); 

        float lxyDs  = std::sqrt(std::pow((sv_x - tv_x),2) + std::pow((sv_y - tv_y),2) ); 
        float lxyzDs = std::sqrt(std::pow((sv_x - tv_x),2) + std::pow((sv_y - tv_y),2) + std::pow((sv_z - tv_z),2) ); 

        float lxyPhi  = std::sqrt(std::pow((tv_x - fv_x),2) + std::pow((tv_y - fv_y),2) ); 
        float lxyzPhi = std::sqrt(std::pow((tv_x - fv_x),2) + std::pow((tv_y - fv_y),2) + std::pow((tv_z - fv_z),2) ); 
       
        bs.addUserFloat("lxy_bs",lxyBs);
        bs.addUserFloat("lxyz_bs",lxyzBs);

        bs.addUserFloat("lxy_ds",lxyDs);
        bs.addUserFloat("lxyz_ds",lxyzDs);

        bs.addUserFloat("lxy_phi",lxyPhi);
        bs.addUserFloat("lxyz_phi",lxyzPhi);

        // dxy(z) is the impact parameter in the xy(z) plane(space), i.e. the distance to the PV
        // TODO: check the errors of dxy and dz      
        // TODO: bestTrack() is not refitted -> bad?
 
        float dxyMu = muPtr->bestTrack()->dxy(pv.position());  
        float dxyMuErr = muPtr->bestTrack()->dxyError(pv.position(),pv.covariance());  
        float dxyMuSig = dxyMu/dxyMuErr;

        float dzMu = muPtr->bestTrack()->dz(pv.position());  
        float dzMuErr = muPtr->bestTrack()->dzError();  
        float dzMuSig = dzMu/dzMuErr ; 

        float dxyPi = piPtr->bestTrack()->dxy(pv.position());  //maybe useful for Ds* vs Ds ? 
        float dxyPiErr = piPtr->bestTrack()->dxyError(pv.position(),pv.covariance());  
        float dxyPiSig = dxyPi/dxyPiErr;

        float dzPi = piPtr->bestTrack()->dz(pv.position());  
        float dzPiErr = piPtr->bestTrack()->dzError();  
        float dzPiSig = dzPi/dzPiErr ; 

        float dxyK1 = k1Ptr->bestTrack()->dxy(pv.position()); //needed ? 
        float dxyK1Err = k1Ptr->bestTrack()->dxyError(pv.position(),pv.covariance());
        float dxyK1Sig = dxyK1/dxyK1Err;

        float dzK1 = k1Ptr->bestTrack()->dz(pv.position());  
        float dzK1Err = k1Ptr->bestTrack()->dzError();  
        float dzK1Sig = dzK1/dzK1Err ; 

        float dxyK2 = k2Ptr->bestTrack()->dxy(pv.position()); //needed ? 
        float dxyK2Err = k2Ptr->bestTrack()->dxyError(pv.position(),pv.covariance());
        float dxyK2Sig = dxyK2/dxyK2Err;

        float dzK2 = k2Ptr->bestTrack()->dz(pv.position());  
        float dzK2Err = k2Ptr->bestTrack()->dzError();  
        float dzK2Sig = dzK2/dzK2Err ; 

        bs.addUserFloat("dxy_mu", dxyMu);
        bs.addUserFloat("dz_mu", dzMu);
        bs.addUserFloat("dxy_mu_err", dxyMuErr);
        bs.addUserFloat("dz_mu_err", dzMuErr);
        bs.addUserFloat("dxy_mu_sig", dxyMuSig);
        bs.addUserFloat("dz_mu_sig", dzMuSig);

        bs.addUserFloat("dxy_pi", dxyPi);
        bs.addUserFloat("dz_pi", dzPi);
        bs.addUserFloat("dxy_pi_err", dxyPiErr);
        bs.addUserFloat("dz_pi_err", dzPiErr);
        bs.addUserFloat("dxy_pi_sig", dxyPiSig);
        bs.addUserFloat("dz_pi_sig", dzPiSig);

        bs.addUserFloat("dxy_k1", dxyK1);
        bs.addUserFloat("dz_k1", dzK1);
        bs.addUserFloat("dxy_k1_err", dxyK1Err);
        bs.addUserFloat("dz_k1_err", dzK1Err);
        bs.addUserFloat("dxy_k1_sig", dxyK1Sig);
        bs.addUserFloat("dz_k1_sig", dzK1Sig);

        bs.addUserFloat("dxy_k2", dxyK2);
        bs.addUserFloat("dz_k2", dzK2);
        bs.addUserFloat("dxy_k2_err", dxyK2Err);
        bs.addUserFloat("dz_k2_err", dzK2Err);
        bs.addUserFloat("dxy_k2_sig", dxyK2Sig);
        bs.addUserFloat("dz_k2_sig", dzK2Sig);

        //std::cout << "7" << std::endl;

        ////////////////////////////////////////// deltaR (defined above) /////////////////////////////////////

        bs.addUserFloat("muonK1dR",muonK1dR);
        bs.addUserFloat("muonK2dR",muonK2dR);
        bs.addUserFloat("muonPidR",muonPidR);

        //////////////////////////////////////////corrected Bs mass ///////////////////////////////////////////
        //build corrected Bs mass accoring to https://journals.aps.org/prd/pdf/10.1103/PhysRevD.100.112006 eq 3.
        double bsMassCorr;
       
        TVector3 bFlightDir;
        bFlightDir.SetXYZ(pv_x - sv_x, pv_y - sv_y, 0.0);
        TLorentzVector dsMuTlv;
        dsMuTlv.SetXYZM(dsMu.px(),dsMu.py(),0.0,dsMu.mass());       
        
        double dsMuPerp = dsMuTlv.Vect().Perp(bFlightDir); // I checked that Perp is doing the right thing (via sin(..))

        bsMassCorr = std::sqrt(std::pow(dsMu.mass(),2) + std::pow(dsMuPerp,2)) + dsMuPerp;
        bs.addUserFloat("bs_mass_corr", bsMassCorr);

        //std::cout << "8" << std::endl;

        ///////////////////////// reconstruct the B momentum with different methods ///////////////////////////////
        

        //Define 4 vectors of fitted resonances
        TLorentzVector fittedPhi;
        TLorentzVector fittedDs;
        TLorentzVector fittedBs;

        //Define placeholder for collinear approx
        TLorentzVector collBs;

        fittedPhi.SetXYZM(phiParams(3), phiParams(4), phiParams(5), phiParams(6));
        fittedDs.SetXYZM(dsParams(3), dsParams(4), dsParams(5), dsParams(6));
        fittedBs.SetXYZM(bsParams(3), bsParams(4), bsParams(5), bsParams(6));
        collBs.SetPtEtaPhiM(bs.p4().pt(), bs.p4().eta(), bs.p4().phi(), bs.p4().mass()); //the old coll approx (before the fit)

        bs.addUserFloat("bs_fitted_pt",fittedBs.Pt());
        bs.addUserFloat("ds_fitted_pt",fittedDs.Pt());
        bs.addUserFloat("phi_fitted_pt",fittedPhi.Pt());

        //Define 4 vectors of refitted final states
        TLorentzVector refittedK1;
        TLorentzVector refittedK2;
        TLorentzVector refittedPi;
        TLorentzVector refittedMu;

        refittedK1.SetXYZM(phiDau1Params(3), phiDau1Params(4), phiDau1Params(5), phiDau1Params(6));
        refittedK2.SetXYZM(phiDau2Params(3), phiDau2Params(4), phiDau2Params(5), phiDau2Params(6));
        refittedPi.SetXYZM(dsDau1Params(3), dsDau1Params(4), dsDau1Params(5), dsDau1Params(6));
        refittedMu.SetXYZM(bsDau1Params(3), bsDau1Params(4), bsDau1Params(5), bsDau1Params(6));
 
        // first lets do again the coll. approx but now after the fit :)

        TLorentzVector refittedCollBs = fittedDs + refittedMu;
        double refittedDsMuMass = refittedCollBs.M();
        refittedCollBs *= bsMass_ / refittedDsMuMass; //scale it
        bs.addUserFloat("bs_pt_coll",refittedCollBs.Pt());

        //std::cout << "9" << std::endl;

        // now lets use the LHCbs method 
        TVector3 bsFlightDir;
        TVector3 beamAxis;
        TVector3 radialAxis;
        dsMuTlv.SetXYZM(dsMu.px(),dsMu.py(),dsMu.pz(),dsMu.mass());       

        //bsFlightDir.SetXYZ(sv_x - pv_x, sv_y - pv_y, sv_z - pv_z);
        bsFlightDir.SetXYZ(bs.userFloat("gen_ds_vx") - pv_x, bs.userFloat("gen_ds_vy") - pv_y , bs.userFloat("gen_ds_vz") - pv_z);

        beamAxis.SetXYZ(0.0,0.0,1.0); // in z direction
        radialAxis.SetXYZ(1.0,0.0,0.0); //in x direction

        //double theta = bsFlightDir.Angle(beamAxis); //angle between beam axis and bs flight dir
        float theta = bsFlightDir.Theta(); //angle between beam axis and bs flight dir

        TLorentzVector refittedDsMu = fittedDs + refittedMu;
        float lhcbPz = dsMuTlv.Pz() *bsMass_ / refittedDsMu.M();
        float lhcbPt = lhcbPz * std::tan(theta); //angle is in radians! std::tan also takes radians :)

        bs.addUserFloat("lhcb_pz",lhcbPz);
        bs.addUserFloat("theta",theta); 

        bsFlightDir.SetZ(0.0); //project on xy plane for phi calculation
        double lhcbPhi = bsFlightDir.Angle(radialAxis); //get the phi angle

        TLorentzVector refittedLhcbBs; 
        double eta = - std::log(std::tan(theta/2));
        refittedLhcbBs.SetPtEtaPhiM(lhcbPt,eta,lhcbPhi,bsMass_); 
        bs.addUserFloat("bs_pt_lhcb",lhcbPt);

        //std::cout << "10" << std::endl;

        //now lets do the reco method
        double neuPar; //neutrino momentum parallel to bs direction
        double bAbs; //absolute 3 momentum of bs
        double bs_pt_reco;

        bsFlightDir.SetZ(sv_z - pv_z); //reset bs flight direction
        // First define the solution for the parallel neutrino momentum
        
        double par1 = refittedDsMu.Vect().Mag();
        double recoAng = refittedDsMu.Vect().Angle(bsFlightDir); 
        double par2 = std::cos(recoAng);
        double par3 = std::sin(recoAng);
        double par4 = refittedDsMu.E(); 

        //solution according to mitternachtsformel 

        double a = 1.0;
        double b = 2*par1*par2 - 2*par4;
        double c = std::pow(par1,2)*std::pow(par2,2) + std::pow(par1*par3,2) - std::pow(par4,2) + std::pow(bsMass_,2);
        double disc = std::pow(b,2) - 4*a*c;      
        if(disc>0) {
          //non complex root 
          neuPar = (-b + std::sqrt(disc)) / (2*a);
          bAbs = par1*par2 + neuPar; 
          TVector3 bsReco = bsFlightDir;
          bsReco *= bAbs; 
          bs_pt_reco = bsReco.Pt();

        }
        else{ 
        bs_pt_reco = std::nan("");
        }

        bs.addUserFloat("bs_pt_reco",bs_pt_reco); 

        //gen level pt
        //bs.addUserFloat("bs_pt_gen", bsPtGen);
        //std::cout << "11" << std::endl;

        ////////////////// Refitted momenta (and masses for consistency, even if constrained /////////////////////////

        bs.addUserFloat("phi_fitted_px", phiParams(3)); 
        bs.addUserFloat("phi_fitted_py", phiParams(4)); 
        bs.addUserFloat("phi_fitted_pz", phiParams(5)); 
        bs.addUserFloat("phi_fitted_m",  phiParams(6)); 

        bs.addUserFloat("k1_refitted_vx",  phiDau1Params(0)); 
        bs.addUserFloat("k1_refitted_vy",  phiDau1Params(1)); 
        bs.addUserFloat("k1_refitted_vz",  phiDau1Params(2)); 
        bs.addUserFloat("k1_refitted_px", phiDau1Params(3)); 
        bs.addUserFloat("k1_refitted_py", phiDau1Params(4)); 
        bs.addUserFloat("k1_refitted_pz", phiDau1Params(5)); 
        bs.addUserFloat("k1_refitted_m",  phiDau1Params(6)); 

        bs.addUserFloat("k2_refitted_vx",  phiDau2Params(0)); 
        bs.addUserFloat("k2_refitted_vy",  phiDau2Params(1)); 
        bs.addUserFloat("k2_refitted_vz",  phiDau2Params(2)); 
        bs.addUserFloat("k2_refitted_px", phiDau2Params(3)); 
        bs.addUserFloat("k2_refitted_py", phiDau2Params(4)); 
        bs.addUserFloat("k2_refitted_pz", phiDau2Params(5)); 
        bs.addUserFloat("k2_refitted_m",  phiDau2Params(6)); 

        bs.addUserFloat("ds_fitted_px", dsParams(3)); 
        bs.addUserFloat("ds_fitted_py", dsParams(4)); 
        bs.addUserFloat("ds_fitted_pz", dsParams(5)); 
        bs.addUserFloat("ds_fitted_m",  dsParams(6)); 

        bs.addUserFloat("pi_refitted_vx",  dsDau1Params(0)); 
        bs.addUserFloat("pi_refitted_vy",  dsDau1Params(1)); 
        bs.addUserFloat("pi_refitted_vz",  dsDau1Params(2)); 
        bs.addUserFloat("pi_refitted_px", dsDau1Params(3)); 
        bs.addUserFloat("pi_refitted_py", dsDau1Params(4)); 
        bs.addUserFloat("pi_refitted_pz", dsDau1Params(5)); 
        bs.addUserFloat("pi_refitted_m",  dsDau1Params(6)); 

        bs.addUserFloat("phi_refitted_vx",  dsDau2Params(0)); 
        bs.addUserFloat("phi_refitted_vy",  dsDau2Params(1)); 
        bs.addUserFloat("phi_refitted_vz",  dsDau2Params(2)); 
        bs.addUserFloat("phi_refitted_px", dsDau2Params(3)); 
        bs.addUserFloat("phi_refitted_py", dsDau2Params(4)); 
        bs.addUserFloat("phi_refitted_pz", dsDau2Params(5)); 
        bs.addUserFloat("phi_refitted_m",  dsDau2Params(6)); 

        bs.addUserFloat("bs_fitted_px", bsParams(3)); 
        bs.addUserFloat("bs_fitted_py", bsParams(4)); 
        bs.addUserFloat("bs_fitted_pz", bsParams(5)); 
        bs.addUserFloat("bs_fitted_m",  bsParams(6)); 

        bs.addUserFloat("mu_refitted_vx",  bsDau1Params(0)); 
        bs.addUserFloat("mu_refitted_vy",  bsDau1Params(1)); 
        bs.addUserFloat("mu_refitted_vz",  bsDau1Params(2)); 
        bs.addUserFloat("mu_refitted_px", bsDau1Params(3)); 
        bs.addUserFloat("mu_refitted_py", bsDau1Params(4)); 
        bs.addUserFloat("mu_refitted_pz", bsDau1Params(5)); 
        bs.addUserFloat("mu_refitted_m",  bsDau1Params(6)); 

        bs.addUserFloat("ds_refitted_vx",  bsDau2Params(0)); 
        bs.addUserFloat("ds_refitted_vy",  bsDau2Params(1)); 
        bs.addUserFloat("ds_refitted_vz",  bsDau2Params(2)); 
        bs.addUserFloat("ds_refitted_px", bsDau2Params(3)); 
        bs.addUserFloat("ds_refitted_py", bsDau2Params(4)); 
        bs.addUserFloat("ds_refitted_pz", bsDau2Params(5)); 
        bs.addUserFloat("ds_refitted_m",  bsDau2Params(6)); 
        //std::cout << "12" << std::endl;

        /////////////////////// helicity angles in all possibe variations /////////////////////////
        // = angle between one of the kaons and the pi in the rest frame of the phi
  
        TLorentzVector fittedCollDs = fittedDs; //need ds twice

        // energy transferred
        TLorentzVector q = fittedBs - fittedDs; //sign is irrelevant, only needed for q2 in next line
        double q2 = q.M2(); 
        TLorentzVector qColl = collBs - fittedCollDs; //sign is irrelevant, only needed for q2 in next line
        double q2Coll = qColl.M2(); 

        //boost vector of the phi rest frame
        TVector3 boostPhi = fittedPhi.BoostVector();
        //boost vector of the Ds rest frame
        TVector3 boostDs = fittedDs.BoostVector();
        //boost vector of the bs rest frame
        TVector3 boostBs = fittedBs.BoostVector();
        //check the collinear approx as well!
        TVector3 boostCollBs = collBs.BoostVector();

        TLorentzVector refittedPi2 = refittedPi; //need pion twice
        TLorentzVector refittedMuColl = refittedMu; //need muon twice

        //TODO: get refitted phi mass from (refittedK1 + refittedK2).Mag() ?

        //boost kaons and pi into the rest frame of the phi
        refittedK1.Boost(-boostPhi);
        refittedK2.Boost(-boostPhi);
        refittedPi.Boost(-boostPhi);

        //boost Ds in the rest frame of the Bs, once via coll approx
        fittedDs.Boost(-boostBs);
        fittedCollDs.Boost(-boostCollBs);

        //boost muon in the rest frame of the W
        //for this we have to define a TLorentzVector for the W
        TLorentzVector w;
        w.SetVectM(-fittedDs.Vect(),std::sqrt(q2)); //the W is back to back with the Ds in the rest frame of the Bs and has mass q2
        TVector3 boostW= w.BoostVector();       
 
        refittedMu.Boost(-boostW);

        TLorentzVector wColl;
        wColl.SetVectM(-fittedCollDs.Vect(),std::sqrt(q2Coll)); //the W is back to back with the Ds in the rest frame of the Bs and has mass q2
        TVector3 boostWColl= wColl.BoostVector();       
 
        refittedMuColl.Boost(-boostWColl);

        //boost phi and pi into rest frame of Ds
        fittedPhi.Boost(-boostDs);
        refittedPi2.Boost(-boostDs);
        //std::cout << "13" << std::endl;

 
        //definition of all helicity angles
        float angPiK1  = refittedK1.Angle(refittedPi.Vect()); 
        float angPiK2  = refittedK2.Angle(refittedPi.Vect()); 
        float angMuW   = refittedMu.Angle(w.Vect()); 
        float angMuWColl   = refittedMuColl.Angle(wColl.Vect()); 
        float angPhiDs = fittedDs.Angle(fittedPhi.Vect()); 
        float angPiDs  = fittedDs.Angle(refittedPi2.Vect()); 

        float cosPiK1  = cos(angPiK1);
        float cosPiK2  = cos(angPiK2);
        float cosMuW   = cos(angMuW);
        float cosMuWColl   = cos(angMuWColl);
        float cosPhiDs = cos(angPhiDs);
        float cosPiDs  = cos(angPiDs);

        //add them to save
        bs.addUserFloat("angPiK1",  angPiK1); 
        bs.addUserFloat("angPiK2",  angPiK2); 
        bs.addUserFloat("angMuW",   angMuW); 
        bs.addUserFloat("angPhiDs", angPhiDs); 
        bs.addUserFloat("angPiDs",  angPiDs); 

        bs.addUserFloat("cosPiK1",  cosPiK1); 
        bs.addUserFloat("cosPiK2",  cosPiK2); 
        bs.addUserFloat("cosMuW",   cosMuW); 
        bs.addUserFloat("cosPhiDs", cosPhiDs); 
        bs.addUserFloat("cosPiDs",  cosPiDs); 



        arrived = 1;
        bs.addUserInt("arrived", arrived);
        //std::cout << "14" << std::endl;

        //append candidate at the end of our return value :)
        //ret_value can be a vector!!
        ret_value->emplace_back(bs);
        //ret_value_gen->emplace_back(gen);

        //std::cout << "15" << std::endl;
        if(arrived >0) break;
        //std::cout << "saved!" << std::endl; 
        } //closing pi loop
      if(arrived >0) break;
      } //closing k2 loop
    if(arrived >0) break;
    } //closing k1 loop
    if(arrived >0) break;

  } //closing mu loop
  //move and store these two new collections in the event 
  //iEvent.put(std::move(ret_value), "bs");
  //evt.put(std::move(dimuon_tt), "KKPiTransientTracks");
  if(arrived >0){
  iEvent.put(std::move(ret_value), "bs");
  }
  else{
  pat::CompositeCandidate bs;
  bs.addUserInt("arrived",arrived);
  ret_value->emplace_back(bs);
  iEvent.put(std::move(ret_value), "bs"); 
  }
}//clsoing event loop

DEFINE_FWK_MODULE(BsToDsPhiKKPiMuBuilder);
