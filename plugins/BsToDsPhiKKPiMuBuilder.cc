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
//for 3DPoint
#include "DataFormats/Math/interface/Point3D.h"
////////////////////////////////////////////////////////////////////////////////////////////
// TODOS:
//
// - move gen matching in separate module? NO
// - remove hardcoded numbers - DONE
// - helicity plane angles    
// - redefine all variables after the fit? save both ? -YES 
// - beautify the bs.addUserFloat (...) - DONE
// - pos. def. cov matrix
// - add counters before every selection  
// - pruned vs packed -> discuss - DONE
// - output tree has now empty entries when there is no trigger/signal -> how to avoid this? DONE (ed filter)
// - adapt for Hb background sample
// - how to save kk same sign pair? --> DONE
// - generally: save only gen matched signals? --> NO
// - what if an event has two signals?  -> SAVE BOTH!
// - divide into submitter chunks DONE
// - save gen information!! DONE 
// - do gen tests, check f.e. refitted p's with gen p's and unfitted p's DONE
// - put hel angle calc. etc into functions! DONE
// - isAncestor and getAncestor are save, they dont modify the Ptr - CHECKED
///////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
// function which checks if mom is ancestor of dau

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
///////////////////////////////////////////////////////////////////////////////////
// function which returns pt eta phi of the ancestor in order to compare ancestors.

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

///////////////////////////////////////////////////////////////////////////////////
// function which returns pointer to ancestor such that we can match by pointer! :)

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

///////////////////////////////////////////////////////////////////////////////////
// function which removes the un-oscillated ancestor of dau
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
///////////////////////////////////////////////////////////////////////////////////
//function which gets the hel angle between the mu and W
float angMuW (TLorentzVector d, TLorentzVector b, TLorentzVector mu){       

  //get q2
  TLorentzVector q = b - d;
  double q2 = q.M2();

  //boost Ds into Bs rest frame
  TVector3 bBoost = b.BoostVector();
  d.Boost(-bBoost);

  //get W via Ds 
  TLorentzVector w;
  w.SetVectM(-d.Vect(),std::sqrt(q2));           
  // boost it back into lab frame
  w.Boost(bBoost);
  // now take the boost vector of w
  TVector3 wBoost = w.BoostVector();
  //boost the muon into the w rest frame
  mu.Boost(-wBoost); 
  //boost the W back into the bs rest frame
  w.Boost(-bBoost);
 
  //now take the angle
  return w.Angle(mu.Vect());

}
///////////////////////////////////////////////////////////////////////////////////
//function which gets the hel angle between dau1 and dau2 in the rest frame of the restFrame particle

float angDoubleDecay (TLorentzVector restFrame, TLorentzVector dau1, TLorentzVector dau2) {
  TVector3 restBoost = restFrame.BoostVector();
  dau1.Boost(-restBoost);
  dau2.Boost(-restBoost);
  return dau1.Angle(dau2.Vect());
}

///////////////////////////////////////////////////////////////////////////////////
//function which gets angle between decay planes

float angPlane (TLorentzVector d, TLorentzVector b, TLorentzVector mu, TLorentzVector pi) {

  //get q2
  TLorentzVector q = b - d;
  double q2 = q.M2();

  //get Ds boost vector in the lab frame (to boost pi)
  TVector3 dBoost = d.BoostVector();

  //boost Ds into Bs rest frame
  TVector3 bBoost = b.BoostVector();
  d.Boost(-bBoost);

  //get W via Ds 
  TLorentzVector w;
  w.SetVectM(-d.Vect(),std::sqrt(q2));           
  // boost it back into lab frame
  w.Boost(bBoost);
  // now take the boost vector of w
  TVector3 wBoost = w.BoostVector();
  //boost the muon into the w rest frame
  mu.Boost(-wBoost);
  //boost the W back into the bs rest frame
  w.Boost(-bBoost);

  //now boost the pi into the ds rest frame, ds is already in Bs rest frame
  pi.Boost(-dBoost);

  // normal vector on lepton-W plane
  TVector3 n1 = mu.Vect().Cross(w.Vect());
  // normal vector on ds-pi plane
  TVector3 n2 = d.Vect().Cross(pi.Vect()); 

  return n1.Angle(n2);
}

///////////////////////////////////////////////////////////////////////////////////
//function which gets angle of kaon in DsPi plane 

float angPlane2 (TLorentzVector d, TLorentzVector b, TLorentzVector k, TLorentzVector pi) {

  //get Ds boost vector in the lab frame (to boost pi)
  TVector3 dBoost = d.BoostVector();

  //now boost the pi and k1 into the ds rest frame
  pi.Boost(-dBoost);
  k.Boost(-dBoost);

  //boost Ds into Bs rest frame
  TVector3 bBoost = b.BoostVector();
  d.Boost(-bBoost);

  // normal vector on Ds - pi plane
  TVector3 n1 = d.Vect().Cross(pi.Vect());
  // normal vector on k1 - pi plane
  TVector3 n2 = k.Vect().Cross(pi.Vect()); 

  return n1.Angle(n2);
}

///////////////////////////////////////////////////////////////////////////////////
//function which returns the phi difference of two phi variables (in cms coordinate system)
double phiDiff(double phi1, double phi2){

  double dPhi = phi1 - phi2;
  double pi = 3.14159265358979323846;
  while (fabs(dPhi) > pi) {
    int sgn = dPhi > 0? 1 : -1;
    dPhi -= sgn*2*pi;
  }
  return dPhi;
}


///////////////////////////////////////////////////////////////////////////////////
// counters for filters to know when we lose how many particles

int nEvents = 0;        // counts the nr of events in total
int nMuons  = 0;        // counts the nr of muons
int nPacked = 0;        // counts the nr of packed candidates in total

int nMuonsPassed = 0;   // how many muons pass (i.e. when we find a pv)

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

  for(size_t trgMuIdx = 0; trgMuIdx < trgMuons->size(); ++trgMuIdx){
    nMuons++;

    //if there is no trg muon, this loop is empty:)
    edm::Ptr<pat::Muon> muPtr(trgMuons, trgMuIdx);

    std::vector<float> dzMuPV;

    //Fix the primary vertex to be the one closest to the trg Muon in dz
    // more accurate for mu than for tau signals

    float dummy = 1.0;
    int goldenIdx = -1;
    reco::Vertex pv;
    for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){

      edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, vtxIdx);
      float minDz = fabs(muPtr->bestTrack()->dz(vtxPtr->position())); 
      if(minDz < dummy){
        dummy = minDz;
        goldenIdx = vtxIdx;
      }

    }
    // take as the primary vertex the one which has minimal dz with the muon 
    //auto dzMuPVMin = std::min_element(std::begin(dzMuPV),std::end(dzMuPV));
    //int pvIdx = std::distance(std::begin(dzMuPV), dzMuPVMin); 

    if (goldenIdx<0) continue;
    if (goldenIdx >= 0){
    pv = primaryVtx->at(goldenIdx);
    }

    //////////////////////////////////////////////////
    // Loop over k1 and select the good tracks      //
    //////////////////////////////////////////////////

    for(size_t k1Idx = 0; k1Idx < pcand->size(); ++k1Idx) {

      nPacked++;

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
        pat::CompositeCandidate kk;

        //PF canidates are always assigned with the pi mass, we have to force the kaon mass
        math::PtEtaPhiMLorentzVector k1P4(k1Ptr->pt(), k1Ptr->eta(), k1Ptr->phi(), K_MASS);
        math::PtEtaPhiMLorentzVector k2P4(k2Ptr->pt(), k2Ptr->eta(), k2Ptr->phi(), K_MASS);
        math::PtEtaPhiMLorentzVector piP4(piPtr->pt(), piPtr->eta(), piPtr->phi(), PI_MASS); //just to be sure lets also force the pi mass
 
        kk.setP4(k1P4 + k2P4);
      
        //only continue when they build a phi resonance, allow 15MeV:
        if (fabs(kk.mass() - phiMass_) > phiMassAllowance_) continue;     

        kk.setCharge(k1Ptr->charge() + k2Ptr->charge());

        //////////////////////////////////////////////////
        // Build Ds resonance                           //
        // ///////////////////////////////////////////////

        pat::CompositeCandidate phiPi;
        phiPi.setP4(kk.p4() + piP4); 

        //only continue when they build a ds resonance, allow 50MeV:
        if (fabs(phiPi.mass() - dsMass_) > dsMassAllowance_) continue;

        phiPi.setCharge(kk.charge() + piPtr->charge());

        //////////////////////////////////////////////////
        // Build Bs resonance                           //
        //////////////////////////////////////////////////

        pat::CompositeCandidate dsMu;
        dsMu.setP4(phiPi.p4() + muPtr->p4()); 

        dsMu.setCharge(phiPi.charge() + muPtr->charge()); //sanity check:shoould be 0

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

               //std::cout << "before getAncestor" << phiFromK1->pt() << std::endl;
               auto dsFromPhi = getAncestor(phiFromK1,431);
              //std::cout << "after getAncestor" << phiFromK1->pt() << std::endl;

               auto dsFromPi  = getAncestor(piReco,431);
               if( (dsFromPhi != dsFromPi) || (dsFromPhi == nullptr) || (dsFromPi == nullptr)) continue; 

               //std::cout << "searching for bs resonance .. " << std::endl;
               auto bsFromDs = getAncestor(dsFromPhi,531);
               auto bsFromMu = getAncestor(muReco,531);
               if( (bsFromDs != bsFromMu) || (bsFromDs == nullptr) || (bsFromMu == nullptr)) continue; 
 
               std::cout << "(px,vx) before oscillations " << bsFromMu->px() << ", " << bsFromMu->vx() << std::endl;

               //remove oscillations
               auto bsFromMuWOOsc = removeOscillations(bsFromMu);
               //bool related = hasAncestor(dsFromK1,bsFromK1);
               std::cout << "(px,vx) without oscillations " << bsFromMuWOOsc->px() << ", " << bsFromMuWOOsc->vx() << std::endl;
              
               nGenMatches++;
               genMatchSuccess = 1;

               if(nGenMatches > 1) continue; //std::cout <<"there is more than one match!!" << std::endl;

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
               genBsTlv.SetXYZM(  bsFromMuWOOsc->px(), bsFromMuWOOsc->py(), bsFromMuWOOsc->pz(), bsMass_); //changed

               genMissTlv = genBsTlv - (genDsTlv + genMuTlv); 
               genQTlv    = genBsTlv - (genDsTlv); 

               float m2_miss_gen = genMissTlv.M2();
               float q2_gen = genQTlv.M2();

               bs.addUserFloat("m2_miss_gen",m2_miss_gen);
               bs.addUserFloat("q2_gen",q2_gen);

               //vertices
               float pv_x_gen = bsFromMuWOOsc->vx();//This is the bs production vertex!
               float pv_y_gen = bsFromMuWOOsc->vy();
               float pv_z_gen = bsFromMuWOOsc->vz();

               float sv_x_gen = dsFromPi->vx(); //This is the ds production vertex!
               float sv_y_gen = dsFromPi->vy();
               float sv_z_gen = dsFromPi->vz();

               float tv_x_gen = phiFromK1->vx(); //This is the phi production vertex!
               float tv_y_gen = phiFromK1->vy();
               float tv_z_gen = phiFromK1->vz();

               float fv_x_gen = k1PtrGen->vx(); //This is the k1 production vertex!
               float fv_y_gen = k1PtrGen->vy();
               float fv_z_gen = k1PtrGen->vz();

               //save the gen info by adding gen candidates of final states
               bs.addUserCand("mu_gen",muPtrGen);
               bs.addUserCand("k1_gen",k1PtrGen);
               bs.addUserCand("k2_gen",k2PtrGen);
               bs.addUserCand("pi_gen",piPtrGen);


               //and gen info from the resonances
               bs.addUserFloat("phi_gen_px"     ,phiFromK1->px());
               bs.addUserFloat("phi_gen_py"     ,phiFromK1->py());
               bs.addUserFloat("phi_gen_pz"     ,phiFromK1->pz());
               bs.addUserFloat("phi_gen_pt"     ,phiFromK1->pt());
               bs.addUserFloat("phi_gen_eta"    ,phiFromK1->eta());
               bs.addUserFloat("phi_gen_phi"    ,phiFromK1->phi());
               bs.addUserFloat("tv_x_gen"       ,tv_x_gen);//This is the phi production vertex!
               bs.addUserFloat("tv_y_gen"       ,tv_y_gen);
               bs.addUserFloat("tv_z_gen"       ,tv_z_gen);
               bs.addUserFloat("phi_gen_charge" ,phiFromK1->charge());
               bs.addUserInt(  "phi_gen_pdgid"  ,phiFromK1->pdgId());

               bs.addUserFloat("ds_gen_px"     ,dsFromPi->px());
               bs.addUserFloat("ds_gen_py"     ,dsFromPi->py());
               bs.addUserFloat("ds_gen_pz"     ,dsFromPi->pz());
               bs.addUserFloat("ds_gen_pt"     ,dsFromPi->pt());
               bs.addUserFloat("ds_gen_eta"    ,dsFromPi->eta());
               bs.addUserFloat("ds_gen_phi"    ,dsFromPi->phi());
               bs.addUserFloat("sv_x_gen"      ,sv_x_gen);//This is the ds production vertex!
               bs.addUserFloat("sv_y_gen"      ,sv_y_gen);
               bs.addUserFloat("sv_z_gen"      ,sv_z_gen);
               bs.addUserFloat("ds_gen_charge" ,dsFromPi->charge());
               bs.addUserInt(  "ds_gen_pdgid"  ,dsFromPi->pdgId());

               bs.addUserFloat("bs_gen_px"     ,bsFromMuWOOsc->px());
               bs.addUserFloat("bs_gen_py"     ,bsFromMuWOOsc->py());
               bs.addUserFloat("bs_gen_pz"     ,bsFromMuWOOsc->pz());
               bs.addUserFloat("bs_gen_pt"     ,bsFromMuWOOsc->pt());
               bs.addUserFloat("bs_gen_eta"    ,bsFromMuWOOsc->eta());
               bs.addUserFloat("bs_gen_phi"    ,bsFromMuWOOsc->phi());
               bs.addUserFloat("pv_x_gen"      ,pv_x_gen); //This is the bs production vertex!
               bs.addUserFloat("pv_y_gen"      ,pv_y_gen);
               bs.addUserFloat("pv_z_gen"      ,pv_z_gen);
               bs.addUserFloat("bs_gen_charge" ,bsFromMuWOOsc->charge());
               bs.addUserInt(  "bs_gen_pdgid"  ,bsFromMuWOOsc->pdgId());

               //lets also store the fourth vertex ( the k production vertex)
               bs.addUserFloat("fv_x_gen"       ,fv_x_gen);//This is the kaon production vertex!
               bs.addUserFloat("fv_y_gen"       ,fv_y_gen);
               bs.addUserFloat("fv_z_gen"       ,fv_z_gen);

               //define vertex variables
 
               float lxyBsGen   = std::sqrt(std::pow((pv_x_gen - sv_x_gen),2) + std::pow((pv_y_gen - sv_y_gen),2) ); 
               float lxyzBsGen  = std::sqrt(std::pow((pv_x_gen - sv_x_gen),2) + std::pow((pv_y_gen - sv_y_gen),2) + std::pow((pv_z_gen - sv_z_gen),2) ); 
       
               float lxyDsGen   = std::sqrt(std::pow((sv_x_gen - tv_x_gen),2) + std::pow((sv_y_gen - tv_y_gen),2) ); 
               float lxyzDsGen  = std::sqrt(std::pow((sv_x_gen - tv_x_gen),2) + std::pow((sv_y_gen - tv_y_gen),2) + std::pow((sv_z_gen - tv_z_gen),2) ); 
       
               float lxyPhiGen  = std::sqrt(std::pow((tv_x_gen - fv_x_gen),2) + std::pow((tv_y_gen - fv_y_gen),2) ); 
               float lxyzPhiGen = std::sqrt(std::pow((tv_x_gen - fv_x_gen),2) + std::pow((tv_y_gen - fv_y_gen),2) + std::pow((tv_z_gen - fv_z_gen),2) ); 

               /* 
               bs.addUserFloat("lxy_bs_gen"   ,lxyBsGen);
               bs.addUserFloat("lxyz_bs_gen"  ,lxyzBsGen);
       
               bs.addUserFloat("lxy_ds_gen"   ,lxyDsGen);
               bs.addUserFloat("lxyz_ds_gen"  ,lxyzDsGen);
       
               bs.addUserFloat("lxy_phi_gen"  ,lxyPhiGen);
               bs.addUserFloat("lxyz_phi_gen" ,lxyzPhiGen);

               math::XYZPoint pvGen(pv_x_gen,pv_y_gen,pv_z_gen); 
              
                
               float dxyMuGen    = muPtrGen->bestTrack()->dxy(pvGen);  
               //float dxyMuErrGen = muPtrGen->bestTrack()->dxyError(pvGen,pvGen.covariance());  
               //float dxyMuSigGen = dxyMuGen/dxyMuErrGen;
       
               float dzMuGen     = muPtrGen->bestTrack()->dz(pvGen);  
               //float dzMuErrGen  = muPtrGen->bestTrack()->dzError();  
               //float dzMuSigGen  = dzMuGen/dzMuErrGen ; 
       
               float dxyPiGen    = piPtrGen->bestTrack()->dxy(pvGen);  //maybe useful for Ds* vs Ds ? 
               //float dxyPiErrGen = piPtrGen->bestTrack()->dxyError(pvGen,pvGen.covariance());  
               //float dxyPiSigGen = dxyPiGen/dxyPiErrGen;
       
               float dzPiGen     = piPtrGen->bestTrack()->dz(pvGen);  
               //float dzPiErrGen  = piPtrGen->bestTrack()->dzError();  
               //float dzPiSigGen  = dzPiGen/dzPiErrGen ; 
       
               float dxyK1Gen    = k1PtrGen->bestTrack()->dxy(pvGen); //needed ? 
               //float dxyK1ErrGen = k1PtrGen->bestTrack()->dxyError(pvGen,pvGen.covariance());
               //float dxyK1SigGen = dxyK1Gen/dxyK1ErrGen;
       
               float dzK1Gen     = k1PtrGen->bestTrack()->dz(pvGen);  
               //float dzK1ErrGen  = k1PtrGen->bestTrack()->dzError();  
               //float dzK1SigGen  = dzK1Gen/dzK1ErrGen ; 
       
               float dxyK2Gen    = k2PtrGen->bestTrack()->dxy(pvGen); //needed ? 
               //float dxyK2ErrGen = k2PtrGen->bestTrack()->dxyError(pvGen,pvGen.covariance());
               //float dxyK2SigGen = dxyK2Gen/dxyK2ErrGen;
       
               float dzK2Gen     = k2PtrGen->bestTrack()->dz(pvGen);  
               //float dzK2ErrGen  = k2PtrGen->bestTrack()->dzError();  
               //float dzK2SigGen  = dzK2Gen/dzK2ErrGen ; 
       
               bs.addUserFloat("dxy_mu_gen",     dxyMuGen);
               bs.addUserFloat("dz_mu_gen",      dzMuGen);
               //bs.addUserFloat("dxy_mu_err_gen", dxyMuErrGen);
               //bs.addUserFloat("dz_mu_err_gen",  dzMuErrGen);
               //bs.addUserFloat("dxy_mu_sig_gen", dxyMuSigGen);
               //bs.addUserFloat("dz_mu_sig_gen",  dzMuSigGen);
       
               bs.addUserFloat("dxy_pi_gen",     dxyPiGen);
               bs.addUserFloat("dz_pi_gen",      dzPiGen);
               //bs.addUserFloat("dxy_pi_err_gen", dxyPiErrGen);
               //bs.addUserFloat("dz_pi_err_gen",  dzPiErrGen);
               //bs.addUserFloat("dxy_pi_sig_gen", dxyPiSigGen);
               //bs.addUserFloat("dz_pi_sig_gen",  dzPiSigGen);
       
               bs.addUserFloat("dxy_k1_gen",     dxyK1Gen);
               bs.addUserFloat("dz_k1_gen",      dzK1Gen);
               //bs.addUserFloat("dxy_k1_err_gen", dxyK1ErrGen);
               //bs.addUserFloat("dz_k1_err_gen",  dzK1ErrGen);
               //bs.addUserFloat("dxy_k1_sig_gen", dxyK1SigGen);
               //bs.addUserFloat("dz_k1_sig_gen",  dzK1SigGen);
       
               bs.addUserFloat("dxy_k2_gen",     dxyK2Gen);
               bs.addUserFloat("dz_k2_gen",      dzK2Gen);
               //bs.addUserFloat("dxy_k2_err_gen", dxyK2ErrGen);
               //bs.addUserFloat("dz_k2_err_gen",  dzK2ErrGen);
               //bs.addUserFloat("dxy_k2_sig_gen", dxyK2SigGen);
               //bs.addUserFloat("dz_k2_sig_gen",  dzK2SigGen);
               */
 
               // now find the signal ID
               // Ds mu   = 0
               // Ds* mu  = 1
               // Ds tau  = 2
               // Ds* tau = 3
               sigId = 0; //default

               // look for tau in signal 
               //std::cout << "before isAncestor" << muPtrGen->pt() << std::endl;
               if(isAncestor(muPtrGen, 15)) sigId = 2;               
               //std::cout << "after isAncestor" << muPtrGen->pt() << std::endl;

               // now look for possible Ds*
               //std::cout << "before isAncestor" << piPtrGen->pt() << std::endl;
               if(isAncestor(piPtrGen, 433)) sigId += 1; 
               //std::cout << "after isAncestor" << piPtrGen->pt() << std::endl;

               //define helicity angles

               //angle between Mu and W
               float angMuWGen = angMuW(genDsTlv,genBsTlv,genMuTlv); 
               bs.addUserFloat("angMuWGen",angMuWGen);
               bs.addUserFloat("cosMuWGen",cos(angMuWGen));

               //angle between k1(k2) and pion in phi rest frame
               float angPiK1Gen  = angDoubleDecay(genPhiTlv, genK1Tlv,  genPiTlv);
               float angPiK2Gen  = angDoubleDecay(genPhiTlv, genK2Tlv,  genPiTlv);
               bs.addUserFloat("angPiK1Gen", angPiK1Gen);
               bs.addUserFloat("angPiK2Gen", angPiK2Gen);
               bs.addUserFloat("cosPiK1Gen", cos(angPiK1Gen));
               bs.addUserFloat("cosPiK2Gen", cos(angPiK2Gen));

               // equivalently, angle between phi(pi) and bs in ds rest frame
               float angPhiDsGen = angDoubleDecay(genDsTlv,  genPhiTlv, genBsTlv);
               float angPiDsGen  = angDoubleDecay(genDsTlv,  genPiTlv,  genBsTlv);
               bs.addUserFloat("angPhiDsGen", angPhiDsGen);
               bs.addUserFloat("angPiDsGen",  angPiDsGen);
               bs.addUserFloat("cosPhiDsGen", cos(angPhiDsGen));
               bs.addUserFloat("cosPiDsGen",  cos(angPiDsGen));

               //plane angle
               float angPlaneBsGen = angPlane(genDsTlv, genBsTlv, genMuTlv, genPiTlv);
               bs.addUserFloat("angPlaneBsGen", angPlaneBsGen);
               bs.addUserFloat("cosPlaneBsGen", cos(angPlaneBsGen));

               float angPlaneDsGen = angPlane2(genDsTlv, genBsTlv, genK1Tlv, genPiTlv);
               bs.addUserFloat("angPlaneDsGen", angPlaneDsGen);
               bs.addUserFloat("cosPlaneDsGen", cos(angPlaneDsGen));

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

        phiToFit.push_back(pFactory.particle(getTransientTrack(k1Ptr->bestTrack()),kMass, chi,ndf,sigma));
        phiToFit.push_back(pFactory.particle(getTransientTrack(k2Ptr->bestTrack()),kMass, chi,ndf,sigma));
        dsToFit.push_back(pFactory.particle(getTransientTrack(piPtr->bestTrack()), piMass,chi,ndf,sigma));
        bsToFit.push_back(pFactory.particle(getTransientTrack(muPtr->bestTrack()), muMass,chi,ndf,sigma));

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
        bs.addUserFloat("mu_mass",muMass_);

        bs.addUserFloat("sig",sigId);

        //add prefit resonances --> are not part of collection and can thus not be
        //automatically access the pt(), .. as I can for k1,k2,pi,mu in the variables_cff.py

        bs.addUserFloat("kk_pt", kk.pt());
        bs.addUserFloat("kk_eta", kk.eta());
        bs.addUserFloat("kk_phi", kk.phi());
        bs.addUserFloat("kk_mass", kk.mass());
        bs.addUserFloat("kk_charge", kk.charge());
        bs.addUserInt("kkCharge",kkCharge); 
        bs.addUserFloat("kk_deltaR", reco::deltaR(*k1Ptr, *k2Ptr));

        bs.addUserFloat("phiPi_pt", phiPi.pt());
        bs.addUserFloat("phiPi_eta", phiPi.eta());
        bs.addUserFloat("phiPi_phi", phiPi.phi());
        bs.addUserFloat("phiPi_mass", phiPi.mass());
        bs.addUserFloat("phiPi_charge", phiPi.charge());
        bs.addUserFloat("phiPi_deltaR", reco::deltaR(kk, *piPtr));

        bs.addUserFloat("dsMu_pt", dsMu.pt());
        bs.addUserFloat("dsMu_eta", dsMu.eta());
        bs.addUserFloat("dsMu_phi", dsMu.phi());
        bs.addUserFloat("dsMu_mass", dsMu.mass()); //we miss the neutrino mass
        bs.addUserFloat("dsMu_charge", dsMu.charge());
        bs.addUserFloat("dsMu_deltaR", reco::deltaR(phiPi, *muPtr));

        //rel charges
        bs.addUserInt("kk_charge",kkCharge); 
        bs.addUserInt("pi_mu_charge",piMuCharge); 


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

        ///////////////////////// reconstruct the B momentum with different methods ///////////////////////////////
        

        //Define 4 vectors of fitted and refitted resonances
        TLorentzVector fittedPhi;
        TLorentzVector fittedDs;
        TLorentzVector fittedBs;

        TLorentzVector refittedPhi;
        TLorentzVector refittedDs;
        // bs does not have a refitted state :)

        fittedPhi.SetXYZM(   phiParams(3),     phiParams(4),     phiParams(5),      phiParams(6));
        fittedDs.SetXYZM(    dsParams(3),      dsParams(4),      dsParams(5),       dsParams(6));
        fittedBs.SetXYZM(    bsParams(3),      bsParams(4),      bsParams(5),       bsParams(6));

        refittedPhi.SetXYZM( dsDau2Params(3),  dsDau2Params(4),  dsDau2Params(5),   dsDau2Params(6));
        refittedDs.SetXYZM(  bsDau2Params(3),  bsDau2Params(4),  bsDau2Params(5),   bsDau2Params(6));

        //Define 4 vectors of refitted final states
        TLorentzVector refittedK1;
        TLorentzVector refittedK2;
        TLorentzVector refittedPi;
        TLorentzVector refittedMu;

        refittedK1.SetXYZM(  phiDau1Params(3), phiDau1Params(4), phiDau1Params(5),  phiDau1Params(6));
        refittedK2.SetXYZM(  phiDau2Params(3), phiDau2Params(4), phiDau2Params(5),  phiDau2Params(6));
        refittedPi.SetXYZM(  dsDau1Params(3),  dsDau1Params(4),  dsDau1Params(5),   dsDau1Params(6));
        refittedMu.SetXYZM(  bsDau1Params(3),  bsDau1Params(4),  bsDau1Params(5),   bsDau1Params(6));

        /////////////////////////
        // Collinear approx.   //
        ///////////////////////// 

        TLorentzVector collBsTlv = fittedDs + refittedMu;

        TLorentzVector collMissTlv; // for m2 miss
        TLorentzVector collQTlv;    // for q2 miss

        double refittedDsMuMass = collBsTlv.M();
        collBsTlv *= bsMass_ / refittedDsMuMass; //scale it

        collMissTlv = collBsTlv - (fittedDs + refittedMu); //bs - ds+mu
        collQTlv = collBsTlv - fittedDs; // bs - ds
        double m2_miss_coll = collMissTlv.M2();
        double q2_coll = collQTlv.M2();

        bs.addUserFloat("bs_px_coll",collBsTlv.Px());
        bs.addUserFloat("bs_py_coll",collBsTlv.Py());
        bs.addUserFloat("bs_pz_coll",collBsTlv.Pz());
        bs.addUserFloat("bs_pt_coll",collBsTlv.Pt());
        bs.addUserFloat("bs_eta_coll",collBsTlv.Eta());
        bs.addUserFloat("bs_phi_coll",collBsTlv.Phi());
        bs.addUserFloat("m2_miss_coll",m2_miss_coll);
        bs.addUserFloat("q2_coll",q2_coll);

        //////////////////////////
        // LHCb method          //
        //////////////////////////

        TVector3 bsFlightDir;
        TVector3 beamAxis;
        TVector3 radialAxis;

        TLorentzVector lhcbQTlv; // for q2
        TLorentzVector lhcbMissTlv; // for m2 miss

        dsMuTlv.SetXYZM(dsMu.px(),dsMu.py(),dsMu.pz(),dsMu.mass());       

        //bsFlightDir.SetXYZ(sv_x - pv_x, sv_y - pv_y, sv_z - pv_z);
        bsFlightDir.SetXYZ(sv_x - pv_x, sv_y - pv_y , sv_z - pv_z);

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

        TLorentzVector lhcbBsTlv; 
        double eta = - std::log(std::tan(theta/2));
        lhcbBsTlv.SetPtEtaPhiM(lhcbPt,eta,lhcbPhi,bsMass_); 
        
        lhcbMissTlv = lhcbBsTlv - refittedDsMu;
        lhcbQTlv = lhcbBsTlv - fittedDs;

        double m2_miss_lhcb = lhcbMissTlv.M2();
        double q2_lhcb = lhcbQTlv.M2();

        bs.addUserFloat("bs_px_lhcb",lhcbBsTlv.Px());
        bs.addUserFloat("bs_py_lhcb",lhcbBsTlv.Py());
        bs.addUserFloat("bs_pz_lhcb",lhcbBsTlv.Pz());
        bs.addUserFloat("bs_pt_lhcb",lhcbBsTlv.Pt());
        bs.addUserFloat("bs_eta_lhcb",lhcbBsTlv.Eta());
        bs.addUserFloat("bs_phi_lhcb",lhcbBsTlv.Phi());

        bs.addUserFloat("m2_miss_lhcb",m2_miss_lhcb);
        bs.addUserFloat("q2_lhcb",q2_lhcb);


        ///////////////////////////////////////////
        // Another Lhcb method                   //
        ///////////////////////////////////////////

        // reset bs flight direction
        bsFlightDir.SetXYZ(sv_x - pv_x, sv_y - pv_y , sv_z - pv_z);

        TVector3 lhcbAltBs;
        TLorentzVector lhcbAltBsTlv;

        lhcbAltBs = bsFlightDir.Unit();  
        lhcbAltBs *= refittedDsMu.Vect().Mag() * bsMass_ / refittedDsMu.M(); 
        lhcbAltBsTlv.SetXYZM(lhcbAltBs.Px(),lhcbAltBs.Py(),lhcbAltBs.Pz(), bsMass_);      

        TLorentzVector lhcbAltMissTlv = lhcbAltBsTlv - refittedDsMu; // bs - ds+mu
        TLorentzVector lhcbAltQTlv = lhcbAltBsTlv - fittedDs; // bs - ds 

        bs.addUserFloat("bs_px_lhcb_alt",lhcbAltBsTlv.Px());
        bs.addUserFloat("bs_py_lhcb_alt",lhcbAltBsTlv.Py());
        bs.addUserFloat("bs_pz_lhcb_alt",lhcbAltBsTlv.Pz());
        bs.addUserFloat("bs_pt_lhcb_alt",lhcbAltBs.Pt());
        bs.addUserFloat("bs_eta_lhcb_alt",lhcbAltBs.Eta());
        bs.addUserFloat("bs_phi_lhcb_alt",lhcbAltBs.Phi());

        bs.addUserFloat("m2_miss_lhcb_alt",lhcbAltMissTlv.M2());
        bs.addUserFloat("q2_lhcb_alt",lhcbAltQTlv.M2());
        
        ///////////////////////////////////////////////////////////
        // RECO METHOD -> 2 Solutions                            //
        // This method is exact in the 1-neutrino final state :) //
        ///////////////////////////////////////////////////////////

        TLorentzVector recoBsTlv1;
        TLorentzVector recoBsTlv2;

        TLorentzVector miss_1; // for m2 miss
        TLorentzVector miss_2; // "
        TLorentzVector Q_1; // for q2
        TLorentzVector Q_2; // "

        double recoNeutPll_1; // neutrino momentum parallel to bs direction
        double recoNeutPll_2; // "

        double recoBsAbs_1; // absolute 3 momentum of bs
        double recoBsAbs_2; // "

        // reset bs flight direction
        bsFlightDir.SetXYZ(sv_x - pv_x, sv_y - pv_y , sv_z - pv_z);
 
        // angle between the bs flight direction and the Dsmu system (visible)
        double recoAngle = refittedDsMu.Angle(bsFlightDir); 
                    
        //define parameters 
        // momentum of DsMu system parallel to the bs
        double recoDsMuPll = std::cos(recoAngle) * refittedDsMu.Vect().Mag();
        // momentum of DsMu system orthogonal to the bs
        double recoDsMuT = std::sin(recoAngle) * refittedDsMu.Vect().Mag();
        // energy of Dsmu system
        double recoDsMuE = refittedDsMu.E(); 
        // cocktail, drops out of equation
        double recoMix = std::pow(bsMass_,2) + std::pow(recoDsMuPll,2) - std::pow(recoDsMuT,2) - std::pow(recoDsMuE,2);

        // define a,b,c, to give to mitternachtsformel
        double a = 4*(std::pow(recoDsMuPll,2) - std::pow(recoDsMuE,2));
        double b = 4*recoDsMuPll*recoMix;
        double c = std::pow(recoMix,2) - 4*std::pow(recoDsMuE,2)*std::pow(recoDsMuT,2);
        // discriminant
        double disc = std::pow(b,2) - 4*a*c;      

        //prepare m2miss and q2 variables
        double m2_miss_reco_1;
        double q2_reco_1;
        double m2_miss_reco_2;
        double q2_reco_2;

        if( disc >= 0) {

          //non complex root -> nice! 
          recoNeutPll_1 = (-b + std::sqrt(disc)) / (2*a);
          recoNeutPll_2 = (-b - std::sqrt(disc)) / (2*a);

          recoBsAbs_1 = recoDsMuPll + recoNeutPll_1; 
          recoBsAbs_2 = recoDsMuPll + recoNeutPll_2; 

          TVector3 recoBs_1 = bsFlightDir.Unit();
          TVector3 recoBs_2 = bsFlightDir.Unit();

          recoBs_1 *= recoBsAbs_1; 
          recoBs_2 *= recoBsAbs_2; 

          recoBsTlv1.SetXYZM(recoBs_1.Px(),recoBs_1.Py(),recoBs_1.Pz(),bsMass_);
          recoBsTlv2.SetXYZM(recoBs_2.Px(),recoBs_2.Py(),recoBs_2.Pz(),bsMass_);

          miss_1 = recoBsTlv1 - refittedDsMu;  
          miss_2 = recoBsTlv2 - refittedDsMu;  
          Q_1 = recoBsTlv1 - fittedDs;  
          Q_2 = recoBsTlv2 - fittedDs;  

          m2_miss_reco_1 = miss_1.M2();
          m2_miss_reco_2 = miss_2.M2();
          q2_reco_1 = Q_1.M2();
          q2_reco_2 = Q_2.M2();


        }
        else{ 
        //complex root, save nans
        recoBsTlv1.SetXYZM(std::nan(""),std::nan(""),std::nan(""),std::nan(""));
        recoBsTlv2.SetXYZM(std::nan(""),std::nan(""),std::nan(""),std::nan(""));

        m2_miss_reco_1 = std::nan("");
        m2_miss_reco_2 = std::nan("");
        q2_reco_1 = std::nan("");
        q2_reco_2 = std::nan("");
        }

        bs.addUserFloat("bs_px_reco_1",recoBsTlv1.Px());
        bs.addUserFloat("bs_py_reco_1",recoBsTlv1.Py());
        bs.addUserFloat("bs_pz_reco_1",recoBsTlv1.Pz());
        bs.addUserFloat("bs_pt_reco_1",recoBsTlv1.Pt()); 
        bs.addUserFloat("bs_eta_reco_1",recoBsTlv1.Eta()); 
        bs.addUserFloat("bs_phi_reco_1",recoBsTlv1.Phi()); 

        bs.addUserFloat("bs_px_reco_2",recoBsTlv2.Px());
        bs.addUserFloat("bs_py_reco_2",recoBsTlv2.Py());
        bs.addUserFloat("bs_pz_reco_2",recoBsTlv2.Pz());
        bs.addUserFloat("bs_pt_reco_2",recoBsTlv2.Pt()); 
        bs.addUserFloat("bs_eta_reco_2",recoBsTlv2.Eta()); 
        bs.addUserFloat("bs_phi_reco_2",recoBsTlv2.Phi()); 

        bs.addUserFloat("m2_miss_reco_1",m2_miss_reco_1); 
        bs.addUserFloat("m2_miss_reco_2",m2_miss_reco_2); 
        bs.addUserFloat("q2_reco_1",q2_reco_1); 
        bs.addUserFloat("q2_reco_2",q2_reco_2); 

        ////////////////// Save momenta (and masses for consistency, even if constrained /////////////////////////

        // Phi fit
        bs.addUserFloat("k1_refitted_px",   refittedK1.Px()); 
        bs.addUserFloat("k1_refitted_py",   refittedK1.Py()); 
        bs.addUserFloat("k1_refitted_pz",   refittedK1.Pz()); 
        bs.addUserFloat("k1_refitted_pt",   refittedK1.Pt()); 
        bs.addUserFloat("k1_refitted_eta",  refittedK1.Eta()); 
        bs.addUserFloat("k1_refitted_phi",  refittedK1.Phi()); 
        bs.addUserFloat("k1_refitted_m",    refittedK1.M()); 

        bs.addUserFloat("k2_refitted_px",   refittedK2.Px()); 
        bs.addUserFloat("k2_refitted_py",   refittedK2.Py()); 
        bs.addUserFloat("k2_refitted_pz",   refittedK2.Pz()); 
        bs.addUserFloat("k2_refitted_pt",   refittedK2.Pt()); 
        bs.addUserFloat("k2_refitted_eta",  refittedK2.Eta()); 
        bs.addUserFloat("k2_refitted_phi",  refittedK2.Phi()); 
        bs.addUserFloat("k2_refitted_m",    refittedK2.M()); 
 
        bs.addUserFloat("phi_fitted_px",    fittedPhi.Px()); 
        bs.addUserFloat("phi_fitted_py",    fittedPhi.Py()); 
        bs.addUserFloat("phi_fitted_pz",    fittedPhi.Pz()); 
        bs.addUserFloat("phi_fitted_pt",    fittedPhi.Pt()); 
        bs.addUserFloat("phi_fitted_eta",   fittedPhi.Eta()); 
        bs.addUserFloat("phi_fitted_phi",   fittedPhi.Phi()); 
        bs.addUserFloat("phi_fitted_m",     fittedPhi.M()); 

        // Ds fit
        bs.addUserFloat("pi_refitted_px",   refittedPi.Px()); 
        bs.addUserFloat("pi_refitted_py",   refittedPi.Py()); 
        bs.addUserFloat("pi_refitted_pz",   refittedPi.Pz()); 
        bs.addUserFloat("pi_refitted_pt",   refittedPi.Pt()); 
        bs.addUserFloat("pi_refitted_eta",  refittedPi.Eta()); 
        bs.addUserFloat("pi_refitted_phi",  refittedPi.Phi()); 
        bs.addUserFloat("pi_refitted_m",    refittedPi.M()); 

        bs.addUserFloat("phi_refitted_px",  refittedPhi.Px()); 
        bs.addUserFloat("phi_refitted_py",  refittedPhi.Py()); 
        bs.addUserFloat("phi_refitted_pz",  refittedPhi.Pz()); 
        bs.addUserFloat("phi_refitted_pt",  refittedPhi.Pt()); 
        bs.addUserFloat("phi_refitted_eta", refittedPhi.Eta()); 
        bs.addUserFloat("phi_refitted_phi", refittedPhi.Phi()); 
        bs.addUserFloat("phi_refitted_m",   refittedPhi.M()); 

        bs.addUserFloat("ds_fitted_px",     fittedDs.Px()); 
        bs.addUserFloat("ds_fitted_py",     fittedDs.Py()); 
        bs.addUserFloat("ds_fitted_pz",     fittedDs.Pz()); 
        bs.addUserFloat("ds_fitted_pt",     fittedDs.Pt()); 
        bs.addUserFloat("ds_fitted_eta",    fittedDs.Eta()); 
        bs.addUserFloat("ds_fitted_phi",    fittedDs.Phi()); 
        bs.addUserFloat("ds_fitted_m",      fittedDs.M()); 

        // bs fit
        bs.addUserFloat("mu_refitted_px",   refittedMu.Px()); 
        bs.addUserFloat("mu_refitted_py",   refittedMu.Py()); 
        bs.addUserFloat("mu_refitted_pz",   refittedMu.Pz()); 
        bs.addUserFloat("mu_refitted_pt",   refittedMu.Pt()); 
        bs.addUserFloat("mu_refitted_eta",  refittedMu.Eta()); 
        bs.addUserFloat("mu_refitted_phi",  refittedMu.Phi()); 
        bs.addUserFloat("mu_refitted_m",    refittedMu.M()); 

        bs.addUserFloat("ds_refitted_px",   refittedDs.Px()); 
        bs.addUserFloat("ds_refitted_py",   refittedDs.Py()); 
        bs.addUserFloat("ds_refitted_pz",   refittedDs.Pz()); 
        bs.addUserFloat("ds_refitted_pt",   refittedDs.Pt()); 
        bs.addUserFloat("ds_refitted_eta",  refittedDs.Eta()); 
        bs.addUserFloat("ds_refitted_phi",  refittedDs.Phi()); 
        bs.addUserFloat("ds_refitted_m",    refittedDs.M()); 

        bs.addUserFloat("bs_fitted_px",     fittedBs.Px()); 
        bs.addUserFloat("bs_fitted_py",     fittedBs.Py()); 
        bs.addUserFloat("bs_fitted_pz",     fittedBs.Pz()); 
        bs.addUserFloat("bs_fitted_pt",     fittedBs.Pt()); 
        bs.addUserFloat("bs_fitted_eta",    fittedBs.Eta()); 
        bs.addUserFloat("bs_fitted_phi",    fittedBs.Phi()); 
        bs.addUserFloat("bs_fitted_m",      fittedBs.M()); 

        /////////////////////// helicity angles in all possibe variations /////////////////////////

        //angle between Mu and W
        float angMuWColl    = angMuW(fittedDs,collBsTlv,   refittedMu);
        float angMuWLhcb    = angMuW(fittedDs,lhcbBsTlv,   refittedMu);
        float angMuWLhcbAlt = angMuW(fittedDs,lhcbAltBsTlv,refittedMu);
        float angMuWReco1   = angMuW(fittedDs,recoBsTlv1,  refittedMu);
        float angMuWReco2   = angMuW(fittedDs,recoBsTlv2,  refittedMu);

        bs.addUserFloat("angMuWColl",       angMuWColl); 
        bs.addUserFloat("angMuWLhcb",       angMuWLhcb); 
        bs.addUserFloat("angMuWLhcbAlt",    angMuWLhcbAlt); 
        bs.addUserFloat("angMuWReco1",      angMuWReco1); 
        bs.addUserFloat("angMuWReco2",      angMuWReco2); 

        bs.addUserFloat("cosMuWColl",       cos(angMuWColl)); 
        bs.addUserFloat("cosMuWLhcb",       cos(angMuWLhcb)); 
        bs.addUserFloat("cosMuWLhcbAlt",    cos(angMuWLhcbAlt)); 
        bs.addUserFloat("cosMuWReco1",      cos(angMuWReco1)); 
        bs.addUserFloat("cosMuWReco2",      cos(angMuWReco2)); 


        //angle between k1(k2) and pion in phi rest frame
        float angPiK1  = angDoubleDecay(fittedPhi,refittedK1, refittedPi);
        float angPiK2  = angDoubleDecay(fittedPhi,refittedK2, refittedPi);
        // equivalently, angle between phi(pi) and bs in ds rest frame
        float angPhiDsColl    = angDoubleDecay(fittedDs, fittedPhi,  collBsTlv);
        float angPhiDsLhcb    = angDoubleDecay(fittedDs, fittedPhi,  lhcbBsTlv);
        float angPhiDsLhcbAlt = angDoubleDecay(fittedDs, fittedPhi,  lhcbAltBsTlv);
        float angPhiDsReco1   = angDoubleDecay(fittedDs, fittedPhi,  recoBsTlv1);
        float angPhiDsReco2   = angDoubleDecay(fittedDs, fittedPhi,  recoBsTlv2);

        float angPiDsColl    = angDoubleDecay(fittedDs, refittedPi,  collBsTlv);
        float angPiDsLhcb    = angDoubleDecay(fittedDs, refittedPi,  lhcbBsTlv);
        float angPiDsLhcbAlt = angDoubleDecay(fittedDs, refittedPi,  lhcbAltBsTlv);
        float angPiDsReco1   = angDoubleDecay(fittedDs, refittedPi,  recoBsTlv1);
        float angPiDsReco2   = angDoubleDecay(fittedDs, refittedPi,  recoBsTlv2);

        bs.addUserFloat("angPiK1",         angPiK1); 
        bs.addUserFloat("angPiK2",         angPiK2); 

        bs.addUserFloat("angPhiDsColl",    angPhiDsColl); 
        bs.addUserFloat("angPhiDsLhcb",    angPhiDsLhcb); 
        bs.addUserFloat("angPhiDsLhcbAlt", angPhiDsLhcbAlt); 
        bs.addUserFloat("angPhiDsReco1",   angPhiDsReco1); 
        bs.addUserFloat("angPhiDsReco2",   angPhiDsReco2); 

        bs.addUserFloat("angPiDsColl",     angPiDsColl); 
        bs.addUserFloat("angPiDsLhcb",     angPiDsLhcb); 
        bs.addUserFloat("angPiDsLhcbAlt",  angPiDsLhcbAlt); 
        bs.addUserFloat("angPiDsReco1",    angPiDsReco1); 
        bs.addUserFloat("angPiDsReco2",    angPiDsReco2); 

        bs.addUserFloat("cosPiK1",         cos(angPiK1)); 
        bs.addUserFloat("cosPiK2",         cos(angPiK2)); 

        bs.addUserFloat("cosPhiDsColl",    cos(angPhiDsColl)); 
        bs.addUserFloat("cosPhiDsLhcb",    cos(angPhiDsLhcb)); 
        bs.addUserFloat("cosPhiDsLhcbAlt", cos(angPhiDsLhcbAlt)); 
        bs.addUserFloat("cosPhiDsReco1",   cos(angPhiDsReco1)); 
        bs.addUserFloat("cosPhiDsReco2",   cos(angPhiDsReco2)); 

        bs.addUserFloat("cosPiDsColl",     cos(angPiDsColl)); 
        bs.addUserFloat("cosPiDsLhcb",     cos(angPiDsLhcb)); 
        bs.addUserFloat("cosPiDsLhcbAlt",  cos(angPiDsLhcbAlt)); 
        bs.addUserFloat("cosPiDsReco1",    cos(angPiDsReco1)); 
        bs.addUserFloat("cosPiDsReco2",    cos(angPiDsReco2)); 

        // helicity plane angle
        float angPlaneBsColl    = angPlane(fittedDs, collBsTlv,    refittedMu, refittedPi);
        float angPlaneBsLhcb    = angPlane(fittedDs, lhcbBsTlv,    refittedMu, refittedPi);
        float angPlaneBsLhcbAlt = angPlane(fittedDs, lhcbAltBsTlv, refittedMu, refittedPi);
        float angPlaneBsReco1   = angPlane(fittedDs, recoBsTlv1,   refittedMu, refittedPi);
        float angPlaneBsReco2   = angPlane(fittedDs, recoBsTlv2,   refittedMu, refittedPi);

        bs.addUserFloat("cosPlaneBsColl",    cos(angPlaneBsColl));
        bs.addUserFloat("cosPlaneBsLhcb",    cos(angPlaneBsLhcb));
        bs.addUserFloat("cosPlaneBsLhcbAlt", cos(angPlaneBsLhcbAlt));
        bs.addUserFloat("cosPlaneBsReco1",   cos(angPlaneBsReco1));
        bs.addUserFloat("cosPlaneBsReco2",   cos(angPlaneBsReco2));

        float angPlaneDsColl    = angPlane2(fittedDs, collBsTlv,    refittedK1, refittedPi);
        float angPlaneDsLhcb    = angPlane2(fittedDs, lhcbBsTlv,    refittedK1, refittedPi);
        float angPlaneDsLhcbAlt = angPlane2(fittedDs, lhcbAltBsTlv, refittedK1, refittedPi);
        float angPlaneDsReco1   = angPlane2(fittedDs, recoBsTlv1,   refittedK1, refittedPi);
        float angPlaneDsReco2   = angPlane2(fittedDs, recoBsTlv2,   refittedK1, refittedPi);

        bs.addUserFloat("cosPlaneDsColl",    cos(angPlaneDsColl));
        bs.addUserFloat("cosPlaneDsLhcb",    cos(angPlaneDsLhcb));
        bs.addUserFloat("cosPlaneDsLhcbAlt", cos(angPlaneDsLhcbAlt));
        bs.addUserFloat("cosPlaneDsReco1",   cos(angPlaneDsReco1));
        bs.addUserFloat("cosPlaneDsReco2",   cos(angPlaneDsReco2));


        /////////////////////// END OF VARIABLE DEFINITION //////////////////////

        arrived = 1;
        bs.addUserInt("arrived", arrived);

        //append candidate at the end of our return value :)
        //ret_value can be a vector!!
        ret_value->emplace_back(bs);
        //ret_value_gen->emplace_back(gen);

        } //closing pi loop
      } //closing k2 loop
    } //closing k1 loop
  } //closing trg muon loop

  if(arrived >0){
  iEvent.put(std::move(ret_value), "bs");
  }
  else{
  pat::CompositeCandidate bs;
  bs.addUserInt("arrived",arrived);
  ret_value->emplace_back(bs);
  iEvent.put(std::move(ret_value), "bs"); 
  }
}//closing event loop

DEFINE_FWK_MODULE(BsToDsPhiKKPiMuBuilder);
