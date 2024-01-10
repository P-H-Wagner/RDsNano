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
#include "TrackingTools/Records/interface/TransientTrackRecord.h" // for the vertex fitting!
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" //for vertex fitting!
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h" //for vertex fitting!
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" //for vertex fitting!
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h" //for vertex fitting!!
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h" //for the vertex fitting!!
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h" //for the vertex fitting!!
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h" //for the vertex fitting!!
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h" //for the vertex fitting!!
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicConstraint.h" //for the vertex fitting!
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"//for the vertex fitting!!
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"//for the vertex fitting!!
#include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h" //for vertex fitting!!
#include "DataFormats/GeometryVector/interface/GlobalPoint.h" //for vertex fitting!!
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
//B field
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

class BsToDsPhiKKPiMuBuilder : public edm::global::EDProducer<> {

public:
  
  //typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  //constructor
  explicit BsToDsPhiKKPiMuBuilder(const edm::ParameterSet&);
  //destructor
  ~BsToDsPhiKKPiMuBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

  //function which produces a transient Track out of a track reference
  //reco::TransientTrack getTransientTrack(const reco::TrackRef& trackRef) {    
  //    reco::TransientTrack transientTrack(trackRef, paramField);
  //    return transientTrack;
  //  }  

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
  const StringCutObjectSelector<pat::PackedCandidate> hadSelection_; // cut on k1
  const double maxdRHadMuon_;
  const double mindRHadMuon_;
  const double maxdzDiffHadMuon_; 
  const double phiMassAllowance_;
  const double dsMassAllowance_;
  const double piMass_;
  const double kMass_;
  const double phiMass_;
  const double dsMass_;
  const double muMass_;
  const double bsMass_;


  //const StringCutObjectSelector<pat::PackedCandidate> k2Selection_; // cut on k2
  //const StringCutObjectSelector<pat::CompositeCandidate> preVtxSelection_; 
  //const StringCutObjectSelector<pat::CompositeCandidate> postVtxSelection_; 
 
  //tokens to access data later
  //edm::Input tag can not be directly initialized inside the construcor! Why did it work fro Trigger.cc??
  //anyway ... 
  const edm::InputTag srcTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> src_;

  //const edm::InputTag ttrackTag;
  //const edm::EDGetTokenT<TransientTrackCollection> ttracksSrc_;

  //const edm::InputTag vertexTag;
  //const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

  //const edm::InputTag beamspotTag;
  //const edm::EDGetTokenT<reco::BeamSpot> beamspot_;

  //for the muons

  const edm::InputTag trgMuonTag;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuons_;

  //const edm::InputTag selMuonTag;
  //const edm::EDGetTokenT<pat::MuonCollection> selMuons_;

  //const edm::InputTag ttrackMuonTag;
  //const edm::EDGetTokenT<TransientTrackCollection> muonTracks_;

  const edm::InputTag primaryVtxTag;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVtx_;
 
  //const edm::ESInputTag ttkTag;
  //const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;

};

//define the constructor
BsToDsPhiKKPiMuBuilder::BsToDsPhiKKPiMuBuilder(const edm::ParameterSet& iConfig):
    // f.e. hadSelection_ = cfg.getPatameter...
    hadSelection_(iConfig.getParameter<std::string>("hadSelection")) ,
    maxdRHadMuon_(iConfig.getParameter<double>("maxdRHadMuon")),
    mindRHadMuon_(iConfig.getParameter<double>("mindRHadMuon")),
    maxdzDiffHadMuon_(iConfig.getParameter<double>("maxdzDiffHadMuon")),
    phiMassAllowance_(iConfig.getParameter<double>("phiMassAllowance")),
    dsMassAllowance_(iConfig.getParameter<double>("dsMassAllowance")),

    piMass_(iConfig.getParameter<double>("piMass")),
    kMass_(iConfig.getParameter<double>("kMass")),
    phiMass_(iConfig.getParameter<double>("phiMass")),
    dsMass_(iConfig.getParameter<double>("dsMass")),
    muMass_(iConfig.getParameter<double>("muMass")),
    bsMass_(iConfig.getParameter<double>("bsMass")),
 
    srcTag(iConfig.getParameter<edm::InputTag>("pfCand")),
    src_(consumes<pat::PackedCandidateCollection>(srcTag)), 

    //ttrackTag(iConfig.getParameter<edm::InputTag>("transientTracksSrc")),
    //ttracksSrc_(consumes<TransientTrackCollection>(ttrackTag)),

    //vertexTag(iConfig.getParameter<edm::InputTag>("vertexCollection")),
    //vertexSrc_(consumes<reco::VertexCollection> (vertexTag)),

    //beamspotTag(iConfig.getParameter<edm::InputTag>("beamSpot")),
    //beamspot_(consumes<reco::BeamSpot>(beamspotTag)),
    //load the muons from Trigger.cc

    trgMuonTag(iConfig.getParameter<edm::InputTag>("muCand")),
    trgMuons_(consumes<pat::MuonCollection>(trgMuonTag)), 

    //selMuonTag(iConfig.getParameter<edm::InputTag>("SelectedMuons")),
    //selMuons_(consumes<pat::MuonCollection>(selMuonTag)), 

    //ttrackMuonTag(iConfig.getParameter<edm::InputTag>("SelectedTransientTracks")),
    //muonTracks_(consumes<TransientTrackCollection>(ttrackMuonTag)), 

    //transienTrackBuilder
    //ttkTag(iConfig.getParameter<edm::ESInputTag>("test")),
    //ttkToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttkTag))

    primaryVtxTag(iConfig.getParameter<edm::InputTag>("pvCand")),
    primaryVtx_(consumes<reco::VertexCollection>(primaryVtxTag)){
       //body of the constructor
       //define edm to be filled collections
       produces<pat::CompositeCandidateCollection>("bs");
       //produces<TransientTrackCollection>("kkTransientTracks");
    }

//check const keywords 
void BsToDsPhiKKPiMuBuilder::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  //input
  edm::Handle<pat::PackedCandidateCollection> pcand;
  iEvent.getByToken(src_, pcand);
  
  //edm::Handle<TransientTrackCollection> ttracks;
  //iEvent.getByToken(ttracksSrc_, ttracks);

  //edm::Handle<reco::BeamSpot> beamspot;
  //iEvent.getByToken(beamspot_, beamspot);

  //edm::Handle<reco::VertexCollection> thePrimaryVerticesHandle;
  //iEvent.getByToken(vertexSrc_, thePrimaryVerticesHandle);

  edm::Handle<pat::MuonCollection> trgMuons;
  iEvent.getByToken(trgMuons_,trgMuons);
 
  //edm::Handle<pat::MuonCollection> selMuons;
  //iEvent.getByToken(selMuons_,selMuons);

  //edm::Handle<pat::MuonCollection> muonTracks;
  //iEvent.getByToken(muonTracks_,muonTracks);

  edm::Handle<reco::VertexCollection> primaryVtx;
  iEvent.getByToken(primaryVtx_,primaryVtx);

  //get the TransientTrackBuilder
  //const TransientTrackBuilder* theB = &iSetup.getData(ttkToken_);


  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  //std::unique_ptr<TransientTrackCollection> kkpi_ttrack(new TransientTrackCollection);

  //over which muons do we need to loop??
  for(size_t muIdx = 0; muIdx < trgMuons->size(); ++muIdx){
  
    //std::cout << "Lets loop over this matched muon!" << std::endl; 
 
    //define a pointer to the muon called mu_ptr
    edm::Ptr<reco::Muon> muPtr(trgMuons, muIdx);
 
    std::vector<float> dzMuPV;

    //Fix the primary vertex to be the one closest to the trg Muon in dz
    // more accurate for mu than for tau signals
    for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){

      // when do we use bestTrack() and when do we use TransientTracks()?
      edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, ++vtxIdx);
      dzMuPV.push_back(muPtr->bestTrack()->dz(vtxPtr->position()));      
    }
    
    // take as the primary vertex the one which has minimal dz with the muon 
    auto dzMuPVMin = std::min_element(std::begin(dzMuPV),std::end(dzMuPV));
    int pvIdx = std::distance(std::begin(dzMuPV), dzMuPVMin); 
    reco::Vertex pv = primaryVtx->at(pvIdx);
 
    //if this event does not have a matching muon-triggermuon pair we dont want it
    //  if (iEvent.id().event() != mu_ptr->userInt("eventNr")) continue;

    //////////////////////////////////////////////////
    // Loop over k1 and select the good tracks      //
    //////////////////////////////////////////////////

    for(size_t k1Idx = 0; k1Idx < pcand->size(); ++k1Idx) {

      //define a pointer to the kaon at position k1Idx
      edm::Ptr<pat::PackedCandidate> k1Ptr(pcand, k1Idx);

      // why is hadSelection_ a function?
      if (!hadSelection_(*k1Ptr)) continue; 

      //the PF algorithm assigns a pdgId hypothesis, generall it distinguishes between:
      // photons, electron/muon, charged hadron, neutral hadrons
      // and we trust the algorithm that when it says its an electron (11) or muon (13), that it is not a kaon or pion

      float muonK1dR = reco::deltaR(*k1Ptr,*muPtr);

      bool k1Sel = (( muonK1dR < maxdRHadMuon_ ) && (reco::deltaR(*k1Ptr, *muPtr) > mindRHadMuon_) && (abs(k1Ptr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_));

      if (!k1Sel) continue;

    //////////////////////////////////////////////////
    // Loop over k2 and select the good tracks      //
    //////////////////////////////////////////////////

      for(size_t k2Idx = k1Idx + 1; k2Idx < pcand->size(); ++k2Idx) {

      // pointer to the second kaon candidate
      edm::Ptr<pat::PackedCandidate> k2Ptr(pcand, k2Idx);
     
      // if this kaon does not pass the selection, jump to the next!
      if(!hadSelection_(*k2Ptr)) continue;

      float muonK2dR = reco::deltaR(*k2Ptr,*muPtr);

      bool k2Sel = (( muonK2dR < maxdRHadMuon_ ) && (reco::deltaR(*k2Ptr, *muPtr) > mindRHadMuon_) && (abs(k2Ptr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_));

      //k1 and k2 must have oppoiste charge
      if ((k1Ptr->charge() * k2Ptr->charge()) > 0) continue; //No!! lets keep everything for background estimation
      int kkCharge = k1Ptr->charge() * k2Ptr->charge();


      if (!k2Sel) continue;

      //////////////////////////////////////////////////
      // Loop over pi and select the good tracks      //
      //////////////////////////////////////////////////

      for(size_t piIdx = 0; piIdx < pcand->size(); ++piIdx) {

        // pointer to the second kaon candidate
        edm::Ptr<pat::PackedCandidate> piPtr(pcand, piIdx);

        //make sure the pion is none of the kaons:
        if((piIdx == k1Idx) || (piIdx == k2Idx)) continue;
     
        // if this pion does not pass the selection, jump to the next!
        if(!hadSelection_(*piPtr)) continue;

        float muonPidR = reco::deltaR(*piPtr,*muPtr);

        bool piSel = ((muonPidR < maxdRHadMuon_) && (reco::deltaR(*piPtr, *muPtr) > mindRHadMuon_) &&(abs(piPtr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_));

        //pi and mu must have opposite charge
        if ((piPtr->charge() * muPtr->charge()) > 0) continue; //No!! Lets keep everyting for background esimtaion
        //but we can save it as mask:
        int piMuCharge = piPtr->charge() * muPtr->charge();

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
        kkPair.addUserFloat("kkPairDeltaR", reco::deltaR(*k1Ptr, *k2Ptr));

        //save index of first and second  kaon
        kkPair.addUserInt("k1Idx", k1Idx );
        kkPair.addUserInt("k2Idx", k2Idx );

        pat::CompositeCandidate ds;
        ds.setP4(kkPair.p4() + piP4); 

        //only continue when they build a ds resonance, allow 50MeV:
        if (fabs(ds.mass() - dsMass_) > dsMassAllowance_) continue;


        ds.setCharge(kkPair.charge() + piPtr->charge());
        ds.addUserFloat("phiPiDeltaR", reco::deltaR(kkPair, *piPtr));

        //save index of kaons and pi 
        ds.addUserInt("k1Idx", k1Idx);
        ds.addUserInt("k2Idx", k2Idx);
        ds.addUserInt("piIdx", piIdx);

        //////////////////////////////////////////////////
        // Build Bs resonance                           //
        //////////////////////////////////////////////////

        pat::CompositeCandidate dsMu;
        dsMu.setP4(ds.p4() + muPtr->p4()); 

        dsMu.setCharge(ds.charge() + muPtr->charge()); //sanity check:shoould be 0
        dsMu.addUserFloat("dsMuDeltaR", reco::deltaR(ds, *muPtr));

        //build bs with collinear approximation, below we fit
        pat::CompositeCandidate bs;
        bs.setP4(dsMu.p4() * bsMass_ / dsMu.mass()); //the bs_mass will thus be fixed at 536688 (peak in the histo)
        bs.setCharge(dsMu.charge());
        bs.addUserFloat("dsMuDeltaR", reco::deltaR(ds, *muPtr));

        //add final states

        bs.addUserCand("k1",k1Ptr);  //be aware that this ptr has the wrong mass, need to assign it in variables_cff.py
        bs.addUserFloat("k1_mass",kMass_);

        bs.addUserCand("k2",k2Ptr);  // "
        bs.addUserFloat("k2_mass",kMass_);

        bs.addUserCand("pi",piPtr);
        bs.addUserFloat("pi_mass",piMass_);

        bs.addUserCand("mu",muPtr);  

        //add resonances --> are not part of collection and can thus not be
        //automatically access the pt(), .. as I can for k1,k2,pi,mu in the variables_cff.py

        bs.addUserFloat("phi_pt", kkPair.pt());
        bs.addUserFloat("phi_eta", kkPair.eta());
        bs.addUserFloat("phi_phi", kkPair.phi());
        bs.addUserFloat("phi_mass", kkPair.mass());
        bs.addUserFloat("phi_charge", kkPair.charge());
        bs.addUserFloat("phi_deltaR", kkPair.userFloat("kkPairDeltaR"));
        bs.addUserInt("kkCharge",kkCharge); 

        bs.addUserFloat("ds_pt", ds.pt());
        bs.addUserFloat("ds_eta", ds.eta());
        bs.addUserFloat("ds_phi", ds.phi());
        bs.addUserFloat("ds_mass", ds.mass());
        bs.addUserFloat("ds_charge", ds.charge());
        bs.addUserFloat("ds_deltaR", ds.userFloat("phiPiDeltaR"));

        bs.addUserFloat("dsMu_pt", dsMu.pt());
        bs.addUserFloat("dsMu_eta", dsMu.eta());
        bs.addUserFloat("dsMu_phi", dsMu.phi());
        bs.addUserFloat("dsMu_mass", dsMu.mass()); //we miss the neutrino mass
        bs.addUserFloat("dsMu_charge", dsMu.charge());
        bs.addUserFloat("dsMu_deltaR", dsMu.userFloat("dsMuDeltaR"));


        //rel charges
        bs.addUserInt("kkCharge",kkCharge); 
        bs.addUserInt("piMuCharge",piMuCharge); 

        // helicity angles in all possibe variations

        //float helAngle_k1 = 
        //float helAngle_k2 = 
        //float helAngle_k+ = 

        //easy fit as a sanity check
        //KinVtxFitter easyFitter(
        //{getTransientTrack(k1Ptr->bestTrack()), getTransientTrack(k2Ptr->bestTrack())},
        //{K_MASS, K_MASS},
        //{0.00005,0.00005} //some small sigma for the lepton mass
        //);
        //if(!easyFitter.success()) continue;
        //std::cout << "phi easy fit vx: " << fitter.fitted_vtx().x() << std::endl;
        //std::cout << "phi easy fit vy: " << fitter.fitted_vtx().y() << std::endl;
        //std::cout << "phi easy fit vz: " << fitter.fitted_vtx().z() << std::endl;

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

        //////////////////////////////////// end of global fitter /////////////////////////////////////

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

        /////////////////////// helicity angles in all possibe variations /////////////////////////
        // = angle between one of the kaons and the pi in the rest frame of the phi
  
        //fitted resonances
        TLorentzVector fittedPhi;
        TLorentzVector fittedDs;
        TLorentzVector fittedBs;
        TLorentzVector collBs;

        fittedPhi.SetXYZM(phiParams(3), phiParams(4), phiParams(5), phiParams(6));
        fittedDs.SetXYZM(dsParams(3), dsParams(4), dsParams(5), dsParams(6));
        fittedBs.SetXYZM(bsParams(3), bsParams(4), bsParams(5), bsParams(6));
        collBs.SetPtEtaPhiM(bs.p4().pt(), bs.p4().eta(), bs.p4().phi(), bs.p4().mass());

        TLorentzVector fittedCollDs = fittedDs; //need ds twice

        // energy transfered
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

        //refitted final states
        TLorentzVector refittedK1;
        TLorentzVector refittedK2;
        TLorentzVector refittedPi;
        TLorentzVector refittedMu;

        refittedK1.SetXYZM(phiDau1Params(3), phiDau1Params(4), phiDau1Params(5), phiDau1Params(6));
        refittedK2.SetXYZM(phiDau2Params(3), phiDau2Params(4), phiDau2Params(5), phiDau2Params(6));
        refittedPi.SetXYZM(dsDau1Params(3), dsDau1Params(4), dsDau1Params(5), dsDau1Params(6));
        refittedMu.SetXYZM(bsDau1Params(3), bsDau1Params(4), bsDau1Params(5), bsDau1Params(6));

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
        w.SetVectM(-fittedDs.Vect(),q2); //the W is back to back with the Ds in the rest frame of the Bs and has mass q2
        TVector3 boostW= w.BoostVector();       
 
        refittedMu.Boost(-boostW);

        TLorentzVector wColl;
        wColl.SetVectM(-fittedCollDs.Vect(),q2Coll); //the W is back to back with the Ds in the rest frame of the Bs and has mass q2
        TVector3 boostWColl= wColl.BoostVector();       
 
        refittedMuColl.Boost(-boostWColl);

        //boost phi and pi into rest frame of Ds
        fittedPhi.Boost(-boostDs);
        refittedPi2.Boost(-boostDs);

 
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





        ret_value->emplace_back(bs);
        //std::cout << "saved!" << std::endl; 
      } //closing k2 loop
    } //closing k1 loop
  } //closing pi loop
  //move and store these two new collections in the event 
  //iEvent.put(std::move(ret_value), "bs");
  //evt.put(std::move(dimuon_tt), "KKPiTransientTracks");
} //closing muon loop
  iEvent.put(std::move(ret_value), "bs");

}
DEFINE_FWK_MODULE(BsToDsPhiKKPiMuBuilder);
