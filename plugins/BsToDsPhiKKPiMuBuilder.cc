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
//#include "KinVtxFitter.h" --> not needed now 
#include "helper.h" // ---> for now!!

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
    hadSelection_(iConfig.getParameter<std::string>("hadSelection")),
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
    for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){

      // when do we use bestTrack() and when do we use TransientTracks()?
      edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, ++vtxIdx);
      dzMuPV.push_back(muPtr->bestTrack()->dz(vtxPtr->position()));      
    }
    
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

      // why is k1Selection_ a function?
      //if (!k1Selection_(*k1Ptr)) continue; 

      //can we put all of the below in the k1Selection?
      //the PF algorithm assigns a pdgId hypothesis, generall it distinguishes between:
      // photons, electron/muon, charged hadron, neutral hadrons
      // and we trust the algorithm that when it says its an electron (11) or muon (13), that it is not a kaon or pion

      bool k1Sel = ((k1Ptr->pdgId() != 11 ) && (k1Ptr->pdgId() != 13) && (k1Ptr->charge() != 0) && (k1Ptr->pt() > 1.0) && (fabs(k1Ptr->eta()) < 2.4) && (k1Ptr->hasTrackDetails()) && (reco::deltaR(*k1Ptr,*muPtr) < 1.2) && (reco::deltaR(*k1Ptr, *muPtr) > 0.005) && ((k1Ptr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < 0.5));

      if (!k1Sel) continue;

    //////////////////////////////////////////////////
    // Loop over k2 and select the good tracks      //
    //////////////////////////////////////////////////

      for(size_t k2Idx = k1Idx + 1; k2Idx < pcand->size(); ++k2Idx) {

      // pointer to the second kaon candidate
      edm::Ptr<pat::PackedCandidate> k2Ptr(pcand, k2Idx);
     
      // if this kaon does not pass the selection, jump to the next!
      //if(!k2Selection_(*k2Ptr)) continue;

      bool k2Sel = ((k2Ptr->pdgId() != 11 ) && (k2Ptr->pdgId() != 13) && (k2Ptr->charge() != 0) && (k2Ptr->pt() > 1.0) && fabs(k2Ptr->eta()) < 2.4 && (k2Ptr->hasTrackDetails()) && (reco::deltaR(*k2Ptr,*muPtr) < 1.2) && (reco::deltaR(*k2Ptr, *muPtr) > 0.005) &&((k2Ptr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < 0.5));

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
        // if(!piSelection_(*piPtr)) continue;

        bool piSel = ((piPtr->pdgId() != 11 ) && (piPtr->pdgId() != 13) && (piPtr->charge() != 0) && (piPtr->pt() > 1.0) && fabs(piPtr->eta()) < 2.4 && (piPtr->hasTrackDetails()) && (reco::deltaR(*piPtr,*muPtr) < 1.2) && (reco::deltaR(*piPtr, *muPtr) > 0.005) &&((piPtr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < 0.5));

        //pi and mu must have opposite charge
        if ((piPtr->charge() * muPtr->charge()) > 0) continue; //No!! Lets keep everyting for background esimtaion
        //but we can save it as mask:
        int piMuCharge = piPtr->charge() * muPtr->charge();

        if (!piSel) continue;

        //////////////////////////////////////////////////
        // Build Phi and Ds resonance if available      //
        //////////////////////////////////////////////////

        //define a composite candidate pair 
        pat::CompositeCandidate kkPair;

        //PF canidates are always assigned with the pi mass, we have to force the kaon mass
        math::PtEtaPhiMLorentzVector k1P4(k1Ptr->pt(), k1Ptr->eta(), k1Ptr->phi(), K_MASS);
        math::PtEtaPhiMLorentzVector k2P4(k2Ptr->pt(), k2Ptr->eta(), k2Ptr->phi(), K_MASS);
        math::PtEtaPhiMLorentzVector piP4(piPtr->pt(), piPtr->eta(), piPtr->phi(), PI_MASS); //just to be sure lets also force the pi mass
 
        kkPair.setP4(k1P4 + k2P4);
      
        //only continue when they build a phi resonance, allow 15MeV:
        if (fabs(kkPair.mass() - 1.019461) > 0.015) continue;
        std::cout << "this is a true phi resonance" << std::endl;

        kkPair.setCharge(k1Ptr->charge() + k2Ptr->charge());
        kkPair.addUserFloat("kkPairDeltaR", reco::deltaR(*k1Ptr, *k2Ptr));

        //save index of first and second  kaon
        kkPair.addUserInt("k1Idx", k1Idx );
        kkPair.addUserInt("k2Idx", k2Idx );

        pat::CompositeCandidate ds;
        ds.setP4(kkPair.p4() + piP4); 

        //only continue when they build a phi resonance, allow 50MeV:
        if (fabs(ds.mass() - 1.9683) > 0.05) continue;

        std::cout << "this is a true ds resonance" << std::endl;

        ds.setCharge(kkPair.charge() + piPtr->charge());
        ds.addUserFloat("phiPiDeltaR", reco::deltaR(kkPair, *piPtr));

        //save index of kaons and pi 
        ds.addUserInt("k1Idx", k1Idx);
        ds.addUserInt("k2Idx", k2Idx);
        ds.addUserInt("piIdx", piIdx);

        // TODO: collinear approx for Bs
        pat::CompositeCandidate dsMu;
        dsMu.setP4(ds.p4() + muPtr->p4()); 

        dsMu.setCharge(ds.charge() + muPtr->charge()); //sanity check:shoould be 0
        dsMu.addUserFloat("dsMuDeltaR", reco::deltaR(ds, *muPtr));

        pat::CompositeCandidate bs;
        bs.setP4(dsMu.p4() * 5.36688 / dsMu.mass());
        bs.setCharge(dsMu.charge());
        bs.addUserFloat("dsMuDeltaR", reco::deltaR(ds, *muPtr));

        //add final states

        bs.addUserCand("k1",k1Ptr);  //be aware that this ptr has the wrong mass, need to assign it in variables_cff.py
        bs.addUserFloat("k1_mass",K_MASS);

        bs.addUserCand("k2",k2Ptr);  // "
        bs.addUserFloat("k2_mass",K_MASS);

        bs.addUserCand("pi",piPtr);
        bs.addUserFloat("pi_mass",PI_MASS);

        bs.addUserCand("mu",muPtr);  

        //add resonances --> a little bit ugly, beautify later

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
        bs.addUserFloat("dsMu_mass", dsMu.mass());
        bs.addUserFloat("dsMu_charge", dsMu.charge());
        bs.addUserFloat("dsMu_deltaR", dsMu.userFloat("dsMuDeltaR"));


        //rel charges
        bs.addUserInt("kkCharge",kkCharge); 
        bs.addUserInt("piMuCharge",piMuCharge); 

        // helicity angles in all possibe variations

        //float helAngle_k1 = 
        //float helAngle_k2 = 
        //float helAngle_k+ = 

        //define a factory
        KinematicParticleFactoryFromTransientTrack pFactory;
        //define the vector for the particles to be fitted
        std::vector<RefCountedKinematicParticle> phiToFit;
        std::vector<RefCountedKinematicParticle> dsToFit;
        std::vector<RefCountedKinematicParticle> bsToFit;

        // add the final states
        ParticleMass piMass = 0.13957039;
        ParticleMass kMass = 0.493677;
        ParticleMass phiMass = 1.019461;
        ParticleMass dsMass = 1.96834;
        ParticleMass muMass = 0.105658;

        float chi = 0.0;
        float ndf = 0.0;
        float sigma = 0.0;
        float phiMassSigma = 0.000005; //discuss
        float dsMassSigma = 0.00005;  //discuss
        float muMassSigma = 0.000005;  //discuss

        phiToFit.push_back(pFactory.particle(getTransientTrack(k1Ptr->bestTrack()),kMass,chi,ndf,sigma));
        phiToFit.push_back(pFactory.particle(getTransientTrack(k2Ptr->bestTrack()),kMass,chi,ndf,sigma));
        dsToFit.push_back(pFactory.particle(getTransientTrack(piPtr->bestTrack()),piMass,chi,ndf,sigma));
        bsToFit.push_back(pFactory.particle(getTransientTrack(muPtr->bestTrack()),muMass,chi,ndf,sigma));

        // constraints
        MultiTrackKinematicConstraint* phiConstr = new TwoTrackMassKinematicConstraint(phiMass);
        MultiTrackKinematicConstraint* dsConstr = new TwoTrackMassKinematicConstraint(dsMass);
        //fitter
        KinematicConstrainedVertexFitter kcvFitter;

        // PHI VERTEX FIT
        RefCountedKinematicTree phiTree = kcvFitter.fit(phiToFit, phiConstr);
        if (!phiTree->isValid()) continue;

        phiTree->movePointerToTheTop();
        RefCountedKinematicParticle phiParticle = phiTree->currentParticle();
        RefCountedKinematicVertex phiDecayVertex = phiTree->currentDecayVertex();
        std::vector<RefCountedKinematicParticle> phiDaus = phiTree->finalStateParticles();

        std::cout << "we should have 2 of them: " << phiDaus.size() << std::endl;

        dsToFit.push_back(phiParticle);


        // DS VERTEX FIT
        RefCountedKinematicTree dsTree = kcvFitter.fit(dsToFit, dsConstr);
        if (!dsTree->isValid()) continue;

        dsTree->movePointerToTheTop();
        RefCountedKinematicParticle dsParticle = dsTree->currentParticle();
        RefCountedKinematicVertex dsDecayVertex = dsTree->currentDecayVertex();
        bsToFit.push_back(dsParticle);

        // BS VERTEX FIT
        RefCountedKinematicTree bsTree = kcvFitter.fit(bsToFit); // no constraint for bs bc missing momentum
        if (!bsTree->isValid()) continue; 

        bsTree->movePointerToTheTop();
        RefCountedKinematicParticle bsParticle = bsTree->currentParticle();
        RefCountedKinematicVertex bsDecayVertex = bsTree->currentDecayVertex();



        //access the parameters of the reconstructed ds (vector of size7, matrix of sie 7x7)
        //defined /RecoVertex/KinematicFitPrimitives/interface/Matrices.h
        AlgebraicVector7  phiParams = phiParticle->currentState().kinematicParameters().vector();
        //AlgebraicMatrix77 phiErr = phiParticle->currentState().kinematicParametersError().matrix();
        AlgebraicVector7  dsParams = dsParticle->currentState().kinematicParameters().vector();
        //AlgebraicMatrix77 dsErr = dsParticle->currentState().kinematicParametersError().matrix();
        AlgebraicVector7  bsParams = bsParticle->currentState().kinematicParameters().vector();
        //AlgebraicMatrix77 bsErr = bsParticle->currentState().kinematicParametersError().matrix();

        //check if vertex position is truly at 0,1,2: result: Yes it is!
        //double xdecay = phiDecayVertex->position().x();
        //std::cout << xdecay << "---" << phiParams(0) << std::endl; 
   
        // save fitter info
        bs.addUserFloat("phi_fitted_x",  phiParams(0)); 
        bs.addUserFloat("phi_fitted_y",  phiParams(1)); 
        bs.addUserFloat("phi_fitted_z",  phiParams(2)); 
        bs.addUserFloat("phi_fitted_px", phiParams(3)); 
        bs.addUserFloat("phi_fitted_py", phiParams(4)); 
        bs.addUserFloat("phi_fitted_pz", phiParams(5)); 
        bs.addUserFloat("phi_fitted_m",  phiParams(6)); 

        bs.addUserFloat("ds_fitted_x",  dsParams(0)); 
        bs.addUserFloat("ds_fitted_y",  dsParams(1)); 
        bs.addUserFloat("ds_fitted_z",  dsParams(2)); 
        bs.addUserFloat("ds_fitted_px", dsParams(3)); 
        bs.addUserFloat("ds_fitted_py", dsParams(4)); 
        bs.addUserFloat("ds_fitted_pz", dsParams(5)); 
        bs.addUserFloat("ds_fitted_m",  dsParams(6)); 

        bs.addUserFloat("bs_fitted_x",  bsParams(0)); 
        bs.addUserFloat("bs_fitted_y",  bsParams(1)); 
        bs.addUserFloat("bs_fitted_z",  bsParams(2)); 
        bs.addUserFloat("bs_fitted_px", bsParams(3)); 
        bs.addUserFloat("bs_fitted_py", bsParams(4)); 
        bs.addUserFloat("bs_fitted_pz", bsParams(5)); 
        bs.addUserFloat("bs_fitted_m",  bsParams(6)); 

        // helicity angles in all possibe variations
        // = angle between one of the kaons and the pi in the rest frame of the phi
 
        //math::PxPyPzMLorentzVector fittedPhiP4(phiParams(3), phiParams(4), phiParams(5), phiParams(6));
        //math::PxPyPzMLorentzVector fittedPhiP4(phiParams(3), phiParams(4), phiParams(5), phiParams(6));
        //math::PxPyPzMLorentzVector fittedPhiP4(phiParams(3), phiParams(4), phiParams(5), phiParams(6));
   
 
        //float helAngle_k1 = 
        //float helAngle_k2 = 
        //float helAngle_k+ = 



        //std::cout << bs.pt() << std::endl;

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
