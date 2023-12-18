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
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"//for the vertex fitting!!
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"//for the vertex fitting!!
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
  //const StringCutObjectSelector<pat::PackedCandidate> k1Selection_; // cut on k1
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
 
  const edm::ESInputTag ttkTag;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
};

//define the constructor
BsToDsPhiKKPiMuBuilder::BsToDsPhiKKPiMuBuilder(const edm::ParameterSet& iConfig):
    // f.e. k1Selection_ = cfg.getPatameter...
    //k1Selection_(iConfig.getParameter<std::string>("k1Selection")) ,
    //k2Selection_(iConfig.getParameter<std::string>("k2Selection")),
    //preVtxSelection_(iConfig.getParameter<std::string>("preVtxSelection")  ),
    //postVtxSelection_(iConfig.getParameter<std::string>("postVtxSelection")),

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

    primaryVtxTag(iConfig.getParameter<edm::InputTag>("pvCand")),
    primaryVtx_(consumes<reco::VertexCollection>(primaryVtxTag)),
  
    //transienTrackBuilder
    ttkTag(iConfig.getParameter<edm::ESInputTag>("test")),
    ttkToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttkTag)){
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
  const TransientTrackBuilder* theB = &iSetup.getData(ttkToken_);


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
        //std::cout << "this is a true phi resonance" << std::endl;

        kkPair.setCharge(k1Ptr->charge() + k2Ptr->charge());
        kkPair.addUserFloat("kkPairDeltaR", reco::deltaR(*k1Ptr, *k2Ptr));

        //save index of first and second  kaon
        kkPair.addUserInt("k1Idx", k1Idx );
        kkPair.addUserInt("k2Idx", k2Idx );

        pat::CompositeCandidate ds;
        ds.setP4(kkPair.p4() + piP4); 

        //only continue when they build a phi resonance, allow 50MeV:
        if (fabs(ds.mass() - 1.9683) > 0.05) continue;

        //std::cout << "this is a true ds resonance" << std::endl;

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
        std::vector<RefCountedKinematicParticle> partToFit;
        // add the final states

        //intermediate step
        const reco::Track* test = piPtr->bestTrack();
        reco::TransientTrack test2 = getTransientTrack(test);
        
        ParticleMass piMass = 0.13957039;
        ParticleMass kMass = 0.493677;
        ParticleMass phiMass = 1.019461;
        ParticleMass dsMass = 1.86965;

        float chi = 0.0;
        float ndf = 0.0;
        float sigma = 0.0;
        partToFit.push_back(pFactory.particle(getTransientTrack(k1Ptr->bestTrack()),kMass,chi,ndf,sigma));
        partToFit.push_back(pFactory.particle(getTransientTrack(k2Ptr->bestTrack()),kMass,chi,ndf,sigma));
        partToFit.push_back(pFactory.particle(getTransientTrack(piPtr->bestTrack()),piMass,chi,ndf,sigma));


        //constraints
        MultiTrackKinematicConstraint* phiConstr = new TwoTrackMassKinematicConstraint(phiMass);
        //create the fitter
        KinematicConstrainedVertexFitter kcvFitter;
        RefCountedKinematicTree myTree = kcvFitter.fit(partToFit, phiConstr);
        
        /*
        //  TODO, where are these defined?
        if( !preVtxSelection_(kkPair)) continue;
        if( !preVtxSelection_(ds)) continue;
        if( !preVtxSelection_(bs)) continue;

      
        // Fit the two kaon tracks
        KinVtxFitter fitterPhi(
          {k1Ptr->bestTrack(), k2Ptr->bestTrack()},
          {k2Ptr->mass(), k1Ptr->mass()},
          {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass ?? where are these defined?
          );
        if( !postVtxSelection_(kkPair) ) continue;//?

        // Fit the two phi and pi track
        KinVtxFitter fitterDs(
          {kkPair.bestTrack(),piPtr->bestTrack()},
          {kkPair.mass(), piPtr->mass()},
          {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass ?? where are these defined?
          );
        if( !postVtxSelection_(dsMu) ) continue;//?

        // Fit the Ds and mu track
        KinVtxFitter fitterDsMu(
          {ds.bestTrack(),muPtr->bestTrack()},
          {ds.mass(), muPtr->mass()},
          {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass ?? where are these defined?
          );

        if( !postVtxSelection_(dsMu) ) continue;//?

        //vertuces must be valid!
        if ((!fitterPhi.fitted_vtx().vertexIsValid()) || (!fitterDs.fitted_vtx().vertexIsValid()) || (!fitterDsMu.fitted_vtx().vertexIsValid()) ) continue;

        // save fitter info
        kkPair.addUserFloat("phiVtx_success", fitterPhi.success()); // float?
        kkPair.addUserFloat("phiVtx_prob", fitterPhi.prob());
        kkPair.addUserFloat("phiVtx_chi2", fitterPhi.chi2());
        
        kkPair.addUserFloat("phiVtx_position", fitterPhi.fitted_vtx().x()); //why x? 
        kkPair.addUserFloat("phiVtx_ndof", fitterPhi.dof()); //flaot?
        
        auto fit_p4 = fitter.fitted_p4();
       
        //syntax explanation:
        // if fitter.success() == true, return the fitter.fitted_candidate().mass(), otherwise return -1
        kkPair.addUserFloat("phi_m"      , fitterPhi.success() ? fitterPhi.fitted_candidate().mass() : -1); 
        kkPair.addUserFloat("phi_merr"   , fitterPhi.success() ? sqrt(fitterPhi.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1);
        
        kkPair.addUserFloat("phiVtx_x"   , fitterPhi.fitted_vtx().x());
        kkPair.addUserFloat("phiVtx_y"   , fitterPhi.fitted_vtx().y());
        kkPair.addUserFloat("phiVtx_z"   , fitterPhi.fitted_vtx().z());
        kkPair.addUserFloat("phiVtx_ex"  , sqrt(fitterPhi.fitted_vtx_uncertainty().cxx()));
        kkPair.addUserFloat("phiVtx_ey"  , sqrt(fitterPhi.fitted_vtx_uncertainty().cyy()));
        kkPair.addUserFloat("phiVtx_ez"  , sqrt(fitterPhi.fitted_vtx_uncertainty().czz()));
        kkPair.addUserFloat("phiVtx_chi2", fitterPhi.chi2());

        kkPair.addUserFloat("phi_pt"  , fit_p4.pt());
        kkPair.addUserFloat("phi_eta" , fit_p4.eta());
        kkPair.addUserFloat("phi_phi" , fit_p4.phi());
        kkPair.addUserFloat("phi_x"   , fit_p4.x());
        kkPair.addUserFloat("phi_y"   , fit_p4.y());
        kkPair.addUserFloat("phi_z"   , fit_p4.z());

        kkPair.addUserFloat("k1_pt"  , fitterPhi.daughter_p4(0).pt());
        kkPair.addUserFloat("k1_eta" , fitterPhi.daughter_p4(0).eta());
        kkPair.addUserFloat("k1_phi" , fitterPhi.daughter_p4(0).phi());
        kkPair.addUserFloat("k1_x"   , fitterPhi.daughter_p4(0).x());
        kkPair.addUserFloat("k1_y"   , fitterPhi.daughter_p4(0).y());
        kkPair.addUserFloat("k1_z"   , fitterPhi.daughter_p4(0).z());
        kkPair.addUserFloat("k2_pt"  , fitterPhi.daughter_p4(1).pt());
        kkPair.addUserFloat("k2_eta" , fitterPhi.daughter_p4(1).eta());
        kkPair.addUserFloat("k2_phi" , fitterPhi.daughter_p4(1).phi());
        kkPair.addUserFloat("k2_x"   , fitterPhi.daughter_p4(1).x());
        kkPair.addUserFloat("k2_y"   , fitterPhi.daughter_p4(1).y());
        kkPair.addUserFloat("k2_z"   , fitterPhi.daughter_p4(1).z());


        kkPair.addUserFloat("phi_err00", fitterPhi.fitted_candidate().kinematicParametersError().matrix()(0,0));
        kkPair.addUserFloat("phi_err11", fitterPhi.fitted_candidate().kinematicParametersError().matrix()(1,1));
        kkPair.addUserFloat("phi_err22", fitterPhi.fitted_candidate().kinematicParametersError().matrix()(2,2));
        kkPair.addUserFloat("phi_err01", fitterPhi.fitted_candidate().kinematicParametersError().matrix()(0,1));
        kkPair.addUserFloat("phi_err02", fitterPhi.fitted_candidate().kinematicParametersError().matrix()(0,2));
        kkPair.addUserFloat("phi_err12", fitterPhi.fitted_candidate().kinematicParametersError().matrix()(1,2));
      
        Ds.addUserFloat("ds_m"      , fitterDs.success() ? fitter.fitted_candidate().mass() : -1); 
        Ds.addUserFloat("ds_merr"   , fitterDs.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1);
        Ds.addUserFloat("dsVtx_x"   , fitterDs.fitted_vtx().x());
        Ds.addUserFloat("dsVtx_y"   , fitterDs.fitted_vtx().y());
        Ds.addUserFloat("dsVtx_z"   , fitterDs.fitted_vtx().z());
        Ds.addUserFloat("dsVtx_ex"  , sqrt(fitterDs.fitted_vtx_uncertainty().cxx()));
        Ds.addUserFloat("dsVtx_ey"  , sqrt(fitterDs.fitted_vtx_uncertainty().cyy()));
        Ds.addUserFloat("dsVtx_ez"  , sqrt(fitterDs.fitted_vtx_uncertainty().czz()));
        Ds.addUserFloat("dsVtx_chi2", fitterDs.chi2());

        Ds.addUserFloat("ds_pt"  , fit_p4.pt());
        Ds.addUserFloat("ds_eta" , fit_p4.eta());
        Ds.addUserFloat("ds_phi" , fit_p4.phi());
        Ds.addUserFloat("ds_x"   , fit_p4.x());
        Ds.addUserFloat("ds_y"   , fit_p4.y());
        Ds.addUserFloat("ds_z"   , fit_p4.z());

        Ds.addUserFloat("phi_pt"  , fitterDs.daughter_p4(0).pt());
        Ds.addUserFloat("phi_eta" , fitterDs.daughter_p4(0).eta());
        Ds.addUserFloat("phi_phi" , fitterDs.daughter_p4(0).phi());
        Ds.addUserFloat("phi_x"   , fitterDs.daughter_p4(0).x());
        Ds.addUserFloat("phi_y"   , fitterDs.daughter_p4(0).y());
        Ds.addUserFloat("phi_z"   , fitterDs.daughter_p4(0).z());
        Ds.addUserFloat("pi_pt"  , fitterDs.daughter_p4(1).pt());
        Ds.addUserFloat("pi_eta" , fitterDs.daughter_p4(1).eta());
        Ds.addUserFloat("pi_phi" , fitterDs.daughter_p4(1).phi());
        Ds.addUserFloat("pi_x"   , fitterDs.daughter_p4(1).x());
        Ds.addUserFloat("pi_y"   , fitterDs.daughter_p4(1).y());
        Ds.addUserFloat("pi_z"   , fitterDs.daughter_p4(1).z());


        Ds.addUserFloat("ds_err00", fitterDs.fitted_candidate().kinematicParametersError().matrix()(0,0));
        Ds.addUserFloat("ds_err11", fitterDs.fitted_candidate().kinematicParametersError().matrix()(1,1));
        Ds.addUserFloat("ds_err22", fitterDs.fitted_candidate().kinematicParametersError().matrix()(2,2));
        Ds.addUserFloat("ds_err01", fitterDs.fitted_candidate().kinematicParametersError().matrix()(0,1));
        Ds.addUserFloat("ds_err02", fitterDs.fitted_candidate().kinematicParametersError().matrix()(0,2));
        Ds.addUserFloat("ds_err12", fitterDs.fitted_candidate().kinematicParametersError().matrix()(1,2));
      
        dsMu.addUserFloat("dsMu_m"      , fitterDs.success() ? fitter.fitted_candidate().mass() : -1); 
        dsMu.addUserFloat("dsMu_merr"   , fitterDs.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1);
        dsMu.addUserFloat("dsMuVtx_x"   , fitterDs.fitted_vtx().x());
        dsMu.addUserFloat("dsMuVtx_y"   , fitterDs.fitted_vtx().y());
        dsMu.addUserFloat("dsMuVtx_z"   , fitterDs.fitted_vtx().z());
        dsMu.addUserFloat("dsMuVtx_ex"  , sqrt(fitterDs.fitted_vtx_uncertainty().cxx()));
        dsMu.addUserFloat("dsMuVtx_ey"  , sqrt(fitterDs.fitted_vtx_uncertainty().cyy()));
        dsMu.addUserFloat("dsMuVtx_ez"  , sqrt(fitterDs.fitted_vtx_uncertainty().czz()));
        dsMu.addUserFloat("dsMuVtx_chi2", fitterDs.chi2());

        dsMu.addUserFloat("dsMu_pt"  , fit_p4.pt());
        dsMu.addUserFloat("dsMu_eta" , fit_p4.eta());
        dsMu.addUserFloat("dsMu_phi" , fit_p4.phi());
        dsMu.addUserFloat("dsMu_x"   , fit_p4.x());
        dsMu.addUserFloat("dsMu_y"   , fit_p4.y());
        dsMu.addUserFloat("dsMu_z"   , fit_p4.z());

        dsMu.addUserFloat("dsMu_pt"  , fitterDs.daughter_p4(0).pt());
        dsMu.addUserFloat("dsMu_eta" , fitterDs.daughter_p4(0).eta());
        dsMu.addUserFloat("dsMu_phi" , fitterDs.daughter_p4(0).phi());
        dsMu.addUserFloat("dsMu_x"   , fitterDs.daughter_p4(0).x());
        dsMu.addUserFloat("dsMu_y"   , fitterDs.daughter_p4(0).y());
        dsMu.addUserFloat("dsMu_z"   , fitterDs.daughter_p4(0).z());
        dsMu.addUserFloat("dsMu_pt"  , fitterDs.daughter_p4(1).pt());
        dsMu.addUserFloat("dsMu_eta" , fitterDs.daughter_p4(1).eta());
        dsMu.addUserFloat("dsMu_phi" , fitterDs.daughter_p4(1).phi());
        dsMu.addUserFloat("dsMu_x"   , fitterDs.daughter_p4(1).x());
        dsMu.addUserFloat("dsMu_y"   , fitterDs.daughter_p4(1).y());
        dsMu.addUserFloat("dsMu_z"   , fitterDs.daughter_p4(1).z());


        dsMu.addUserFloat("dsMu_err00", fitterDs.fitted_candidate().kinematicParametersError().matrix()(0,0));
        dsMu.addUserFloat("dsMu_err11", fitterDs.fitted_candidate().kinematicParametersError().matrix()(1,1));
        dsMu.addUserFloat("dsMu_err22", fitterDs.fitted_candidate().kinematicParametersError().matrix()(2,2));
        dsMu.addUserFloat("dsMu_err01", fitterDs.fitted_candidate().kinematicParametersError().matrix()(0,1));
        dsMu.addUserFloat("dsMu_err02", fitterDs.fitted_candidate().kinematicParametersError().matrix()(0,2));
        dsMu.addUserFloat("dsMu_err12", fitterDs.fitted_candidate().kinematicParametersError().matrix()(1,2));

        //calculate the decay lenght of the phi in the transverse plane 
         
        const reco::GlobalPoint& phiPosi =  fitterPhi.fitted_vtx().position();
        const reco::GlobalPoint& dsPosi  =  fitterDs.fitted_vtx().position(); 

        reco::GlobalPoint lxyzPosi = phiPosi - dsPosi;
        reco::GlobalPoint lxyPosi = (lxyzPosi.x(),lxyzPosi.y(),0.0);

        kkPair.addUserFloat("phiLxyz",lxyzPosi.mag());
        kkPair.addUserFloat("phiLxy",lxyPosi.mag());

        // if needed, add here more stuff
        // cut on the SV info
        // const reco::TransientTrack& fitted_candidate_ttrk()
        if(!fitter.fitted_candidate_ttrk().isValid()) continue;
        const reco::TransientTrack& dimuonTT = fitter.fitted_candidate_ttrk();
        int pvIdx = getPVIdx(vertices, dimuonTT);
        
  
        dimuon_tt->emplace_back(fitter.fitted_candidate_ttrk());
        //ret_value->push_back(muon_pair);
        */

        //kkPair.addUserFloat("sumpt"  , k1NewPtr->pt() + k2NewPtr->pt());

        //std::cout << iEvent.id().event() << std::endl;
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
