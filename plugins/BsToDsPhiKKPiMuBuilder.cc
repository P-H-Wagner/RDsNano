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
#include <iostream>
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
// photon and taus
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"

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

// vtx probability
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

// for lxy
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

////////////////////////////////////////////////////////////////////////////////////////////
// TODOS:
//
// - move gen matching in separate module?             - DO THIS AT SOME POINT
// - remove hardcoded numbers                          - DONE
// - helicity plane angles                             - DONE 
// - redefine all variables after the fit? save both ? - YES 
// - beautify the bs.addUserFloat (...)                - DONE
// - pos. def. cov matrix                              - DONE
// - add counters before every selection               - DONE 
// - pruned vs packed -> discuss                       - DONE
// - output tree has now empty entries when there is no trigger/signal -> DONE (ed filter)
// - adapt for Hb background sample                    - DONE
// - how to save kk same sign pair?                    - DONE
// - generally: save only gen matched signals?         - NO
// - what if an event has two signals?                 - SAVE BOTH!
// - divide into submitter chunks                      - DONE
// - save gen information!!                            - DONE 
// - do gen tests, check f.e. refitted p's with gen p's and unfitted p's - DONE
// - put hel angle calc. etc into functions!           - DONE
// - isAncestor and getAncestor are save, they dont modify the Ptr - CHECKED
// - add Ds boost as dicrimanting variable             - DONE
// - add all impact parameters correctly and use refitted tracks - DONE
// - add mu isolation, e* miss and pt miss, photon energy - DONE 
// - add and save Vcb Variables
// - Oscillate back also all b hadrons, pick first beauty ancestor
///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
// counters for filters to know when we lose how many particles

int nMuons  = 0;        // counts the nr of muons
int nTracks = 0;        // counts the nr of tracks in total
int nPv     = 0;        // counts the nr of tracks in total

int muSelCounter  = 0;   // how many muons pass (i.e. when we find a pv)

int k1Sel1Counter = 0;   
int k1Sel2Counter = 0;

int k2Sel1Counter = 0;
int k2Sel2Counter = 0;

int piSel1Counter = 0;
int piSel2Counter = 0;

int gSelCounter   = 0;   

int nKKPiMu       = 0;
int nPhiMassCut   = 0;
int nDsMassCut    = 0;
int nBsMassCut    = 0;

int nPhiFit       = 0;
int nDsFit        = 0;
int nBsFit        = 0;

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
  virtual void endJob() override; // NEW!

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
  const StringCutObjectSelector<pat::PackedCandidate> gSelection_; // cut on hadrons
  const StringCutObjectSelector<reco::GenParticle> hadSelectionGen_; // cut on gen hadrons for test with pruned only!

  const double maxdRHadMuon_;
  const double mindRHadMuon_;
  const double maxdRPhotonDs_;
  const double maxdzDiffHadMuon_; 
  const double maxdxyHadPv_; 
  const double phiMassAllowance_;
  const double dsMassAllowance_;
  const double dsStarMassAllowance_;
  const double drMatchGen_;
  const double maxBsMass_;
  const double piMass_;
  const double piMassSigma_;
  const double kMass_;
  const double kMassSigma_;
  const double phiMass_;
  const bool   constrainPhiMass_;
  const double minPhiVtxProb_;
  const double dsMass_;
  const bool   constrainDsMass_;
  const double minDsVtxProb_;
  const double dsStarMass_;
  const double muMass_;
  const double muMassSigma_;
  const double bsMass_;
  const double isoCone_;
  //tokens to access data later
  //edm::Input tag can not be directly initialized inside the construcor! Why did it work fro Trigger.cc??
  //anyway ... 

  const edm::InputTag beamSpotTag;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;

  const edm::InputTag srcTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> src_;

  //for the muons

  const edm::InputTag trgMuonTag;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuons_;

  // vertices
  const edm::InputTag primaryVtxTag;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVtx_;

  // tracks for isolation
  const edm::InputTag isoTracksTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isoTracks_;

  // lost tracks for isolation
  const edm::InputTag tracksLostTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksLost_;
 
  //gen for gen-matching
  const edm::InputTag prunedGenTag; //pruned is a compressed packed format
  const edm::EDGetTokenT<reco::GenParticleCollection> prunedGen_;

  const edm::InputTag packedGenTag; //packed contains much more info->most likely not needed!
  const edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGen_;

  const edm::InputTag photonTag; 
  const edm::EDGetTokenT<reco::PhotonCoreCollection> photon_;

  const edm::InputTag unpackedTracksTag;
  const edm::EDGetTokenT<reco::TrackCollection> unpackedTracks_;

  //const edm::InputTag refitVtxTag;
  //const edm::EDGetTokenT<reco::VertexCollection> refitVtx_;


};

//define the constructor
BsToDsPhiKKPiMuBuilder::BsToDsPhiKKPiMuBuilder(const edm::ParameterSet& iConfig):
    // f.e. hadSelection_ = cfg.getPatameter...
    hadSelection_(iConfig.getParameter<std::string>("hadSelection")),
    gSelection_(iConfig.getParameter<std::string>("gSelection")),
    hadSelectionGen_(iConfig.getParameter<std::string>("hadSelectionGen")),
    maxdRHadMuon_(iConfig.getParameter<double>("maxdRHadMuon")),
    mindRHadMuon_(iConfig.getParameter<double>("mindRHadMuon")),
    maxdRPhotonDs_(iConfig.getParameter<double>("maxdRPhotonDs")),
    maxdzDiffHadMuon_(iConfig.getParameter<double>("maxdzDiffHadMuon")),
    maxdxyHadPv_(iConfig.getParameter<double>("maxdxyHadPv")),
    phiMassAllowance_(iConfig.getParameter<double>("phiMassAllowance")),
    dsMassAllowance_(iConfig.getParameter<double>("dsMassAllowance")),
    dsStarMassAllowance_(iConfig.getParameter<double>("dsStarMassAllowance")),
    drMatchGen_(iConfig.getParameter<double>("drMatchGen")),
    maxBsMass_(iConfig.getParameter<double>("maxBsMass")),

    piMass_(iConfig.getParameter<double>("piMass")),
    piMassSigma_(iConfig.getParameter<double>("piMassSigma")),
    kMass_(iConfig.getParameter<double>("kMass")),
    kMassSigma_(iConfig.getParameter<double>("kMassSigma")),
    phiMass_(iConfig.getParameter<double>("phiMass")),
    constrainPhiMass_(iConfig.getParameter<bool>("constrainPhiMass")),
    minPhiVtxProb_(iConfig.getParameter<double>("minPhiVtxProb")),
    dsMass_(iConfig.getParameter<double>("dsMass")),
    constrainDsMass_(iConfig.getParameter<bool>("constrainDsMass")),
    minDsVtxProb_(iConfig.getParameter<double>("minDsVtxProb")),
    dsStarMass_(iConfig.getParameter<double>("dsStarMass")),
    muMass_(iConfig.getParameter<double>("muMass")),
    muMassSigma_(iConfig.getParameter<double>("muMassSigma")),
    bsMass_(iConfig.getParameter<double>("bsMass")),
    isoCone_(iConfig.getParameter<double>("isoCone")),

    beamSpotTag(iConfig.getParameter<edm::InputTag>("offBeamSpot")),
    beamSpot_(consumes<reco::BeamSpot>(beamSpotTag)),

    srcTag(iConfig.getParameter<edm::InputTag>("pfCand")),
    src_(consumes<pat::PackedCandidateCollection>(srcTag)), 

    trgMuonTag(iConfig.getParameter<edm::InputTag>("muCand")),
    trgMuons_(consumes<pat::MuonCollection>(trgMuonTag)), 

    primaryVtxTag(iConfig.getParameter<edm::InputTag>("pvCand")),
    primaryVtx_(consumes<reco::VertexCollection>(primaryVtxTag)),

    isoTracksTag(iConfig.getParameter<edm::InputTag>("tracks")),
    isoTracks_(consumes<pat::PackedCandidateCollection>(isoTracksTag)),

    tracksLostTag(iConfig.getParameter<edm::InputTag>("lostTracks")),
    tracksLost_(consumes<pat::PackedCandidateCollection>(tracksLostTag)),

    prunedGenTag(iConfig.getParameter<edm::InputTag>("prunedCand")),
    prunedGen_(consumes<reco::GenParticleCollection>(prunedGenTag)),
    packedGenTag(iConfig.getParameter<edm::InputTag>("packedCand")),
    packedGen_(consumes<pat::PackedGenParticleCollection>(packedGenTag)),
 
    photonTag(iConfig.getParameter<edm::InputTag>("photonCand")),
    photon_(consumes<reco::PhotonCoreCollection>(photonTag)),

    unpackedTracksTag(iConfig.getParameter<edm::InputTag>("unpackedTracksCand")),
    unpackedTracks_(consumes<reco::TrackCollection>(unpackedTracksTag))
 
    //refitVtxTag(iConfig.getParameter<edm::InputTag>("refitVtxCand")),
    //refitVtx_(consumes<reco::VertexCollection>(refitVtxTag))
 
   {
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
 
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpot_, beamSpotHandle);
 
  edm::Handle<pat::MuonCollection> trgMuons;
  iEvent.getByToken(trgMuons_,trgMuons);
 
  edm::Handle<reco::VertexCollection> primaryVtx;
  iEvent.getByToken(primaryVtx_,primaryVtx);

  edm::Handle<pat::PackedCandidateCollection> isoTracks;
  iEvent.getByToken(isoTracks_,isoTracks);

  edm::Handle<pat::PackedCandidateCollection> tracksLost;
  iEvent.getByToken(tracksLost_,tracksLost);

  edm::Handle<reco::GenParticleCollection> prunedGen;
  iEvent.getByToken(prunedGen_,prunedGen);

  edm::Handle<pat::PackedGenParticleCollection> packedGen;
  iEvent.getByToken(packedGen_,packedGen);

  edm::Handle<reco::PhotonCoreCollection> photon;
  iEvent.getByToken(photon_,photon);

  edm::Handle<reco::TrackCollection> unpackedTracks;
  iEvent.getByToken(unpackedTracks_, unpackedTracks);

  //edm::Handle<reco::VertexCollection> refitVtx;
  //iEvent.getByToken(refitVtx_, refitVtx);


  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> bsCandidates(new pat::CompositeCandidateCollection());
  //std::unique_ptr<pat::CompositeCandidateCollection> ret_value_gen(new pat::CompositeCandidateCollection());
  //std::unique_ptr<TransientTrackCollection> kkpi_ttrack(new TransientTrackCollection);

  //std::cout << "---------------- NEW EVENT ---------------" << std::endl;

  pat::PackedCandidateCollection mergedTracks;
  //append lost tracks to packed tracks
  for (auto& normCand : *pcand) {
    mergedTracks.push_back(normCand);
  }
  for (auto& lostCand : *tracksLost) {
    mergedTracks.push_back(lostCand);
  }



  //////////////////////////////////////////////////////
  // Match the trigger muon with a muon from the      //
  // packed pat collection. Reminder, we often have   //
  // only one trigger muon, but sometimes two! :)     //
  // --> checked, and its valid to have more than one //
  //////////////////////////////////////////////////////

  for(size_t trgMuIdx = 0; trgMuIdx < trgMuons->size(); ++trgMuIdx){

    nMuons++;
    //std::cout << "muon loop" << std::endl; 
    //if there is no trg muon, this loop is empty:)
    edm::Ptr<pat::Muon> muPtr(trgMuons, trgMuIdx);

    std::vector<float> dzMuPV;

    // For the first candidate selection, we pick the PV 
    // to be the one closest to the trg Muon in dz.
    // The final PV is calculated later.
    // More accurate for mu than for tau signals

    float dummy = 1.0;
    int pvIdx_dz= -1;
    reco::Vertex pv_dz;

    for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){
      nPv++;
      edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, vtxIdx);
      float minDz = fabs(muPtr->bestTrack()->dz(vtxPtr->position())); 
      if(minDz < dummy){
        dummy = minDz;
        pvIdx_dz = vtxIdx;
      }

    }

    if (pvIdx_dz <0) continue;
    muSelCounter++;
    if (pvIdx_dz >= 0){
    pv_dz = primaryVtx->at(pvIdx_dz);
    }

   


    //if (true ) std::cout << "found mu: " << muPtr->pt() << std::endl;
    //////////////////////////////////////////////////
    // Loop over k1 and select the good tracks      //
    //////////////////////////////////////////////////

    for(size_t k1Idx = 0; k1Idx < pcand->size() + tracksLost->size() ; ++k1Idx) {

      nTracks++;

      //define a pointer to the kaon at position k1Idx
      edm::Ptr<pat::PackedCandidate> k1Ptr;
      if (k1Idx < pcand->size()) k1Ptr = edm::Ptr<pat::PackedCandidate>(pcand, k1Idx); //normal tracks
      else k1Ptr = edm::Ptr<pat::PackedCandidate>(tracksLost, k1Idx - pcand->size());  //lost tracks


      //if (k1Ptr->pt() > 2.05273 & k1Ptr->pt() < 2.05274) //std::cout << "possible k1: " << k1Ptr->pt() << std::endl;

      if (!hadSelection_(*k1Ptr)) continue; 
      k1Sel1Counter++;

      //std::cout << "passed had selection" << k1Ptr->pt() << std::endl;
      //the PF algorithm assigns a pdgId hypothesis, generall it distinguishes between:
      // photons, electron/muon, charged hadron, neutral hadrons
      // and we trust the algorithm that when it says its an electron (11) or muon (13), that it is not a kaon or pion

      float muonK1dR = reco::deltaR(*k1Ptr,*muPtr);

      bool k1Sel = (( muonK1dR < maxdRHadMuon_ ) && 
      (reco::deltaR(*k1Ptr, *muPtr) > mindRHadMuon_) && 
      (abs(k1Ptr->bestTrack()->dz(pv_dz.position()) - muPtr->bestTrack()->dz(pv_dz.position())) < maxdzDiffHadMuon_)) &&
      (abs(k1Ptr->bestTrack()->dxy(pv_dz.position()) < maxdxyHadPv_ )) ;

      if (!k1Sel) continue;
      //if (true) std::cout << " found k1: " << k1Ptr->pt() << std::endl;
      k1Sel2Counter++;
      //////////////////////////////////////////////////
      // Loop over k2 and select the good tracks      //
      //////////////////////////////////////////////////

      for(size_t k2Idx = k1Idx + 1; k2Idx < pcand->size()+ tracksLost->size() ; ++k2Idx) {

      //make sure k2 is not k1
      if (k2Idx == k1Idx) continue;

      edm::Ptr<pat::PackedCandidate> k2Ptr;
      if (k2Idx < pcand->size()) k2Ptr = edm::Ptr<pat::PackedCandidate>(pcand, k2Idx); //normal tracks
      else k2Ptr = edm::Ptr<pat::PackedCandidate>(tracksLost, k2Idx - pcand->size());  //lost tracks
    
      // if this kaon does not pass the selection, jump to the next!
      if(!hadSelection_(*k2Ptr)) continue;
      k2Sel1Counter++;

      float muonK2dR = reco::deltaR(*k2Ptr,*muPtr);

      bool k2Sel = (( muonK2dR < maxdRHadMuon_ ) && 
      (reco::deltaR(*k2Ptr, *muPtr) > mindRHadMuon_) && 
      (abs(k2Ptr->bestTrack()->dz(pv_dz.position()) - muPtr->bestTrack()->dz(pv_dz.position())) < maxdzDiffHadMuon_)) &&
      (abs(k2Ptr->bestTrack()->dxy(pv_dz.position()) < maxdxyHadPv_ )) ;

      //std::cout << "passed had selection" << k1Ptr->pt() << std::endl;
      //k1 and k2 must have oppoiste charge -> only for signal tests, later we keep everything
      int kkCharge = k1Ptr->charge() * k2Ptr->charge();

      //if (kkCharge > 0) continue; //To be commented out 


      if (!k2Sel) continue;
      //if (true) std::cout << " --- found k2: " << k2Ptr->pt() << std::endl;
      k2Sel2Counter++;

      //////////////////////////////////////////////////
      // Loop over pi and select the good tracks      //
      //////////////////////////////////////////////////

      for(size_t piIdx = 0; piIdx < pcand->size()+ tracksLost->size() ; ++piIdx) {

        //make sure the pion is none of the kaons:
        if((piIdx == k1Idx) || (piIdx == k2Idx)) continue;

        edm::Ptr<pat::PackedCandidate> piPtr;
        if (piIdx < pcand->size()) piPtr = edm::Ptr<pat::PackedCandidate>(pcand, piIdx); //normal tracks
        else piPtr = edm::Ptr<pat::PackedCandidate>(tracksLost, piIdx - pcand->size());  //lost tracks

        // if this pion does not pass the selection, jump to the next!
        if(!hadSelection_(*piPtr)) continue;
        piSel1Counter++;

        float muonPidR = reco::deltaR(*piPtr,*muPtr);

        bool piSel = ((muonPidR < maxdRHadMuon_) && 
        (reco::deltaR(*piPtr, *muPtr) > mindRHadMuon_) &&
        (abs(piPtr->bestTrack()->dz(pv_dz.position()) - muPtr->bestTrack()->dz(pv_dz.position())) < maxdzDiffHadMuon_)) &&
        (abs(piPtr->bestTrack()->dxy(pv_dz.position()) < maxdxyHadPv_ )) ;

        //std::cout << "passed had selection" << k1Ptr->pt() << std::endl;
        //pi and mu must have opposite charge -> only for signal tests, later we keep everything
        int piMuCharge = piPtr->charge() * muPtr->charge();
        
        //if (piMuCharge > 0) continue; //To be commented out
        //std::cout << piPtr->pt() <<"and" <<piPtr->pdgId() << std::endl;
        if (!piSel) continue;
        piSel2Counter++;
	

        //std::cout << " ----- found pi: " << piPtr->pt() << std::endl;
        //////////////////////////////////////////////////
        // Build Phi resonance                          //
        //////////////////////////////////////////////////

        nKKPiMu++; // found candidate

        //define a composite candidate pair 
        pat::CompositeCandidate kk;

        //PF canidates are always assigned with the pi mass, we have to force the kaon mass
        math::PtEtaPhiMLorentzVector k1P4(k1Ptr->pt(), k1Ptr->eta(), k1Ptr->phi(), K_MASS);
        math::PtEtaPhiMLorentzVector k2P4(k2Ptr->pt(), k2Ptr->eta(), k2Ptr->phi(), K_MASS);
        math::PtEtaPhiMLorentzVector piP4(piPtr->pt(), piPtr->eta(), piPtr->phi(), PI_MASS); //just to be sure lets also force the pi mass
 
        kk.setP4(k1P4 + k2P4);
      
        //std::cout << "found kkpi candidate with pt: mu k1 k2 pi and mass (KK) "  << std::endl; 
        //std::cout << muPtr->pt() << std::endl;
        //std::cout << k1Ptr->pt() << std::endl;
        //std::cout << k2Ptr->pt() << std::endl;
        //std::cout << piPtr->pt() << std::endl;
        //std::cout << kk.mass() << std::endl;
 

        //only continue when they build a phi resonance, allow 15MeV:
        if (fabs(kk.mass() - phiMass_) > phiMassAllowance_) continue;     
        //std::cout << "we passed the phi resonance" << std::endl; 
        kk.setCharge(k1Ptr->charge() + k2Ptr->charge());
        nPhiMassCut++;

        //////////////////////////////////////////////////
        // Build Ds resonance                           //
        // ///////////////////////////////////////////////

        pat::CompositeCandidate phiPi;
        phiPi.setP4(kk.p4() + piP4); 

        //only continue when they build a ds resonance, allow 50MeV:
        if (fabs(phiPi.mass() - dsMass_) > dsMassAllowance_) continue;

        //std::cout << "we passed the ds resonance" << std::endl; 
        phiPi.setCharge(kk.charge() + piPtr->charge());
        //std::cout << "found ds resonance" << std::endl;
        nDsMassCut++;



        //////////////////////////////////////////////////
        // Build Bs resonance                           //
        //////////////////////////////////////////////////

        pat::CompositeCandidate dsMu;
        dsMu.setP4(phiPi.p4() + muPtr->p4()); 
        dsMu.setCharge(phiPi.charge() + muPtr->charge()); //sanity check:shoould be 0

        if(dsMu.mass() > maxBsMass_) continue;
        nBsMassCut++;

        //std::cout << "we passed the bs cut" << std::endl; 
        //std::cout << "and have ds mass"<< dsMu.mass() << std::endl; 
        //build bs with collinear approximation
        pat::CompositeCandidate bs;

        bs.setP4(dsMu.p4() * bsMass_ / dsMu.mass()); //the bs_mass will thus be fixed at 536688 (peak in the histo)
        bs.setCharge(dsMu.charge());

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
        ParticleMass piMass  = piMass_;
        ParticleMass kMass   = kMass_;
        ParticleMass phiMass = phiMass_;
        ParticleMass dsMass  = dsMass_;
        ParticleMass muMass  = muMass_;

        
        float ndf = 0.0;
        float chi = 0.0;
        float kMassSigma = kMassSigma_ ;
        float piMassSigma = piMassSigma_;
        float muMassSigma = muMassSigma_;

        // fix the tracks to have pos def covariance matrix, taken from:
        // https://github.com/CMSKStarMuMu/miniB0KstarMuMu/blob/master/miniKstarMuMu/plugins/miniKstarMuMu.cc#L1611-L1678

        const reco::Track k1Track = *k1Ptr->bestTrack();
        const reco::Track k2Track = *k2Ptr->bestTrack();
        const reco::Track piTrack = *piPtr->bestTrack();
        const reco::Track muTrack = *muPtr->bestTrack();

        reco::Track k1TrackCorr   = correctCovMat(&k1Track, 1e-8);
        reco::Track k2TrackCorr   = correctCovMat(&k2Track, 1e-8);
        reco::Track piTrackCorr   = correctCovMat(&piTrack, 1e-8);
        reco::Track muTrackCorr   = correctCovMat(&muTrack, 1e-8);

        reco::TransientTrack ttK1 = ttBuilder->build(k1TrackCorr);
        reco::TransientTrack ttK2 = ttBuilder->build(k2TrackCorr);

        phiToFit.push_back( pFactory.particle(ttK1, kMass,  chi, ndf, kMassSigma ));
        phiToFit.push_back( pFactory.particle(ttK2, kMass,  chi, ndf, kMassSigma ));
        dsToFit.push_back(  pFactory.particle(getTransientTrack(piTrackCorr), piMass, chi, ndf, piMassSigma));
        bsToFit.push_back(  pFactory.particle(getTransientTrack(muTrackCorr), muMass, chi, ndf, muMassSigma));

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


        // PHI VERTEX FIT

        // baby phi fit
        //KinVtxFitter easyFitter(
        //{ttK1, ttK2, getTransientTrack(piTrackCorr)},
        //{kMass, kMass, piMass},
        //{kMassSigma,kMassSigma,piMassSigma}
        //);
        //if(!easyFitter.success()) continue;
        //std::cout << "ds easy fit chi2: " << easyFitter.chi2() << std::endl;
        //std::cout << "ds easy fit ndof: " << easyFitter.dof() << std::endl;
        //std::cout << "ds easy fit prob: " << easyFitter.prob() << std::endl;
        //std::cout << "ds easy fit vx " << easyFitter.fitted_vtx().x() << std::endl;
        


        RefCountedKinematicTree phiTree  = vertexFit(phiToFit, phiMass, constrainPhiMass_);

        if (!phiTree->isValid() || phiTree->isEmpty() || !phiTree->isConsistent()) continue; //check if fit result is valid

        //access the fitted resonance and vertex 
        phiTree->movePointerToTheTop();
        RefCountedKinematicParticle phiParticle = phiTree->currentParticle();

        //get vtx chi2 and ndof

        auto phiVtx = phiTree->currentDecayVertex();
        if (!phiVtx->vertexIsValid() || !phiParticle->currentState().isValid() ) continue; //check if fit result is valid

        float phiVtxChi2    = phiVtx->chiSquared();
        //std::cout << "phi chi2:" << phiVtxChi2 << std::endl;
        if (phiVtxChi2 < 0) continue;

        float phiVtxNDof    = phiVtx->degreesOfFreedom();
        float phiVtxRedChi2 = phiVtxChi2 / phiVtxNDof; 
        float phiVtxProb    = ChiSquaredProbability(phiVtxChi2, phiVtxNDof); 
        //std::cout << "phi ndof:" << phiVtxNDof << std::endl;
        //std::cout << "phi prob:" << phiVtxProb << std::endl;
        if (phiVtxProb < 0.01) continue;
        //std::cout << "we passed the phi vtx fit prob" << std::endl; 
        //std::cout << "phi prob:" << ChiSquaredProbability(phiVtxChi2,phiVtxNDof) << std::endl;
        
        nPhiFit++;

        // access refitted children
        phiTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle phiDau1 = phiTree->currentParticle();
        phiTree->movePointerToTheNextChild();
        RefCountedKinematicParticle phiDau2 = phiTree->currentParticle();

        // get the vectors full of fit information (vertex and momenta and mass)
        AlgebraicVector7 phiParams     = phiParticle->currentState().kinematicParameters().vector();     
        AlgebraicVector7 phiDau1Params = phiDau1->currentState().kinematicParameters().vector();     
        AlgebraicVector7 phiDau2Params = phiDau2->currentState().kinematicParameters().vector();
     
        //std::cout << "phi fit vx " << phiParams(0) << std::endl;
        // add the phi to the list of particles (pi) to fit the ds 
        dsToFit.push_back(phiParticle);

        // DS VERTEX FIT

        RefCountedKinematicTree dsTree = vertexFit(dsToFit, dsMass, constrainDsMass_);
        if (!dsTree->isValid() || dsTree->isEmpty() ) continue; //check if fit result is valid

        // access the fitted resonance and the refitted children
        dsTree->movePointerToTheTop();
        RefCountedKinematicParticle dsParticle = dsTree->currentParticle();

        // get vtx chi2 and ndof
        RefCountedKinematicVertex dsVtx = dsTree->currentDecayVertex(); //compare to the access via AlgebraicVector7
        if (!dsVtx->vertexIsValid()) continue; //check if fit result is valid

        float dsVtxChi2    = dsVtx->chiSquared();
        //std::cout << "chi2 of ds fit is:" << dsVtxChi2 << std::endl;

        //std::cout << muPtr->pt() << std::endl;
        //std::cout << k1Ptr->pt() << std::endl;
        //std::cout << k2Ptr->pt() << std::endl;
        //std::cout << piPtr->pt() << std::endl;

        if (dsVtxChi2 < 0) continue;
        float dsVtxNDof    = dsVtx->degreesOfFreedom();
        float dsVtxRedChi2 = dsVtxChi2 / dsVtxNDof; 
        float dsVtxProb    = ChiSquaredProbability(dsVtxChi2, dsVtxNDof); 
        //std::cout << "ds chi2:" << dsVtxChi2 << std::endl;
        //std::cout << "ds ndof:" << dsVtxNDof << std::endl;
        //std::cout << "ds prob:" << dsVtxProb << std::endl;

        if (dsVtxProb < 0.01) continue;
        //std::cout << "ds chi2:" << dsVtxChi2 << std::endl;
        //std::cout << "ds prob:" << ChiSquaredProbability(dsVtxChi2,dsVtxNDof) << std::endl;

        nDsFit++;

        //std::cout << "we passed the ds vtx fit prob" << std::endl; 
        dsTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle dsDau1 = dsTree->currentParticle();
        dsTree->movePointerToTheNextChild();
        RefCountedKinematicParticle dsDau2 = dsTree->currentParticle();

        // get the vectors full of fit information (vertex and momenta and mass)
        AlgebraicVector7 dsParams     = dsParticle->currentState().kinematicParameters().vector();     
        AlgebraicVector7 dsDau1Params = dsDau1->currentState().kinematicParameters().vector();     
        AlgebraicVector7 dsDau2Params = dsDau2->currentState().kinematicParameters().vector();

        // add the ds to the list of particles to fit the bs
        bsToFit.push_back(dsParticle);
        
        // BS VERTEX FIT
        KinematicConstrainedVertexFitter bsFitter;

        RefCountedKinematicTree bsTree = bsFitter.fit(bsToFit); // no constraint for bs because missing momentum
        if (!bsTree->isValid() || bsTree->isEmpty() ) continue; //check if fit result is valid 
        //std::cout << "we passed the bs tree" << std::endl;
        // access the fitted resonance and the refitted children
        bsTree->movePointerToTheTop();
        RefCountedKinematicParticle bsParticle = bsTree->currentParticle();

        // get vtx chi2 and ndof
        RefCountedKinematicVertex bsVtx = bsTree->currentDecayVertex();
        if (!bsVtx->vertexIsValid()) continue; //check if fit result is valid

        //std::cout << "bs vtx valid " << std::endl;
        float bsVtxChi2    = bsVtx->chiSquared();
        //std::cout << "bs chi2:" << bsVtxChi2 << std::endl;
        if (bsVtxChi2 < 0) continue;
        float bsVtxNDof    = bsVtx->degreesOfFreedom();
        float bsVtxRedChi2 = bsVtxChi2 / bsVtxNDof; 
        float bsVtxProb    = ChiSquaredProbability(bsVtxChi2, bsVtxNDof); 
 
        nBsFit++;

        //std::cout << "bs chi2 > 0 " << std::endl;
        //std::cout << "bs chi2:" << bsVtxChi2 << std::endl;
        //std::cout << "bs ndof:" << bsVtxNDof << std::endl;
        //std::cout << "bs prob:" << bsVtxProb << std::endl;
        //std::cout << "bs prob:" << ChiSquaredProbability(bsVtxChi2,bsVtxNDof) << std::endl;

        bsTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle bsDau1 = bsTree->currentParticle();
        bsTree->movePointerToTheNextChild();
        RefCountedKinematicParticle bsDau2 = bsTree->currentParticle();

        // get the vectors full of fit information (vertex and momenta and mass)
        AlgebraicVector7 bsParams = bsParticle->currentState().kinematicParameters().vector();     
        AlgebraicVector7 bsDau1Params = bsDau1->currentState().kinematicParameters().vector();     
        AlgebraicVector7 bsDau2Params = bsDau2->currentState().kinematicParameters().vector();

        //std::cout << "we passed all the vtx fit" << std::endl;

        //////////////////////////////////// end of global fitter /////////////////////////////////////


        ///////////////////////////////////////////////////
        // Match final state track with AOD tracks       //
        ///////////////////////////////////////////////////
  
        
        int foundTrackMu    = 0;
        int foundTrackK1    = 0;
        int foundTrackK2    = 0;
        int foundTrackPi    = 0;
  
        int trkIdxMu        = -1;
        int trkIdxK1        = -1;
        int trkIdxK2        = -1;
        int trkIdxPi        = -1;

        float minDeltaDrMu  = 0;
        float minDeltaDrK1  = 0;
        float minDeltaDrK2  = 0;
        float minDeltaDrPi  = 0;
       
        for(size_t trkIdx = 0; trkIdx < unpackedTracks->size(); ++trkIdx) {

          //edm::Ptr<reco::Track> trkPtr;
          //trkPtr = edm::Ptr<reco::Track>(unpackedTracks, trkIdx); // aod tracks 

          edm::Ref<reco::TrackCollection> trkRef(unpackedTracks, trkIdx);
          const reco::Track& track = *trkRef; // Access the track object

          /* 
          double drMu = reco::deltaR(*trkPtr, *muPtr); 
          double drK1 = reco::deltaR(*trkPtr, *k1Ptr); 
          double drK2 = reco::deltaR(*trkPtr, *k2Ptr); 
          double drPi = reco::deltaR(*trkPtr, *piPtr); 

          if (((drMu < minDeltaDrMu) || (foundTrackMu == 0)) && ( muPtr->charge() * trkPtr->charge() > 0)){
            
            foundTrackMu = 1;
            minDeltaDrMu = drMu;
            trkIdxMu     = trkIdx;

          }
          if (((drK1 < minDeltaDrK1) || (foundTrackK1 == 0)) && ( k1Ptr->charge() * trkPtr->charge() > 0)){
            
            foundTrackK1 = 1;
            minDeltaDrK1 = drK1;
            trkIdxK1     = trkIdx;

          }
          if (((drK2 < minDeltaDrK2) || (foundTrackK2 == 0)) && ( k2Ptr->charge() * trkPtr->charge() > 0)){
            
            foundTrackK2 = 1;
            minDeltaDrK2 = drK2;
            trkIdxK2     = trkIdx;

          }
          if (((drPi < minDeltaDrPi) || (foundTrackPi == 0)) && ( piPtr->charge() * trkPtr->charge() > 0)){
            
            foundTrackPi = 1;
            minDeltaDrPi = drPi;
            trkIdxPi     = trkIdx;

          }

          */
        }

        //std::set<int> signalTrkIdx = {trkIdxMu, trkIdxK1, trkIdxK2, trkIdxPi};

        //edm::Ptr<reco::Track> aodTrkMu;
        //edm::Ptr<reco::Track> aodTrkK1;
        //edm::Ptr<reco::Track> aodTrkK2;
        //edm::Ptr<reco::Track> aodTrkPi;

        //aodTrkMu = edm::Ptr<reco::Track>(unpackedTracks, trkIdxMu); 
        //aodTrkK1 = edm::Ptr<reco::Track>(unpackedTracks, trkIdxK1); 
        //aodTrkK2 = edm::Ptr<reco::Track>(unpackedTracks, trkIdxK2); 
        //aodTrkPi = edm::Ptr<reco::Track>(unpackedTracks, trkIdxPi); 


        // now fill a vector of aod tracks, except for the signal tracks
        // want to remove signal tracks from vtx fit

        //std::vector<reco::Track> toFit;

        //for(size_t trkIdx = 0; trkIdx < unpackedTracks->size(); ++trkIdx) {

        //  //skip the signal tracks 
        //  if (std::find(signalTrkIdx.begin(), signalTrkIdx.end(), trkIdx) != signalTrkIdx.end()) continue; 
  
        //  //only add other tracks for the fit
        //  edm::Ptr<reco::Track> trkPtr;
        //  trkPtr = edm::Ptr<reco::Track>(unpackedTracks, trkIdx); // aod tracks 
        //  toFit.push_back(*trkPtr);

        //}


        //////////////////////////////////////////////////
        // Look for possible photons (g for gamma)      //
        //////////////////////////////////////////////////

        int foundPhoton       = 0;
        float minDeltaDsStar  = 0;
        int photonIdx         = -1;

        float dr_photon_ds = -9999; 
        float dsPhoton_m   = -9999;
        float dsPhotonMu_m = -9999;
        float photon_pt    = -9999; 
        float photon_eta   = -9999; 
        float photon_phi   = -9999; 
        int   photon_pdgid = -9999; 

        TLorentzVector photonTlv(std::nan(""),std::nan(""),std::nan(""),std::nan(""));
        edm::Ptr<pat::PackedCandidate> phtPtr;


        int photonMultiplicity = 0;

        for (size_t gIdx = 0; gIdx < pcand->size() + tracksLost->size() ; ++gIdx){

          edm::Ptr<pat::PackedCandidate> gPtr;
          if (gIdx < pcand->size()) gPtr = edm::Ptr<pat::PackedCandidate>(pcand, gIdx); //normal tracks
          else gPtr = edm::Ptr<pat::PackedCandidate>(tracksLost, gIdx - pcand->size());  //lost tracks 
  
          // if the photon does not pass the selection, jump to the next!
          if(!gSelection_(*gPtr)) continue;
        
          photonMultiplicity++; 
          gSelCounter++;
          
          // define photon momentum
          //math::PtEtaPhiMLorentzVector photonP4(gPtr->pt(), gPtr->eta(), gPtr->phi(), 0.0);
          pat::CompositeCandidate dsStar;

          dsStar.setP4(gPtr->p4() + phiPi.p4());
 
          float dr           = reco::deltaR(*gPtr, phiPi);
          float mDsStar      = dsStar.mass();
          float deltaDsStar  = abs(mDsStar - DSSTAR_MASS); 
 
          if (((dr < minDeltaDsStar) || (foundPhoton == 0))  && ( deltaDsStar < dsStarMassAllowance_)){

            foundPhoton    = 1;
            minDeltaDsStar = dr; 
            photonIdx      = gIdx;        
            phtPtr         = gPtr;
            dr_photon_ds   = dr;       
            dsPhoton_m     = mDsStar;
 
          }
        }   

        if (photonIdx != -1){

          photon_pt    = phtPtr->pt();     
          photon_eta   = phtPtr->eta();     
          photon_phi   = phtPtr->phi();     
          photon_pdgid = phtPtr->pdgId();     

          pat::CompositeCandidate dsStarMu;
          dsStarMu.setP4(phtPtr->p4() + phiPi.p4() + muPtr->p4());
          dsPhotonMu_m   = dsStarMu.mass();

          photonTlv.SetXYZM(phtPtr->px(), phtPtr->py(), phtPtr->pz(), 0.0);

        }

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

        // refitted ds+mu 
        TLorentzVector refittedDsMu = refittedDs + refittedMu ;


        // set sv,tv, and fv from fit
        float sv_x = bsParams(0);
        float sv_y = bsParams(1);
        float sv_z = bsParams(2);

        float tv_x = dsParams(0);
        float tv_y = dsParams(1);
        float tv_z = dsParams(2);

        float fv_x = phiParams(0);
        float fv_y = phiParams(1);
        float fv_z = phiParams(2);

        ///////////////////////////////////////////////////////////////
        // Redefinition of PV, this is more excact! (Around 5%)      //
        // Use this definition for all the variables depending on PV //
        ///////////////////////////////////////////////////////////////

        float dummy2   = 1.0;
        int pvIdx      = -1;
        reco::Vertex pv;
 
        for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){
          edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, vtxIdx);

          float ip3d = getIP3D(vtxPtr, sv_x, sv_y, sv_z, refittedDs, refittedMu);

          if(ip3d < dummy2){
            dummy2 = ip3d;
            pvIdx  = vtxIdx;
           }

         }

        if (pvIdx <  0) continue;
        if (pvIdx >= 0){
          pv = primaryVtx->at(pvIdx);
        }


        // get the x,y position of the pv by evaluating the beamspot density
        // at pv_z

        if (!beamSpotHandle.isValid()) continue;
        const reco::BeamSpot& beamSpot = *beamSpotHandle;

        float pv_z = pv.z();
        float pv_x = beamSpot.x(pv_z);
        float pv_y = beamSpot.y(pv_z);
       

        //////////////////////////////////////////////

        //2d cosine

        ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> kkPi_xy;
        ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> ds_xy;
        kkPi_xy.SetXYZ(phiPi.px(), phiPi.py(),0.0);
        ds_xy.SetXYZ(tv_x - beamSpot.x(pv_z), tv_y - beamSpot.y(pv_z) ,0.0);

        if (!(ds_xy.R() > 0.0)) continue;
         
        float ds_vtx_cosine = kkPi_xy.Dot(ds_xy) / (kkPi_xy.R() * ds_xy.R());
        //std::cout << "ds vtc cosine" << ds_vtx_cosine << std::endl;

        if (ds_vtx_cosine < 0.8) continue;

        //std::cout << "survived cos 2d cut" << std::endl;

        ///////////////////////////////////////////////////////////////////////////////////////////
        // store basic variables                                                                 //
        ///////////////////////////////////////////////////////////////////////////////////////////

        //save the indices of the final states in the pruned Collection
        bs.addUserInt("k1_idx",k1Idx);
        bs.addUserInt("k2_idx",k2Idx);
        bs.addUserInt("pi_idx",piIdx);

        //bs.addUserInt("mu_idx",muIdx); //always 0 :)

        //add final states as Candidates
        bs.addUserCand("k1",k1Ptr);  //be aware that this ptr has the wrong mass, need to assign it in variables_cff.py
        bs.addUserFloat("k1_m",kMass_);

        bs.addUserCand("k2",k2Ptr);  // "
        bs.addUserFloat("k2_m",kMass_);

        bs.addUserCand("pi",piPtr);
        bs.addUserFloat("pi_m",piMass_);

        bs.addUserCand("mu",muPtr); //muon mass is correct -> check
        bs.addUserFloat("mu_m",muMass_);

        TLorentzVector muTlv;
        muTlv.SetXYZM(muPtr->px(), muPtr->py(), muPtr->pz(), muMass_);

        // add information about trg muon
        bs.addUserInt("mu_is_tracker",     muPtr->isTrackerMuon());
        bs.addUserInt("mu_is_pf",          muPtr->isPFMuon());
        bs.addUserInt("mu_is_global",      muPtr->isGlobalMuon());
        bs.addUserInt("mu_is_standalone",  muPtr->isStandAloneMuon());

        bs.addUserFloat("photon_pt",    photon_pt);
        bs.addUserFloat("photon_eta",   photon_eta);
        bs.addUserFloat("photon_phi",   photon_phi);
        bs.addUserInt("photon_pdgid", photon_pdgid);
        bs.addUserFloat("dr_photon_ds", dr_photon_ds);
        bs.addUserFloat("dsPhoton_m", dsPhoton_m);
        bs.addUserFloat("dsPhotonMu_m", dsPhotonMu_m);
        bs.addUserInt("foundPhoton", foundPhoton);
        bs.addUserInt("photonMultiplicity", photonMultiplicity);

        //add prefit resonances --> are not part of collection and can thus not be
        //automatically access the pt(), .. as I can for k1,k2,pi,mu in the variables_cff.py

        bs.addUserFloat("kk_pt", kk.pt());
        bs.addUserFloat("kk_eta", kk.eta());
        bs.addUserFloat("kk_phi", kk.phi());
        bs.addUserFloat("kk_m", kk.mass());
        bs.addUserFloat("kk_charge", kk.charge());
        bs.addUserFloat("kk_deltaR", reco::deltaR(*k1Ptr, *k2Ptr));

        bs.addUserFloat("phiPi_pt", phiPi.pt());
        bs.addUserFloat("phiPi_eta", phiPi.eta());
        bs.addUserFloat("phiPi_phi", phiPi.phi());
        bs.addUserFloat("phiPi_m", phiPi.mass());
        bs.addUserFloat("phiPi_charge", phiPi.charge());
        bs.addUserFloat("phiPi_deltaR", reco::deltaR(kk, *piPtr));

        bs.addUserFloat("dsMu_pt", dsMu.pt());
        bs.addUserFloat("dsMu_eta", dsMu.eta());
        bs.addUserFloat("dsMu_phi", dsMu.phi());
        bs.addUserFloat("dsMu_m", dsMu.mass()); //we miss the neutrino contribution
        bs.addUserFloat("dsMu_charge", dsMu.charge());
        bs.addUserFloat("dsMu_deltaR", reco::deltaR(phiPi, *muPtr));

        //rel charges (not additive!)
        bs.addUserFloat("pi_mu_charge",piMuCharge); 
        bs.addUserFloat("k_k_charge"  ,kkCharge); 


        ///////////////////////////////////////////////////////////////////////////////////////////
        // store  ID's                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////

        // for ID numbering check: https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonReco/interface/Muon.h

        bs.addUserInt("mu_id_loose",               muPtr->passed(0)); // 
        bs.addUserInt("mu_id_medium",              muPtr->passed(1)); // 
        bs.addUserInt("mu_id_medium_prompt",       muPtr->passed(2)); // 
        bs.addUserInt("mu_id_tight",               muPtr->passed(3)); // 
        bs.addUserInt("mu_id_global_high_pt",      muPtr->passed(4)); // 
        bs.addUserInt("mu_id_trk_high_pt",         muPtr->passed(5)); // 
        bs.addUserInt("mu_pf_iso_very_loose",      muPtr->passed(6)); // 
        bs.addUserInt("mu_pf_iso_loose",           muPtr->passed(7)); // 
        bs.addUserInt("mu_pf_iso_medium",          muPtr->passed(8)); // 
        bs.addUserInt("mu_pf_iso_tight",           muPtr->passed(9)); // 
        bs.addUserInt("mu_pf_iso_very_tight",      muPtr->passed(10)); // 

        bs.addUserInt("mu_tk_iso_loose",           muPtr->passed(11)); // 
        bs.addUserInt("mu_tk_iso_tight",           muPtr->passed(12)); // 
        bs.addUserInt("mu_id_soft",                muPtr->passed(13)); // 
        bs.addUserInt("mu_id_soft_mva",            muPtr->passed(14)); // 
        bs.addUserInt("mu_mva_loose",              muPtr->passed(15)); // 
        bs.addUserInt("mu_mva_medium",             muPtr->passed(16)); // 
        bs.addUserInt("mu_mva_tight",              muPtr->passed(17)); // 
        bs.addUserInt("mu_mini_iso_loose",         muPtr->passed(18)); // 
        bs.addUserInt("mu_mini_iso_medium",        muPtr->passed(19)); // 
        bs.addUserInt("mu_mini_iso_tight",         muPtr->passed(20)); // 

        bs.addUserInt("mu_mini_iso_very_tight",    muPtr->passed(21)); // 
        bs.addUserInt("mu_trigger_id_loose",       muPtr->passed(22)); // 
        bs.addUserInt("mu_in_time_muon",           muPtr->passed(23)); // 
        bs.addUserInt("mu_pf_iso_very_very_tight", muPtr->passed(24)); // 
        bs.addUserInt("mu_multi_iso_loose",        muPtr->passed(25)); // 
        bs.addUserInt("mu_multi_iso_medium",       muPtr->passed(26)); // 
        bs.addUserInt("mu_puppi_iso_loose",        muPtr->passed(27)); // 
        bs.addUserInt("mu_puppi_iso_medium",       muPtr->passed(28)); // 
        bs.addUserInt("mu_puppi_iso_tight",        muPtr->passed(29)); // 
        bs.addUserInt("mu_mva_v_tight",            muPtr->passed(30)); // 

        bs.addUserInt("mu_mva_vv_tight",           muPtr->passed(31)); // 
        bs.addUserInt("mu_low_pt_mva_loose",       muPtr->passed(32)); // 
        bs.addUserInt("mu_low_pt_mva_medium",      muPtr->passed(33)); // 
        bs.addUserInt("mu_mv_id_wp_medium",        muPtr->passed(34)); // 
        bs.addUserInt("mu_mv_id_wp_tight",         muPtr->passed(35)); // 


        //AlgebraicMatrix77 phiErr = phiParticle->currentState().kinematicParametersError().matrix();
        //AlgebraicMatrix77 dsErr = dsParticle->currentState().kinematicParametersError().matrix();
        //AlgebraicMatrix77 bsErr = bsParticle->currentState().kinematicParametersError().matrix();
        //check if vertex position is truly at 0,1,2: result: Yes it is!
        //double xdecay = phiDecayVertex->position().x();

        // save fitter info
        // Remark: fitted describes the first fit, refitted describes the refit due to daughters.
        // TODO: In the end its sequential, and not fitting all at the same time.. maybe there is a possibility? -> keep looking around
        // TODO: Double check if the order of daughters is kept according to the **ToFit vectors -> include sanity checks

        ///////////////////////////////////////////////////////////////////////////////////////////
        // store  Vertices                                                                       //
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Remark: currentDecayVertex() returns the same as the ALgebraic Vector components 0,1,2--> checked :)
        // Remark: the daughters of a fitted vertex return the same vertex via Algebraic vector 0,1,2 like the
        // fitted mother --> checked 
   

        // Beam spot
        bs.addUserFloat("bs_x0",         beamSpot.x0());
        bs.addUserFloat("bs_y0",         beamSpot.y0());
        bs.addUserFloat("bs_z0",         beamSpot.z0());

        // primary vertex ( = Bs production vertex)

        bs.addUserFloat("pv_x",       pv_x); 
        bs.addUserFloat("pv_y",       pv_y); 
        bs.addUserFloat("pv_z",       pv_z); 
        bs.addUserFloat("pv_chi2",    pv.chi2()); // can access this directlyfor the reco::Vertex class
        bs.addUserFloat("pv_ndof",    pv.ndof());
        bs.addUserFloat("pv_redchi2", pv.normalizedChi2());
        bs.addUserFloat("pv_prob",    ChiSquaredProbability(pv.chi2(), pv.ndof()));
        bs.addUserInt  ("pv_idx",     pvIdx);

        // save also the x,y components without beamspot calculation 
        bs.addUserFloat("pv_x_wout_bs",  pv.x()); 
        bs.addUserFloat("pv_y_wout_bs",  pv.y()); 

        // The first definition using minimal dz between pv and muon
        bs.addUserFloat("pv_dz_x",  pv_dz.x()); 
        bs.addUserFloat("pv_dz_y",  pv_dz.y()); 
        bs.addUserFloat("pv_dz_z",  pv_dz.z()); 

        // secondary vertex ( = Bs decay vertex)  

        bs.addUserFloat("sv_x",       sv_x); 
        bs.addUserFloat("sv_y",       sv_y); 
        bs.addUserFloat("sv_z",       sv_z); 
        bs.addUserFloat("sv_chi2",    bsVtxChi2);
        bs.addUserFloat("sv_ndof",    bsVtxNDof);
        bs.addUserFloat("sv_redchi2", bsVtxRedChi2);
        bs.addUserFloat("sv_prob",    bsVtxProb);

        // tertiary vertex ( = Ds decay vertex) 

        bs.addUserFloat("tv_x",       tv_x); 
        bs.addUserFloat("tv_y",       tv_y); 
        bs.addUserFloat("tv_z",       tv_z); 
        bs.addUserFloat("tv_chi2",    dsVtxChi2);
        bs.addUserFloat("tv_ndof",    dsVtxNDof);
        bs.addUserFloat("tv_redchi2", dsVtxRedChi2);
        bs.addUserFloat("tv_prob",    dsVtxProb);

        // fourth vertex ( = Phi Decay Vertex)

        bs.addUserFloat("fv_x",       fv_x); 
        bs.addUserFloat("fv_y",       fv_y); 
        bs.addUserFloat("fv_z",       fv_z); 
        bs.addUserFloat("fv_chi2",    phiVtxChi2);
        bs.addUserFloat("fv_ndof",    phiVtxNDof);
        bs.addUserFloat("fv_redchi2", phiVtxRedChi2);
        bs.addUserFloat("fv_prob",    phiVtxProb);

        // opening angle between bs - ds vtx directioini and kkpi flight direction
        bs.addUserFloat("ds_vtx_cosine", ds_vtx_cosine);


        ///////////////////////////////////////////////////////////////////////////////////////////
        // store  Lxy, Dxy                                                                       //
        ///////////////////////////////////////////////////////////////////////////////////////////

        // lxy(z) is the flight distance in the xy(z) plane(space)

        Measurement1D lxyDsVect       = lxyz(bsVtx->vertexState(), dsVtx->vertexState(),  true);        
        Measurement1D lxyzDsVect      = lxyz(bsVtx->vertexState(), dsVtx->vertexState(),  false);        

        Measurement1D lxyPhiVect      = lxyz(dsVtx->vertexState(), phiVtx->vertexState(), true);        
        Measurement1D lxyzPhiVect     = lxyz(dsVtx->vertexState(), phiVtx->vertexState(), false);        

        Measurement1D lxyBsVect       = VertexDistanceXY().distance(  pv, dsVtx->vertexState());
        Measurement1D lxyzBsVect      = VertexDistance3D().distance( pv, dsVtx->vertexState()); 

        bs.addUserFloat("lxy_bs",     lxyBsVect.value());
        bs.addUserFloat("lxy_bs_err", lxyBsVect.error());
        bs.addUserFloat("lxy_bs_sig", lxyBsVect.significance());

        bs.addUserFloat("lxyz_bs",    lxyzBsVect.value());
        bs.addUserFloat("lxyz_bs_err", lxyzBsVect.error());
        bs.addUserFloat("lxyz_bs_sig", lxyzBsVect.significance());

        bs.addUserFloat("lxy_ds",     lxyDsVect.value());
        bs.addUserFloat("lxy_ds_err", lxyDsVect.error());
        bs.addUserFloat("lxy_ds_sig", lxyDsVect.significance());

        bs.addUserFloat("lxyz_ds",    lxyzDsVect.value());
        bs.addUserFloat("lxyz_ds_err", lxyzDsVect.error());
        bs.addUserFloat("lxyz_ds_sig", lxyzDsVect.significance());

        bs.addUserFloat("lxy_phi",    lxyPhiVect.value());
        bs.addUserFloat("lxy_phi_err", lxyPhiVect.error());
        bs.addUserFloat("lxy_phi_sig", lxyPhiVect.significance());

        bs.addUserFloat("lxyz_phi",   lxyzPhiVect.value());
        bs.addUserFloat("lxyz_phi_err", lxyzPhiVect.error());
        bs.addUserFloat("lxyz_phi_sig", lxyzPhiVect.significance());

        // dxy(z) is the impact parameter in the xy(z) plane(space), i.e. the distance to the PV
        // TODO: check the errors of dxy and dz      
        // TODO: bestTrack() is not refitted -> bad? #can be changed but i think its good! bc the refitted track
        // may be wrongly crorected when its coming from a tau, nbc we fit it with the ds directly to the bs
 
        float dxyMu    = muPtr->bestTrack()->dxy(pv.position());  
        float dxyMuErr = muPtr->bestTrack()->dxyError(pv.position(),pv.covariance());  
        float dxyMuSig = dxyMu/dxyMuErr;

        float dzMu     = muPtr->bestTrack()->dz(pv.position());  
        float dzMuErr  = muPtr->bestTrack()->dzError();  
        float dzMuSig  = dzMu/dzMuErr ; 

        float dxyPi    = piPtr->bestTrack()->dxy(pv.position());  //maybe useful for Ds* vs Ds ? 
        float dxyPiErr = piPtr->bestTrack()->dxyError(pv.position(),pv.covariance());  
        float dxyPiSig = dxyPi/dxyPiErr;

        float dzPi     = piPtr->bestTrack()->dz(pv.position());  
        float dzPiErr  = piPtr->bestTrack()->dzError();  
        float dzPiSig  = dzPi/dzPiErr ; 

        float dxyK1    = k1Ptr->bestTrack()->dxy(pv.position()); //needed ? 
        float dxyK1Err = k1Ptr->bestTrack()->dxyError(pv.position(),pv.covariance());
        float dxyK1Sig = dxyK1/dxyK1Err;

        float dzK1     = k1Ptr->bestTrack()->dz(pv.position());  
        float dzK1Err  = k1Ptr->bestTrack()->dzError();  
        float dzK1Sig  = dzK1/dzK1Err ; 

        float dxyK2    = k2Ptr->bestTrack()->dxy(pv.position()); //needed ? 
        float dxyK2Err = k2Ptr->bestTrack()->dxyError(pv.position(),pv.covariance());
        float dxyK2Sig = dxyK2/dxyK2Err;

        float dzK2     = k2Ptr->bestTrack()->dz(pv.position());  
        float dzK2Err  = k2Ptr->bestTrack()->dzError();  
        float dzK2Sig  = dzK2/dzK2Err ; 

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

        //////////////////////////perp Ds+Mu and /corrected Bs mass ///////////////////////////////////////////
        //accoring to https://journals.aps.org/prd/pdf/10.1103/PhysRevD.100.112006 eq 3.

        TVector3 bFlightDir;
        bFlightDir.SetXYZ(pv_x - sv_x, pv_y - sv_y, pv_z - sv_z);

        double dsPerp     = refittedDs.  Vect().Perp(bFlightDir);
        double dsMuPerp   = refittedDsMu.Vect().Perp(bFlightDir);
        double bsMassCorr = std::sqrt(std::pow(refittedDsMu.M(),2) + std::pow(dsMuPerp,2)) + dsMuPerp;

        bs.addUserFloat("bs_m_corr",  bsMassCorr);
        bs.addUserFloat("ds_perp",    dsPerp);
        bs.addUserFloat("ds_mu_perp", dsMuPerp);


        ///////////////////////// reconstruct the B momentum with different methods ///////////////////////////////


        // For the ds mass
        TLorentzVector refittedPhiPi;
        refittedPhiPi = refittedPhi + refittedPi; 

        /////////////////////////
        // Collinear approx.   //
        ///////////////////////// 

        TLorentzVector collBsTlv = collMethod(refittedDsMu, bsMass_);

        TLorentzVector collMissTlv; // for m2 miss
        TLorentzVector collQTlv;    // for q2 miss

        collMissTlv = collBsTlv - refittedDsMu; // bs - ds+mu
        collQTlv    = collBsTlv - refittedDs;   // bs - ds

        double m2_miss_coll = collMissTlv.M2();
        double pt_miss_coll = collMissTlv.Pt();
        double q2_coll      = collQTlv.M2();
        float  e_star_coll  = getEStar(collBsTlv,refittedMu); 

        bs.addUserFloat("bs_px_coll",collBsTlv.Px());
        bs.addUserFloat("bs_py_coll",collBsTlv.Py());
        bs.addUserFloat("bs_pz_coll",collBsTlv.Pz());
        bs.addUserFloat("bs_pt_coll",collBsTlv.Pt());
        bs.addUserFloat("bs_eta_coll",collBsTlv.Eta());
        bs.addUserFloat("bs_phi_coll",collBsTlv.Phi());

        bs.addUserFloat("m2_miss_coll",m2_miss_coll);
        bs.addUserFloat("pt_miss_coll",pt_miss_coll);
        bs.addUserFloat("q2_coll",q2_coll);
        bs.addUserFloat("e_star_coll",e_star_coll);

        bs.addUserFloat("b_boost_coll",collBsTlv.BoostVector().Mag());
        bs.addUserFloat("b_boost_coll_pt",collBsTlv.BoostVector().Pt());
        bs.addUserFloat("b_boost_coll_eta",collBsTlv.BoostVector().Eta());
        bs.addUserFloat("b_boost_coll_phi",collBsTlv.BoostVector().Phi());

        //////////////////////////
        // LHCb method          //
        //////////////////////////

        TLorentzVector lhcbBsTlv = lhcbMethod(refittedDsMu, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_);      
        TLorentzVector lhcbMissTlv = lhcbBsTlv - refittedDsMu;
        TLorentzVector lhcbQTlv = lhcbBsTlv - fittedDs;

        float m2_miss_lhcb = lhcbMissTlv.M2();
        float pt_miss_lhcb = lhcbMissTlv.Pt();
        float q2_lhcb      = lhcbQTlv.M2();
        float e_star_lhcb  = getEStar(lhcbBsTlv,refittedMu); 

        bs.addUserFloat("bs_px_lhcb",lhcbBsTlv.Px());
        bs.addUserFloat("bs_py_lhcb",lhcbBsTlv.Py());
        bs.addUserFloat("bs_pz_lhcb",lhcbBsTlv.Pz());
        bs.addUserFloat("bs_pt_lhcb",lhcbBsTlv.Pt());
        bs.addUserFloat("bs_eta_lhcb",lhcbBsTlv.Eta());
        bs.addUserFloat("bs_phi_lhcb",lhcbBsTlv.Phi());

        bs.addUserFloat("m2_miss_lhcb",m2_miss_lhcb);
        bs.addUserFloat("pt_miss_lhcb",pt_miss_lhcb);
        bs.addUserFloat("q2_lhcb",q2_lhcb);
        bs.addUserFloat("e_star_lhcb",e_star_lhcb);

        bs.addUserFloat("b_boost_lhcb",lhcbBsTlv.BoostVector().Mag());
        bs.addUserFloat("b_boost_lhcb_pt",lhcbBsTlv.BoostVector().Pt());
        bs.addUserFloat("b_boost_lhcb_eta",lhcbBsTlv.BoostVector().Eta());
        bs.addUserFloat("b_boost_lhcb_phi",lhcbBsTlv.BoostVector().Phi());


        ///////////////////////////////////////////
        // Another Lhcb method                   //
        ///////////////////////////////////////////

        TLorentzVector lhcbAltBsTlv = lhcbAltMethod(refittedDsMu, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_);      

        TLorentzVector lhcbAltMissTlv = lhcbAltBsTlv - refittedDsMu; // bs - ds+mu
        TLorentzVector lhcbAltQTlv = lhcbAltBsTlv - fittedDs; // bs - ds 

        float e_star_lhcb_alt  = getEStar(lhcbAltBsTlv,refittedMu); 

        bs.addUserFloat("bs_px_lhcb_alt",lhcbAltBsTlv.Px());
        bs.addUserFloat("bs_py_lhcb_alt",lhcbAltBsTlv.Py());
        bs.addUserFloat("bs_pz_lhcb_alt",lhcbAltBsTlv.Pz());
        bs.addUserFloat("bs_pt_lhcb_alt",lhcbAltBsTlv.Pt());
        bs.addUserFloat("bs_eta_lhcb_alt",lhcbAltBsTlv.Eta());
        bs.addUserFloat("bs_phi_lhcb_alt",lhcbAltBsTlv.Phi());

        bs.addUserFloat("m2_miss_lhcb_alt",lhcbAltMissTlv.M2());
        bs.addUserFloat("pt_miss_lhcb_alt",lhcbAltMissTlv.Pt());
        bs.addUserFloat("q2_lhcb_alt",lhcbAltQTlv.M2());
        bs.addUserFloat("e_star_lhcb_alt",e_star_lhcb_alt);

        bs.addUserFloat("b_boost_lhcb_alt",lhcbAltBsTlv.BoostVector().Mag());
        bs.addUserFloat("b_boost_lhcb_alt_pt",lhcbAltBsTlv.BoostVector().Pt());
        bs.addUserFloat("b_boost_lhcb_alt_eta",lhcbAltBsTlv.BoostVector().Eta());
        bs.addUserFloat("b_boost_lhcb_alt_phi",lhcbAltBsTlv.BoostVector().Phi());

        ///////////////////////////////////////////////////////////
        // RECO METHOD -> 2 Solutions                            //
        // This method is exact in the 1-neutrino final state :) //
        ///////////////////////////////////////////////////////////
        //std::cout << "now i calc the reco method:" << std::endl;
        std::tuple<std::vector<TLorentzVector>,float> recoResult = recoMethod(refittedDsMu, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_); 

        //std::tuple<std::vector<TLorentzVector>,bool> recoResult = recoMethod(fittedDs + muTlv, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_);

        std::vector<TLorentzVector> recosBs = std::get<0>(recoResult);
        float discNegativity               = std::get<1>(recoResult);

        int discIsNegative = 0;
        if (discNegativity > 0) discIsNegative = 1;

        TLorentzVector recoBsTlv1 = recosBs.at(0);
        TLorentzVector recoBsTlv2 = recosBs.at(1);

        TLorentzVector miss_1 = recoBsTlv1 - refittedDsMu;  
        TLorentzVector miss_2 = recoBsTlv2 - refittedDsMu;  
        TLorentzVector Q_1 = recoBsTlv1 - fittedDs;  
        TLorentzVector Q_2 = recoBsTlv2 - fittedDs;  

        float m2_miss_reco_1 = miss_1.M2();
        float m2_miss_reco_2 = miss_2.M2();
        float pt_miss_reco_1 = miss_1.Pt();
        float pt_miss_reco_2 = miss_2.Pt();

        float q2_reco_1 = Q_1.M2();
        float q2_reco_2 = Q_2.M2();
        float e_star_reco_1  = getEStar(recoBsTlv1,refittedMu); 
        float e_star_reco_2  = getEStar(recoBsTlv2,refittedMu); 


        bs.addUserInt("disc_is_negative",discIsNegative);
        bs.addUserFloat("disc_negativity",discNegativity);

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
        bs.addUserFloat("pt_miss_reco_1",pt_miss_reco_1); 
        bs.addUserFloat("pt_miss_reco_2",pt_miss_reco_2); 


        bs.addUserFloat("q2_reco_1",q2_reco_1); 
        bs.addUserFloat("q2_reco_2",q2_reco_2); 
        bs.addUserFloat("e_star_reco_1",e_star_reco_1);
        bs.addUserFloat("e_star_reco_2",e_star_reco_2);

        bs.addUserFloat("b_boost_reco_1",recoBsTlv1.BoostVector().Mag());
        bs.addUserFloat("b_boost_reco_1_pt",recoBsTlv1.BoostVector().Pt());
        bs.addUserFloat("b_boost_reco_1_eta",recoBsTlv1.BoostVector().Eta());
        bs.addUserFloat("b_boost_reco_1_phi",recoBsTlv1.BoostVector().Phi());


        bs.addUserFloat("b_boost_reco_2",recoBsTlv2.BoostVector().Mag());
        bs.addUserFloat("b_boost_reco_2_pt",recoBsTlv2.BoostVector().Pt());
        bs.addUserFloat("b_boost_reco_2_eta",recoBsTlv2.BoostVector().Eta());
        bs.addUserFloat("b_boost_reco_2_phi",recoBsTlv2.BoostVector().Phi());


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
        bs.addUserFloat("ds_fitted_boost",  fittedDs.BoostVector().Mag()); 

        bs.addUserFloat("phiPi_refitted_px",     refittedPhiPi.Px()); 
        bs.addUserFloat("phiPi_refitted_py",     refittedPhiPi.Py()); 
        bs.addUserFloat("phiPi_refitted_pz",     refittedPhiPi.Pz()); 
        bs.addUserFloat("phiPi_refitted_pt",     refittedPhiPi.Pt()); 
        bs.addUserFloat("phiPi_refitted_eta",    refittedPhiPi.Eta()); 
        bs.addUserFloat("phiPi_refitted_phi",    refittedPhiPi.Phi()); 
        bs.addUserFloat("phiPi_refitted_m",      refittedPhiPi.M()); 
        bs.addUserFloat("phiPi_refitted_boost",  refittedPhiPi.BoostVector().Mag()); 

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
        bs.addUserFloat("ds_refitted_boost",refittedDs.BoostVector().Mag()); 

        bs.addUserFloat("bs_fitted_px",     fittedBs.Px()); 
        bs.addUserFloat("bs_fitted_py",     fittedBs.Py()); 
        bs.addUserFloat("bs_fitted_pz",     fittedBs.Pz()); 
        bs.addUserFloat("bs_fitted_pt",     fittedBs.Pt()); 
        bs.addUserFloat("bs_fitted_eta",    fittedBs.Eta()); 
        bs.addUserFloat("bs_fitted_phi",    fittedBs.Phi()); 
        bs.addUserFloat("bs_fitted_m",      fittedBs.M()); 

        /////////////////////// helicity angles in all possibe variations /////////////////////////

        //angle between Mu and W
        float angMuWColl    = angMuW(refittedDs,collBsTlv,   refittedMu);
        float angMuWLhcb    = angMuW(refittedDs,lhcbBsTlv,   refittedMu);
        float angMuWLhcbAlt = angMuW(refittedDs,lhcbAltBsTlv,refittedMu);
        float angMuWReco1   = angMuW(refittedDs,recoBsTlv1,  refittedMu);
        float angMuWReco2   = angMuW(refittedDs,recoBsTlv2,  refittedMu);

        bs.addUserFloat("cosMuWColl",       cos(angMuWColl)); 
        bs.addUserFloat("cosMuWLhcb",       cos(angMuWLhcb)); 
        bs.addUserFloat("cosMuWLhcbAlt",    cos(angMuWLhcbAlt)); 
        bs.addUserFloat("cosMuWReco1",      cos(angMuWReco1)); 
        bs.addUserFloat("cosMuWReco2",      cos(angMuWReco2)); 


        //angle between k1(k2) and pion in phi rest frame
        float angPiK1  = angDoubleDecay(refittedPhi,refittedK1, refittedPi);
        float angPiK2  = angDoubleDecay(refittedPhi,refittedK2, refittedPi);
        // equivalently, angle between phi(pi) and bs in ds rest frame
        float angPhiDsColl    = angDsPi(refittedDs, refittedPhi,  collBsTlv);
        float angPhiDsLhcb    = angDsPi(refittedDs, refittedPhi,  lhcbBsTlv);
        float angPhiDsLhcbAlt = angDsPi(refittedDs, refittedPhi,  lhcbAltBsTlv);
        float angPhiDsReco1   = angDsPi(refittedDs, refittedPhi,  recoBsTlv1);
        float angPhiDsReco2   = angDsPi(refittedDs, refittedPhi,  recoBsTlv2);

        float angPiDsColl    = angDsPi(refittedDs, refittedPi,  collBsTlv);
        float angPiDsLhcb    = angDsPi(refittedDs, refittedPi,  lhcbBsTlv);
        float angPiDsLhcbAlt = angDsPi(refittedDs, refittedPi,  lhcbAltBsTlv);
        float angPiDsReco1   = angDsPi(refittedDs, refittedPi,  recoBsTlv1);
        float angPiDsReco2   = angDsPi(refittedDs, refittedPi,  recoBsTlv2);
        

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
        float angPlaneBsColl    = angPlane(refittedDs, collBsTlv,    refittedMu, refittedPi);
        float angPlaneBsLhcb    = angPlane(refittedDs, lhcbBsTlv,    refittedMu, refittedPi);
        float angPlaneBsLhcbAlt = angPlane(refittedDs, lhcbAltBsTlv, refittedMu, refittedPi);
        float angPlaneBsReco1   = angPlane(refittedDs, recoBsTlv1,   refittedMu, refittedPi);
        float angPlaneBsReco2   = angPlane(refittedDs, recoBsTlv2,   refittedMu, refittedPi);

        bs.addUserFloat("cosPlaneBsColl",    cos(angPlaneBsColl));
        bs.addUserFloat("cosPlaneBsLhcb",    cos(angPlaneBsLhcb));
        bs.addUserFloat("cosPlaneBsLhcbAlt", cos(angPlaneBsLhcbAlt));
        bs.addUserFloat("cosPlaneBsReco1",   cos(angPlaneBsReco1));
        bs.addUserFloat("cosPlaneBsReco2",   cos(angPlaneBsReco2));

        float angPlaneDsColl    = angPlane2(refittedDs, collBsTlv,    refittedK1, refittedPi);
        float angPlaneDsLhcb    = angPlane2(refittedDs, lhcbBsTlv,    refittedK1, refittedPi);
        float angPlaneDsLhcbAlt = angPlane2(refittedDs, lhcbAltBsTlv, refittedK1, refittedPi);
        float angPlaneDsReco1   = angPlane2(refittedDs, recoBsTlv1,   refittedK1, refittedPi);
        float angPlaneDsReco2   = angPlane2(refittedDs, recoBsTlv2,   refittedK1, refittedPi);

        bs.addUserFloat("cosPlaneDsColl",    cos(angPlaneDsColl));
        bs.addUserFloat("cosPlaneDsLhcb",    cos(angPlaneDsLhcb));
        bs.addUserFloat("cosPlaneDsLhcbAlt", cos(angPlaneDsLhcbAlt));
        bs.addUserFloat("cosPlaneDsReco1",   cos(angPlaneDsReco1));
        bs.addUserFloat("cosPlaneDsReco2",   cos(angPlaneDsReco2));


        // muon isolation (simply looping and summing track pt's is too simple!)
        // idea: charged part: take all charged hadron pts and subtract kkpi cand
        //       neutral part: take neutral hadronand photon ET and subtract contributions
        //                     from PU vertices (assumed to be about ~50%, see paper)

        auto muIso03 = muPtr->pfIsolationR03(); //dR = 0.3
        auto muIso04 = muPtr->pfIsolationR04(); //dR = 0.4

        
        float iso03KKPi = 0.0;
        float iso04KKPi = 0.0;

        // subtract kkpi candidate pt if in cone   

        if( (reco::deltaR(*k1Ptr,*muPtr) < 0.3) && ( abs(muPtr->bestTrack()->dz(pv.position()) - k1Ptr->dz(pv.position()) ) < 0.2 ) ) iso03KKPi += k1Ptr->pt();
        if( (reco::deltaR(*k2Ptr,*muPtr) < 0.3) && ( abs(muPtr->bestTrack()->dz(pv.position()) - k2Ptr->dz(pv.position()) ) < 0.2 ) ) iso03KKPi += k2Ptr->pt();
        if( (reco::deltaR(*piPtr,*muPtr) < 0.3) && ( abs(muPtr->bestTrack()->dz(pv.position()) - piPtr->dz(pv.position()) ) < 0.2 ) ) iso03KKPi += piPtr->pt();

        if( (reco::deltaR(*k1Ptr,*muPtr) < 0.4) && ( abs(muPtr->bestTrack()->dz(pv.position()) - k1Ptr->dz(pv.position()) ) < 0.2 ) ) iso04KKPi += k1Ptr->pt();
        if( (reco::deltaR(*k2Ptr,*muPtr) < 0.4) && ( abs(muPtr->bestTrack()->dz(pv.position()) - k2Ptr->dz(pv.position()) ) < 0.2 ) ) iso04KKPi += k2Ptr->pt();
        if( (reco::deltaR(*piPtr,*muPtr) < 0.4) && ( abs(muPtr->bestTrack()->dz(pv.position()) - piPtr->dz(pv.position()) ) < 0.2 ) ) iso04KKPi += piPtr->pt();


        // factor 0.5 from https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=MUO-16-001&tp=an&id=1779&ancode=MUO-16-001

        //float zero = 0.0; // avoid problems with std::max function

        float iso03charged    = std::max(0.0f, muIso03.sumChargedHadronPt - iso03KKPi);  
        float iso04charged    = std::max(0.0f, muIso04.sumChargedHadronPt - iso04KKPi);  

        float iso03neutral    = std::max(0.0, muIso03.sumNeutralHadronEt + muIso03.sumPhotonEt - 0.5 * muIso03.sumPUPt);  
        float iso04neutral    = std::max(0.0, muIso04.sumNeutralHadronEt + muIso04.sumPhotonEt - 0.5 * muIso04.sumPUPt);  

        //std::cout << "neutral part:" << iso03charged<< std::endl;
        //std::cout << "neutral part:" << iso03neutral<< std::endl;

        float iso03 = iso03charged + iso03neutral;
        float iso04 = iso04charged + iso04neutral;

        float relIso03 = iso03 / muPtr->pt();
        float relIso04 = iso04 / muPtr->pt();
        float refittedRelIso03 = iso03 / refittedMu.Pt();
        float refittedRelIso04 = iso04 / refittedMu.Pt();

        bs.addUserFloat("mu_iso_03", iso03);
        bs.addUserFloat("mu_iso_04", iso04);
        bs.addUserFloat("mu_rel_iso_03", relIso03);
        bs.addUserFloat("mu_rel_iso_04", relIso04);
        bs.addUserFloat("mu_rel_iso_03_refitted", refittedRelIso03);
        bs.addUserFloat("mu_rel_iso_04_refitted", refittedRelIso04);
        //std::cout << "mu iso 03: " << relIso03 << std::endl;

        // e gamma
        float e_gamma = getEGamma(refittedDs, dsMass_,dsStarMass_);
        bs.addUserFloat("e_gamma", e_gamma);

        // add kappa
        float kappa = getKappa(refittedDs, muTlv);
        bs.addUserFloat("kappa", kappa);



        ////////////////////////////////////////////////
        // Redo all the methods if we have a photon!  //
        ////////////////////////////////////////////////
        double dsPerpPhoton                 = -9999;
        double dsMuPerpPhoton               = -9999;
        double bsMassCorrPhoton             = -9999;

        double m2_miss_coll_photon          = -9999;
        double pt_miss_coll_photon          = -9999;
        double q2_coll_photon               = -9999;
        float  e_star_coll_photon           = -9999;

        double bs_px_coll_photon            = -9999;
        double bs_py_coll_photon            = -9999;
        double bs_pz_coll_photon            = -9999;
        double bs_pt_coll_photon            = -9999;
        double bs_eta_coll_photon           = -9999;
        double bs_phi_coll_photon           = -9999;

        double b_boost_coll_photon         = -9999;
        double b_boost_coll_pt_photon      = -9999;
        double b_boost_coll_eta_photon     = -9999;
        double b_boost_coll_phi_photon     = -9999;

        double m2_miss_lhcb_alt_photon      = -9999;
        double pt_miss_lhcb_alt_photon      = -9999; 
        double q2_lhcb_alt_photon           = -9999; 
        double e_star_lhcb_alt_photon       = -9999;  

        double bs_px_lhcb_alt_photon        = -9999;
        double bs_py_lhcb_alt_photon        = -9999;
        double bs_pz_lhcb_alt_photon        = -9999;
        double bs_pt_lhcb_alt_photon        = -9999;
        double bs_eta_lhcb_alt_photon       = -9999;
        double bs_phi_lhcb_alt_photon       = -9999;

        double b_boost_lhcb_alt_photon     = -9999;
        double b_boost_lhcb_alt_pt_photon  = -9999;
        double b_boost_lhcb_alt_eta_photon = -9999;
        double b_boost_lhcb_alt_phi_photon = -9999;

        double m2_miss_reco_1_photon        = -9999;
        double pt_miss_reco_1_photon        = -9999; 
        double q2_reco_1_photon             = -9999; 
        double e_star_reco_1_photon         = -9999;  

        double bs_px_reco_1_photon          = -9999;
        double bs_py_reco_1_photon          = -9999;
        double bs_pz_reco_1_photon          = -9999;
        double bs_pt_reco_1_photon          = -9999;
        double bs_eta_reco_1_photon         = -9999;
        double bs_phi_reco_1_photon         = -9999;

        double b_boost_reco_1_photon       = -9999;
        double b_boost_reco_1_pt_photon    = -9999;
        double b_boost_reco_1_eta_photon   = -9999;
        double b_boost_reco_1_phi_photon   = -9999;

        double m2_miss_reco_2_photon        = -9999;
        double pt_miss_reco_2_photon        = -9999; 
        double q2_reco_2_photon             = -9999; 
        double e_star_reco_2_photon         = -9999;  

        double bs_px_reco_2_photon          = -9999;
        double bs_py_reco_2_photon          = -9999;
        double bs_pz_reco_2_photon          = -9999;
        double bs_pt_reco_2_photon          = -9999;
        double bs_eta_reco_2_photon         = -9999;
        double bs_phi_reco_2_photon         = -9999;

        double b_boost_reco_2_photon       = -9999;
        double b_boost_reco_2_pt_photon    = -9999;
        double b_boost_reco_2_eta_photon   = -9999;
        double b_boost_reco_2_phi_photon   = -9999;

        double discNegativityPhoton         = -9999;
        double discIsNegativePhoton         = -9999;

        //angle between Mu and W
        float cosMuWCollPhoton             = -9999;
        float cosMuWLhcbAltPhoton          = -9999;
        float cosMuWReco1Photon            = -9999;
        float cosMuWReco2Photon            = -9999;

        float cosPhiDsCollPhoton           = -9999; 
        float cosPhiDsLhcbAltPhoton        = -9999; 
        float cosPhiDsReco1Photon          = -9999; 
        float cosPhiDsReco2Photon          = -9999; 

        float cosPiDsCollPhoton            = -9999; 
        float cosPiDsLhcbAltPhoton         = -9999; 
        float cosPiDsReco1Photon           = -9999; 
        float cosPiDsReco2Photon           = -9999; 

        // helicity plane angle
        float cosPlaneBsCollPhoton         = -9999;
        float cosPlaneBsLhcbAltPhoton      = -9999;
        float cosPlaneBsReco1Photon        = -9999;
        float cosPlaneBsReco2Photon        = -9999;

        float cosPlaneDsCollPhoton         = -9999;
        float cosPlaneDsLhcbAltPhoton      = -9999;
        float cosPlaneDsReco1Photon        = -9999;
        float cosPlaneDsReco2Photon        = -9999;



        if (foundPhoton){

          TLorentzVector refittedDsMuPhoton = refittedDsMu + photonTlv;
          TLorentzVector refittedDsPhoton   = refittedDs   + photonTlv;

          /////////////////////////////
          // Ds + mu perp variable   //
          /////////////////////////////
         
          double dsPerpPhoton   = refittedDsPhoton.  Vect().Perp(bFlightDir);
          double dsMuPerpPhoton = refittedDsMuPhoton.Vect().Perp(bFlightDir);
  
          bsMassCorrPhoton = std::sqrt(std::pow(refittedDsMuPhoton.M(),2) + std::pow(dsMuPerpPhoton,2)) + dsMuPerpPhoton;
  
          bs.addUserFloat("bs_m_corr_photon",  bsMassCorrPhoton);
          bs.addUserFloat("ds_perp_photon",    dsPerpPhoton);
          bs.addUserFloat("ds_mu_perp_photon", dsMuPerpPhoton);

          /////////////////////////////
          // Collinear approximation //
          /////////////////////////////

  
          TLorentzVector collBsTlvPhoton = collMethod(refittedDsMuPhoton, bsMass_);
  
          TLorentzVector collMissTlvPhoton; // for m2 miss
          TLorentzVector collQTlvPhoton;    // for q2 miss
  
          collMissTlvPhoton = collBsTlvPhoton - refittedDsMuPhoton; // bs - ds+mu
          collQTlvPhoton    = collBsTlvPhoton - refittedDsPhoton;   // bs - ds
  
          m2_miss_coll_photon         = collMissTlvPhoton.M2();
          pt_miss_coll_photon         = collMissTlvPhoton.Pt();
          q2_coll_photon              = collQTlvPhoton.M2();
          e_star_coll_photon          = getEStar(collBsTlvPhoton,refittedMu); 

          bs_px_coll_photon           = collBsTlvPhoton.Px(); 
          bs_py_coll_photon           = collBsTlvPhoton.Py(); 
          bs_pz_coll_photon           = collBsTlvPhoton.Pz(); 
          bs_pt_coll_photon           = collBsTlvPhoton.Pt(); 
          bs_eta_coll_photon          = collBsTlvPhoton.Eta(); 
          bs_phi_coll_photon          = collBsTlvPhoton.Phi(); 

          b_boost_coll_photon        = collBsTlvPhoton.BoostVector().Mag(); 
          b_boost_coll_pt_photon     = collBsTlvPhoton.BoostVector().Pt(); 
          b_boost_coll_eta_photon    = collBsTlvPhoton.BoostVector().Eta(); 
          b_boost_coll_phi_photon    = collBsTlvPhoton.BoostVector().Phi(); 
 
          /////////////////////////////
          // Lhcb method             //
          /////////////////////////////

          TLorentzVector lhcbAltBsTlvPhoton   = lhcbAltMethod(refittedDsMuPhoton, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_);      
  
          TLorentzVector lhcbAltMissTlvPhoton = lhcbAltBsTlv - refittedDsMuPhoton; // bs - ds+mu
          TLorentzVector lhcbAltQTlvPhoton    = lhcbAltBsTlv - refittedDsPhoton;   // bs - ds 
  
          m2_miss_lhcb_alt_photon      = lhcbAltMissTlvPhoton.M2();
          pt_miss_lhcb_alt_photon      = lhcbAltMissTlvPhoton.Pt();
          q2_lhcb_alt_photon           = lhcbAltQTlvPhoton.M2();
          e_star_lhcb_alt_photon       = getEStar(lhcbAltBsTlvPhoton,refittedMu); 

          bs_px_lhcb_alt_photon        = lhcbAltBsTlvPhoton.Px(); 
          bs_py_lhcb_alt_photon        = lhcbAltBsTlvPhoton.Py(); 
          bs_pz_lhcb_alt_photon        = lhcbAltBsTlvPhoton.Pz(); 
          bs_pt_lhcb_alt_photon        = lhcbAltBsTlvPhoton.Pt(); 
          bs_eta_lhcb_alt_photon       = lhcbAltBsTlvPhoton.Eta(); 
          bs_phi_lhcb_alt_photon       = lhcbAltBsTlvPhoton.Phi(); 

          b_boost_lhcb_alt_photon     = lhcbAltBsTlvPhoton.BoostVector().Mag(); 
          b_boost_lhcb_alt_pt_photon  = lhcbAltBsTlvPhoton.BoostVector().Pt(); 
          b_boost_lhcb_alt_eta_photon = lhcbAltBsTlvPhoton.BoostVector().Eta(); 
          b_boost_lhcb_alt_phi_photon = lhcbAltBsTlvPhoton.BoostVector().Phi(); 

          /////////////////////////////
          // Reco method             //
          /////////////////////////////

          std::tuple<std::vector<TLorentzVector>,float> recoResultPhoton = recoMethod(refittedDsMuPhoton, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_); 
  
          std::vector<TLorentzVector> recosBsPhoton = std::get<0>(recoResultPhoton);
          discNegativityPhoton                      = std::get<1>(recoResultPhoton);
  
          discIsNegativePhoton = 0;
          if (discNegativityPhoton > 0) discIsNegativePhoton = 1;
  
          TLorentzVector recoBsTlv1Photon = recosBsPhoton.at(0);
          TLorentzVector recoBsTlv2Photon = recosBsPhoton.at(1);
  
          TLorentzVector miss_1Photon     = recoBsTlv1Photon - refittedDsMuPhoton;  
          TLorentzVector miss_2Photon     = recoBsTlv2Photon - refittedDsMuPhoton;  
          TLorentzVector Q_1Photon        = recoBsTlv1Photon - refittedDsPhoton;  
          TLorentzVector Q_2Photon        = recoBsTlv2Photon - refittedDsPhoton;  
  
          m2_miss_reco_1_photon           = miss_1Photon.M2();
          pt_miss_reco_1_photon           = miss_1Photon.Pt();
          q2_reco_1_photon                = Q_1Photon.M2();
          e_star_reco_1_photon            = getEStar(recoBsTlv1Photon,refittedMu); 
 
          bs_px_reco_1_photon             = recoBsTlv1Photon.Px(); 
          bs_py_reco_1_photon             = recoBsTlv1Photon.Py(); 
          bs_pz_reco_1_photon             = recoBsTlv1Photon.Pz(); 
          bs_pt_reco_1_photon             = recoBsTlv1Photon.Pt(); 
          bs_eta_reco_1_photon            = recoBsTlv1Photon.Eta(); 
          bs_phi_reco_1_photon            = recoBsTlv1Photon.Phi(); 

          b_boost_reco_1_photon          = recoBsTlv1Photon.BoostVector().Mag(); 
          b_boost_reco_1_pt_photon       = recoBsTlv1Photon.BoostVector().Pt(); 
          b_boost_reco_1_eta_photon      = recoBsTlv1Photon.BoostVector().Eta(); 
          b_boost_reco_1_phi_photon      = recoBsTlv1Photon.BoostVector().Phi(); 

          m2_miss_reco_2_photon           = miss_2Photon.M2();
          pt_miss_reco_2_photon           = miss_2Photon.Pt();
          q2_reco_2_photon                = Q_2Photon.M2();
          e_star_reco_2_photon            = getEStar(recoBsTlv2Photon,refittedMu); 

          bs_px_reco_2_photon             = recoBsTlv2Photon.Px(); 
          bs_py_reco_2_photon             = recoBsTlv2Photon.Py(); 
          bs_pz_reco_2_photon             = recoBsTlv2Photon.Pz(); 
          bs_pt_reco_2_photon             = recoBsTlv2Photon.Pt(); 
          bs_eta_reco_2_photon            = recoBsTlv2Photon.Eta(); 
          bs_phi_reco_2_photon            = recoBsTlv2Photon.Phi(); 

          b_boost_reco_2_photon          = recoBsTlv2Photon.BoostVector().Mag(); 
          b_boost_reco_2_pt_photon       = recoBsTlv2Photon.BoostVector().Pt(); 
          b_boost_reco_2_eta_photon      = recoBsTlv2Photon.BoostVector().Eta(); 
          b_boost_reco_2_phi_photon      = recoBsTlv2Photon.BoostVector().Phi(); 

          /////////////////////////////
          // Hel angles              // 
          /////////////////////////////

          double angMuWCollPhoton         = angMuW(refittedDsPhoton, collBsTlvPhoton,   refittedMu);
          double angMuWLhcbAltPhoton      = angMuW(refittedDsPhoton, lhcbAltBsTlvPhoton,refittedMu);
          double angMuWReco1Photon        = angMuW(refittedDsPhoton, recoBsTlv1Photon,  refittedMu);
          double angMuWReco2Photon        = angMuW(refittedDsPhoton, recoBsTlv2Photon,  refittedMu);

          // equivalently, angle between phi(pi) and bs in ds rest frame
          double angPhiDsCollPhoton       = angDsPi(refittedDsPhoton, refittedPhi,  collBsTlvPhoton);
          double angPhiDsLhcbAltPhoton    = angDsPi(refittedDsPhoton, refittedPhi,  lhcbAltBsTlvPhoton);
          double angPhiDsReco1Photon      = angDsPi(refittedDsPhoton, refittedPhi,  recoBsTlv1Photon);
          double angPhiDsReco2Photon      = angDsPi(refittedDsPhoton, refittedPhi,  recoBsTlv2Photon);
  
          double angPiDsCollPhoton        = angDsPi(refittedDsPhoton, refittedPi,  collBsTlvPhoton);
          double angPiDsLhcbAltPhoton     = angDsPi(refittedDsPhoton, refittedPi,  lhcbAltBsTlvPhoton);
          double angPiDsReco1Photon       = angDsPi(refittedDsPhoton, refittedPi,  recoBsTlv1Photon);
          double angPiDsReco2Photon       = angDsPi(refittedDsPhoton, refittedPi,  recoBsTlv2Photon);


          cosMuWCollPhoton         =  cos(angMuWCollPhoton        ); 
          cosMuWLhcbAltPhoton      =  cos(angMuWLhcbAltPhoton     ); 
          cosMuWReco1Photon        =  cos(angMuWReco1Photon       ); 
          cosMuWReco2Photon        =  cos(angMuWReco2Photon       ); 

          cosPhiDsCollPhoton       =  cos(angPhiDsCollPhoton      ); 
          cosPhiDsLhcbAltPhoton    =  cos(angPhiDsLhcbAltPhoton   ); 
          cosPhiDsReco1Photon      =  cos(angPhiDsReco1Photon     ); 
          cosPhiDsReco2Photon      =  cos(angPhiDsReco2Photon     ); 

          cosPiDsCollPhoton        =  cos(angPiDsCollPhoton       ); 
          cosPiDsLhcbAltPhoton     =  cos(angPiDsLhcbAltPhoton    ); 
          cosPiDsReco1Photon       =  cos(angPiDsReco1Photon      ); 
          cosPiDsReco2Photon       =  cos(angPiDsReco2Photon      ); 
          
          // helicity plane angle
          float angPlaneBsCollPhoton    = angPlane(refittedDsPhoton, collBsTlvPhoton,    refittedMu, refittedPi);
          float angPlaneBsLhcbAltPhoton = angPlane(refittedDsPhoton, lhcbAltBsTlvPhoton, refittedMu, refittedPi);
          float angPlaneBsReco1Photon   = angPlane(refittedDsPhoton, recoBsTlv1Photon,   refittedMu, refittedPi);
          float angPlaneBsReco2Photon   = angPlane(refittedDsPhoton, recoBsTlv2Photon,   refittedMu, refittedPi);
  
          float angPlaneDsCollPhoton    = angPlane2(refittedDsPhoton, collBsTlvPhoton,    refittedK1, refittedPi);
          float angPlaneDsLhcbAltPhoton = angPlane2(refittedDsPhoton, lhcbAltBsTlvPhoton, refittedK1, refittedPi);
          float angPlaneDsReco1Photon   = angPlane2(refittedDsPhoton, recoBsTlv1Photon,   refittedK1, refittedPi);
          float angPlaneDsReco2Photon   = angPlane2(refittedDsPhoton, recoBsTlv2Photon,   refittedK1, refittedPi);
 
          cosPlaneBsCollPhoton    =      cos( angPlaneBsCollPhoton   );  
          cosPlaneBsLhcbAltPhoton =      cos( angPlaneBsLhcbAltPhoton);
          cosPlaneBsReco1Photon   =      cos( angPlaneBsReco1Photon  );
          cosPlaneBsReco2Photon   =      cos( angPlaneBsReco2Photon  );

          cosPlaneDsCollPhoton    =      cos( angPlaneDsCollPhoton   );
          cosPlaneDsLhcbAltPhoton =      cos( angPlaneDsLhcbAltPhoton);
          cosPlaneDsReco1Photon   =      cos( angPlaneDsReco1Photon  );
          cosPlaneDsReco2Photon   =      cos( angPlaneDsReco2Photon  );
          
           

        }

        bs.addUserFloat("bs_m_corr_photon",  bsMassCorrPhoton);
        bs.addUserFloat("ds_perp_photon",    dsPerpPhoton);
        bs.addUserFloat("ds_mu_perp_photon", dsMuPerpPhoton);

        bs.addUserFloat("bs_px_coll_photon",        bs_px_coll_photon);
        bs.addUserFloat("bs_py_coll_photon",        bs_py_coll_photon);
        bs.addUserFloat("bs_pz_coll_photon",        bs_pz_coll_photon);
        bs.addUserFloat("bs_pt_coll_photon",        bs_pt_coll_photon);
        bs.addUserFloat("bs_eta_coll_photon",       bs_eta_coll_photon);
        bs.addUserFloat("bs_phi_coll_photon",       bs_phi_coll_photon);

        bs.addUserFloat("m2_miss_coll_photon",      m2_miss_coll_photon);
        bs.addUserFloat("pt_miss_coll_photon",      pt_miss_coll_photon);
        bs.addUserFloat("q2_coll_photon",           q2_coll_photon);
        bs.addUserFloat("e_star_coll_photon",       e_star_coll_photon);

        bs.addUserFloat("b_boost_coll_photon",      b_boost_coll_photon);
        bs.addUserFloat("b_boost_coll_pt_photon",   b_boost_coll_pt_photon);
        bs.addUserFloat("b_boost_coll_eta_photon",  b_boost_coll_eta_photon);
        bs.addUserFloat("b_boost_coll_phi_photon",  b_boost_coll_phi_photon);

        bs.addUserFloat("bs_px_lhcb_alt_photon",        bs_px_lhcb_alt_photon);
        bs.addUserFloat("bs_py_lhcb_alt_photon",        bs_py_lhcb_alt_photon);
        bs.addUserFloat("bs_pz_lhcb_alt_photon",        bs_pz_lhcb_alt_photon);
        bs.addUserFloat("bs_pt_lhcb_alt_photon",        bs_pt_lhcb_alt_photon);
        bs.addUserFloat("bs_eta_lhcb_alt_photon",       bs_eta_lhcb_alt_photon);
        bs.addUserFloat("bs_phi_lhcb_alt_photon",       bs_phi_lhcb_alt_photon);

        bs.addUserFloat("m2_miss_lhcb_alt_photon",      m2_miss_lhcb_alt_photon);
        bs.addUserFloat("pt_miss_lhcb_alt_photon",      pt_miss_lhcb_alt_photon);
        bs.addUserFloat("q2_lhcb_alt_photon",           q2_lhcb_alt_photon);
        bs.addUserFloat("e_star_lhcb_alt_photon",       e_star_lhcb_alt_photon);

        bs.addUserFloat("b_boost_lhcb_alt_photon",      b_boost_lhcb_alt_photon);
        bs.addUserFloat("b_boost_lhcb_alt_pt_photon",   b_boost_lhcb_alt_pt_photon);
        bs.addUserFloat("b_boost_lhcb_alt_eta_photon",  b_boost_lhcb_alt_eta_photon);
        bs.addUserFloat("b_boost_lhcb_alt_phi_photon",  b_boost_lhcb_alt_phi_photon);

        bs.addUserFloat("bs_px_reco_1_photon",      bs_px_reco_1_photon);
        bs.addUserFloat("bs_py_reco_1_photon",      bs_py_reco_1_photon);
        bs.addUserFloat("bs_pz_reco_1_photon",      bs_pz_reco_1_photon);
        bs.addUserFloat("bs_pt_reco_1_photon",      bs_pt_reco_1_photon);
        bs.addUserFloat("bs_eta_reco_1_photon",     bs_eta_reco_1_photon);
        bs.addUserFloat("bs_phi_reco_1_photon",     bs_phi_reco_1_photon);

        bs.addUserFloat("m2_miss_reco_1_photon",    m2_miss_reco_1_photon);
        bs.addUserFloat("pt_miss_reco_1_photon",    pt_miss_reco_1_photon);
        bs.addUserFloat("q2_reco_1_photon",         q2_reco_1_photon);
        bs.addUserFloat("e_star_reco_1_photon",     e_star_reco_1_photon);

        bs.addUserFloat("b_boost_reco_1_photon",     b_boost_reco_1_photon);
        bs.addUserFloat("b_boost_reco_1_pt_photon",  b_boost_reco_1_pt_photon);
        bs.addUserFloat("b_boost_reco_1_eta_photon", b_boost_reco_1_eta_photon);
        bs.addUserFloat("b_boost_reco_1_phi_photon", b_boost_reco_1_phi_photon);

        bs.addUserFloat("bs_px_reco_2_photon",      bs_px_reco_2_photon);
        bs.addUserFloat("bs_py_reco_2_photon",      bs_py_reco_2_photon);
        bs.addUserFloat("bs_pz_reco_2_photon",      bs_pz_reco_2_photon);
        bs.addUserFloat("bs_pt_reco_2_photon",      bs_pt_reco_2_photon);
        bs.addUserFloat("bs_eta_reco_2_photon",     bs_eta_reco_2_photon);
        bs.addUserFloat("bs_phi_reco_2_photon",     bs_phi_reco_2_photon);

        bs.addUserFloat("m2_miss_reco_2_photon",    m2_miss_reco_2_photon);
        bs.addUserFloat("pt_miss_reco_2_photon",    pt_miss_reco_2_photon);
        bs.addUserFloat("q2_reco_2_photon",         q2_reco_2_photon);
        bs.addUserFloat("e_star_reco_2_photon",     e_star_reco_2_photon);

        bs.addUserFloat("b_boost_reco_2_photon",     b_boost_reco_2_photon);
        bs.addUserFloat("b_boost_reco_2_pt_photon",  b_boost_reco_2_pt_photon);
        bs.addUserFloat("b_boost_reco_2_eta_photon", b_boost_reco_2_eta_photon);
        bs.addUserFloat("b_boost_reco_2_phi_photon", b_boost_reco_2_phi_photon);

        bs.addUserInt("disc_is_negative_photon",discIsNegativePhoton);
        bs.addUserFloat("disc_negativity_photon",discNegativityPhoton);

        bs.addUserFloat("cosMuWColl_photon",      cosMuWCollPhoton); 
        bs.addUserFloat("cosMuWLhcbAlt_photon",   cosMuWLhcbAltPhoton); 
        bs.addUserFloat("cosMuWReco1_photon",     cosMuWReco1Photon); 
        bs.addUserFloat("cosMuWReco2_photon",     cosMuWReco2Photon); 

        bs.addUserFloat("cosPhiDsColl_photon",    cosPhiDsCollPhoton); 
        bs.addUserFloat("cosPhiDsLhcbAlt_photon", cosPhiDsLhcbAltPhoton); 
        bs.addUserFloat("cosPhiDsReco1_photon",   cosPhiDsReco1Photon);
        bs.addUserFloat("cosPhiDsReco2_photon",   cosPhiDsReco2Photon); 

        bs.addUserFloat("cosPiDsColl_photon",     cosPiDsCollPhoton); 
        bs.addUserFloat("cosPiDsLhcbAlt_photon",  cosPiDsLhcbAltPhoton); 
        bs.addUserFloat("cosPiDsReco1_photon",    cosPiDsReco1Photon); 
        bs.addUserFloat("cosPiDsReco2_photon",    cosPiDsReco2Photon); 

        bs.addUserFloat("cosPlaneBsColl_photon",    cosPlaneBsCollPhoton);
        bs.addUserFloat("cosPlaneBsLhcbAlt_photon", cosPlaneBsLhcbAltPhoton);
        bs.addUserFloat("cosPlaneBsReco1_photon",   cosPlaneBsReco1Photon);
        bs.addUserFloat("cosPlaneBsReco2_photon",   cosPlaneBsReco2Photon);

        bs.addUserFloat("cosPlaneDsColl_photon",    cosPlaneDsCollPhoton);
        bs.addUserFloat("cosPlaneDsLhcbAlt_photon", cosPlaneDsLhcbAltPhoton);
        bs.addUserFloat("cosPlaneDsReco1_photon",   cosPlaneDsReco1Photon);
        bs.addUserFloat("cosPlaneDsReco2_photon",   cosPlaneDsReco2Photon);



        /////////////////////// END OF VARIABLE DEFINITION //////////////////////
        //std::cout << "mu pt is" << muPtr->pt() << std::endl;
        //std::cout << "saving.."<< std::endl; 
        //arrived = 1;
        //bs.addUserInt("arrived", arrived);
        //arrived = -1;
        //append candidate at the end of our return value :)
        //ret_value can be a vector!!
        bsCandidates->emplace_back(bs);
        //ret_value_gen->emplace_back(gen);

        } //closing pi loop
      } //closing k2 loop
    } //closing k1 loop
  } //closing trg muon loop
  
  iEvent.put(std::move(bsCandidates), "bs");
  /*
  if(arrived >0){
  //std::cout << "arrived: " << arrived << std::endl;
  iEvent.put(std::move(bsCandidates), "bs");
  }
  else{

  //std::cout << "arrived: " << arrived << std::endl;
  pat::CompositeCandidate bs;
  bs.addUserInt("arrived",arrived);
  bsCandidates->emplace_back(bs);
  iEvent.put(std::move(bsCandidates), "bs"); 
  }
  */
}//closing event loop

void BsToDsPhiKKPiMuBuilder::endJob(){
// Printouts:

std::cout << "\n--------- Bs BUIDLER MODULE ----------\n" << std::endl;
std::cout << "#Muons  in file                                    : " << nMuons  << std::endl;
std::cout << "#Tracks in file                                    : " << nTracks << std::endl;
std::cout << "#Pv     in file                                    : " << nPv     << std::endl;
std::cout << "\n#Muons for which we found primary vertex         : " << muSelCounter  << std::endl;
std::cout << "#Kaon 1 which passed the hadronic selection        : " << k1Sel1Counter << std::endl;
std::cout << "#Kaon 1 which passed the angular  selection        : " << k1Sel2Counter << std::endl;
std::cout << "#Kaon 2 which passed the hadronic selection        : " << k2Sel1Counter << std::endl;
std::cout << "#Kaon 2 which passed the angular  selection        : " << k2Sel2Counter << std::endl;
std::cout << "#Pions  which passed the hadronic selection        : " << piSel1Counter << std::endl;
std::cout << "#Pions  which passed the angular  selection        : " << piSel2Counter << std::endl;

std::cout << "\n#KKPiMu combinations:" << nKKPiMu << std::endl;
std::cout << "#KKPiMu combinations which passed the Phi Mass cut : " << nPhiMassCut << std::endl;
std::cout << "#KKPiMu combinations which passed the Ds Mass cut  : " << nDsMassCut << std::endl;
std::cout << "#KKPiMu combinations which passed the Bs Mass cut  : " << nBsMassCut << std::endl;
std::cout << "#KKPiMu combinations which passed the Phi fit      : " << nPhiFit << std::endl;
std::cout << "#KKPiMu combinations which passed the Ds  fit      : " << nDsFit << std::endl;
std::cout << "#KKPiMu combinations which passed the Bs  fit      : " << nBsFit << std::endl;
}


DEFINE_FWK_MODULE(BsToDsPhiKKPiMuBuilder);
