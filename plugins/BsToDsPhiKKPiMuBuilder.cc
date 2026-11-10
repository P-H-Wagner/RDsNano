// Standards
#include <vector>
#include <memory>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <utility>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

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
#include "KinVtxFitter.h"   //--> not needed now 

// photons 
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"

// Muon, tracks and vertices
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" 
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h" 

// Helper functions and ROOT 
#include "trackCorrection.h"
#include "helper.h"         // helper functions
#include "PVRefitter.h"     // refit the primary vertex 
#include "TLorentzVector.h" 
#include "TVector3.h" 
#include "TMatrixDSym.h"
#include "TVectorD.h"

// for gen matching
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h" 
// B field
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
// for 3DPoint
#include "DataFormats/Math/interface/Point3D.h"

// vtx probability and lxy
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"

///////////////////////////////////////////////////////////////////////////////////////////////
// TODOS:                                                                                    // 
//                                                                                           //
// - move gen matching in separate module?             - DONE                                //
// - remove hardcoded numbers                          - DONE                                //
// - helicity plane angles                             - DONE                                //
// - redefine all variables after the fit? save both ? - YES                                 //
// - beautify the bs.addUserFloat (...)                - DONE                                //
// - pos. def. cov matrix                              - DONE                                //
// - add counters before every selection               - DONE                                //
// - pruned vs packed -> discuss                       - DONE                                //
// - output tree has now empty entries when there is no trigger/signal -> DONE (ed filter)   //
// - adapt for Hb background sample                    - DONE                                //
// - how to save kk same sign pair?                    - DONE                                //
// - generally: save only gen matched signals?         - NO                                  //
// - what if an event has two signals?                 - SAVE BOTH!                          //
// - divide into submitter chunks                      - DONE                                //
// - save gen information!!                            - DONE                                // 
// - do gen tests, check f.e. refitted p's with gen p's and unfitted p's - DONE              //
// - put hel angle calc. etc into functions!           - DONE                                //
// - isAncestor and getAncestor are save, they dont modify the Ptr - CHECKED                 //
// - add Ds boost as dicrimanting variable             - DONE                                //
// - add all impact parameters correctly and use refitted tracks - DONE                      //
// - add mu isolation, e* miss and pt miss, photon energy - DONE                             //
// - add and save Vcb Variables                        - DONE                                //
// - Oscillate back also all b hadrons, pick first beauty ancestor - DONE                    //
// - muon isolation, ds isolation (also wrt SV)        - DONE                                //       
// - muon impact parameter                             - DONE                                //
// - move photon tagging to separate module                                                  //
// - Question: Muon is from pat::Muons --> how to remove it from pv cands? To check!         //
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// counters for filters to know when we lose how many particles                              //
///////////////////////////////////////////////////////////////////////////////////////////////

//this is not thread save (use one atomic)
int  nMuons  = 0;        
int  nTracks = 0;        
int  nPv     = 0;        
int  muSelCounter  = 0;   
int  k1Sel1Counter = 0;   
int  k2Sel1Counter = 0;
int  piSel1Counter = 0;
int  gSelCounter   = 0;   
int  nKKPiMu       = 0;
int  nPhiMassCut   = 0;
int  nDsMassCut    = 0;
int  nBsMassCut    = 0;
int  nPhiFit       = 0;
int  nDsFit        = 0;
int  nBsFit        = 0;
int nAod           = 0;
int nVtxCut        = 0;
int nCosCut        = 0;
int nSurvived      = 0;

///////////////////////////////////////////////////////////////////////////////////////////////
// Start of class definition                                                                 //
///////////////////////////////////////////////////////////////////////////////////////////////

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
  //reco::TransientTrack getTransientTrack(const reco::Track track) const {    
  //    reco::TransientTrack transientTrack(track, paramField);

  //    return transientTrack;
  //}


private:
  
  //Bfield
  //OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");

  //defintion of these constants can be found in python/BsToDsPhiKKPiMu_cff.py

  //cuts 
  const StringCutObjectSelector<pat::PackedCandidate> hadSelection_;     // cut on hadrons
  const StringCutObjectSelector<pat::PackedCandidate> gSelection_;       // cut on photons (gamma)
  const StringCutObjectSelector<reco::GenParticle>    hadSelectionGen_;  // cut on gen hadrons for test with pruned only!

  //variables
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

  //beamspot
  const edm::InputTag beamSpotTag;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;

  //packed candidates
  const edm::InputTag srcTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> src_;

  //for the muons
  const edm::InputTag trgMuonTag;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuons_;

  //primary vertices
  const edm::InputTag primaryVtxTag;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVtx_;

  //tracks for isolation
  const edm::InputTag isoTracksTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isoTracks_;

  //lost tracks for isolation
  const edm::InputTag tracksLostTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksLost_;
 
  //gen for gen-matching
  const edm::InputTag prunedGenTag; //pruned is a compressed packed format
  const edm::EDGetTokenT<reco::GenParticleCollection> prunedGen_;

  const edm::InputTag packedGenTag; //packed contains much more info->most likely not needed!
  const edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGen_;
 
  //photons
  const edm::InputTag photonTag; 
  const edm::EDGetTokenT<reco::PhotonCoreCollection> photon_;

  //unpacked tracks for pv refit
  const edm::InputTag unpackedTracksTag;
  const edm::EDGetTokenT<reco::TrackCollection> unpackedTracks_;

  //pv with  bs info
  const edm::InputTag vtxWithBsTag;
  const edm::EDGetTokenT<reco::VertexCollection> vtxWithBs_;

  // track scale correction
  const std::string  trkCalibFile_;
  std::unique_ptr<const rdscalib::TrackCalib> trkCalib;
  std::unique_ptr<const rdscalib::TrackCalib> muCalib;

  // for the mass error propagation (load json file and access a and c)
  boost::property_tree::ptree jsonCalib;
  double a2[3];
  double c2;


};

//define the constructor
BsToDsPhiKKPiMuBuilder::BsToDsPhiKKPiMuBuilder(const edm::ParameterSet& iConfig):

    hadSelection_       (iConfig.getParameter<std::string>("hadSelection"       )),
    gSelection_         (iConfig.getParameter<std::string>("gSelection"         )),
    hadSelectionGen_    (iConfig.getParameter<std::string>("hadSelectionGen"    )),
    maxdRHadMuon_       (iConfig.getParameter<double>     ("maxdRHadMuon"       )),
    mindRHadMuon_       (iConfig.getParameter<double>     ("mindRHadMuon"       )),
    maxdRPhotonDs_      (iConfig.getParameter<double>     ("maxdRPhotonDs"      )),
    maxdzDiffHadMuon_   (iConfig.getParameter<double>     ("maxdzDiffHadMuon"   )),
    maxdxyHadPv_        (iConfig.getParameter<double>     ("maxdxyHadPv"        )),
    phiMassAllowance_   (iConfig.getParameter<double>     ("phiMassAllowance"   )),
    dsMassAllowance_    (iConfig.getParameter<double>     ("dsMassAllowance"    )),
    dsStarMassAllowance_(iConfig.getParameter<double>     ("dsStarMassAllowance")),
    drMatchGen_         (iConfig.getParameter<double>     ("drMatchGen"         )),
    maxBsMass_          (iConfig.getParameter<double>     ("maxBsMass"          )),

    piMass_             (iConfig.getParameter<double>     ("piMass"             )),
    piMassSigma_        (iConfig.getParameter<double>     ("piMassSigma"        )),
    kMass_              (iConfig.getParameter<double>     ("kMass"              )),
    kMassSigma_         (iConfig.getParameter<double>     ("kMassSigma"         )),
    phiMass_            (iConfig.getParameter<double>     ("phiMass"            )),
    constrainPhiMass_   (iConfig.getParameter<bool>       ("constrainPhiMass"   )),
    minPhiVtxProb_      (iConfig.getParameter<double>     ("minPhiVtxProb"      )),
    dsMass_             (iConfig.getParameter<double>     ("dsMass"             )),
    constrainDsMass_    (iConfig.getParameter<bool>       ("constrainDsMass"    )),
    minDsVtxProb_       (iConfig.getParameter<double>     ("minDsVtxProb"       )),
    dsStarMass_         (iConfig.getParameter<double>     ("dsStarMass"         )),
    muMass_             (iConfig.getParameter<double>     ("muMass"             )),
    muMassSigma_        (iConfig.getParameter<double>     ("muMassSigma"        )),
    bsMass_             (iConfig.getParameter<double>     ("bsMass"             )),
    isoCone_            (iConfig.getParameter<double>     ("isoCone"            )),

    beamSpotTag(iConfig.getParameter<edm::InputTag>       ("offBeamSpot"        )),
    beamSpot_(consumes<reco::BeamSpot>                    (beamSpotTag          )),

    srcTag(iConfig.getParameter<edm::InputTag>            ("pfCand"             )),
    src_(consumes<pat::PackedCandidateCollection>         (srcTag               )), 

    trgMuonTag(iConfig.getParameter<edm::InputTag>        ("muCand"             )),
    trgMuons_(consumes<pat::MuonCollection>               (trgMuonTag           )), 

    primaryVtxTag(iConfig.getParameter<edm::InputTag>     ("pvCand"             )),
    primaryVtx_(consumes<reco::VertexCollection>          (primaryVtxTag        )),

    isoTracksTag(iConfig.getParameter<edm::InputTag>      ("tracks"             )),
    isoTracks_(consumes<pat::PackedCandidateCollection>   (isoTracksTag         )),

    tracksLostTag(iConfig.getParameter<edm::InputTag>     ("lostTracks"         )),
    tracksLost_(consumes<pat::PackedCandidateCollection>  (tracksLostTag        )),

    prunedGenTag(iConfig.getParameter<edm::InputTag>      ("prunedCand"         )),
    prunedGen_(consumes<reco::GenParticleCollection>      (prunedGenTag         )),
    packedGenTag(iConfig.getParameter<edm::InputTag>      ("packedCand"         )),
    packedGen_(consumes<pat::PackedGenParticleCollection> (packedGenTag         )),
 
    photonTag(iConfig.getParameter<edm::InputTag>         ("photonCand"         )),
    photon_(consumes<reco::PhotonCoreCollection>          (photonTag            )),

    unpackedTracksTag(iConfig.getParameter<edm::InputTag> ("unpackedTracksCand" )),
    unpackedTracks_(consumes<reco::TrackCollection>       (unpackedTracksTag    )),
 
    vtxWithBsTag(iConfig.getParameter<edm::InputTag>      ("vtxWithBsCand"      )),
    vtxWithBs_(consumes<reco::VertexCollection>           (vtxWithBsTag         )),
    trkCalibFile_       (iConfig.getParameter<std::string>("trkCalibFile"       ))
   {
    //for the track scale correction
    const std::string calibPath = edm::FileInPath(trkCalibFile_).fullPath();
    trkCalib = std::make_unique<const rdscalib::TrackCalib>(calibPath);
    muCalib  = std::make_unique<const rdscalib::TrackCalib>(calibPath);

    //for the mass error propagation
    jsonCalib = loadTrackCalibJSON(edm::FileInPath(trkCalibFile_).fullPath());
    int i = 0;
    for (const auto& entry : jsonCalib.get_child("per_bin.a2")) a2[i++] = entry.second.get_value<double>();
    c2 = jsonCalib.get_child("per_bin.c2").begin()->second.get_value<double>();

    // output collection
    produces<pat::CompositeCandidateCollection>        ("bs"                  );
   }

///////////////////////////////////////////////////////////////////////////////////////////////
// Start of Event loop                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////

void BsToDsPhiKKPiMuBuilder::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  //std::cout << "---------------- NEW EVENT ---------------" << std::endl;

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

  edm::Handle<reco::VertexCollection> vtxWithBs;
  iEvent.getByToken(vtxWithBs_, vtxWithBs);

  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> bsCandidates(new pat::CompositeCandidateCollection());

  // for printing
  //std::cout << std::fixed << std::setprecision(10);

  const auto& id    = iEvent.id();


  // consider also lost tracks!
  const size_t nPf = pcand->size() + tracksLost->size();

  // helper function 
  auto pfAt = [&](size_t pfIdx) -> edm::Ptr<pat::PackedCandidate> {

    if (pfIdx < pcand->size()) return edm::Ptr<pat::PackedCandidate>(pcand, pfIdx); //normal tracks
    else                       return edm::Ptr<pat::PackedCandidate>(tracksLost, pfIdx - pcand->size()); //lost tracks

  };


  ///////////////////////////////////////////////////
  // Check if we have a valid beamspot             //
  ///////////////////////////////////////////////////

  if (!beamSpotHandle.isValid()) {
    iEvent.put(std::make_unique<pat::CompositeCandidateCollection>(), "bs");
    return;
  }
  const reco::BeamSpot& beamSpot = *beamSpotHandle;


  //////////////////////////
  // Prepare the PV refit //
  //////////////////////////

  // we refit all vertices from the vtxWithBs collection
  PVRefitter fitter;

  std::vector<reco::Vertex> refitAll;                          // refit with all tracks
  std::vector<std::vector<edm::Ptr<reco::Track>>> vtxTracks;   // tracks of each vertex
  
  for (size_t v = 0; v < vtxWithBs->size(); ++v) {
    const reco::Vertex& vtx = vtxWithBs->at(v);
  
    std::vector<edm::Ptr<reco::Track>> ptrs;
    std::vector<reco::Track> tracks;
    for (size_t t = 0; t < vtx.tracksSize(); ++t) {
      edm::Ptr<reco::Track> p = vtx.trackRefAt(t).castTo<edm::Ptr<reco::Track>>();
      ptrs.push_back(p);
      tracks.push_back(*p);
    }
    vtxTracks.push_back(ptrs);
    refitAll.push_back(fitter.fit(tracks, beamSpot));
  }


  //////////////////////////////
  // collect all GOOD photons //
  //////////////////////////////

  std::vector<size_t> goodPhotonsIdx;
  std::vector<edm::Ptr<pat::PackedCandidate>> goodPhotons;

  for (size_t gIdx = 0; gIdx < nPf ; ++gIdx){

    edm::Ptr<pat::PackedCandidate> gPtr;
    if (gIdx < pcand->size()) gPtr = edm::Ptr<pat::PackedCandidate>(pcand, gIdx); //normal tracks
    else gPtr = edm::Ptr<pat::PackedCandidate>(tracksLost, gIdx - pcand->size());  //lost tracks 
  
    // if the photon does not pass the selection, jump to the next!
    if(!gSelection_(*gPtr)) continue;

    goodPhotonsIdx.push_back(gIdx);
    goodPhotons   .push_back(gPtr);
  }

  //////////////////////////////
  // prepare AOD matching     //
  //////////////////////////////

  //this is a lambda function, takes my reference, returns int
  //matches reco to AOD in dr and charge
  auto matchAodTrack = [&](double eta, double phi, int charge) -> int {

  int    bestIdx  = -1;
  double mindR    = 1e9;
  //loop over AOD and match
  for (size_t t = 0; t < unpackedTracks->size(); ++t) {

    const reco::Track& trk = (*unpackedTracks)[t];

    //match in charge and dR
    if (trk.charge() != charge) continue;
    const double dr = reco::deltaR(trk.eta(), trk.phi(), eta, phi);

    if (dr < mindR) { 
      mindR = dr; 
      bestIdx = (int)t; 
    }
   }

   return bestIdx;
  };

  //prepare a vector of size nPf (normal + lost), filled with -2
  std::vector<int> aodMatchPf(nPf, -2);

  // another lambda explicitly to match on pf collection
  auto matchAodTrackPf = [&](size_t idx) -> int {
    if (aodMatchPf[idx] == -2) {
      const auto p = pfAt(idx);
      aodMatchPf[idx] = matchAodTrack(p->eta(), p->phi(), p->charge());
    }
    return aodMatchPf[idx];
  };

  //////////////////////////////////////////////////////
  // Acces the trigger muon inside the                //
  // packed pat collection. Reminder, we often have   //
  // only one trigger muon, but sometimes two! :)     //
  // --> checked, and its valid to have more than one //
  //////////////////////////////////////////////////////

  for(size_t trgMuIdx = 0; trgMuIdx < trgMuons->size(); ++trgMuIdx){

    nMuons++;

    edm::Ptr<pat::Muon> muPtr(trgMuons, trgMuIdx);

    std::vector<float> dzMuPV;

    // For the first candidate selection, we pick the PV to be the one closest to the trg Muon in dz.
    // The final PV is calculated later.

    float minDz  = 0.0;
    int pvIdx_dz = -1;
    reco::Vertex pv_dz;
    const reco::Track muTrack = *muPtr->bestTrack();

    for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){
      nPv++;
      edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, vtxIdx);
      float dz  = fabs(muTrack.dz (vtxPtr->position())); 
      if((dz < minDz) || (pvIdx_dz < 0)){
        minDz    = dz;
        pvIdx_dz = vtxIdx;
      }
    }

    float minDz_w_bs  = 0.0;
    int pvIdx_dz_w_bs = -1;
    reco::Vertex pv_dz_w_bs;

    for(size_t vtxIdx = 0; vtxIdx < vtxWithBs->size(); ++vtxIdx){
      nPv++;
      edm::Ptr<reco::Vertex> vtxPtr(vtxWithBs, vtxIdx);
      float dz  = fabs(muTrack.dz (vtxPtr->position())); 
      if((dz < minDz_w_bs) || (pvIdx_dz_w_bs < 0)){
        minDz_w_bs    = dz;
        pvIdx_dz_w_bs = vtxIdx;
      }
    }

    if (pvIdx_dz <0) continue;
    else pv_dz = primaryVtx->at(pvIdx_dz);

    if (pvIdx_dz_w_bs <0) continue;
    else pv_dz_w_bs = vtxWithBs->at(pvIdx_dz_w_bs);


    /////////////////////////////
    // get matching aod track  //
    /////////////////////////////

    const int trkIdxMu = matchAodTrack(muTrack.eta(), muTrack.phi(), muPtr->charge());
    if (trkIdxMu < 0) continue; 


    //////////////////////////////////////
    // correct momentum scale and smear //
    //////////////////////////////////////

    double f_mu = -1.0; 
    math::PtEtaPhiMLorentzVector muVec(muPtr->pt(), muPtr->eta(), muPtr->phi(), muMass_);
    if (!iEvent.isRealData()) {

      const double z_mu = muCalib ->gauss(id.run(), id.luminosityBlock(), static_cast<int64_t>(id.event()), muPtr->userInt("muIdx")); //take the index in the reco::muon collection
      f_mu              = muCalib ->factor(muPtr->p(), muPtr->pt(), muPtr->eta(), muPtr->charge(), z_mu, "nominal", rdscalib::Kind::Muon);
      muVec = math::PtEtaPhiMLorentzVector(muVec.Pt() * f_mu, muVec.Eta(), muVec.Phi(), muMass_);

    }

    muSelCounter++;

    ///////////////////////////////////////////////////
    // Before we start the kkpi event loop,          //
    // define already collections of                 //
    // charged/neutral pf Idx for the muon isolation //
    ///////////////////////////////////////////////////

    std::set<int> neutralPdgId  = {22,  111, 130, 310, 2112 };
    std::set<int> chargedPdgId  = {211, 321, 11,  13,  2212, 999211};
  
    std::vector<int> chargedPfIdx;
    std::vector<int> neutralPfIdx;
    std::vector<edm::Ptr<pat::PackedCandidate>>  chargedPfPtr; 
    std::vector<edm::Ptr<pat::PackedCandidate>>  neutralPfPtr; 
    std::vector<reco::TransientTrack> chargedTT;
 
    for(size_t pfIdx = 0; pfIdx < nPf ; ++pfIdx) {
  
      //define a pointer to the pf cand at position pfIdx
      edm::Ptr<pat::PackedCandidate> pfPtr = pfAt(pfIdx);
    
      float deltaRMuTrk = reco::deltaR(*pfPtr, *muPtr);
      if      ( neutralPdgId.count(abs(pfPtr->pdgId())) && (deltaRMuTrk < 0.4) && (deltaRMuTrk > 0.01  ) && (pfPtr->pt() > 0.5 )) {
        neutralPfIdx.push_back(pfIdx); 
        neutralPfPtr.push_back(pfPtr); 
      } 
      else if ( chargedPdgId.count(abs(pfPtr->pdgId())) && (deltaRMuTrk < 0.4) && (deltaRMuTrk > 0.0001)                        ) {

        if ((abs(pfPtr->pdgId()) == 11) || (abs(pfPtr->pdgId()) == 13)) continue;
        if (!pfPtr->hasTrackDetails())  continue;

        chargedPfIdx.push_back(pfIdx); 
        chargedPfPtr.push_back(pfPtr); 
        chargedTT.push_back(ttBuilder->build(pfPtr->bestTrack()));
      }    
    }

    /////////////////////////////////
    // colelct all GOOD hadrons :) //
    /////////////////////////////////

    std::vector<edm::Ptr<pat::PackedCandidate>> goodHads;
    std::vector<size_t> goodHadsIdx; //keep the PF index

    for(size_t hadIdx = 0; hadIdx < nPf ; ++hadIdx) {

      nTracks++;
      //define a pointer to the kaon at position hadIdx
      edm::Ptr<pat::PackedCandidate> hadPtr;
      if (hadIdx < pcand->size()) hadPtr = edm::Ptr<pat::PackedCandidate>(pcand, hadIdx); //normal tracks
      else hadPtr = edm::Ptr<pat::PackedCandidate>(tracksLost, hadIdx - pcand->size());  //lost tracks
    
      if (!hadSelection_(*hadPtr)) continue; 
    
      float muonHdR = reco::deltaR(*hadPtr,*muPtr);

      bool hSel = (( muonHdR < maxdRHadMuon_ )                                                               && 
      (reco::deltaR(*hadPtr, *muPtr) > mindRHadMuon_)                                                        && 
      (fabs(hadPtr->bestTrack()->dz(pv_dz.position()) - muTrack.dz(pv_dz.position())) < maxdzDiffHadMuon_)) &&
      (fabs(hadPtr->bestTrack()->dxy(pv_dz.position())) < maxdxyHadPv_ ) ;
 
      if (!hSel) continue;

      //append to good hadrons
      goodHads   .push_back(hadPtr);
      goodHadsIdx.push_back(hadIdx);

    }

    ///////////////////////////////////
    // Loop over k1 good tracks      //
    ///////////////////////////////////

    for(size_t a = 0; a < goodHads.size() ; ++a) {
      k1Sel1Counter++;

      ///////////////////////////////////
      // Loop over k2 good tracks      //
      ///////////////////////////////////

      for(size_t b = a + 1; b < goodHads.size() ; ++b) {

        //make sure k2 is not k1
        if (b == a) continue;
        k2Sel1Counter++;

        edm::Ptr<pat::PackedCandidate> k1Ptr = goodHads[a];
        edm::Ptr<pat::PackedCandidate> k2Ptr = goodHads[b];
        size_t k1Idx = goodHadsIdx[a];
        size_t k2Idx = goodHadsIdx[b];
  
        if (k1Ptr->pt() < k2Ptr->pt()) {
          std::swap(k1Ptr, k2Ptr);
          std::swap(k1Idx, k2Idx);
        }

        float muonK1dR = reco::deltaR(*k1Ptr,*muPtr);
        float muonK2dR = reco::deltaR(*k2Ptr,*muPtr);

        //std::cout << " found k1: " << k1Ptr->pt() << std::endl;
        //std::cout << " --- found k2: " << k2Ptr->pt() << std::endl;

        //////////////////////////////////////////////////
        // Build Phi resonance                          //
        //////////////////////////////////////////////////

        math::PtEtaPhiMLorentzVector k1Vec(k1Ptr->pt(), k1Ptr->eta(), k1Ptr->phi(), kMass_ );
        math::PtEtaPhiMLorentzVector k2Vec(k2Ptr->pt(), k2Ptr->eta(), k2Ptr->phi(), kMass_ );

        // correct momentum scale and smear
        double f_k1 = -1.0; 
        double f_k2 = -1.0; 

        if (!iEvent.isRealData()) {

          const double z_k1 = trkCalib->gauss(id.run(), id.luminosityBlock(), static_cast<int64_t>(id.event()), k1Idx);
          const double z_k2 = trkCalib->gauss(id.run(), id.luminosityBlock(), static_cast<int64_t>(id.event()), k2Idx);

          f_k1 = trkCalib->factor(k1Ptr->p(), k1Ptr->pt(), k1Ptr->eta(), k1Ptr->charge(), z_k1);
          f_k2 = trkCalib->factor(k2Ptr->p(), k2Ptr->pt(), k2Ptr->eta(), k2Ptr->charge(), z_k2);

          k1Vec = math::PtEtaPhiMLorentzVector(k1Vec.Pt() * f_k1, k1Vec.Eta(), k1Vec.Phi(), kMass_);
          k2Vec = math::PtEtaPhiMLorentzVector(k2Vec.Pt() * f_k2, k2Vec.Eta(), k2Vec.Phi(), kMass_);
        }
 

        //define a composite candidate pair 
        math::PtEtaPhiMLorentzVector kkVec = k1Vec + k2Vec;

        //only continue when they build a phi resonance (see cfg for mass window):
        if (fabs(kkVec.M() - phiMass_) > phiMassAllowance_) continue;     
        
        int kkCharge = k1Ptr->charge() + k2Ptr->charge();
        nPhiMassCut++;

        //////////////////////////////////////////////////
        // Loop over pi and select the good tracks      //
        //////////////////////////////////////////////////

        for(size_t c = 0; c < goodHads.size() ; ++c) {

          //make sure the pion is none of the kaons:
          if((c == a) || (c == b)) continue;

          piSel1Counter++;
          edm::Ptr<pat::PackedCandidate> piPtr = goodHads[c];
          const size_t piIdx = goodHadsIdx[c];

          float muonPidR = reco::deltaR(*piPtr,*muPtr);
          int piMuCharge = piPtr->charge() * muPtr->charge();
          //std::cout << " ----- found pi: " << piPtr->pt() << std::endl;

          // correct momentum scale and smear 
          //std::cout << "After pi: pt=" << piVec.Pt() << ", eta=" << piVec.Eta() << ", phi=" << piVec.Phi() << ", px=" << piVec.Px() << ", py=" << piVec.Py() << ", pz=" << piVec.Pz() << std::endl;
          math::PtEtaPhiMLorentzVector piVec(piPtr->pt(), piPtr->eta(), piPtr->phi(), piMass_); //just to be sure lets also force the pi mass

          double f_pi = -1.0; 
 
          if (!iEvent.isRealData()) {

            const double z_pi = trkCalib->gauss(id.run(), id.luminosityBlock(), static_cast<int64_t>(id.event()), piIdx);
            f_pi = trkCalib->factor(piPtr->p(), piPtr->pt(), piPtr->eta(), piPtr->charge(), z_pi);
            piVec = math::PtEtaPhiMLorentzVector(piVec.Pt() * f_pi, piVec.Eta(), piVec.Phi(), piMass_);
          
          }


          // check if any momentum correction is ill defined
          if (!iEvent.isRealData()) {
            if (!(f_mu > 0) || !(f_k1 > 0) || !(f_k2 > 0) || !(f_pi > 0) || !std::isfinite(f_mu + f_k1 + f_k2 + f_pi)) {
              //std::cout << "bad calib factor: " << f_mu << " " << f_k1 << " " << f_k2 << " " << f_pi << " event " << iEvent.id().event() << std::endl;
              continue;
            }
          }

          if (!std::isfinite(
              muVec.Px() + muVec.Py() + muVec.Pz() +
              k1Vec.Px() + k1Vec.Py() + k1Vec.Pz() +
              k2Vec.Px() + k2Vec.Py() + k2Vec.Pz() +
              piVec.Px() + piVec.Py() + piVec.Pz())) continue;



          //std::cout << "After pi: pt=" << piVec.Pt() << " k1: pt " << k1Vec.Pt() << "k2: pt" << k2Vec.Pt() << std::endl;  
          nKKPiMu++; // found candidate

          //////////////////////////////////////////////////
          // Build Ds resonance                           //
          // ///////////////////////////////////////////////

          math::PtEtaPhiMLorentzVector phiPiVec = kkVec + piVec;
          
          //only continue when they build a ds resonance, (see cfg for mass window):

          // CHANGE ME BACK!!
          //if (fabs(phiPi.mass() - dsMass_) > dsMassAllowance_) continue;
          if (phiPiVec.M() - dsMass_ < -0.220) continue;
          if (phiPiVec.M() - dsMass_ >  0.182) continue;

          //phiPi.setCharge(kk.charge() + piPtr->charge());
          int phiPiCharge = k1Ptr->charge() + k2Ptr->charge() + piPtr->charge();
          nDsMassCut++;

          //////////////////////////////////////////////////
          // Build Bs resonance                           //
          //////////////////////////////////////////////////

          math::PtEtaPhiMLorentzVector dsMuVec = phiPiVec + muVec;

          if(dsMuVec.M() > maxBsMass_) continue;
          nBsMassCut++;

          int dsMuCharge = k1Ptr->charge() + k2Ptr->charge() + piPtr->charge() + muPtr->charge();

          //build bs with collinear approximation
          pat::CompositeCandidate bs;

          bs.setP4(dsMuVec * bsMass_ / dsMuVec.M()); //the bs_mass will thus be fixed at 5.36688 GeV (peak in the histo)
          bs.setCharge(dsMuCharge);

          ///////////////////////////////////////////////////
          // Now we do a proper fit                        //
          ///////////////////////////////////////////////////

          //define a factory
          KinematicParticleFactoryFromTransientTrack pFactory;

          //define the vector for the particles to be fitted
          std::vector<RefCountedKinematicParticle> phiToFit;
          std::vector<RefCountedKinematicParticle> dsToFit;
          std::vector<RefCountedKinematicParticle> bsToFit;

          // add masses
          ParticleMass piMass  = piMass_ ;
          ParticleMass kMass   = kMass_  ;
          ParticleMass phiMass = phiMass_;
          ParticleMass dsMass  = dsMass_ ;
          ParticleMass muMass  = muMass_ ;

          
          float ndf = 0.0;
          float chi = 0.0;
          float kMassSigma  = kMassSigma_ ;
          float piMassSigma = piMassSigma_;
          float muMassSigma = muMassSigma_;

          // fix the tracks to have pos def covariance matrix, taken from:
          // https://github.com/CMSKStarMuMu/miniB0KstarMuMu/blob/master/miniKstarMuMu/plugins/miniKstarMuMu.cc#L1611-L1678

          const reco::Track k1Track = *k1Ptr->bestTrack();
          const reco::Track k2Track = *k2Ptr->bestTrack();
          const reco::Track piTrack = *piPtr->bestTrack();

          //////////////////////////////////////////////
          // Build track for corrected momentum scale //
          //////////////////////////////////////////////

          reco::Track::Vector trackVectorMu(muVec.Px(), muVec.Py(), muVec.Pz());
          reco::Track::Vector trackVectorK1(k1Vec.Px(), k1Vec.Py(), k1Vec.Pz());
          reco::Track::Vector trackVectorK2(k2Vec.Px(), k2Vec.Py(), k2Vec.Pz());
          reco::Track::Vector trackVectorPi(piVec.Px(), piVec.Py(), piVec.Pz());


          //std::cout << "k1 track pt before: " << k1Track.pt() << std::endl;


          reco::Track newTrackMu ( muTrack.chi2(), muTrack.ndof(), muTrack.referencePoint(),
                                   trackVectorMu, //thats the only thing that changed!
                                   muTrack.charge(), muTrack.covariance() , muTrack.algo());

          reco::Track newTrackK1 ( k1Track.chi2(), k1Track.ndof(), k1Track.referencePoint(),
                                   trackVectorK1, //thats the only thing that changed!
                                   k1Track.charge(), k1Track.covariance() , k1Track.algo());

          reco::Track newTrackK2 ( k2Track.chi2(), k2Track.ndof(), k2Track.referencePoint(),
                                   trackVectorK2, //thats the only thing that changed!
                                   k2Track.charge(), k2Track.covariance() , k2Track.algo());

          reco::Track newTrackPi ( piTrack.chi2(), piTrack.ndof(), piTrack.referencePoint(),
                                   trackVectorPi, //thats the only thing that changed!
                                   piTrack.charge(), piTrack.covariance() , piTrack.algo());

          //std::cout << "k1 track pt: " << newTrackK1.pt() << std::endl;

          //////////////////////////////////////////////////
          // Fix covariane matrix to be positive definite //
          //////////////////////////////////////////////////

          //those tracks have the new covariance matrix assigned :) 
          reco::Track muTrackCorr   = correctCovMat(&newTrackMu, 1e-8);
          reco::Track k1TrackCorr   = correctCovMat(&newTrackK1, 1e-8);
          reco::Track k2TrackCorr   = correctCovMat(&newTrackK2, 1e-8);
          reco::Track piTrackCorr   = correctCovMat(&newTrackPi, 1e-8);

          reco::TransientTrack ttMu = ttBuilder->build(muTrackCorr);
          reco::TransientTrack ttK1 = ttBuilder->build(k1TrackCorr);
          reco::TransientTrack ttK2 = ttBuilder->build(k2TrackCorr);
          reco::TransientTrack ttPi = ttBuilder->build(piTrackCorr);

          phiToFit.push_back(  pFactory.particle(ttK1, kMass,  chi, ndf, kMassSigma ));
          phiToFit.push_back(  pFactory.particle(ttK2, kMass,  chi, ndf, kMassSigma ));
          dsToFit.push_back (  pFactory.particle(ttPi, piMass, chi, ndf, piMassSigma));
          bsToFit.push_back (  pFactory.particle(ttMu, muMass, chi, ndf, muMassSigma));

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

          ///////////////////////////////////////////////////
          // PHI VTX FIT                                   //
          ///////////////////////////////////////////////////
           
          RefCountedKinematicTree phiTree;
          try { phiTree = vertexFit(phiToFit, phiMass, constrainPhiMass_); }
          catch (const std::exception&) { continue; } // degenerate input: drop the candidate, not the job

          if (!phiTree->isValid() || phiTree->isEmpty() || !phiTree->isConsistent()) continue; //check if fit result is valid

          //access the fitted resonance  
          phiTree->movePointerToTheTop();
          RefCountedKinematicParticle phiParticle = phiTree->currentParticle();

          //get vtx chi2 and ndof
          auto phiVtx = phiTree->currentDecayVertex();
          if (!phiVtx->vertexIsValid() || !phiParticle->currentState().isValid() ) continue; //check if fit result is valid

          float phiVtxChi2    = phiVtx->chiSquared();
          if (phiVtxChi2 < 0) continue;

          float phiVtxNDof    = phiVtx->degreesOfFreedom();
          float phiVtxRedChi2 = phiVtxChi2 / phiVtxNDof; 
          float phiVtxProb    = ChiSquaredProbability(phiVtxChi2, phiVtxNDof); 
          if (phiVtxProb < minPhiVtxProb_ ) continue;

          //std::cout << "we passed the phi vtx fit prob" << std::endl; 
          
          nPhiFit++;

          // access refitted children
          phiTree->movePointerToTheFirstChild();
          RefCountedKinematicParticle phiDau1 = phiTree->currentParticle();
          phiTree->movePointerToTheNextChild();
          RefCountedKinematicParticle phiDau2 = phiTree->currentParticle();

          // get the vectors full of fit information (vertex and momenta and mass)
          AlgebraicVector7 phiParams     = phiParticle->currentState().kinematicParameters().vector();     
          //AlgebraicVector7 phiDau1Params = phiDau1->currentState().kinematicParameters().vector();     
          //AlgebraicVector7 phiDau2Params = phiDau2->currentState().kinematicParameters().vector();
     
          // add the phi to the list of particles (pi) to fit the ds 
          dsToFit.push_back(phiParticle);

          ///////////////////////////////////////////////////
          // DS VTX FIT                                    //
          ///////////////////////////////////////////////////

          RefCountedKinematicTree dsTree;
          try { dsTree = vertexFit(dsToFit, dsMass, constrainDsMass_); }
          catch (const std::exception&) { continue; } // degenerate input: drop the candidate, not the job
 
          if (!dsTree->isValid() || dsTree->isEmpty() ) continue; //check if fit result is valid

          // access the fitted resonance 
          dsTree->movePointerToTheTop();
          RefCountedKinematicParticle dsParticle = dsTree->currentParticle();

          // get vtx chi2 and ndof
          RefCountedKinematicVertex dsVtx = dsTree->currentDecayVertex(); //compare to the access via AlgebraicVector7
          if (!dsVtx->vertexIsValid()) continue; //check if fit result is valid

          float dsVtxChi2    = dsVtx->chiSquared();
          if (dsVtxChi2 < 0) continue;

          float dsVtxNDof    = dsVtx->degreesOfFreedom();
          float dsVtxRedChi2 = dsVtxChi2 / dsVtxNDof; 
          float dsVtxProb    = ChiSquaredProbability(dsVtxChi2, dsVtxNDof); 
          if (dsVtxProb < minDsVtxProb_) continue;

          //std::cout << "we passed the ds vtx fit prob" << std::endl; 

          nDsFit++;

          // access refitted children
          dsTree->movePointerToTheFirstChild();
          RefCountedKinematicParticle dsDau1 = dsTree->currentParticle();
          dsTree->movePointerToTheNextChild();
          RefCountedKinematicParticle dsDau2 = dsTree->currentParticle();

          // get the vectors full of fit information (vertex and momenta and mass)
          AlgebraicVector7 dsParams     = dsParticle->currentState().kinematicParameters().vector();     
          //AlgebraicVector7 dsDau1Params = dsDau1->currentState().kinematicParameters().vector();     
          //AlgebraicVector7 dsDau2Params = dsDau2->currentState().kinematicParameters().vector();

          // add the ds to the list of particles to fit the bs
          bsToFit.push_back(dsParticle);
          
          ///////////////////////////////////////////////////
          // BS VTX FIT                                    //
          ///////////////////////////////////////////////////

          RefCountedKinematicTree bsTree;
          try { bsTree = vertexFit(bsToFit, 0., false); } // no constraint for bs because missing momentum
          catch (const std::exception&) { continue; }   // degenerate input: drop the candidate, not the job

          if (!bsTree->isValid() || bsTree->isEmpty() ) continue; //check if fit result is valid 

          // access the fitted resonance 
          bsTree->movePointerToTheTop();
          RefCountedKinematicParticle bsParticle = bsTree->currentParticle();

          // get vtx chi2 and ndof
          RefCountedKinematicVertex bsVtx = bsTree->currentDecayVertex();
          if (!bsVtx->vertexIsValid()) continue; //check if fit result is valid

          float bsVtxChi2    = bsVtx->chiSquared();
          if (bsVtxChi2 < 0) continue;

          float bsVtxNDof    = bsVtx->degreesOfFreedom();
          float bsVtxRedChi2 = bsVtxChi2 / bsVtxNDof; 
          float bsVtxProb    = ChiSquaredProbability(bsVtxChi2, bsVtxNDof); 
 
          nBsFit++;

          // access refitted children
          bsTree->movePointerToTheFirstChild();
          RefCountedKinematicParticle bsDau1 = bsTree->currentParticle();
          bsTree->movePointerToTheNextChild();
          RefCountedKinematicParticle bsDau2 = bsTree->currentParticle();

          // get the vectors full of fit information (vertex and momenta and mass)
          AlgebraicVector7 bsParams = bsParticle->currentState().kinematicParameters().vector();     
          //AlgebraicVector7 bsDau1Params = bsDau1->currentState().kinematicParameters().vector();     
          //AlgebraicVector7 bsDau2Params = bsDau2->currentState().kinematicParameters().vector();

          //std::cout << "we passed all the vtx fit" << std::endl;

          //////////////////////////////////// end of fitter /////////////////////////////////////

          ///////////////////////////////////////////////////
          // Calculate the mass resolution in a custom way //
          ///////////////////////////////////////////////////


          std::array<std::array<double, 3>, 3> covK1 = scaleAndCalibrateCovMat(&newTrackK1, a2, c2 ,f_k1);
          std::array<std::array<double, 3>, 3> covK2 = scaleAndCalibrateCovMat(&newTrackK2, a2, c2 ,f_k2);
          std::array<std::array<double, 3>, 3> covPi = scaleAndCalibrateCovMat(&newTrackPi, a2, c2 ,f_pi);

          //stack tracks and covariance matrices
          std::vector<std::array<std::array<double, 3>, 3>> covs;
          std::vector<reco::Track> tracks;
          std::vector<double> masses;
 
          covs.push_back(covK1); tracks.push_back(newTrackK1); masses.push_back(kMass_ );
          covs.push_back(covK2); tracks.push_back(newTrackK2); masses.push_back(kMass_ );
          covs.push_back(covPi); tracks.push_back(newTrackPi); masses.push_back(piMass_);


          //std::cout << "phiPi_m: " << phiPiVec.mass() << std::endl;     
          double dsErrCustom =  getMassSigma(&tracks, &covs, &masses);
          //std::cout << "mass error: " << dsErrCustom << std::endl;     
 

          // save the fitted resonances        
          TLorentzVector fittedPhi;
          TLorentzVector fittedDs;
          TLorentzVector fittedBs;

          fittedPhi.SetXYZM(   phiParams(3),     phiParams(4),     phiParams(5),      phiParams(6));
          fittedDs.SetXYZM(    dsParams(3),      dsParams(4),      dsParams(5),       dsParams(6));
          fittedBs.SetXYZM(    bsParams(3),      bsParams(4),      bsParams(5),       bsParams(6));


          ///////////////////////////////////////////////////
          // Match final state track with AOD tracks       //
          // We do this to perform a PV refit. The latter  //
          // requires more general tracks than in miniAOD  //
          ///////////////////////////////////////////////////
          
          // match them by delta R + charge. Matching by pointer is not possible, since the collections
          // are truly different! 
          const int trkIdxK1 = matchAodTrackPf(k1Idx);
          const int trkIdxK2 = matchAodTrackPf(k2Idx);
          const int trkIdxPi = matchAodTrackPf(piIdx);
       
          //match found? 
          if (trkIdxK1 < 0 || trkIdxK2 < 0 || trkIdxPi < 0) continue; 
          
          // all four AOD tracks must be distinct
          if (trkIdxK1 == trkIdxK2 || trkIdxK1 == trkIdxPi || trkIdxK2 == trkIdxPi ||
              trkIdxMu == trkIdxK1 || trkIdxMu == trkIdxK2 || trkIdxMu == trkIdxPi) continue;
       
          nAod++;
 
          edm::Ptr<reco::Track> trkMu(unpackedTracks, trkIdxMu);
          edm::Ptr<reco::Track> trkK1(unpackedTracks, trkIdxK1);
          edm::Ptr<reco::Track> trkK2(unpackedTracks, trkIdxK2);
          edm::Ptr<reco::Track> trkPi(unpackedTracks, trkIdxPi);

          //prepare the loop to remove signal tracks from refitted vertices 
          bool matchedMu = false;
          bool matchedK1 = false;
          bool matchedK2 = false;
          bool matchedPi = false;

          //copy the refitted Vtx
          std::vector<reco::Vertex> refitVtx = refitAll;   // copy: all vertices, signal tracks still inside

          for (size_t v = 0; v < vtxTracks.size(); ++v) {
          
            bool containsSignal = false;
            std::vector<reco::Track> toFit;
       
            //loop over tracks p associated to that vertex 
            for (const auto& p : vtxTracks[v]) {
              if (p == trkMu) { matchedMu = true; containsSignal = true; continue; }
              if (p == trkK1) { matchedK1 = true; containsSignal = true; continue; }
              if (p == trkK2) { matchedK2 = true; containsSignal = true; continue; }
              if (p == trkPi) { matchedPi = true; containsSignal = true; continue; }
              //only non-signal tracks want to be refitted
              toFit.push_back(*p);
            }
          
            // if the vertex indeed was associated with signal tracks, we refit
            if (containsSignal) refitVtx[v] = fitter.fit(toFit, beamSpot);   // overwrite only this one
          }
       
 
          //////////////////////////////////////////////////
          // Look for possible photons (g for gamma)      //
          //////////////////////////////////////////////////

          int foundPhoton       = 0;
          float minDeltaDsStar  = 0;
          int photonIdx         = -1;

          float dr_photon_ds = 0.0; 
          float dsPhoton_m   = 0.0;
          float dsPhotonMu_m = 0.0;
          float photon_pt    = 0.0; 
          float photon_eta   = 0.0; 
          float photon_phi   = 0.0; 
          int   photon_pdgid = 0.0; 

          math::PtEtaPhiMLorentzVector photonVec(std::nan(""),std::nan(""),std::nan(""),std::nan(""));
          edm::Ptr<pat::PackedCandidate> phtPtr;


          int photonMultiplicity = 0;


          for (size_t g = 0; g < goodPhotons.size(); ++g) {

            photonMultiplicity++; 
            gSelCounter++;
          
            const auto&  gPtr = goodPhotons[g];      
            const size_t gIdx = goodPhotonsIdx[g];   
  
            // define photon momentum
            math::PtEtaPhiMLorentzVector gVec(gPtr->pt(), gPtr->eta(), gPtr->phi(), 0.0);
            // define Ds* Momentum
            math::PtEtaPhiMLorentzVector dsStarVec   = phiPiVec + gVec;

            float mDsStar      = dsStarVec  .M();
            float dr           = reco::deltaR(gVec, phiPiVec);
            float deltaDsStar  = fabs(mDsStar - dsStarMass_); 
 
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

            // define photon Momentum
            photonVec = math::PtEtaPhiMLorentzVector(phtPtr->pt(), phtPtr->eta(), phtPtr->phi(), 0.0); 

            // define Ds* + Mu Momentum
            math::PtEtaPhiMLorentzVector dsStarMuVec = dsMuVec + photonVec;

            dsPhotonMu_m   = dsStarMuVec.M();

          }

          //get reco::Vertices from fit
          reco::Vertex recoSv = getRecoVertex(bsVtx );
          reco::Vertex recoTv = getRecoVertex(dsVtx );
          reco::Vertex recoFv = getRecoVertex(phiVtx);

          // set sv,tv, and fv from fit
          float sv_x = recoSv.position().x(); //bsParams(0);
          float sv_y = recoSv.position().y(); //bsParams(1);
          float sv_z = recoSv.position().z(); //bsParams(2);

          float tv_x = recoTv.position().x(); //dsParams(0);
          float tv_y = recoTv.position().y(); //dsParams(1);
          float tv_z = recoTv.position().z(); //dsParams(2);

          float fv_x = recoFv.position().x(); //phiParams(0);
          float fv_y = recoFv.position().y(); //phiParams(1);
          float fv_z = recoFv.position().z(); //phiParams(2);

          ///////////////////////////////////////////////////////////////
          // Redefinition of PV, this with ip3d is more excact!        //
          // (Around 5%). Use refitted Vtx without signal tracks       //
          // Use this definition for all the variables depending on PV //
          ///////////////////////////////////////////////////////////////

          float minIp3d_general    = 0.0;
          int pvIdx_general        = -1;
          reco::Vertex pv_general;

          for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){
            edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, vtxIdx);

            float ip3d = getIP3D(vtxPtr, sv_x, sv_y, sv_z, phiPiVec, muVec);

            if((ip3d < minIp3d_general) || (pvIdx_general < 0)){
              minIp3d_general = ip3d;
              pvIdx_general  = vtxIdx;
             }

           }

          if (pvIdx_general <  0) continue;
          if (pvIdx_general >= 0){
            pv_general = primaryVtx->at(pvIdx_general);
          }
 
    

          // get the x,y position of the pv by evaluating the beamspot density
          // at pv_z

          float pv_general_z = pv_general.z();
          float pv_general_x = beamSpot.x(pv_general_z);
          float pv_general_y = beamSpot.y(pv_general_z);

       
          //now use refitted vertices!

          float minIp3d= 0.0;
          int pvIdx= -1;
          reco::Vertex pv;
 
          for(size_t vtxIdx = 0; vtxIdx < refitVtx.size(); ++vtxIdx){

            const reco::Vertex& vtx = refitVtx[vtxIdx];

            float ip3d = getIP3D(vtx, sv_x, sv_y, sv_z, phiPiVec, muVec);

            if((ip3d < minIp3d) || (pvIdx < 0)){
              minIp3d= ip3d;
              pvIdx= vtxIdx;
             }

           }

          if (pvIdx<  0) continue;
          if (pvIdx>= 0){
            pv= refitVtx.at(pvIdx);
          }


          float pv_z = pv.z();
          float pv_x = pv.x(); 
          float pv_y = pv.y(); 

          float pv_x_w_bs = beamSpot.x(pv_z);
          float pv_y_w_bs = beamSpot.y(pv_z);

          //////////////////////////////////////////////

          //2d cosine

          ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> kkPi_xy;
          ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> kkPi_xyz;

          ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> ds_xy;
          ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> ds_xy_pv;

          ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> ds_xyz;
          ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> ds_xyz_pv;

          kkPi_xy   .SetXYZ(phiPiVec.px()  , phiPiVec.py()  , 0.0);
          kkPi_xyz  .SetXYZ(phiPiVec.px()  , phiPiVec.py()  , phiPiVec.pz());

          ds_xy     .SetXYZ(tv_x - sv_x , tv_y - sv_y , 0.0);
          ds_xy_pv  .SetXYZ(tv_x - pv_x , tv_y - pv_y , 0.0);

          ds_xyz_pv .SetXYZ(tv_x - pv_x , tv_y - pv_y , tv_z - pv_z);
          ds_xyz    .SetXYZ(tv_x - sv_x , tv_y - sv_y , tv_z - sv_z);

          // check if any vertex position is nan/inf
          if (!std::isfinite(sv_x) || !std::isfinite(sv_y) || !std::isfinite(sv_z) ||
              !std::isfinite(tv_x) || !std::isfinite(tv_y) || !std::isfinite(tv_z) ||
              !std::isfinite(pv_x) || !std::isfinite(pv_y) || !std::isfinite(pv_z)) continue;        


          nVtxCut++;

          double den = 0.0;
          
          den = kkPi_xy.R() * ds_xy.R();
          if (!(den > 1e-12)) continue;
          float ds_vtx_cosine_xy = kkPi_xy.Dot(ds_xy) / den;
          
          den = kkPi_xy.R() * ds_xy_pv.R();
          if (!(den > 1e-12)) continue;
          float ds_vtx_cosine_xy_pv = kkPi_xy.Dot(ds_xy_pv) / den;
          
          den = kkPi_xyz.R() * ds_xyz.R();
          if (!(den > 1e-12)) continue;
          float ds_vtx_cosine_xyz = kkPi_xyz.Dot(ds_xyz) / den;
          
          den = kkPi_xyz.R() * ds_xyz_pv.R();
          if (!(den > 1e-12)) continue;
          float ds_vtx_cosine_xyz_pv = kkPi_xyz.Dot(ds_xyz_pv) / den;
 
          // no abs here, since we really want the directions to be in the SAME direction, not antiparallel
          if (ds_vtx_cosine_xyz_pv < 0.8) continue;

          nCosCut++;

          ///////////////////////////////////////////////////////////////////////////////////////////
          // store basic variables                                                                 //
          ///////////////////////////////////////////////////////////////////////////////////////////

          //save the indices of the final states in the pruned Collection
          bs.addUserInt("k1_idx",k1Idx);
          bs.addUserInt("k2_idx",k2Idx);
          bs.addUserInt("pi_idx",piIdx);

          // still need to add the candidates (uncorrected) to match in dR in genMatcher
          bs.addUserCand("k1",k1Ptr);
          bs.addUserCand("k2",k2Ptr);
          bs.addUserCand("pi",piPtr);
          bs.addUserCand("mu",muPtr); 

          //add final states 
          bs.addUserFloat("mu_px"    , muVec.px()     );
          bs.addUserFloat("mu_py"    , muVec.py()     );
          bs.addUserFloat("mu_pz"    , muVec.pz()     );
          bs.addUserFloat("mu_pt"    , muVec.pt()     );
          bs.addUserFloat("mu_eta"   , muVec.eta()    );
          bs.addUserFloat("mu_phi"   , muVec.phi()    );
          bs.addUserFloat("mu_m"     , muVec.mass()   );
          bs.addUserFloat("mu_charge", muPtr->charge());
          bs.addUserInt  ("mu_pdgId" , muPtr->pdgId() );
          
          // Kaon 1
          bs.addUserFloat("k1_px"    , k1Vec.px()     );
          bs.addUserFloat("k1_py"    , k1Vec.py()     );
          bs.addUserFloat("k1_pz"    , k1Vec.pz()     );
          bs.addUserFloat("k1_pt"    , k1Vec.pt()     );
          bs.addUserFloat("k1_eta"   , k1Vec.eta()    );
          bs.addUserFloat("k1_phi"   , k1Vec.phi()    );
          bs.addUserFloat("k1_m"     , k1Vec.mass()   );
          bs.addUserFloat("k1_charge", k1Ptr->charge());
          bs.addUserInt  ("k1_pdgId" , k1Ptr->pdgId() );
          
          // Kaon 2
          bs.addUserFloat("k2_px"    , k2Vec.px()     );
          bs.addUserFloat("k2_py"    , k2Vec.py()     );
          bs.addUserFloat("k2_pz"    , k2Vec.pz()     );
          bs.addUserFloat("k2_pt"    , k2Vec.pt()     );
          bs.addUserFloat("k2_eta"   , k2Vec.eta()    );
          bs.addUserFloat("k2_phi"   , k2Vec.phi()    );
          bs.addUserFloat("k2_m"     , k2Vec.mass()   );
          bs.addUserFloat("k2_charge", k2Ptr->charge());
          bs.addUserInt  ("k2_pdgId" , k2Ptr->pdgId() );
          
          // Pion
          bs.addUserFloat("pi_px"    , piVec.px()     );
          bs.addUserFloat("pi_py"    , piVec.py()     );
          bs.addUserFloat("pi_pz"    , piVec.pz()     );
          bs.addUserFloat("pi_pt"    , piVec.pt()     );
          bs.addUserFloat("pi_eta"   , piVec.eta()    );
          bs.addUserFloat("pi_phi"   , piVec.phi()    );
          bs.addUserFloat("pi_m"     , piVec.mass()   );
          bs.addUserFloat("pi_charge", piPtr->charge());
          bs.addUserInt  ("pi_pdgId" , piPtr->pdgId() );

          // track covariance matrix
          bs.addUserFloat("mu_cov_00", muTrackCorr.covariance(0,0));
          bs.addUserFloat("mu_cov_01", muTrackCorr.covariance(0,1));
          bs.addUserFloat("mu_cov_02", muTrackCorr.covariance(0,2));
          bs.addUserFloat("mu_cov_03", muTrackCorr.covariance(0,3));
          bs.addUserFloat("mu_cov_04", muTrackCorr.covariance(0,4));
          bs.addUserFloat("mu_cov_11", muTrackCorr.covariance(1,1));
          bs.addUserFloat("mu_cov_12", muTrackCorr.covariance(1,2));
          bs.addUserFloat("mu_cov_13", muTrackCorr.covariance(1,3));
          bs.addUserFloat("mu_cov_14", muTrackCorr.covariance(1,4));
          bs.addUserFloat("mu_cov_22", muTrackCorr.covariance(2,2));
          bs.addUserFloat("mu_cov_23", muTrackCorr.covariance(2,3));
          bs.addUserFloat("mu_cov_24", muTrackCorr.covariance(2,4));
          bs.addUserFloat("mu_cov_33", muTrackCorr.covariance(3,3));
          bs.addUserFloat("mu_cov_34", muTrackCorr.covariance(3,4));
          bs.addUserFloat("mu_cov_44", muTrackCorr.covariance(4,4));
          
          bs.addUserFloat("k1_cov_00", k1TrackCorr.covariance(0,0));
          bs.addUserFloat("k1_cov_01", k1TrackCorr.covariance(0,1));
          bs.addUserFloat("k1_cov_02", k1TrackCorr.covariance(0,2));
          bs.addUserFloat("k1_cov_03", k1TrackCorr.covariance(0,3));
          bs.addUserFloat("k1_cov_04", k1TrackCorr.covariance(0,4));
          bs.addUserFloat("k1_cov_11", k1TrackCorr.covariance(1,1));
          bs.addUserFloat("k1_cov_12", k1TrackCorr.covariance(1,2));
          bs.addUserFloat("k1_cov_13", k1TrackCorr.covariance(1,3));
          bs.addUserFloat("k1_cov_14", k1TrackCorr.covariance(1,4));
          bs.addUserFloat("k1_cov_22", k1TrackCorr.covariance(2,2));
          bs.addUserFloat("k1_cov_23", k1TrackCorr.covariance(2,3));
          bs.addUserFloat("k1_cov_24", k1TrackCorr.covariance(2,4));
          bs.addUserFloat("k1_cov_33", k1TrackCorr.covariance(3,3));
          bs.addUserFloat("k1_cov_34", k1TrackCorr.covariance(3,4));
          bs.addUserFloat("k1_cov_44", k1TrackCorr.covariance(4,4));
          
          bs.addUserFloat("k2_cov_00", k2TrackCorr.covariance(0,0));
          bs.addUserFloat("k2_cov_01", k2TrackCorr.covariance(0,1));
          bs.addUserFloat("k2_cov_02", k2TrackCorr.covariance(0,2));
          bs.addUserFloat("k2_cov_03", k2TrackCorr.covariance(0,3));
          bs.addUserFloat("k2_cov_04", k2TrackCorr.covariance(0,4));
          bs.addUserFloat("k2_cov_11", k2TrackCorr.covariance(1,1));
          bs.addUserFloat("k2_cov_12", k2TrackCorr.covariance(1,2));
          bs.addUserFloat("k2_cov_13", k2TrackCorr.covariance(1,3));
          bs.addUserFloat("k2_cov_14", k2TrackCorr.covariance(1,4));
          bs.addUserFloat("k2_cov_22", k2TrackCorr.covariance(2,2));
          bs.addUserFloat("k2_cov_23", k2TrackCorr.covariance(2,3));
          bs.addUserFloat("k2_cov_24", k2TrackCorr.covariance(2,4));
          bs.addUserFloat("k2_cov_33", k2TrackCorr.covariance(3,3));
          bs.addUserFloat("k2_cov_34", k2TrackCorr.covariance(3,4));
          bs.addUserFloat("k2_cov_44", k2TrackCorr.covariance(4,4));
          
          bs.addUserFloat("pi_cov_00", piTrackCorr.covariance(0,0));
          bs.addUserFloat("pi_cov_01", piTrackCorr.covariance(0,1));
          bs.addUserFloat("pi_cov_02", piTrackCorr.covariance(0,2));
          bs.addUserFloat("pi_cov_03", piTrackCorr.covariance(0,3));
          bs.addUserFloat("pi_cov_04", piTrackCorr.covariance(0,4));
          bs.addUserFloat("pi_cov_11", piTrackCorr.covariance(1,1));
          bs.addUserFloat("pi_cov_12", piTrackCorr.covariance(1,2));
          bs.addUserFloat("pi_cov_13", piTrackCorr.covariance(1,3));
          bs.addUserFloat("pi_cov_14", piTrackCorr.covariance(1,4));
          bs.addUserFloat("pi_cov_22", piTrackCorr.covariance(2,2));
          bs.addUserFloat("pi_cov_23", piTrackCorr.covariance(2,3));
          bs.addUserFloat("pi_cov_24", piTrackCorr.covariance(2,4));
          bs.addUserFloat("pi_cov_33", piTrackCorr.covariance(3,3));
          bs.addUserFloat("pi_cov_34", piTrackCorr.covariance(3,4));
          bs.addUserFloat("pi_cov_44", piTrackCorr.covariance(4,4));

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
          bs.addUserInt("found_photon", foundPhoton);
          bs.addUserInt("photon_multiplicity", photonMultiplicity);

          //add prefit resonances --> are not part of collection and can thus not be
          //automatically access the pt(), .. as I can for k1,k2,pi,mu in the variables_cff.py

          bs.addUserFloat("kk_pt"    , kkVec.pt()  );
          bs.addUserFloat("kk_eta"   , kkVec.eta() );
          bs.addUserFloat("kk_phi"   , kkVec.phi() );
          bs.addUserFloat("kk_m"     , kkVec.mass());
          bs.addUserFloat("kk_charge", kkCharge    );
          bs.addUserFloat("kk_deltaR", reco::deltaR(k1Vec, k2Vec));

          bs.addUserFloat("phiPi_pt"    , phiPiVec.pt()  );
          bs.addUserFloat("phiPi_eta"   , phiPiVec.eta() );
          bs.addUserFloat("phiPi_phi"   , phiPiVec.phi() );
          bs.addUserFloat("phiPi_m"     , phiPiVec.mass());
          bs.addUserFloat("phiPi_charge", phiPiCharge    );
          bs.addUserFloat("phiPi_deltaR", reco::deltaR(kkVec, piVec));

          bs.addUserFloat("dsMu_pt"    , dsMuVec.pt()  );
          bs.addUserFloat("dsMu_eta"   , dsMuVec.eta() );
          bs.addUserFloat("dsMu_phi"   , dsMuVec.phi() );
          bs.addUserFloat("dsMu_m"     , dsMuVec.mass()); //we miss the neutrino contribution
          bs.addUserFloat("dsMu_charge", dsMuCharge    );
          bs.addUserFloat("dsMu_deltaR", reco::deltaR(phiPiVec, muVec));

          //rel charges (not additive!)
          bs.addUserFloat("pi_mu_charge",piMuCharge); 

          bs.addUserInt("matchedMu",matchedMu); 
          bs.addUserInt("matchedK1",matchedK1); 
          bs.addUserInt("matchedK2",matchedK2); 
          bs.addUserInt("matchedPi",matchedPi); 

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


          const float phiErr = std::sqrt(std::max(0., phiParticle->currentState().kinematicParametersError().matrix()(6,6)));
          const float dsErr  = std::sqrt(std::max(0., dsParticle->currentState() .kinematicParametersError().matrix()(6,6))); 
          const float bsErr  = std::sqrt(std::max(0., bsParticle->currentState() .kinematicParametersError().matrix()(6,6)));
          
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
          bs.addUserFloat("pv_general_x",       pv_general_x); 
          bs.addUserFloat("pv_general_y",       pv_general_y); 
          bs.addUserFloat("pv_general_z",       pv_general_z); 

          bs.addUserFloat("pv_x",       pv_x); 
          bs.addUserFloat("pv_y",       pv_y); 
          bs.addUserFloat("pv_z",       pv_z); 

          bs.addUserFloat("pv_x_w_bs",       pv_x_w_bs); 
          bs.addUserFloat("pv_y_w_bs",       pv_y_w_bs); 

          bs.addUserFloat("pv_chi2",    pv.chi2()); // can access this directlyfor the reco::Vertex class
          bs.addUserFloat("pv_ndof",    pv.ndof());
          bs.addUserFloat("pv_redchi2", pv.normalizedChi2());
          bs.addUserFloat("pv_prob",    ChiSquaredProbability(pv.chi2(), pv.ndof()));
          bs.addUserInt  ("pv_idx",     pvIdx);

          // save also the x,y components without beamspot calculation 
          bs.addUserFloat("pv_x_wout_bs",  pv_general.x()); 
          bs.addUserFloat("pv_y_wout_bs",  pv_general.y()); 

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
          bs.addUserFloat("ds_vtx_cosine_xy"    , ds_vtx_cosine_xy);
          bs.addUserFloat("ds_vtx_cosine_xy_pv" , ds_vtx_cosine_xy_pv);

          bs.addUserFloat("ds_vtx_cosine_xyz"   , ds_vtx_cosine_xyz);
          bs.addUserFloat("ds_vtx_cosine_xyz_pv", ds_vtx_cosine_xyz_pv);


          ///////////////////////////////////////////////////////////////////////////////////////////
          // store  Lxy, Dxy                                                                       //
          ///////////////////////////////////////////////////////////////////////////////////////////

          // lxy(z) is the flight distance in the xy(z) plane(space)

          Measurement1D lxyDsVect       = VertexDistanceXY().distance(recoSv, recoTv ); //dsVtx->vertexState());
          Measurement1D lxyzDsVect      = VertexDistance3D().distance(recoSv, recoTv );

          Measurement1D lxyPhiVect      = VertexDistanceXY().distance(recoTv, recoFv );        
          Measurement1D lxyzPhiVect     = VertexDistance3D().distance(recoTv, recoFv );        

          Measurement1D lxyBsVect       = VertexDistanceXY().distance( pv, recoSv    );
          Measurement1D lxyzBsVect      = VertexDistance3D().distance( pv, recoSv    ); 

          bs.addUserFloat("lxy_bs",       lxyBsVect.value());
          bs.addUserFloat("lxy_bs_err",   lxyBsVect.error());
          bs.addUserFloat("lxy_bs_sig",   lxyBsVect.significance());

          bs.addUserFloat("lxyz_bs",      lxyzBsVect.value());
          bs.addUserFloat("lxyz_bs_err",  lxyzBsVect.error());
          bs.addUserFloat("lxyz_bs_sig",  lxyzBsVect.significance());

          bs.addUserFloat("lxy_ds",       lxyDsVect.value());
          bs.addUserFloat("lxy_ds_err",   lxyDsVect.error());
          bs.addUserFloat("lxy_ds_sig",   lxyDsVect.significance());

          bs.addUserFloat("lxyz_ds",      lxyzDsVect.value());
          bs.addUserFloat("lxyz_ds_err",  lxyzDsVect.error());
          bs.addUserFloat("lxyz_ds_sig",  lxyzDsVect.significance());

          bs.addUserFloat("lxy_phi",      lxyPhiVect.value());
          bs.addUserFloat("lxy_phi_err",  lxyPhiVect.error());
          bs.addUserFloat("lxy_phi_sig",  lxyPhiVect.significance());

          bs.addUserFloat("lxyz_phi",     lxyzPhiVect.value());
          bs.addUserFloat("lxyz_phi_err", lxyzPhiVect.error());
          bs.addUserFloat("lxyz_phi_sig", lxyzPhiVect.significance());

          // dxy(z) is the impact parameter in the xy(z) plane(space), i.e. the distance to the PV

          auto direction_ds        =  GlobalVector(phiPiVec.px(),  phiPiVec.py(),  phiPiVec.pz()  ).unit();
          auto direction_bs        =  GlobalVector(sv_x - pv_x, sv_y - pv_y, sv_z - pv_z ).unit();

          // First use signedDecayLength3D
          std::pair<bool, Measurement1D> signed_decay_ip3DMu_ds_sv_pair  = IPTools::signedDecayLength3D(ttMu, direction_ds, recoSv);
          float signed_decay_ip3DMu_ds_sv       = signed_decay_ip3DMu_ds_sv_pair.second.value();
          float signed_decay_ip3DMu_ds_err_sv   = signed_decay_ip3DMu_ds_sv_pair.second.error();
          float signed_decay_ip3DMu_ds_sig_sv   = signed_decay_ip3DMu_ds_sv_pair.second.significance();

          std::pair<bool, Measurement1D> signed_decay_ip3DMu_bs_sv_pair  = IPTools::signedDecayLength3D(ttMu, direction_bs, recoSv);
          float signed_decay_ip3DMu_bs_sv       = signed_decay_ip3DMu_bs_sv_pair.second.value();
          float signed_decay_ip3DMu_bs_err_sv   = signed_decay_ip3DMu_bs_sv_pair.second.error();
          float signed_decay_ip3DMu_bs_sig_sv   = signed_decay_ip3DMu_bs_sv_pair.second.significance();

          // Then use signedImpactParameter3D
          std::pair<bool, Measurement1D> signed_impact_ip3DMu_ds_sv_pair  = IPTools::signedImpactParameter3D(ttMu, direction_ds, recoSv);
          float signed_impact_ip3DMu_ds_sv       = signed_impact_ip3DMu_ds_sv_pair.second.value();
          float signed_impact_ip3DMu_ds_err_sv   = signed_impact_ip3DMu_ds_sv_pair.second.error();
          float signed_impact_ip3DMu_ds_sig_sv   = signed_impact_ip3DMu_ds_sv_pair.second.significance();

          std::pair<bool, Measurement1D> signed_impact_ip3DMu_bs_sv_pair  = IPTools::signedImpactParameter3D(ttMu, direction_bs, recoSv);
          float signed_impact_ip3DMu_bs_sv       = signed_impact_ip3DMu_bs_sv_pair.second.value();
          float signed_impact_ip3DMu_bs_err_sv   = signed_impact_ip3DMu_bs_sv_pair.second.error();
          float signed_impact_ip3DMu_bs_sig_sv   = signed_impact_ip3DMu_bs_sv_pair.second.significance();

          // use btv 
          SignedImpactParameter3D signedObject;

          std::pair<bool, Measurement1D> btv_ip3DMu_ds_sv_pair = signedObject.apply(ttMu, direction_ds, recoSv);
          float btv_ip3DMu_ds_sv       = btv_ip3DMu_ds_sv_pair.second.value();
          float btv_ip3DMu_ds_err_sv   = btv_ip3DMu_ds_sv_pair.second.error();
          float btv_ip3DMu_ds_sig_sv   = btv_ip3DMu_ds_sv_pair.second.significance();

          std::pair<bool, Measurement1D> btv_ip3DMu_bs_sv_pair = signedObject.apply(ttMu, direction_bs, recoSv);
          float btv_ip3DMu_bs_sv       = btv_ip3DMu_bs_sv_pair.second.value();
          float btv_ip3DMu_bs_err_sv   = btv_ip3DMu_bs_sv_pair.second.error();
          float btv_ip3DMu_bs_sig_sv   = btv_ip3DMu_bs_sv_pair.second.significance();


 
          float dxyMu_pv           = muTrackCorr.dxy(pv.position());  
          float dxyMuErr_pv        = muTrackCorr.dxyError(pv.position(),pv.covariance());  
          float dxyMuSig_pv        = dxyMu_pv/dxyMuErr_pv;

          float dxyMu_pv_dz        = muTrackCorr.dxy(pv_dz.position());  
          float dxyMuErr_pv_dz     = muTrackCorr.dxyError(pv_dz.position(),pv_dz.covariance());  
          float dxyMuSig_pv_dz     = dxyMu_pv_dz/dxyMuErr_pv_dz;

          float dxyMu_pv_dz_w_bs        = muTrackCorr.dxy(pv_dz_w_bs.position());  
          float dxyMuErr_pv_dz_w_bs     = muTrackCorr.dxyError(pv_dz_w_bs.position(),pv_dz_w_bs.covariance());  
          float dxyMuSig_pv_dz_w_bs     = dxyMu_pv_dz_w_bs/dxyMuErr_pv_dz_w_bs;

          //std::cout << "pv general position x is: " << pv_general.position().x() << std::endl;
          //std::cout << "pv general _x is: : " << pv_general_x << std::endl;

          float dxyMu_pv_general        = muTrackCorr.dxy(pv_general.position());  
          float dxyMuErr_pv_general     = muTrackCorr.dxyError(pv_general.position(),pv_general.covariance());  
          float dxyMuSig_pv_general     = dxyMu_pv_general/dxyMuErr_pv_general;


          reco::Vertex::CovarianceMatrix bsCovFixed;
          
          auto bsCov = beamSpot.covariance();
          
          bsCovFixed(0,0) = bsCov(0,0);
          bsCovFixed(0,1) = bsCov(0,1);
          bsCovFixed(0,2) = bsCov(0,2);
          bsCovFixed(1,0) = bsCov(1,0);
          bsCovFixed(1,1) = bsCov(1,1);
          bsCovFixed(1,2) = bsCov(1,2);
          bsCovFixed(2,0) = bsCov(2,0);
          bsCovFixed(2,1) = bsCov(2,1);
          bsCovFixed(2,2) = bsCov(2,2);


          float dxyMu_beamSpot        = muTrackCorr.dxy(beamSpot.position());  
          float dxyMuErr_beamSpot     = muTrackCorr.dxyError(beamSpot.position(),bsCovFixed);  
          float dxyMuSig_beamSpot     = dxyMu_beamSpot/dxyMuErr_beamSpot;


          float dzMu_pv            = muTrackCorr.dz(pv.position());  
          float dzMuErr_pv         = muTrackCorr.dzError();  
          float dzMuSig_pv         = dzMu_pv/dzMuErr_pv ; 

          float dxyMu_sv           = muTrackCorr.dxy(recoSv.position());  
          float dxyMuErr_sv        = muTrackCorr.dxyError(recoSv.position(),recoSv.covariance());  
          float dxyMuSig_sv        = dxyMu_sv/dxyMuErr_sv;

          float dzMu_sv            = muTrackCorr.dz(recoSv.position());  
          float dzMuErr_sv         = muTrackCorr.dzError();  
          float dzMuSig_sv         = dzMu_sv/dzMuErr_sv ; 

          float dxyPi_pv           = piTrackCorr.dxy(pv.position());  //maybe useful for Ds* vs Ds ? 
          float dxyPiErr_pv        = piTrackCorr.dxyError(pv.position(),pv.covariance());  
          float dxyPiSig_pv        = dxyPi_pv/dxyPiErr_pv;

          float dzPi_pv            = piTrackCorr.dz(pv.position());  
          float dzPiErr_pv         = piTrackCorr.dzError();  
          float dzPiSig_pv         = dzPi_pv/dzPiErr_pv ; 

          float dxyPi_sv           = piTrackCorr.dxy(recoSv.position());  //maybe useful for Ds* vs Ds ? 
          float dxyPiErr_sv        = piTrackCorr.dxyError(recoSv.position(),recoSv.covariance());  
          float dxyPiSig_sv        = dxyPi_sv/dxyPiErr_sv;

          float dzPi_sv            = piTrackCorr.dz(recoSv.position());  
          float dzPiErr_sv         = piTrackCorr.dzError();  
          float dzPiSig_sv         = dzPi_sv/dzPiErr_sv ; 

          bs.addUserFloat("signed_decay_ip3d_mu_ds_sv",     signed_decay_ip3DMu_ds_sv    );
          bs.addUserFloat("signed_decay_ip3d_mu_ds_err_sv", signed_decay_ip3DMu_ds_err_sv);
          bs.addUserFloat("signed_decay_ip3d_mu_ds_sig_sv", signed_decay_ip3DMu_ds_sig_sv);

          bs.addUserFloat("signed_decay_ip3d_mu_bs_sv",     signed_decay_ip3DMu_bs_sv    );
          bs.addUserFloat("signed_decay_ip3d_mu_bs_err_sv", signed_decay_ip3DMu_bs_err_sv);
          bs.addUserFloat("signed_decay_ip3d_mu_bs_sig_sv", signed_decay_ip3DMu_bs_sig_sv);

          bs.addUserFloat("signed_impact_ip3d_mu_ds_sv",     signed_impact_ip3DMu_ds_sv    );
          bs.addUserFloat("signed_impact_ip3d_mu_ds_err_sv", signed_impact_ip3DMu_ds_err_sv);
          bs.addUserFloat("signed_impact_ip3d_mu_ds_sig_sv", signed_impact_ip3DMu_ds_sig_sv);

          bs.addUserFloat("signed_impact_ip3d_mu_bs_sv",     signed_impact_ip3DMu_bs_sv    );
          bs.addUserFloat("signed_impact_ip3d_mu_bs_err_sv", signed_impact_ip3DMu_bs_err_sv);
          bs.addUserFloat("signed_impact_ip3d_mu_bs_sig_sv", signed_impact_ip3DMu_bs_sig_sv);

          bs.addUserFloat("btv_ip3d_mu_ds_sv",     btv_ip3DMu_ds_sv    );
          bs.addUserFloat("btv_ip3d_mu_ds_err_sv", btv_ip3DMu_ds_err_sv);
          bs.addUserFloat("btv_ip3d_mu_ds_sig_sv", btv_ip3DMu_ds_sig_sv);

          bs.addUserFloat("btv_ip3d_mu_bs_sv",     btv_ip3DMu_bs_sv    );
          bs.addUserFloat("btv_ip3d_mu_bs_err_sv", btv_ip3DMu_bs_err_sv);
          bs.addUserFloat("btv_ip3d_mu_bs_sig_sv", btv_ip3DMu_bs_sig_sv);



          bs.addUserFloat("dxy_mu_sv",     dxyMu_sv);
          bs.addUserFloat("dxy_mu_err_sv", dxyMuErr_sv);
          bs.addUserFloat("dxy_mu_sig_sv", dxyMuSig_sv);

          bs.addUserFloat("dz_mu_sv",      dzMu_sv);
          bs.addUserFloat("dz_mu_err_sv",  dzMuErr_sv);
          bs.addUserFloat("dz_mu_sig_sv",  dzMuSig_sv);

          bs.addUserFloat("dxy_mu_pv",     dxyMu_pv);
          bs.addUserFloat("dxy_mu_err_pv", dxyMuErr_pv);
          bs.addUserFloat("dxy_mu_sig_pv", dxyMuSig_pv);

          bs.addUserFloat("dxy_mu_beamSpot",     dxyMu_beamSpot);
          bs.addUserFloat("dxy_mu_err_beamSpot", dxyMuErr_beamSpot);
          bs.addUserFloat("dxy_mu_sig_beamSpot", dxyMuSig_beamSpot);

          bs.addUserFloat("dxy_mu_pv_dz",     dxyMu_pv_dz);
          bs.addUserFloat("dxy_mu_err_pv_dz", dxyMuErr_pv_dz);
          bs.addUserFloat("dxy_mu_sig_pv_dz", dxyMuSig_pv_dz);

          bs.addUserFloat("dxy_mu_pv_dz_w_bs",     dxyMu_pv_dz_w_bs);
          bs.addUserFloat("dxy_mu_err_pv_dz_w_bs", dxyMuErr_pv_dz_w_bs);
          bs.addUserFloat("dxy_mu_sig_pv_dz_w_bs", dxyMuSig_pv_dz_w_bs);


          bs.addUserFloat("dxy_mu_pv_general",     dxyMu_pv_general);
          bs.addUserFloat("dxy_mu_err_pv_general", dxyMuErr_pv_general);
          bs.addUserFloat("dxy_mu_sig_pv_general", dxyMuSig_pv_general);


          bs.addUserFloat("dz_mu_pv",      dzMu_pv);
          bs.addUserFloat("dz_mu_err_pv",  dzMuErr_pv);
          bs.addUserFloat("dz_mu_sig_pv",  dzMuSig_pv);

          bs.addUserFloat("dxy_pi_pv",     dxyPi_pv);
          bs.addUserFloat("dxy_pi_err_pv", dxyPiErr_pv);
          bs.addUserFloat("dxy_pi_sig_pv", dxyPiSig_pv);

          bs.addUserFloat("dz_pi_pv",      dzPi_pv);
          bs.addUserFloat("dz_pi_err_pv",  dzPiErr_pv);
          bs.addUserFloat("dz_pi_sig_pv",  dzPiSig_pv);


          bs.addUserFloat("dxy_pi_sv",     dxyPi_sv);
          bs.addUserFloat("dxy_pi_err_sv", dxyPiErr_sv);
          bs.addUserFloat("dxy_pi_sig_sv", dxyPiSig_sv);

          bs.addUserFloat("dz_pi_sv",      dzPi_sv);
          bs.addUserFloat("dz_pi_err_sv",  dzPiErr_sv);
          bs.addUserFloat("dz_pi_sig_sv",  dzPiSig_sv);


          ////////////////////////////////////////// deltaR (defined above) /////////////////////////////////////

          bs.addUserFloat("muonK1dR",muonK1dR);
          bs.addUserFloat("muonK2dR",muonK2dR);
          bs.addUserFloat("muonPidR",muonPidR);

          //////////////////////////perp Ds+Mu and /corrected Bs mass ///////////////////////////////////////////
          //accoring to https://journals.aps.org/prd/pdf/10.1103/PhysRevD.100.112006 eq 3.

          TVector3 bFlightDir;
          bFlightDir.SetXYZ(pv_x - sv_x, pv_y - sv_y, pv_z - sv_z);

          double dsPerp      = phiPiVec.Vect().Cross(bFlightDir).R() / bFlightDir.Mag(); //phiPiVec.Vect().Perp(bFlightDir);
          double dsMuPerp    = dsMuVec .Vect().Cross(bFlightDir).R() / bFlightDir.Mag(); //dsMuVec .Vect().Perp(bFlightDir);
          double bsMassCorr  = std::sqrt(std::pow(dsMuVec.M(),2) + std::pow(dsMuPerp,2)) + dsMuPerp;

          bs.addUserFloat("bs_m_corr",  bsMassCorr);
          bs.addUserFloat("ds_perp",    dsPerp);
          bs.addUserFloat("ds_mu_perp", dsMuPerp);


          ///////////////////////// reconstruct the B momentum with different methods ///////////////////////////////

          /////////////////////////
          // Collinear approx.   //
          ///////////////////////// 

          math::PtEtaPhiMLorentzVector collBsTlv = collMethod(dsMuVec, bsMass_);

          math::PtEtaPhiMLorentzVector collMissTlv = collBsTlv - dsMuVec ; // bs - ds+mu
          math::PtEtaPhiMLorentzVector collQTlv    = collBsTlv - phiPiVec; // bs - ds

          double m2_miss_coll = collMissTlv.M2();
          double pt_miss_coll = collMissTlv.Pt();
          double q2_coll      = collQTlv.M2();
          float  e_star_coll  = getEStar(collBsTlv, muVec); 

          bs.addUserFloat("bs_px_coll" ,collBsTlv.Px());
          bs.addUserFloat("bs_py_coll" ,collBsTlv.Py());
          bs.addUserFloat("bs_pz_coll" ,collBsTlv.Pz());
          bs.addUserFloat("bs_pt_coll" ,collBsTlv.Pt());
          bs.addUserFloat("bs_eta_coll",collBsTlv.Eta());
          bs.addUserFloat("bs_phi_coll",collBsTlv.Phi());

          bs.addUserFloat("m2_miss_coll",m2_miss_coll);
          bs.addUserFloat("pt_miss_coll",pt_miss_coll);
          bs.addUserFloat("q2_coll"     ,q2_coll);
          bs.addUserFloat("e_star_coll" ,e_star_coll);

          bs.addUserFloat("b_boost_coll"    ,-collBsTlv.BoostToCM().R());
          bs.addUserFloat("b_boost_coll_x"  ,-collBsTlv.BoostToCM().X());
          bs.addUserFloat("b_boost_coll_y"  ,-collBsTlv.BoostToCM().Y());
          bs.addUserFloat("b_boost_coll_z"  ,-collBsTlv.BoostToCM().Z());

          //////////////////////////
          // LHCb method          //
          //////////////////////////

          math::PtEtaPhiMLorentzVector lhcbBsTlv   = lhcbMethod(dsMuVec, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_);      

          math::PtEtaPhiMLorentzVector lhcbMissTlv = lhcbBsTlv - dsMuVec ;
          math::PtEtaPhiMLorentzVector lhcbQTlv    = lhcbBsTlv - phiPiVec;

          float m2_miss_lhcb = lhcbMissTlv.M2();
          float pt_miss_lhcb = lhcbMissTlv.Pt();
          float q2_lhcb      = lhcbQTlv.M2();
          float e_star_lhcb  = getEStar(lhcbBsTlv, muVec); 

          bs.addUserFloat("bs_px_lhcb" ,lhcbBsTlv.Px());
          bs.addUserFloat("bs_py_lhcb" ,lhcbBsTlv.Py());
          bs.addUserFloat("bs_pz_lhcb" ,lhcbBsTlv.Pz());
          bs.addUserFloat("bs_pt_lhcb" ,lhcbBsTlv.Pt());
          bs.addUserFloat("bs_eta_lhcb",lhcbBsTlv.Eta());
          bs.addUserFloat("bs_phi_lhcb",lhcbBsTlv.Phi());

          bs.addUserFloat("m2_miss_lhcb",m2_miss_lhcb);
          bs.addUserFloat("pt_miss_lhcb",pt_miss_lhcb);
          bs.addUserFloat("q2_lhcb"     ,q2_lhcb     );
          bs.addUserFloat("e_star_lhcb" ,e_star_lhcb );

          bs.addUserFloat("b_boost_lhcb"    ,-lhcbBsTlv.BoostToCM().R());
          bs.addUserFloat("b_boost_lhcb_x"  ,-lhcbBsTlv.BoostToCM().X());
          bs.addUserFloat("b_boost_lhcb_y"  ,-lhcbBsTlv.BoostToCM().Y());
          bs.addUserFloat("b_boost_lhcb_z"  ,-lhcbBsTlv.BoostToCM().Z());


          ///////////////////////////////////////////
          // Another Lhcb method                   //
          ///////////////////////////////////////////

          math::PtEtaPhiMLorentzVector lhcbAltBsTlv   = lhcbAltMethod(dsMuVec, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_);      

          math::PtEtaPhiMLorentzVector lhcbAltMissTlv = lhcbAltBsTlv - dsMuVec ; // bs - ds+mu
          math::PtEtaPhiMLorentzVector lhcbAltQTlv    = lhcbAltBsTlv - phiPiVec; // bs - ds 

          float e_star_lhcb_alt  = getEStar(lhcbAltBsTlv, muVec); 

          bs.addUserFloat("bs_px_lhcb_alt", lhcbAltBsTlv.Px());
          bs.addUserFloat("bs_py_lhcb_alt", lhcbAltBsTlv.Py());
          bs.addUserFloat("bs_pz_lhcb_alt", lhcbAltBsTlv.Pz());
          bs.addUserFloat("bs_pt_lhcb_alt", lhcbAltBsTlv.Pt());
          bs.addUserFloat("bs_eta_lhcb_alt",lhcbAltBsTlv.Eta());
          bs.addUserFloat("bs_phi_lhcb_alt",lhcbAltBsTlv.Phi());

          bs.addUserFloat("m2_miss_lhcb_alt",lhcbAltMissTlv.M2());
          bs.addUserFloat("pt_miss_lhcb_alt",lhcbAltMissTlv.Pt());
          bs.addUserFloat("q2_lhcb_alt"     ,lhcbAltQTlv.M2())   ;
          bs.addUserFloat("e_star_lhcb_alt" ,e_star_lhcb_alt)    ;

          bs.addUserFloat("b_boost_lhcb_alt"    , -lhcbAltBsTlv.BoostToCM().R());
          bs.addUserFloat("b_boost_lhcb_alt_x"  , -lhcbAltBsTlv.BoostToCM().X());
          bs.addUserFloat("b_boost_lhcb_alt_y"  , -lhcbAltBsTlv.BoostToCM().Y());
          bs.addUserFloat("b_boost_lhcb_alt_z"  , -lhcbAltBsTlv.BoostToCM().Z());

          ///////////////////////////////////////////////////////////
          // RECO METHOD -> 2 Solutions                            //
          // This method is exact in the 1-neutrino final state :) //
          ///////////////////////////////////////////////////////////
          //std::cout << "now i calc the reco method:" << std::endl;
          std::tuple<std::vector<math::PtEtaPhiMLorentzVector>,float> recoResult = recoMethod(dsMuVec, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_); 

          std::vector<math::PtEtaPhiMLorentzVector> recosBs = std::get<0>(recoResult);
          float discNegativity                              = std::get<1>(recoResult);

          int discIsNegative = 0;
          if (discNegativity > 0) discIsNegative = 1;

          math::PtEtaPhiMLorentzVector recoBsTlv1 = recosBs.at(0);
          math::PtEtaPhiMLorentzVector recoBsTlv2 = recosBs.at(1);

          math::PtEtaPhiMLorentzVector miss_1 = recoBsTlv1 - dsMuVec;  
          math::PtEtaPhiMLorentzVector miss_2 = recoBsTlv2 - dsMuVec;  
          math::PtEtaPhiMLorentzVector Q_1    = recoBsTlv1 - phiPiVec;  
          math::PtEtaPhiMLorentzVector Q_2    = recoBsTlv2 - phiPiVec;  

          float m2_miss_reco_1 = miss_1.M2();
          float m2_miss_reco_2 = miss_2.M2();
          float pt_miss_reco_1 = miss_1.Pt();
          float pt_miss_reco_2 = miss_2.Pt();

          float q2_reco_1 = Q_1.M2();
          float q2_reco_2 = Q_2.M2();
          float e_star_reco_1  = getEStar(recoBsTlv1, muVec); 
          float e_star_reco_2  = getEStar(recoBsTlv2, muVec); 

          bs.addUserInt("disc_is_negative" ,discIsNegative);
          bs.addUserFloat("disc_negativity",discNegativity);

          bs.addUserFloat("bs_px_reco_1" ,recoBsTlv1.Px());
          bs.addUserFloat("bs_py_reco_1" ,recoBsTlv1.Py());
          bs.addUserFloat("bs_pz_reco_1" ,recoBsTlv1.Pz());
          bs.addUserFloat("bs_pt_reco_1" ,recoBsTlv1.Pt()); 
          bs.addUserFloat("bs_eta_reco_1",recoBsTlv1.Eta()); 
          bs.addUserFloat("bs_phi_reco_1",recoBsTlv1.Phi()); 

          bs.addUserFloat("bs_px_reco_2" ,recoBsTlv2.Px());
          bs.addUserFloat("bs_py_reco_2" ,recoBsTlv2.Py());
          bs.addUserFloat("bs_pz_reco_2" ,recoBsTlv2.Pz());
          bs.addUserFloat("bs_pt_reco_2" ,recoBsTlv2.Pt()); 
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

          bs.addUserFloat("b_boost_reco_1"    ,-recoBsTlv1.BoostToCM().R());
          bs.addUserFloat("b_boost_reco_1_x"  ,-recoBsTlv1.BoostToCM().X());
          bs.addUserFloat("b_boost_reco_1_y"  ,-recoBsTlv1.BoostToCM().Y());
          bs.addUserFloat("b_boost_reco_1_z"  ,-recoBsTlv1.BoostToCM().Z());


          bs.addUserFloat("b_boost_reco_2"    ,-recoBsTlv2.BoostToCM().R());
          bs.addUserFloat("b_boost_reco_2_x"  ,-recoBsTlv2.BoostToCM().X());
          bs.addUserFloat("b_boost_reco_2_y"  ,-recoBsTlv2.BoostToCM().Y());
          bs.addUserFloat("b_boost_reco_2_z"  ,-recoBsTlv2.BoostToCM().Z());


          ////////////////// Save momenta (and masses for consistency, even if constrained /////////////////////////

          // Phi fit
          bs.addUserFloat("phi_fitted_px",    fittedPhi.Px()); 
          bs.addUserFloat("phi_fitted_py",    fittedPhi.Py()); 
          bs.addUserFloat("phi_fitted_pz",    fittedPhi.Pz()); 
          bs.addUserFloat("phi_fitted_pt",    fittedPhi.Pt()); 
          bs.addUserFloat("phi_fitted_eta",   fittedPhi.Eta()); 
          bs.addUserFloat("phi_fitted_phi",   fittedPhi.Phi()); 
          bs.addUserFloat("phi_fitted_m",     fittedPhi.M()); 
          bs.addUserFloat("phi_fitted_m_err", phiErr); 

          // bs fit
          bs.addUserFloat("ds_fitted_px",     fittedDs.Px()); 
          bs.addUserFloat("ds_fitted_py",     fittedDs.Py()); 
          bs.addUserFloat("ds_fitted_pz",     fittedDs.Pz()); 
          bs.addUserFloat("ds_fitted_pt",     fittedDs.Pt()); 
          bs.addUserFloat("ds_fitted_eta",    fittedDs.Eta()); 
          bs.addUserFloat("ds_fitted_phi",    fittedDs.Phi()); 
          bs.addUserFloat("ds_fitted_m",      fittedDs.M()); 
          bs.addUserFloat("ds_fitted_m_err",  dsErr); 
          bs.addUserFloat("ds_custom_m_err",  dsErrCustom); 

          // bs fit
          bs.addUserFloat("bs_fitted_px",     fittedBs.Px()); 
          bs.addUserFloat("bs_fitted_py",     fittedBs.Py()); 
          bs.addUserFloat("bs_fitted_pz",     fittedBs.Pz()); 
          bs.addUserFloat("bs_fitted_pt",     fittedBs.Pt()); 
          bs.addUserFloat("bs_fitted_eta",    fittedBs.Eta()); 
          bs.addUserFloat("bs_fitted_phi",    fittedBs.Phi()); 
          bs.addUserFloat("bs_fitted_m",      fittedBs.M()); 
          bs.addUserFloat("bs_fitted_m_err",  bsErr); 

          /////////////////////// helicity angles in all possibe variations /////////////////////////

          //angle between Mu and W
          float angMuWColl    = angMuW(phiPiVec, collBsTlv,   muVec);
          float angMuWLhcb    = angMuW(phiPiVec, lhcbBsTlv,   muVec);
          float angMuWLhcbAlt = angMuW(phiPiVec, lhcbAltBsTlv,muVec);
          float angMuWReco1   = angMuW(phiPiVec, recoBsTlv1,  muVec);
          float angMuWReco2   = angMuW(phiPiVec, recoBsTlv2,  muVec);

          bs.addUserFloat("cosMuWColl",       cos(angMuWColl)); 
          bs.addUserFloat("cosMuWLhcb",       cos(angMuWLhcb)); 
          bs.addUserFloat("cosMuWLhcbAlt",    cos(angMuWLhcbAlt)); 
          bs.addUserFloat("cosMuWReco1",      cos(angMuWReco1)); 
          bs.addUserFloat("cosMuWReco2",      cos(angMuWReco2)); 


          //angle between k1(k2) and pion in phi rest frame
          float angPiK1  = angDoubleDecay(kkVec, k1Vec, piVec);
          float angPiK2  = angDoubleDecay(kkVec, k2Vec, piVec);
          // equivalently, angle between phi(pi) and bs in ds rest frame
          float angPhiDsColl    = angDsPi( phiPiVec, kkVec,  collBsTlv);
          float angPhiDsLhcb    = angDsPi( phiPiVec, kkVec,  lhcbBsTlv);
          float angPhiDsLhcbAlt = angDsPi( phiPiVec, kkVec,  lhcbAltBsTlv);
          float angPhiDsReco1   = angDsPi( phiPiVec, kkVec,  recoBsTlv1);
          float angPhiDsReco2   = angDsPi( phiPiVec, kkVec,  recoBsTlv2);

          float angPiDsColl    = angDsPi ( phiPiVec, piVec,  collBsTlv);
          float angPiDsLhcb    = angDsPi ( phiPiVec, piVec,  lhcbBsTlv);
          float angPiDsLhcbAlt = angDsPi ( phiPiVec, piVec,  lhcbAltBsTlv);
          float angPiDsReco1   = angDsPi ( phiPiVec, piVec,  recoBsTlv1);
          float angPiDsReco2   = angDsPi ( phiPiVec, piVec,  recoBsTlv2);
          

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
          float angPlaneBsColl    = angPlane( phiPiVec, collBsTlv,    muVec, piVec);
          float angPlaneBsLhcb    = angPlane( phiPiVec, lhcbBsTlv,    muVec, piVec);
          float angPlaneBsLhcbAlt = angPlane( phiPiVec, lhcbAltBsTlv, muVec, piVec);
          float angPlaneBsReco1   = angPlane( phiPiVec, recoBsTlv1,   muVec, piVec);
          float angPlaneBsReco2   = angPlane( phiPiVec, recoBsTlv2,   muVec, piVec);

          bs.addUserFloat("cosPlaneBsColl",    cos(angPlaneBsColl));
          bs.addUserFloat("cosPlaneBsLhcb",    cos(angPlaneBsLhcb));
          bs.addUserFloat("cosPlaneBsLhcbAlt", cos(angPlaneBsLhcbAlt));
          bs.addUserFloat("cosPlaneBsReco1",   cos(angPlaneBsReco1));
          bs.addUserFloat("cosPlaneBsReco2",   cos(angPlaneBsReco2));

          float angPlaneDsColl    = angPlane2(phiPiVec, collBsTlv,    k1Vec, piVec);
          float angPlaneDsLhcb    = angPlane2(phiPiVec, lhcbBsTlv,    k1Vec, piVec);
          float angPlaneDsLhcbAlt = angPlane2(phiPiVec, lhcbAltBsTlv, k1Vec, piVec);
          float angPlaneDsReco1   = angPlane2(phiPiVec, recoBsTlv1,   k1Vec, piVec);
          float angPlaneDsReco2   = angPlane2(phiPiVec, recoBsTlv2,   k1Vec, piVec);

          bs.addUserFloat("cosPlaneDsColl",    cos(angPlaneDsColl));
          bs.addUserFloat("cosPlaneDsLhcb",    cos(angPlaneDsLhcb));
          bs.addUserFloat("cosPlaneDsLhcbAlt", cos(angPlaneDsLhcbAlt));
          bs.addUserFloat("cosPlaneDsReco1",   cos(angPlaneDsReco1));
          bs.addUserFloat("cosPlaneDsReco2",   cos(angPlaneDsReco2));

          ///////////////////////////////////////////////////////////////////////////////////////////
          // Calculate muon ISOLATION                                                              //
          //                                                                                       //
          // idea: charged part: take all charged hadron pts and subtract kkpi cand                //
          //       neutral part: take neutral hadronand photon ET and subtract contributions       //
          //                     from PU vertices (assumed to be about ~50%, see paper)            //
          ///////////////////////////////////////////////////////////////////////////////////////////

          // STANDARD: Using PF-built-In Isolation function

          // pfIsolationR0{X} takes tracks associated with 1st vtx in offlineSlimmedPrimaryVertices
          auto muIso03 = muPtr->pfIsolationR03(); //dR = 0.3
          auto muIso04 = muPtr->pfIsolationR04(); //dR = 0.4
          
          float iso03KKPi = 0.0;
          float iso04KKPi = 0.0;

          int idxMinVtxDz_pi = -1;
          int idxMinVtxDz_k1 = -1;
          int idxMinVtxDz_k2 = -1;

          float minVtxDz_pi = 0;
          float minVtxDz_k1 = 0;
          float minVtxDz_k2 = 0;

          // for standard isolation, check if vertex associated to pf candidates is
          // the first in the offlineSlimmedPrimaryVertexCollection
          for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){

            edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, vtxIdx);


            // and minimize the dz           
            float vtxDz_pi = fabs( piPtr->vz() - vtxPtr->z() );  
            float vtxDz_k1 = fabs( k1Ptr->vz() - vtxPtr->z() );  
            float vtxDz_k2 = fabs( k2Ptr->vz() - vtxPtr->z() );  

            if ( ( vtxDz_pi < minVtxDz_pi) || (idxMinVtxDz_pi < 0) ) {

              minVtxDz_pi    = vtxDz_pi;
              idxMinVtxDz_pi = vtxIdx;         

            } 
            if ( ( vtxDz_k1 < minVtxDz_k1) || (idxMinVtxDz_k1 < 0) ) {

              minVtxDz_k1    = vtxDz_k1;
              idxMinVtxDz_k1 = vtxIdx;         

            } 
            if ( ( vtxDz_k2 < minVtxDz_k2) || (idxMinVtxDz_k2 < 0) ) {

              minVtxDz_k2    = vtxDz_k2;
              idxMinVtxDz_k2 = vtxIdx;         

            } 

          }

          float pu03KKPi = 0; 
          float pu04KKPi = 0;

          // add KKPi to isolation or to pileup, depending if they come from the first vertex or not.
 
          if     ( (reco::deltaR(*k1Ptr,*muPtr) < 0.3) && (idxMinVtxDz_k1 == 0 ))  iso03KKPi += k1Vec.Pt();
          else if( (reco::deltaR(*k1Ptr,*muPtr) < 0.3) && (idxMinVtxDz_k1 != 0 ))  pu03KKPi  += k1Vec.Pt();

          if     ( (reco::deltaR(*k2Ptr,*muPtr) < 0.3) && (idxMinVtxDz_k2 == 0 ))  iso03KKPi += k2Vec.Pt();
          else if( (reco::deltaR(*k2Ptr,*muPtr) < 0.3) && (idxMinVtxDz_k2 != 0 ))  pu03KKPi  += k2Vec.Pt();

          if     ( (reco::deltaR(*piPtr,*muPtr) < 0.3) && (idxMinVtxDz_pi == 0 ))  iso03KKPi += piVec.Pt();
          else if( (reco::deltaR(*piPtr,*muPtr) < 0.3) && (idxMinVtxDz_pi != 0 ))  pu03KKPi  += piVec.Pt();

          if     ( (reco::deltaR(*k1Ptr,*muPtr) < 0.4) && (idxMinVtxDz_k1 == 0 ))  iso04KKPi += k1Vec.Pt();
          else if( (reco::deltaR(*k1Ptr,*muPtr) < 0.4) && (idxMinVtxDz_k1 != 0 ))  pu04KKPi  += k1Vec.Pt();

          if     ( (reco::deltaR(*k2Ptr,*muPtr) < 0.4) && (idxMinVtxDz_k2 == 0 ))  iso04KKPi += k2Vec.Pt();
          else if( (reco::deltaR(*k2Ptr,*muPtr) < 0.4) && (idxMinVtxDz_k2 != 0 ))  pu04KKPi  += k2Vec.Pt();

          if     ( (reco::deltaR(*piPtr,*muPtr) < 0.4) && (idxMinVtxDz_pi == 0 ))  iso04KKPi += piVec.Pt();
          else if( (reco::deltaR(*piPtr,*muPtr) < 0.4) && (idxMinVtxDz_pi != 0 ))  pu04KKPi  += piVec.Pt();


          // factor 0.5 from https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=MUO-16-001&tp=an&id=1779&ancode=MUO-16-001
          //float zero = 0.0; // avoid problems with std::max function

          float iso03charged    = std::max(0.0f, muIso03.sumChargedHadronPt - iso03KKPi);  
          float iso04charged    = std::max(0.0f, muIso04.sumChargedHadronPt - iso04KKPi);  

          float iso03pileup     = muIso03.sumPUPt - pu03KKPi;
          float iso04pileup     = muIso04.sumPUPt - pu04KKPi;

          float iso03neutral    = std::max(0.0, muIso03.sumNeutralHadronEt + muIso03.sumPhotonEt - 0.5 * iso03pileup);  
          float iso04neutral    = std::max(0.0, muIso04.sumNeutralHadronEt + muIso04.sumPhotonEt - 0.5 * iso04pileup);  

          float iso03 = iso03charged + iso03neutral;
          float iso04 = iso04charged + iso04neutral;

          float relIso03 = iso03 / muVec.Pt();
          float relIso04 = iso04 / muVec.Pt();

          bs.addUserFloat("iso_03", iso03);
          bs.addUserFloat("iso_04", iso04);
          bs.addUserFloat("rel_iso_03", relIso03);
          bs.addUserFloat("rel_iso_04", relIso04);

          ///////////////////////////////////////////////////////////////////////////////////////////////////////
          // ADVANCED: Calculate isolation using tracks corresponding to our custom vertices, not  the default //
          // ( default is the 1st vtx in offlineSlimmedPrimaryVertices, the latter corresponds to highest      //
          // sumPT, not the use-case for B-Physics :)                                                          //
          //                                                                                                   // 
          // Implemented for our custom PV, SV and TV                                                          //
          ///////////////////////////////////////////////////////////////////////////////////////////////////////

          // holds charged tracks 
          float iso03_ch_pv = 0.; 
          float iso04_ch_pv = 0.; 
          float iso03_ch_sv = 0.; 
          float iso04_ch_sv = 0.; 
          float iso03_ch_tv = 0.; 
          float iso04_ch_tv = 0.; 

          // holds neutral photon tracks
          float iso03_np_pv = 0.; 
          float iso04_np_pv = 0.; 
          float iso03_np_sv = 0.; 
          float iso04_np_sv = 0.; 
          float iso03_np_tv = 0.; 
          float iso04_np_tv = 0.; 

          // holds neutral hadron tracks
          float iso03_nh_pv = 0.; 
          float iso04_nh_pv = 0.; 
          float iso03_nh_sv = 0.; 
          float iso04_nh_sv = 0.; 
          float iso03_nh_tv = 0.; 
          float iso04_nh_tv = 0.; 

          // holds pile up tracks
          float iso03_pu_pv = 0.; 
          float iso04_pu_pv = 0.; 
          float iso03_pu_sv = 0.; 
          float iso04_pu_sv = 0.; 
          float iso03_pu_tv = 0.; 
          float iso04_pu_tv = 0.; 

          // Fill charged isolations
          for (size_t k = 0; k < chargedPfIdx.size(); ++k) {

            const size_t idx = chargedPfIdx[k];
            // skip kkpi signal candidate 
            if ((idx == piIdx) || (idx == k1Idx) || (idx == k2Idx)) continue;

            const auto& chPtr = chargedPfPtr[k];
            const auto& tt    = chargedTT[k];

            //////////////////////////////////////////////////////////////
            // Correct the momentum scale and smear of the PF candidate //
            //////////////////////////////////////////////////////////////

            //charged PF candidates have the pion mass per default
            math::PtEtaPhiMLorentzVector pfVec(chPtr->pt(), chPtr->eta(), chPtr->phi(), chPtr->mass());


            if (!iEvent.isRealData()) {

            const auto& id    = iEvent.id();
            const double z_PF = trkCalib->gauss(id.run(), id.luminosityBlock(), static_cast<int64_t>(id.event()), idx);
            const double f_PF = trkCalib->factor(chPtr->p(), chPtr->pt(), chPtr->eta(), chPtr->charge(), z_PF);

            pfVec = math::PtEtaPhiMLorentzVector(pfVec.Pt() * f_PF, pfVec.Eta(), pfVec.Phi(), pfVec.M());

            }


            ////////////////////////////////////////////////////////
            // For each charged particle, get the assigned vertex // 
            ////////////////////////////////////////////////////////

            float minVtxDz    = 0.0;
            int   idxMinVtxDz = -1;
              
            for(size_t vtxIdx = 0; vtxIdx < refitVtx.size(); ++vtxIdx){

              // loop over all refit vertices 
              const reco::Vertex& vtx = refitVtx[vtxIdx];

              // and minimize the dz           
              float vtxDz = fabs( chPtr->vz() - vtx.z() );  
              if ( ( vtxDz < minVtxDz) || (idxMinVtxDz < 0) ) {
 
                minVtxDz    = vtxDz;
                idxMinVtxDz = vtxIdx;         

              } 

            }
           
            ///////////////////////////////////////
            // Fill for pv as reference vertex   //
            ///////////////////////////////////////

            float deltaR_mu_chTrk = reco::deltaR(*muPtr, *chPtr);

            // if closest refit vtx is our custom pv, count it to isolation!
            // remove electrons and muons (11 and 13) as in cmssw isolation definition
            if  ((idxMinVtxDz == pvIdx) && (abs(chPtr->pdgId()) != 11) && (abs(chPtr->pdgId()) != 13)) {
              iso04_ch_pv += pfVec.Pt();
              if (deltaR_mu_chTrk < 0.3) iso03_ch_pv += pfVec.Pt();
            } 

            // displaced track, count it to pile up
            else if ((idxMinVtxDz != pvIdx) && ( deltaR_mu_chTrk > 0.01 ) && (pfVec.Pt() > 0.5)) {
              iso04_pu_pv += pfVec.Pt();
              if (deltaR_mu_chTrk < 0.3) iso03_pu_pv += pfVec.Pt();
            }

            ///////////////////////////////////////
            // Fill for sv as reference vertex   //
            ///////////////////////////////////////

            std::pair<bool, Measurement1D> ip3DSv = IPTools::signedDecayLength3D(tt, direction_bs, recoSv);

            if ( ip3DSv.first && (std::fabs(ip3DSv.second.value()) < 0.01) && (abs(chPtr->pdgId()) != 11) && (abs(chPtr->pdgId()) != 13) ){
              iso04_ch_sv += pfVec.Pt();
              if (deltaR_mu_chTrk < 0.3) iso03_ch_sv += pfVec.Pt();
            }

            // displaced track, count it to pile up (why no deltaR_mu_chTrk > 0.01 here ? )
            else if (( ip3DSv.first && (std::fabs(ip3DSv.second.value()) >= 0.01)) && (pfVec.Pt() > 0.5)) {
              iso04_pu_sv += pfVec.Pt();
              if (deltaR_mu_chTrk < 0.3) iso03_pu_sv += pfVec.Pt();
            }

            ///////////////////////////////////////
            // Fill for tv as reference vertex   // 
            ///////////////////////////////////////

            std::pair<bool, Measurement1D> ip3DTv = IPTools::signedDecayLength3D(tt, direction_ds, recoTv);

            if ( ip3DTv.first && (std::fabs(ip3DTv.second.value()) < 0.01) && (abs(chPtr->pdgId()) != 11) && (abs(chPtr->pdgId()) != 13) ){
              iso04_ch_tv += pfVec.Pt();
              if (deltaR_mu_chTrk < 0.3) iso03_ch_tv += pfVec.Pt();
            }

            // displaced track, count it to pile up (why no deltaR_mu_chTrk > 0.01 here ? )
            else if ( ip3DTv.first && (std::fabs(ip3DTv.second.value()) >= 0.01) && (pfVec.Pt() > 0.5)) {
              iso04_pu_tv += pfVec.Pt();
              if (deltaR_mu_chTrk < 0.3) iso03_pu_tv += pfVec.Pt();
            }


          } 


          // Fill neutral isolations
          for (const auto& nPtr: neutralPfPtr) {

            if  (nPtr -> pdgId() == 22) iso04_np_pv += nPtr->pt();
            else                        iso04_nh_pv += nPtr->pt(); 

            float deltaR_mu_nTrk = reco::deltaR(*muPtr, *nPtr);

            if      ((deltaR_mu_nTrk < 0.3) && (nPtr->pdgId() == 22)) iso03_np_pv += nPtr->pt();
            else if ((deltaR_mu_nTrk < 0.3) && (nPtr->pdgId() != 22)) iso03_nh_pv += nPtr->pt();

          } 

          iso03_np_sv = iso03_np_pv;
          iso04_np_sv = iso04_np_pv;
          iso03_nh_sv = iso03_nh_pv;
          iso04_nh_sv = iso04_nh_pv;

          iso03_np_tv = iso03_np_pv;
          iso04_np_tv = iso04_np_pv;
          iso03_nh_tv = iso03_nh_pv;
          iso04_nh_tv = iso04_nh_pv;


          // 0.5 is an empirical factor. 50 % of the charged pile up can be mapepd to neutral pile up
          // This factor also compensates to have flat efficiency when having higher pile up
          // It is subtracted from all the neutral traacks (which can not be associated with a vertex and thus
          // we pick up also pile up neutrals when summing iso03(4)_nh_refit_pv, iso03(4)_np_refit_pv
          float iso03_pv = std::max(0.0, iso03_ch_pv + iso03_nh_pv + iso03_np_pv - 0.5*iso03_pu_pv);
          float iso04_pv = std::max(0.0, iso04_ch_pv + iso04_nh_pv + iso04_np_pv - 0.5*iso04_pu_pv);

          float iso03_sv = std::max(0.0, iso03_ch_sv + iso03_nh_sv + iso03_np_sv - 0.5*iso03_pu_sv);
          float iso04_sv = std::max(0.0, iso04_ch_sv + iso04_nh_sv + iso04_np_sv - 0.5*iso04_pu_sv);

          float iso03_tv = std::max(0.0, iso03_ch_tv + iso03_nh_tv + iso03_np_tv - 0.5*iso03_pu_tv);
          float iso04_tv = std::max(0.0, iso04_ch_tv + iso04_nh_tv + iso04_np_tv - 0.5*iso04_pu_tv);

          bs.addUserFloat("iso_03_pv", iso03_pv);
          bs.addUserFloat("iso_04_pv", iso04_pv);
          bs.addUserFloat("iso_03_sv", iso03_sv);
          bs.addUserFloat("iso_04_sv", iso04_sv);
          bs.addUserFloat("iso_03_tv", iso03_tv);
          bs.addUserFloat("iso_04_tv", iso04_tv);

          bs.addUserFloat("rel_iso_03_pv", iso03_pv / muVec.Pt());
          bs.addUserFloat("rel_iso_04_pv", iso04_pv / muVec.Pt());
          bs.addUserFloat("rel_iso_03_sv", iso03_sv / muVec.Pt());
          bs.addUserFloat("rel_iso_04_sv", iso04_sv / muVec.Pt());
          bs.addUserFloat("rel_iso_03_tv", iso03_tv / muVec.Pt());
          bs.addUserFloat("rel_iso_04_tv", iso04_tv / muVec.Pt());

          ///////////////////////////////////////////////////////////////////////////////////////////
          // Calculate Ds ISOLATION                                                                //
          //                                                                                       //
          // idea: follow same procedure as for muon isolation. Thus, the only difference is that  //
          //       we take canidates in a deltaR cone (0.3 or 0.4) around the ds (kkpi) candidate  //
          //                                                                                       //
          //       charged part: take all charged hadron pts and subtract kkpi cand                //
          //       neutral part: take neutral hadron and photon ET and subtract contributions      //
          //                     from PU vertices (assumed to be about ~50%, see paper)            //
          ///////////////////////////////////////////////////////////////////////////////////////////


          pat::CompositeCandidate ds;
          math::PtEtaPhiMLorentzVector dsP4(phiPiVec.Pt(), phiPiVec.Eta(), phiPiVec.Phi(), dsMass_);
          ds.setP4(dsP4);
          // select particles in cone 
      
          std::vector<int> chargedPf_Ds_Idx;
          std::vector<int> neutralPf_Ds_Idx;
          std::vector<edm::Ptr<pat::PackedCandidate>>  chargedPf_Ds_Ptr; 
          std::vector<edm::Ptr<pat::PackedCandidate>>  neutralPf_Ds_Ptr; 
          std::vector<reco::TransientTrack> chargedTT_Ds;
 
          for(size_t pfIdx = 0; pfIdx < pcand->size() + tracksLost->size() ; ++pfIdx) {
      
            //define a pointer to the pf cand at position pfIdx
            edm::Ptr<pat::PackedCandidate> pfPtr;
            if (pfIdx < pcand->size()) pfPtr = edm::Ptr<pat::PackedCandidate>(pcand, pfIdx); //normal tracks
            else pfPtr = edm::Ptr<pat::PackedCandidate>(tracksLost, pfIdx - pcand->size()); //lost tracks

            float deltaRDsTrk = reco::deltaR(*pfPtr, phiPiVec);
            if      ( neutralPdgId.count(abs(pfPtr->pdgId())) && (deltaRDsTrk < 0.4) && (deltaRDsTrk > 0.01  ) && (pfPtr->pt() > 0.5 )) {
              neutralPf_Ds_Idx.push_back(pfIdx); 
              neutralPf_Ds_Ptr.push_back(pfPtr); 
            }
            else if ( chargedPdgId.count(abs(pfPtr->pdgId())) && (deltaRDsTrk < 0.4) && (deltaRDsTrk > 0.0001)                        ) {

              if ((abs(pfPtr->pdgId()) == 11) || (abs(pfPtr->pdgId()) == 13)) continue;
              if (!pfPtr->hasTrackDetails())  continue;

              chargedPf_Ds_Idx.push_back(pfIdx); 
              chargedPf_Ds_Ptr.push_back(pfPtr); 
              chargedTT_Ds    .push_back(ttBuilder->build(pfPtr->bestTrack()));
            }    
          }

          // holds charged tracks 
          float iso03_ds_ch_sv = 0.; 
          float iso04_ds_ch_sv = 0.; 

          // holds neutral photon tracks
          float iso03_ds_np_sv = 0.; 
          float iso04_ds_np_sv = 0.; 

          // holds neutral hadron tracks
          float iso03_ds_nh_sv = 0.; 
          float iso04_ds_nh_sv = 0.; 

          // holds pile up tracks
          float iso03_ds_pu_sv = 0.; 
          float iso04_ds_pu_sv = 0.; 

          // Fill charged isolations
          for (size_t k = 0; k < chargedPf_Ds_Idx.size(); ++k) {

            // skip kkpi signal candidate 
            const size_t idx = chargedPf_Ds_Idx[k];
            if ((idx == piIdx) || (idx == k1Idx) || (idx == k2Idx)) continue;

            const auto& chPtr = chargedPf_Ds_Ptr[k];
            const auto& tt    = chargedTT_Ds[k];

            //////////////////////////////////////////////////////////////
            // Correct the momentum scale and smear of the PF candidate //
            //////////////////////////////////////////////////////////////

            //charged PF candidates have the pion mass per default
            math::PtEtaPhiMLorentzVector pfVec(chPtr->pt(), chPtr->eta(), chPtr->phi(), chPtr->mass());

            if (!iEvent.isRealData()) {

            const auto& id    = iEvent.id();
            const double z_PF = trkCalib->gauss(id.run(), id.luminosityBlock(), static_cast<int64_t>(id.event()), idx);
            const double f_PF = trkCalib->factor(chPtr->p(), chPtr->pt(), chPtr->eta(), chPtr->charge(), z_PF);

            pfVec = math::PtEtaPhiMLorentzVector(pfVec.Pt() * f_PF, pfVec.Eta(), pfVec.Phi(), pfVec.M());

            }


            // for each charged particle, check if the assigned vertex is the one we pick up
            float minVtxDz    = 0.0;
            int   idxMinVtxDz = -1;
              
            for(size_t vtxIdx = 0; vtxIdx < refitVtx.size(); ++vtxIdx){

              // loop over all refit vertices 
              const reco::Vertex& vtx = refitVtx[vtxIdx];

              // and minimize the dz           
              float vtxDz = fabs( chPtr->vz() - vtx.z() );  
              if ( ( vtxDz < minVtxDz) || (idxMinVtxDz < 0) ) {
 
                minVtxDz    = vtxDz;
                idxMinVtxDz = vtxIdx;         

              } 

            }
           
            ///////////////////////////////////////
            // Fill for sv as reference vertex   //
            ///////////////////////////////////////

            float deltaR_ds_chTrk = reco::deltaR(phiPiVec, *chPtr);

            std::pair<bool, Measurement1D> ip3DSv = IPTools::signedDecayLength3D(tt, direction_bs, recoSv);

            if ( ip3DSv.first && (std::fabs(ip3DSv.second.value()) < 0.01) && (abs(chPtr->pdgId()) != 11) && (abs(chPtr->pdgId()) != 13) ){
              iso04_ds_ch_sv += pfVec.Pt();
              if (deltaR_ds_chTrk < 0.3) iso03_ds_ch_sv += pfVec.Pt();
            }

            // displaced track, count it to pile up (why no deltaR_ds_chTrk > 0.01 here ? )
            else if (( ip3DSv.first && (std::fabs(ip3DSv.second.value()) >= 0.01)) && (pfVec.Pt() > 0.5)) {
              iso04_ds_pu_sv += pfVec.Pt();
              if (deltaR_ds_chTrk < 0.3) iso03_ds_pu_sv += pfVec.Pt();
            }

          } 


          // Fill neutral isolations
          for (const auto& nPtr: neutralPf_Ds_Ptr) {

            if  (nPtr -> pdgId() == 22) iso04_ds_np_sv += nPtr->pt();
            else                        iso04_ds_nh_sv += nPtr->pt(); 

            //float deltaR_mu_nTrk = reco::deltaR(phiPi, *nPtr);
            float deltaR_mu_nTrk = reco::deltaR(ds, *nPtr);

            if      ((deltaR_mu_nTrk < 0.3) && (nPtr->pdgId() == 22)) iso03_ds_np_sv += nPtr->pt();
            else if ((deltaR_mu_nTrk < 0.3) && (nPtr->pdgId() != 22)) iso03_ds_nh_sv += nPtr->pt();

          } 

          float iso03_ds_sv = std::max(0.0, iso03_ds_ch_sv + iso03_ds_nh_sv + iso03_ds_np_sv - 0.5*iso03_ds_pu_sv);
          float iso04_ds_sv = std::max(0.0, iso04_ds_ch_sv + iso04_ds_nh_sv + iso04_ds_np_sv - 0.5*iso04_ds_pu_sv);

          bs.addUserFloat("iso_03_ds_sv", iso03_ds_sv);
          bs.addUserFloat("iso_04_ds_sv", iso04_ds_sv);

          bs.addUserFloat("rel_iso_03_ds_sv", iso03_ds_sv / phiPiVec.Pt());
          bs.addUserFloat("rel_iso_04_ds_sv", iso04_ds_sv / phiPiVec.Pt());

          // e gamma
          float e_gamma = getEGamma(phiPiVec, dsMass_,dsStarMass_);
          bs.addUserFloat("e_gamma", e_gamma);

          // add kappa
          float kappa = getKappa(phiPiVec, muVec);
          bs.addUserFloat("kappa", kappa);

          ////////////////////////////////////////////////
          // Redo all the methods if we have a photon!  //
          ////////////////////////////////////////////////
          double dsPerpPhoton                 = dsPerp            ;
          double dsMuPerpPhoton               = dsMuPerp          ;
          double bsMassCorrPhoton             = bsMassCorr        ;
                                                                        
          double m2_miss_coll_photon          = m2_miss_coll;
          double pt_miss_coll_photon          = pt_miss_coll;
          double q2_coll_photon               = q2_coll     ;
          float  e_star_coll_photon           = e_star_coll ;
                                                                        
          double bs_px_coll_photon            = collBsTlv.Px() ;
          double bs_py_coll_photon            = collBsTlv.Py() ;
          double bs_pz_coll_photon            = collBsTlv.Pz() ;
          double bs_pt_coll_photon            = collBsTlv.Pt() ;
          double bs_eta_coll_photon           = collBsTlv.Eta();
          double bs_phi_coll_photon           = collBsTlv.Phi();
                                                                       
          double b_boost_coll_photon          = -collBsTlv.BoostToCM().R();
          double b_boost_coll_x_photon        = -collBsTlv.BoostToCM().X();
          double b_boost_coll_y_photon        = -collBsTlv.BoostToCM().Y();
          double b_boost_coll_z_photon        = -collBsTlv.BoostToCM().Z();
                                                                       
          double m2_miss_lhcb_alt_photon      = lhcbAltMissTlv.M2();
          double pt_miss_lhcb_alt_photon      = lhcbAltMissTlv.Pt(); 
          double q2_lhcb_alt_photon           = lhcbAltQTlv.M2()   ; 
          double e_star_lhcb_alt_photon       = e_star_lhcb_alt    ;  

          double bs_px_lhcb_alt_photon        = lhcbAltBsTlv.Px() ;
          double bs_py_lhcb_alt_photon        = lhcbAltBsTlv.Py() ;
          double bs_pz_lhcb_alt_photon        = lhcbAltBsTlv.Pz() ;
          double bs_pt_lhcb_alt_photon        = lhcbAltBsTlv.Pt() ;
          double bs_eta_lhcb_alt_photon       = lhcbAltBsTlv.Eta();
          double bs_phi_lhcb_alt_photon       = lhcbAltBsTlv.Phi();

          double b_boost_lhcb_alt_photon      = -lhcbAltBsTlv.BoostToCM().R();
          double b_boost_lhcb_alt_x_photon    = -lhcbAltBsTlv.BoostToCM().X();
          double b_boost_lhcb_alt_y_photon    = -lhcbAltBsTlv.BoostToCM().Y();
          double b_boost_lhcb_alt_z_photon    = -lhcbAltBsTlv.BoostToCM().Z();

          double m2_miss_reco_1_photon        = m2_miss_reco_1        ;
          double pt_miss_reco_1_photon        = pt_miss_reco_1        ; 
          double q2_reco_1_photon             = q2_reco_1             ; 
          double e_star_reco_1_photon         = e_star_reco_1         ;  
                                                                             
          double bs_px_reco_1_photon          = recoBsTlv1.Px() ;
          double bs_py_reco_1_photon          = recoBsTlv1.Py() ;
          double bs_pz_reco_1_photon          = recoBsTlv1.Pz() ;
          double bs_pt_reco_1_photon          = recoBsTlv1.Pt() ;
          double bs_eta_reco_1_photon         = recoBsTlv1.Eta();
          double bs_phi_reco_1_photon         = recoBsTlv1.Phi();

          double b_boost_reco_1_photon        = -recoBsTlv1.BoostToCM().R();
          double b_boost_reco_1_x_photon      = -recoBsTlv1.BoostToCM().X();
          double b_boost_reco_1_y_photon      = -recoBsTlv1.BoostToCM().Y();
          double b_boost_reco_1_z_photon      = -recoBsTlv1.BoostToCM().Z();

          double m2_miss_reco_2_photon        = m2_miss_reco_2      ;
          double pt_miss_reco_2_photon        = pt_miss_reco_2      ; 
          double q2_reco_2_photon             = q2_reco_2           ; 
          double e_star_reco_2_photon         = e_star_reco_2       ;  
                                                                           
          double bs_px_reco_2_photon          = recoBsTlv2.Px() ;
          double bs_py_reco_2_photon          = recoBsTlv2.Py() ;
          double bs_pz_reco_2_photon          = recoBsTlv2.Pz() ;
          double bs_pt_reco_2_photon          = recoBsTlv2.Pt() ;
          double bs_eta_reco_2_photon         = recoBsTlv2.Eta();
          double bs_phi_reco_2_photon         = recoBsTlv2.Phi();
 
          double b_boost_reco_2_photon        = -recoBsTlv2.BoostToCM().R();
          double b_boost_reco_2_x_photon      = -recoBsTlv2.BoostToCM().X();
          double b_boost_reco_2_y_photon      = -recoBsTlv2.BoostToCM().Y();
          double b_boost_reco_2_z_photon      = -recoBsTlv2.BoostToCM().Z();

          double discNegativityPhoton         = discNegativity         ;
          double discIsNegativePhoton         = discIsNegative         ;

          //angle between Mu and W
          float cosMuWCollPhoton              = cos(angMuWColl)       ;
          float cosMuWLhcbAltPhoton           = cos(angMuWLhcbAlt)    ;
          float cosMuWReco1Photon             = cos(angMuWReco1)      ;
          float cosMuWReco2Photon             = cos(angMuWReco2)      ;


                                                                        
          float cosPhiDsCollPhoton            =  cos(angPhiDsColl)       ; 
          float cosPhiDsLhcbAltPhoton         =  cos(angPhiDsLhcbAlt)    ; 
          float cosPhiDsReco1Photon           =  cos(angPhiDsReco1)      ; 
          float cosPhiDsReco2Photon           =  cos(angPhiDsReco2)      ; 
                                                                        
          float cosPiDsCollPhoton             = cos(angPiDsColl)       ; 
          float cosPiDsLhcbAltPhoton          = cos(angPiDsLhcbAlt)    ; 
          float cosPiDsReco1Photon            = cos(angPiDsReco1)      ; 
          float cosPiDsReco2Photon            = cos(angPiDsReco2)      ; 

          // helicity plane angle
          float cosPlaneBsCollPhoton          = cos(angPlaneBsColl)    ;
          float cosPlaneBsLhcbAltPhoton       = cos(angPlaneBsLhcbAlt) ;
          float cosPlaneBsReco1Photon         = cos(angPlaneBsReco1)   ;
          float cosPlaneBsReco2Photon         = cos(angPlaneBsReco2)   ;
                                                                        
          float cosPlaneDsCollPhoton          = cos(angPlaneDsColl)    ;
          float cosPlaneDsLhcbAltPhoton       = cos(angPlaneDsLhcbAlt) ;
          float cosPlaneDsReco1Photon         = cos(angPlaneDsReco1)   ;
          float cosPlaneDsReco2Photon         = cos(angPlaneDsReco2)   ;
                                                                              
          float iso03Photon                   =  iso03                  ;
          float iso04Photon                   =  iso04                  ;
          float relIso03Photon                =  relIso03               ;
          float relIso04Photon                =  relIso04               ;
                                                                              
          float iso03PhotonPv                 =  iso03_pv;
          float iso04PhotonPv                 =  iso04_pv;
          float relIso03PhotonPv              =  iso03_pv / muVec.Pt()    ;
          float relIso04PhotonPv              =  iso04_pv / muVec.Pt()    ;
                                                                              
          float iso03PhotonSv                 =  iso03_sv;
          float iso04PhotonSv                 =  iso04_sv;
          float relIso03PhotonSv              =  iso03_sv / muVec.Pt()    ;
          float relIso04PhotonSv              =  iso04_sv / muVec.Pt()    ;
                                                 
          float iso03PhotonTv                 =  iso03_tv;
          float iso04PhotonTv                 =  iso04_tv;
          float relIso03PhotonTv              =  iso03_tv / muVec.Pt()    ;
          float relIso04PhotonTv              =  iso04_tv / muVec.Pt()    ;
                                                                              
          float iso03DsPhotonSv               = iso03_ds_sv ;
          float iso04DsPhotonSv               = iso04_ds_sv ;
          float relIso03DsPhotonSv            = iso03_ds_sv / phiPiVec.pt() ;
          float relIso04DsPhotonSv            = iso04_ds_sv / phiPiVec.pt() ;

          if (foundPhoton){

            math::PtEtaPhiMLorentzVector dsMuVecPhoton = dsMuVec + photonVec;
            math::PtEtaPhiMLorentzVector phiPiVecPhoton   = phiPiVec   + photonVec;

            /////////////////////////////
            // Ds + mu perp variable   //
            /////////////////////////////


            dsPerpPhoton   = phiPiVecPhoton.Vect().Cross(bFlightDir).R() / bFlightDir.Mag(); //phiPiVecPhoton.  Vect().Perp(bFlightDir);
            dsMuPerpPhoton = dsMuVecPhoton .Vect().Cross(bFlightDir).R() / bFlightDir.Mag();//dsMuVecPhoton.Vect().Perp(bFlightDir);
            bsMassCorrPhoton = std::sqrt(std::pow(dsMuVecPhoton.M(),2) + std::pow(dsMuPerpPhoton,2)) + dsMuPerpPhoton;
  
            /////////////////////////////
            // Collinear approximation //
            /////////////////////////////

  
            math::PtEtaPhiMLorentzVector collBsTlvPhoton = collMethod(dsMuVecPhoton, bsMass_);
  
            math::PtEtaPhiMLorentzVector collMissTlvPhoton; // for m2 miss
            math::PtEtaPhiMLorentzVector collQTlvPhoton;    // for q2 miss
  
            collMissTlvPhoton = collBsTlvPhoton - dsMuVecPhoton; // bs - ds+mu
            collQTlvPhoton    = collBsTlvPhoton - phiPiVecPhoton;   // bs - ds
  
            m2_miss_coll_photon         = collMissTlvPhoton.M2();
            pt_miss_coll_photon         = collMissTlvPhoton.Pt();
            q2_coll_photon              = collQTlvPhoton.M2();
            e_star_coll_photon          = getEStar(collBsTlvPhoton,muVec); 

            bs_px_coll_photon           = collBsTlvPhoton.Px(); 
            bs_py_coll_photon           = collBsTlvPhoton.Py(); 
            bs_pz_coll_photon           = collBsTlvPhoton.Pz(); 
            bs_pt_coll_photon           = collBsTlvPhoton.Pt(); 
            bs_eta_coll_photon          = collBsTlvPhoton.Eta(); 
            bs_phi_coll_photon          = collBsTlvPhoton.Phi(); 

            b_boost_coll_photon         = -collBsTlvPhoton.BoostToCM().R(); 
            b_boost_coll_x_photon       = -collBsTlvPhoton.BoostToCM().X(); 
            b_boost_coll_y_photon       = -collBsTlvPhoton.BoostToCM().Y(); 
            b_boost_coll_z_photon       = -collBsTlvPhoton.BoostToCM().Z(); 
 
            /////////////////////////////
            // Lhcb method             //
            /////////////////////////////

            math::PtEtaPhiMLorentzVector lhcbAltBsTlvPhoton   = lhcbAltMethod(dsMuVecPhoton, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_);      
  
            math::PtEtaPhiMLorentzVector lhcbAltMissTlvPhoton = lhcbAltBsTlvPhoton - dsMuVecPhoton; // bs - ds+mu
            math::PtEtaPhiMLorentzVector lhcbAltQTlvPhoton    = lhcbAltBsTlvPhoton - phiPiVecPhoton;   // bs - ds 
  
            m2_miss_lhcb_alt_photon      = lhcbAltMissTlvPhoton.M2();
            pt_miss_lhcb_alt_photon      = lhcbAltMissTlvPhoton.Pt();
            q2_lhcb_alt_photon           = lhcbAltQTlvPhoton.M2();
            e_star_lhcb_alt_photon       = getEStar(lhcbAltBsTlvPhoton,muVec); 

            bs_px_lhcb_alt_photon        = lhcbAltBsTlvPhoton.Px(); 
            bs_py_lhcb_alt_photon        = lhcbAltBsTlvPhoton.Py(); 
            bs_pz_lhcb_alt_photon        = lhcbAltBsTlvPhoton.Pz(); 
            bs_pt_lhcb_alt_photon        = lhcbAltBsTlvPhoton.Pt(); 
            bs_eta_lhcb_alt_photon       = lhcbAltBsTlvPhoton.Eta(); 
            bs_phi_lhcb_alt_photon       = lhcbAltBsTlvPhoton.Phi(); 

            b_boost_lhcb_alt_photon     = -lhcbAltBsTlvPhoton.BoostToCM().R(); 
            b_boost_lhcb_alt_x_photon   = -lhcbAltBsTlvPhoton.BoostToCM().X(); 
            b_boost_lhcb_alt_y_photon   = -lhcbAltBsTlvPhoton.BoostToCM().Y(); 
            b_boost_lhcb_alt_z_photon   = -lhcbAltBsTlvPhoton.BoostToCM().Z(); 

            /////////////////////////////
            // Reco method             //
            /////////////////////////////

            std::tuple<std::vector<math::PtEtaPhiMLorentzVector>,float> recoResultPhoton = recoMethod(dsMuVecPhoton, pv_x, pv_y, pv_z, sv_x, sv_y, sv_z, bsMass_); 
  
            std::vector<math::PtEtaPhiMLorentzVector> recosBsPhoton = std::get<0>(recoResultPhoton);
            discNegativityPhoton                      = std::get<1>(recoResultPhoton);
  
            discIsNegativePhoton = 0;
            if (discNegativityPhoton > 0) discIsNegativePhoton = 1;
  
            math::PtEtaPhiMLorentzVector recoBsTlv1Photon = recosBsPhoton.at(0);
            math::PtEtaPhiMLorentzVector recoBsTlv2Photon = recosBsPhoton.at(1);
  
            math::PtEtaPhiMLorentzVector miss_1Photon     = recoBsTlv1Photon - dsMuVecPhoton;  
            math::PtEtaPhiMLorentzVector miss_2Photon     = recoBsTlv2Photon - dsMuVecPhoton;  
            math::PtEtaPhiMLorentzVector Q_1Photon        = recoBsTlv1Photon - phiPiVecPhoton;  
            math::PtEtaPhiMLorentzVector Q_2Photon        = recoBsTlv2Photon - phiPiVecPhoton;  
  
            m2_miss_reco_1_photon           = miss_1Photon.M2();
            pt_miss_reco_1_photon           = miss_1Photon.Pt();
            q2_reco_1_photon                = Q_1Photon.M2();
            e_star_reco_1_photon            = getEStar(recoBsTlv1Photon,muVec); 
 
            bs_px_reco_1_photon             = recoBsTlv1Photon.Px(); 
            bs_py_reco_1_photon             = recoBsTlv1Photon.Py(); 
            bs_pz_reco_1_photon             = recoBsTlv1Photon.Pz(); 
            bs_pt_reco_1_photon             = recoBsTlv1Photon.Pt(); 
            bs_eta_reco_1_photon            = recoBsTlv1Photon.Eta(); 
            bs_phi_reco_1_photon            = recoBsTlv1Photon.Phi(); 

            b_boost_reco_1_photon           = -recoBsTlv1Photon.BoostToCM().R(); 
            b_boost_reco_1_x_photon         = -recoBsTlv1Photon.BoostToCM().X(); 
            b_boost_reco_1_y_photon         = -recoBsTlv1Photon.BoostToCM().Y(); 
            b_boost_reco_1_z_photon         = -recoBsTlv1Photon.BoostToCM().Z(); 

            m2_miss_reco_2_photon           = miss_2Photon.M2();
            pt_miss_reco_2_photon           = miss_2Photon.Pt();
            q2_reco_2_photon                = Q_2Photon.M2();
            e_star_reco_2_photon            = getEStar(recoBsTlv2Photon,muVec); 

            bs_px_reco_2_photon             = recoBsTlv2Photon.Px(); 
            bs_py_reco_2_photon             = recoBsTlv2Photon.Py(); 
            bs_pz_reco_2_photon             = recoBsTlv2Photon.Pz(); 
            bs_pt_reco_2_photon             = recoBsTlv2Photon.Pt(); 
            bs_eta_reco_2_photon            = recoBsTlv2Photon.Eta(); 
            bs_phi_reco_2_photon            = recoBsTlv2Photon.Phi(); 

            b_boost_reco_2_photon          = -recoBsTlv2Photon.BoostToCM().R(); 
            b_boost_reco_2_x_photon        = -recoBsTlv2Photon.BoostToCM().X(); 
            b_boost_reco_2_y_photon        = -recoBsTlv2Photon.BoostToCM().Y(); 
            b_boost_reco_2_z_photon        = -recoBsTlv2Photon.BoostToCM().Z(); 

            /////////////////////////////
            // Hel angles              // 
            /////////////////////////////

            double angMuWCollPhoton         = angMuW(phiPiVecPhoton, collBsTlvPhoton,   muVec);
            double angMuWLhcbAltPhoton      = angMuW(phiPiVecPhoton, lhcbAltBsTlvPhoton,muVec);
            double angMuWReco1Photon        = angMuW(phiPiVecPhoton, recoBsTlv1Photon,  muVec);
            double angMuWReco2Photon        = angMuW(phiPiVecPhoton, recoBsTlv2Photon,  muVec);

            // equivalently, angle between phi(pi) and bs in ds rest frame
            double angPhiDsCollPhoton       = angDsPi(phiPiVecPhoton, kkVec,  collBsTlvPhoton);
            double angPhiDsLhcbAltPhoton    = angDsPi(phiPiVecPhoton, kkVec,  lhcbAltBsTlvPhoton);
            double angPhiDsReco1Photon      = angDsPi(phiPiVecPhoton, kkVec,  recoBsTlv1Photon);
            double angPhiDsReco2Photon      = angDsPi(phiPiVecPhoton, kkVec,  recoBsTlv2Photon);
  
            double angPiDsCollPhoton        = angDsPi(phiPiVecPhoton, piVec,  collBsTlvPhoton);
            double angPiDsLhcbAltPhoton     = angDsPi(phiPiVecPhoton, piVec,  lhcbAltBsTlvPhoton);
            double angPiDsReco1Photon       = angDsPi(phiPiVecPhoton, piVec,  recoBsTlv1Photon);
            double angPiDsReco2Photon       = angDsPi(phiPiVecPhoton, piVec,  recoBsTlv2Photon);


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
            float angPlaneBsCollPhoton    = angPlane(phiPiVecPhoton, collBsTlvPhoton,    muVec, piVec);
            float angPlaneBsLhcbAltPhoton = angPlane(phiPiVecPhoton, lhcbAltBsTlvPhoton, muVec, piVec);
            float angPlaneBsReco1Photon   = angPlane(phiPiVecPhoton, recoBsTlv1Photon,   muVec, piVec);
            float angPlaneBsReco2Photon   = angPlane(phiPiVecPhoton, recoBsTlv2Photon,   muVec, piVec);
  
            float angPlaneDsCollPhoton    = angPlane2(phiPiVecPhoton, collBsTlvPhoton,    k1Vec, piVec);
            float angPlaneDsLhcbAltPhoton = angPlane2(phiPiVecPhoton, lhcbAltBsTlvPhoton, k1Vec, piVec);
            float angPlaneDsReco1Photon   = angPlane2(phiPiVecPhoton, recoBsTlv1Photon,   k1Vec, piVec);
            float angPlaneDsReco2Photon   = angPlane2(phiPiVecPhoton, recoBsTlv2Photon,   k1Vec, piVec);
 
            cosPlaneBsCollPhoton    =      cos( angPlaneBsCollPhoton   );  
            cosPlaneBsLhcbAltPhoton =      cos( angPlaneBsLhcbAltPhoton);
            cosPlaneBsReco1Photon   =      cos( angPlaneBsReco1Photon  );
            cosPlaneBsReco2Photon   =      cos( angPlaneBsReco2Photon  );

            cosPlaneDsCollPhoton    =      cos( angPlaneDsCollPhoton   );
            cosPlaneDsLhcbAltPhoton =      cos( angPlaneDsLhcbAltPhoton);
            cosPlaneDsReco1Photon   =      cos( angPlaneDsReco1Photon  );
            cosPlaneDsReco2Photon   =      cos( angPlaneDsReco2Photon  );
            
            /////////////////////////////
            // Isolations              //
            /////////////////////////////

            //first the classic isolation. Check if photon is in cone 0.3/0.4 around the mu
            float deltaR_ph_mu = reco::deltaR(*phtPtr, *muPtr  );
            float deltaR_ph_ds = reco::deltaR(*phtPtr, phiPiVec);

            iso03Photon   = iso03   ; 
            iso03PhotonPv = iso03_pv; 
            iso03PhotonSv = iso03_sv; 
            iso03PhotonTv = iso03_tv; 

            iso04Photon   = iso04   ; 
            iso04PhotonPv = iso04_pv; 
            iso04PhotonSv = iso04_sv; 
            iso04PhotonTv = iso04_tv; 

            iso03DsPhotonSv = iso03_ds_sv;
            iso04DsPhotonSv = iso04_ds_sv;


            if (deltaR_ph_mu < 0.3) { 

              iso03Photon   -= phtPtr->pt();
              iso03PhotonPv -= phtPtr->pt();
              iso03PhotonSv -= phtPtr->pt();
              iso03PhotonTv -= phtPtr->pt();

              iso04Photon   -= phtPtr->pt();
              iso04PhotonPv -= phtPtr->pt();
              iso04PhotonSv -= phtPtr->pt();
              iso04PhotonTv -= phtPtr->pt();

            }
            else if (deltaR_ph_mu < 0.4) {

              iso04Photon   -= phtPtr->pt();
              iso04PhotonPv -= phtPtr->pt();
              iso04PhotonSv -= phtPtr->pt();
              iso04PhotonTv -= phtPtr->pt();

            } 
 
            if (deltaR_ph_ds < 0.3) { 

              iso03DsPhotonSv -= phtPtr->pt();
              iso04DsPhotonSv -= phtPtr->pt();

            }

            else if (deltaR_ph_ds < 0.4) {

              iso04DsPhotonSv -= phtPtr->pt();

            }

            relIso03Photon     = iso03Photon      / muVec.Pt(); 
            relIso03PhotonPv   = iso03PhotonPv    / muVec.Pt(); 
            relIso03PhotonSv   = iso03PhotonSv    / muVec.Pt(); 
            relIso03PhotonTv   = iso03PhotonTv    / muVec.Pt(); 

            relIso04Photon     = iso04Photon      / muVec.Pt(); 
            relIso04PhotonPv   = iso04PhotonPv    / muVec.Pt(); 
            relIso04PhotonSv   = iso04PhotonSv    / muVec.Pt(); 
            relIso04PhotonTv   = iso04PhotonTv    / muVec.Pt(); 

            relIso03DsPhotonSv = iso03DsPhotonSv  / phiPiVec.Pt(); 
            relIso04DsPhotonSv = iso04DsPhotonSv  / phiPiVec.Pt(); 
 
  
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
          bs.addUserFloat("b_boost_coll_x_photon",   b_boost_coll_x_photon);
          bs.addUserFloat("b_boost_coll_y_photon",  b_boost_coll_y_photon);
          bs.addUserFloat("b_boost_coll_z_photon",  b_boost_coll_z_photon);

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
          bs.addUserFloat("b_boost_lhcb_alt_x_photon",   b_boost_lhcb_alt_x_photon);
          bs.addUserFloat("b_boost_lhcb_alt_y_photon",  b_boost_lhcb_alt_y_photon);
          bs.addUserFloat("b_boost_lhcb_alt_z_photon",  b_boost_lhcb_alt_z_photon);

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
          bs.addUserFloat("b_boost_reco_1_x_photon",  b_boost_reco_1_x_photon);
          bs.addUserFloat("b_boost_reco_1_y_photon", b_boost_reco_1_y_photon);
          bs.addUserFloat("b_boost_reco_1_z_photon", b_boost_reco_1_z_photon);

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
          bs.addUserFloat("b_boost_reco_2_x_photon",  b_boost_reco_2_x_photon);
          bs.addUserFloat("b_boost_reco_2_y_photon", b_boost_reco_2_y_photon);
          bs.addUserFloat("b_boost_reco_2_z_photon", b_boost_reco_2_z_photon);

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

          bs.addUserFloat("iso_03_photon"                ,  iso03Photon                 ) ;
          bs.addUserFloat("iso_04_photon"                ,  iso04Photon                 ) ;
          bs.addUserFloat("rel_iso_03_photon"            ,  relIso03Photon              ) ;
          bs.addUserFloat("rel_iso_04_photon"            ,  relIso04Photon              ) ;

          bs.addUserFloat("iso_03_pv_photon"             ,  iso03PhotonPv            ) ;
          bs.addUserFloat("iso_04_pv_photon"             ,  iso04PhotonPv            ) ;
          bs.addUserFloat("rel_iso_03_pv_photon"         ,  relIso03PhotonPv         ) ;
          bs.addUserFloat("rel_iso_04_pv_photon"         ,  relIso04PhotonPv         ) ;

          bs.addUserFloat("iso_03_sv_photon"             ,  iso03PhotonSv               ) ;
          bs.addUserFloat("iso_04_sv_photon"             ,  iso04PhotonSv               ) ;
          bs.addUserFloat("rel_iso_03_sv_photon"         ,  relIso03PhotonSv            ) ;
          bs.addUserFloat("rel_iso_04_sv_photon"         ,  relIso04PhotonSv            ) ;

          bs.addUserFloat("iso_03_tv_photon"             ,  iso03PhotonTv               ) ;
          bs.addUserFloat("iso_04_tv_photon"             ,  iso04PhotonTv               ) ;
          bs.addUserFloat("rel_iso_03_tv_photon"         ,  relIso03PhotonTv            ) ;
          bs.addUserFloat("rel_iso_04_tv_photon"         ,  relIso04PhotonTv            ) ;

          bs.addUserFloat("iso_03_ds_sv_photon"             ,  iso03DsPhotonSv             ) ;
          bs.addUserFloat("iso_04_ds_sv_photon"             ,  iso04DsPhotonSv             ) ;
          bs.addUserFloat("rel_iso_03_ds_sv_photon"         ,  relIso03DsPhotonSv          ) ;
          bs.addUserFloat("rel_iso_04_ds_sv_photon"         ,  relIso04DsPhotonSv          ) ;

  	  bs.addUserInt("muIdx",                   muPtr->userInt("muIdx"          ));
   	  bs.addUserInt("trgObjIdx",               muPtr->userInt("trgObjIdx"      ));

          bs.addUserInt("mu7_ip4",                 muPtr->userInt("mu7_ip4"       ));
          bs.addUserInt("prescale_mu7_ip4",        muPtr->userInt("prescale_mu7_ip4"   ));

          bs.addUserInt("mu8_ip3",                 muPtr->userInt("mu8_ip3"       ));
          bs.addUserInt("prescale_mu8_ip3",        muPtr->userInt("prescale_mu8_ip3"   ));

          bs.addUserInt("mu8p5_ip3p5",             muPtr->userInt("mu8p5_ip3p5"       ));
          bs.addUserInt("prescale_mu8p5_ip3p5",    muPtr->userInt("prescale_mu8p5_ip3p5"   ));


          bs.addUserInt("mu8_ip5",                 muPtr->userInt("mu8_ip5"       ));
          bs.addUserInt("prescale_mu8_ip5",        muPtr->userInt("prescale_mu8_ip5"   ));

          bs.addUserInt("mu8_ip6",                 muPtr->userInt("mu8_ip6"       ));
          bs.addUserInt("prescale_mu8_ip6",        muPtr->userInt("prescale_mu8_ip6"   ));

          bs.addUserInt("mu9_ip4",                 muPtr->userInt("mu9_ip4"       ));
          bs.addUserInt("prescale_mu9_ip4"   ,     muPtr->userInt("prescale_mu9_ip4"   ));

          bs.addUserInt("mu9_ip5",                 muPtr->userInt("mu9_ip5"       ));
          bs.addUserInt("prescale_mu9_ip5"   ,     muPtr->userInt("prescale_mu9_ip5"   ));

          bs.addUserInt("mu9_ip6",                 muPtr->userInt("mu9_ip6"       ));
          bs.addUserInt("prescale_mu9_ip6"   ,     muPtr->userInt("prescale_mu9_ip6"   ));

          bs.addUserInt("mu10p5_ip3p5",             muPtr->userInt("mu10p5_ip3p5"       ));
          bs.addUserInt("prescale_mu10p5_ip3p5",    muPtr->userInt("prescale_mu10p5_ip3p5"   ));

          bs.addUserInt("mu12_ip6",                 muPtr->userInt("mu12_ip6"       ));
          bs.addUserInt("prescale_mu12_ip6"   ,     muPtr->userInt("prescale_mu12_ip6"   ));

          bs.addUserInt("filter_mu7_4" ,             muPtr->userInt("filter_mu7_4")); 
          bs.addUserInt("filter_mu8_3" ,             muPtr->userInt("filter_mu8_3")); 
          bs.addUserInt("filter_mu8p5_3p5" ,         muPtr->userInt("filter_mu8p5_3p5")); 
          bs.addUserInt("filter_mu8_5" ,             muPtr->userInt("filter_mu8_5")); 
          bs.addUserInt("filter_mu8_6" ,             muPtr->userInt("filter_mu8_6")); 
          bs.addUserInt("filter_mu9_4" ,             muPtr->userInt("filter_mu9_4")); 
          bs.addUserInt("filter_mu9_5" ,             muPtr->userInt("filter_mu9_5")); 
          bs.addUserInt("filter_mu9_6" ,             muPtr->userInt("filter_mu9_6")); 
          bs.addUserInt("filter_mu10p5_3p5" ,        muPtr->userInt("filter_mu10p5_3p5")); 
          bs.addUserInt("filter_mu12_6" ,            muPtr->userInt("filter_mu12_6")); 


          /////////////////////// END OF VARIABLE DEFINITION //////////////////////
          //std::cout << "Found RECO candidate!"<< std::endl; 
          //append candidate at the end of our return value :)
          bsCandidates->emplace_back(bs);

          nSurvived++;

        } //closing pi loop
      } //closing k2 loop
    } //closing k1 loop
  } //closing trg muon loop
  
  iEvent.put(std::move(bsCandidates), "bs");

}//closing event loop

void BsToDsPhiKKPiMuBuilder::endJob(){
// Printouts:

std::cout << "\n--------- Bs BUIDLER MODULE ----------\n" << std::endl;
std::cout << "#Muons  in file                                    : " << nMuons  << std::endl;
std::cout << "#Tracks in file                                    : " << nTracks << std::endl;
std::cout << "#Pv     in file                                    : " << nPv     << std::endl;
std::cout << "\n#Muons for which we found primary vertex         : " << muSelCounter  << std::endl;
std::cout << "#Kaon 1 which passed the hadronic selection        : " << k1Sel1Counter << std::endl;
std::cout << "#Kaon 2 which passed the hadronic selection        : " << k2Sel1Counter << std::endl;
std::cout << "#Pions  which passed the hadronic selection        : " << piSel1Counter << std::endl;

std::cout << "#KK combinations which passed the Phi Mass cut     : " << nPhiMassCut << std::endl;
std::cout << "\n#KKPi combinations before the Ds mass cut        : " << nKKPiMu << std::endl;
std::cout << "#KKPiMu combinations which passed the Ds Mass cut  : " << nDsMassCut << std::endl;
std::cout << "#KKPiMu combinations which passed the Bs Mass cut  : " << nBsMassCut << std::endl;
std::cout << "#KKPiMu combinations which passed the Phi fit      : " << nPhiFit << std::endl;
std::cout << "#KKPiMu combinations which passed the Ds  fit      : " << nDsFit  << std::endl;
std::cout << "#KKPiMu combinations which passed the Bs  fit      : " << nBsFit  << std::endl;
std::cout << "#KKPiMu combinations which have been AOD matched   : " << nAod      << std::endl;
std::cout << "#KKPiMu combinations which have valid vertices     : " << nVtxCut   << std::endl;
std::cout << "#KKPiMu combinations which survive pointing angle  : " << nCosCut   << std::endl;
std::cout << "#KKPiMu combinations to be saved                   : " << nSurvived << std::endl;
}


DEFINE_FWK_MODULE(BsToDsPhiKKPiMuBuilder);
