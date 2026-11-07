//#include "RDs-tools-private/RDsNano/src/helper.h" 
//#include "RDs-tools-private/RDsNano/src/classes.h" 

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/BasicTransientTrack.h"

//B field 
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
// for the cut as string
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <TLorentzVector.h>
#include "helper.h" 

using namespace std;
int nEvents        = 0;
int passesHLT      = 0;
int passesHLTFilter= 0;
int patTriggerCand = 0;
int trgCand        = 0;

class Trigger : public edm::one::EDProducer<> {

public:

    //define the transienttrackcollection type, does not exist!
    //typedef std::vector<reco::TransientTrack> TransientTrackCollection;
    
    //constructor, takes a reference to the ParameterSet 
    explicit Trigger(const edm::ParameterSet&);
    //destructor
    ~Trigger() override {};
    virtual void endJob() override; // NEW! 

private:
   
    OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");

    // muon selection
    const StringCutObjectSelector<pat::Muon> muSelection_;
 
    //physics heart ---- check the constant override stuff!!!
    void produce(edm::Event&, const edm::EventSetup&) override;

    //define tokens to access data later
    const edm::InputTag muonTag;
    const edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;

    const edm::InputTag triggerBitTag;
    const edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

    const edm::InputTag triggerObjectTag;
    const edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;

    const edm::InputTag triggerPrescaleTag;
    const edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

    const edm::InputTag vertexSrcTag;
    const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

    //the maximal dR you allow between pat muon and trigger muon candidate
    const string trgFilterLabelMu7_4_ ;
    const string trgFilterLabelMu8_3_ ;
    const string trgFilterLabelMu8p5_3p5_ ;
    const string trgFilterLabelMu8_5_ ;
    const string trgFilterLabelMu8_6_ ;
    const string trgFilterLabelMu9_4_ ;
    const string trgFilterLabelMu9_5_ ;
    const string trgFilterLabelMu9_6_ ;
    const string trgFilterLabelMu10p5_3p5_ ;
    const string trgFilterLabelMu12_6_;

    const string hlt_7_4_ ;
    const string hlt_8_3_ ;
    const string hlt_8p5_3p5_ ;
    const string hlt_8_5_ ;
    const string hlt_8_6_ ;
    const string hlt_9_4_ ;
    const string hlt_9_5_ ;
    const string hlt_9_6_ ;
    const string hlt_10p5_3p5_ ;
    const string hlt_12_6_;
 
    const double maxdR_; 

};

//constructor definition (outsde class)
Trigger::Trigger(const edm::ParameterSet& iConfig):
 
  //for the muons
  muSelection_(iConfig.getParameter<std::string>                     ("muSelection")), 
  muonTag(iConfig.getParameter<edm::InputTag>                        ("muonCollection")),
  muonSrc_(consumes<std::vector<pat::Muon>>                          (muonTag)),
  //for trigger info
  triggerBitTag(iConfig.getParameter<edm::InputTag>                  ("trgResultsCollection")),
  triggerBits_(consumes<edm::TriggerResults>                         (triggerBitTag)),

  triggerObjectTag(iConfig.getParameter<edm::InputTag>               ("trgObjectsCollection")),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerObjectTag)),

  triggerPrescaleTag(iConfig.getParameter<edm::InputTag>             ("trgPrescaleCollection")),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>            (triggerPrescaleTag)),
  //vertex info
  vertexSrcTag(iConfig.getParameter<edm::InputTag>                   ("vtxCollection")),
  vertexSrc_(consumes<reco::VertexCollection>(vertexSrcTag)), 

  //parameters
  trgFilterLabelMu7_4_ (iConfig.getParameter<string>                    ("trgFilterLabelMu7_4" )),
  trgFilterLabelMu8_3_ (iConfig.getParameter<string>                    ("trgFilterLabelMu8_3" )),
  trgFilterLabelMu8p5_3p5_ (iConfig.getParameter<string>                ("trgFilterLabelMu8p5_3p5" )),
  trgFilterLabelMu8_5_ (iConfig.getParameter<string>                    ("trgFilterLabelMu8_5" )),
  trgFilterLabelMu8_6_ (iConfig.getParameter<string>                    ("trgFilterLabelMu8_6" )),
  trgFilterLabelMu9_4_ (iConfig.getParameter<string>                    ("trgFilterLabelMu9_4" )),
  trgFilterLabelMu9_5_ (iConfig.getParameter<string>                    ("trgFilterLabelMu9_5" )),
  trgFilterLabelMu9_6_ (iConfig.getParameter<string>                    ("trgFilterLabelMu9_6" )),
  trgFilterLabelMu10p5_3p5_ (iConfig.getParameter<string>               ("trgFilterLabelMu10p5_3p5" )),
  trgFilterLabelMu12_6_(iConfig.getParameter<string>                    ("trgFilterLabelMu12_6")),

  hlt_7_4_     (iConfig.getParameter<string>                           ("hlt_7_4")),
  hlt_8_3_     (iConfig.getParameter<string>                           ("hlt_8_3")),
  hlt_8p5_3p5_ (iConfig.getParameter<string>                           ("hlt_8p5_3p5")),
  hlt_8_5_     (iConfig.getParameter<string>                           ("hlt_8_5")),
  hlt_8_6_     (iConfig.getParameter<string>                           ("hlt_8_6")),
  hlt_9_4_     (iConfig.getParameter<string>                           ("hlt_9_4")),
  hlt_9_5_     (iConfig.getParameter<string>                           ("hlt_9_5")),
  hlt_9_6_     (iConfig.getParameter<string>                           ("hlt_9_6")),
  hlt_10p5_3p5_(iConfig.getParameter<string>                           ("hlt_10p5_3p5")),
  hlt_12_6_    (iConfig.getParameter<string>                           ("hlt_12_6")),

  maxdR_(iConfig.getParameter<double>                                ("maxdR_matching"))

{
  // produces 
  produces<pat::MuonCollection>                                      ("trgMuons");  
}



void Trigger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
 
  //Define handles
  //std::cout << "New event!" << std::endl;

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexSrc_, vertexHandle);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);


  // for every event, get the name list of the triggers
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

      
  unsigned int trgBitsMax = triggerBits->size();

  unsigned int index_7_4       = trgBitsMax; 
  unsigned int index_8_3       = trgBitsMax; 
  unsigned int index_8p5_3p5   = trgBitsMax; 
  unsigned int index_8_5       = trgBitsMax; 
  unsigned int index_8_6       = trgBitsMax; 
  unsigned int index_9_4       = trgBitsMax; 
  unsigned int index_9_5       = trgBitsMax; 
  unsigned int index_9_6       = trgBitsMax; 
  unsigned int index_10p5_3p5  = trgBitsMax; 
  unsigned int index_12_6      = trgBitsMax; 

  bool pass_7_4         = false; 
  bool pass_8_3         = false; 
  bool pass_8p5_3p5     = false; 
  bool pass_8_5         = false; 
  bool pass_8_6         = false; 
  bool pass_9_4         = false; 
  bool pass_9_5         = false; 
  bool pass_9_6         = false; 
  bool pass_10p5_3p5    = false; 
  bool pass_12_6        = false; 

  int prescale_7_4      = 0;
  int prescale_8_3      = 0;
  int prescale_8p5_3p5  = 0;
  int prescale_8_5      = 0;
  int prescale_8_6      = 0;
  int prescale_9_4      = 0;
  int prescale_9_5      = 0;
  int prescale_9_6      = 0;
  int prescale_10p5_3p5 = 0;
  int prescale_12_6     = 0;


  ///////////// Debugging:
  for (unsigned int i = 0; i < names.size(); ++i) {

    std::string trgName = names.triggerName(i);

    if (trgName.find(hlt_7_4_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_7_4 = names.triggerIndex(trgName);

      if ((index_7_4 < trgBitsMax) && (triggerBits->accept(index_7_4))){
        //std::cout << "index is: "  << index_7_4 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_7_4) << std::endl;
        pass_7_4 = true;
        prescale_7_4  = triggerPrescales->getPrescaleForIndex (index_7_4 );
        prescale_7_4 *= 5;
      }

    }
    if (trgName.find(hlt_8_3_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_8_3 = names.triggerIndex(trgName);
      if ((index_8_3 < trgBitsMax) && (triggerBits->accept(index_8_3))){
        //std::cout << "index is: "  << index_8_3 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_8_3) << std::endl;
        pass_8_3      = true;
        prescale_8_3  = triggerPrescales->getPrescaleForIndex (index_8_3 );
        prescale_8_3 *= 5;
      }

    }
    if (trgName.find(hlt_8_5_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_8_5 = names.triggerIndex(trgName);
      if ((index_8_5 < trgBitsMax) && (triggerBits->accept(index_8_5))){
        //std::cout << "index is: "  << index_8_5 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_8_5) << std::endl;
        pass_8_5 = true;
        prescale_8_5  = triggerPrescales->getPrescaleForIndex (index_8_5 );
        prescale_8_5 *= 5;
      }

    }
    if (trgName.find(hlt_8_6_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_8_6 = names.triggerIndex(trgName);
      if ((index_8_6 < trgBitsMax) && (triggerBits->accept(index_8_6))){
        //std::cout << "index is: "  << index_8_6 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_8_6) << std::endl;
        pass_8_6 = true;
        prescale_8_6  = triggerPrescales->getPrescaleForIndex (index_8_6 );
        prescale_8_6 *= 5;
      }

    }


    if (trgName.find(hlt_8p5_3p5_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_8p5_3p5 = names.triggerIndex(trgName);
      if ((index_8p5_3p5 < trgBitsMax) && (triggerBits->accept(index_8p5_3p5))){
        //std::cout << "index is: "  << index_8p5_3p5 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_8p5_3p5) << std::endl;
        pass_8p5_3p5 = true;
        prescale_8p5_3p5  = triggerPrescales->getPrescaleForIndex (index_8p5_3p5 );
        prescale_8p5_3p5 *= 5;
      }

    }


    if (trgName.find(hlt_9_4_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_9_4 = names.triggerIndex(trgName);
      if ((index_9_4 < trgBitsMax) && (triggerBits->accept(index_9_4))){
        //std::cout << "index is: "  << index_9_4 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_9_4) << std::endl;
        pass_9_4 = true;
        prescale_9_4  = triggerPrescales->getPrescaleForIndex (index_9_4 );
        prescale_9_4 *= 5;
      }

    }
    if (trgName.find(hlt_9_5_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_9_5 = names.triggerIndex(trgName);
      if ((index_9_5 < trgBitsMax) && (triggerBits->accept(index_9_5))){
        //std::cout << "index is: "  << index_9_5 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_9_5) << std::endl;
        pass_9_5 = true;
        prescale_9_5  = triggerPrescales->getPrescaleForIndex (index_9_5 );
        prescale_9_5 *= 5;
      }

    }
    if (trgName.find(hlt_9_6_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_9_6 = names.triggerIndex(trgName);
      if ((index_9_6 < trgBitsMax) && (triggerBits->accept(index_9_6))){
        //std::cout << "index is: "  << index_9_6 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_9_6) << std::endl;
        pass_9_6 = true;
        prescale_9_6  = triggerPrescales->getPrescaleForIndex (index_9_6 );
        prescale_9_6 *= 5;
      }

    }

    if (trgName.find(hlt_10p5_3p5_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_10p5_3p5 = names.triggerIndex(trgName);
      if ((index_10p5_3p5 < trgBitsMax) && (triggerBits->accept(index_10p5_3p5))){
        //std::cout << "index is: "  << index_10p5_3p5 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_10p5_3p5) << std::endl;
        pass_10p5_3p5 = true;
        prescale_10p5_3p5  = triggerPrescales->getPrescaleForIndex (index_10p5_3p5 );
        prescale_10p5_3p5 *= 5;
      }

    }




    if (trgName.find(hlt_12_6_) != std::string::npos){

      //std::cout << trgName << std::endl;
      index_12_6 = names.triggerIndex(trgName);
      if ((index_12_6 < trgBitsMax) && (triggerBits->accept(index_12_6))){
        //std::cout << "index is: "  << index_12_6 << std::endl;
        //std::cout << "accepted ? " << triggerBits->accept(index_12_6) << std::endl;
        pass_12_6 = true;
        prescale_12_6  = triggerPrescales->getPrescaleForIndex (index_12_6 );
        prescale_12_6 *= 5;
      }

    }


  }
  /////////////

  nEvents++; 

  //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  //pat muons contain more info than reco muons, f.e. trigger info! 
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonSrc_, muons);

  // to save
  std::unique_ptr<pat::MuonCollection> trgMuons(new pat::MuonCollection);
 
 
  //std::cout << "before passing " << muons->size() << std::endl;
  //only continue when we the event passes one of the triggers
  if (pass_7_4 || pass_8_3 || pass_8_5 || pass_8_6 || pass_9_4 || pass_9_5 || pass_9_6 || pass_12_6 || pass_8p5_3p5 || pass_10p5_3p5) {

    //std::cout << "passed :) " << muons->size() << std::endl;
    passesHLT++;

    ///////////////////////////////////////////
    // Loop over trigger objects and unpack  //
    ///////////////////////////////////////////

    std::vector<bool> passTrigger_Vec          ;
    std::vector<bool> failsSelection_Vec       ;
    std::vector<bool> filterLabelMu7_4_Vec     ;
    std::vector<bool> filterLabelMu8_3_Vec     ;
    std::vector<bool> filterLabelMu8p5_3p5_Vec ;
    std::vector<bool> filterLabelMu8_5_Vec     ;
    std::vector<bool> filterLabelMu8_6_Vec     ;
    std::vector<bool> filterLabelMu9_4_Vec     ;
    std::vector<bool> filterLabelMu9_5_Vec     ;
    std::vector<bool> filterLabelMu9_6_Vec     ;
    std::vector<bool> filterLabelMu10p5_3p5_Vec;
    std::vector<bool> filterLabelMu12_6_Vec    ;


    for(unsigned int objIdx=0; objIdx < triggerObjects->size(); ++objIdx){

      const pat::TriggerObjectStandAlone& trgObjRef = (*triggerObjects)[objIdx];

      //check pt, eta, and if trigger object is of type muon
      const bool failsSelection = 
          (std::fabs(trgObjRef.eta()) > 1.6) || 
          (trgObjRef.pt() < 6.0)             || 
          (!trgObjRef.hasTriggerObjectType(trigger::TriggerMuon));

      failsSelection_Vec.push_back(failsSelection);

      // Objects failing the selection are never matched to a muon (see the matching loop below),
      // so there is no need to copy them and unpack their filter labels (the expensive part).
      if (failsSelection) {
        passTrigger_Vec          .push_back(false);
        filterLabelMu7_4_Vec     .push_back(false);
        filterLabelMu8_3_Vec     .push_back(false);
        filterLabelMu8p5_3p5_Vec .push_back(false);
        filterLabelMu8_5_Vec     .push_back(false);
        filterLabelMu8_6_Vec     .push_back(false);
        filterLabelMu9_4_Vec     .push_back(false);
        filterLabelMu9_5_Vec     .push_back(false);
        filterLabelMu9_6_Vec     .push_back(false);
        filterLabelMu10p5_3p5_Vec.push_back(false);
        filterLabelMu12_6_Vec    .push_back(false);
        continue;
      }

      // copy needed: unpacking modifies the object
      pat::TriggerObjectStandAlone trgObj = trgObjRef;

      //unpack trigger labels
      trgObj.unpackFilterLabels(iEvent, *triggerBits);

      int filterLabelMu7_4        = trgObj.hasFilterLabel(trgFilterLabelMu7_4_)     ;
      int filterLabelMu8_3        = trgObj.hasFilterLabel(trgFilterLabelMu8_3_)     ;
      int filterLabelMu8p5_3p5    = trgObj.hasFilterLabel(trgFilterLabelMu8p5_3p5_) ;
      int filterLabelMu8_5        = trgObj.hasFilterLabel(trgFilterLabelMu8_5_)     ;
      int filterLabelMu8_6        = trgObj.hasFilterLabel(trgFilterLabelMu8_6_)     ;
      int filterLabelMu9_4        = trgObj.hasFilterLabel(trgFilterLabelMu9_4_)     ;
      int filterLabelMu9_5        = trgObj.hasFilterLabel(trgFilterLabelMu9_5_)     ;
      int filterLabelMu9_6        = trgObj.hasFilterLabel(trgFilterLabelMu9_6_)     ;
      int filterLabelMu10p5_3p5   = trgObj.hasFilterLabel(trgFilterLabelMu10p5_3p5_);
      int filterLabelMu12_6       = trgObj.hasFilterLabel(trgFilterLabelMu12_6_)    ;
  
      // check if the trigger object was actually firing our trigger
      const bool passTrigger =
          (pass_7_4       && filterLabelMu7_4)      ||
          (pass_8_3       && filterLabelMu8_3)      ||
          (pass_8p5_3p5   && filterLabelMu8p5_3p5)  ||
          (pass_8_5       && filterLabelMu8_5)      ||
          (pass_8_6       && filterLabelMu8_6)      ||
          (pass_9_4       && filterLabelMu9_4)      ||
          (pass_9_5       && filterLabelMu9_5)      ||
          (pass_9_6       && filterLabelMu9_6)      ||
          (pass_10p5_3p5  && filterLabelMu10p5_3p5) ||
          (pass_12_6      && filterLabelMu12_6)     ;

      passTrigger_Vec          .push_back(passTrigger)           ;      
      filterLabelMu7_4_Vec     .push_back(filterLabelMu7_4)     ;
      filterLabelMu8_3_Vec     .push_back(filterLabelMu8_3)     ;
      filterLabelMu8p5_3p5_Vec .push_back(filterLabelMu8p5_3p5) ;
      filterLabelMu8_5_Vec     .push_back(filterLabelMu8_5)     ;
      filterLabelMu8_6_Vec     .push_back(filterLabelMu8_6)     ;
      filterLabelMu9_4_Vec     .push_back(filterLabelMu9_4)     ;
      filterLabelMu9_5_Vec     .push_back(filterLabelMu9_5)     ;
      filterLabelMu9_6_Vec     .push_back(filterLabelMu9_6)     ;
      filterLabelMu10p5_3p5_Vec.push_back(filterLabelMu10p5_3p5);
      filterLabelMu12_6_Vec    .push_back(filterLabelMu12_6)    ;


      if (passTrigger) passesHLTFilter++ ;
    }


    // count the nr of pat::muons which fired the trigger -> reset to 0 in every event
    int counter = 0;

    ///////////////////////////////
    // Loop over muons and match //
    ///////////////////////////////

    //std::cout << "#muons in this event " << muons->size() << std::endl;

    for (unsigned int muIdx=0; muIdx<muons->size(); ++muIdx){
      //std::cout << iEvent.id().event() << std::endl;

      const pat::Muon& muon = (*muons)[muIdx];    
  
      // muon cuts
      if (!muSelection_(muon)) continue; 
      //std::cout<<"found pat muon with pt:"<< muon.pt() << std::endl;
 
      // validd track ?
      if (muon.bestTrack() == nullptr) continue;
      const reco::TransientTrack ttrackTrgMuon(*(muon.bestTrack()), paramField);  
      if (!ttrackTrgMuon.isValid()) continue;

      // this is a possible trigger cand 
      patTriggerCand++;
  
      // initialize start values
      float drMuonTrgObj   = 0.0;
      int   trgObjIdx      = -1;


      //////////////////////////////////////////////
      // Match muon candidates to trigger objects //
      //////////////////////////////////////////////

      for(unsigned int objIdx=0; objIdx < triggerObjects->size(); ++objIdx){

        if (!passTrigger_Vec.at(objIdx)) continue;
        if (failsSelection_Vec.at(objIdx)) continue;

        const pat::TriggerObjectStandAlone& trgObj = (*triggerObjects)[objIdx]; // no copy needed here

        if (trgObj.charge() * muon.charge() < 0) continue; // exclude opposite charge

 
        float dr = reco::deltaR(trgObj, muon);

        if( ((dr < drMuonTrgObj) || (trgObjIdx < 0))  && (dr < maxdR_)){
              //modify dR and the reco and triggering muon index
  	      drMuonTrgObj = dr;
  	      trgObjIdx    = objIdx;
              //std::cout << "pt: " << trgObj.pt() << "eta: " << trgObj.eta() << std::endl;
              //std::cout << "ID: " << trgObj.hasTriggerObjectType(trigger::TriggerMuon) << std::endl;
        }
  
      } // closing loop over trg muons
  
      //save pat muon if we found a matching trigger object 
      if(trgObjIdx != -1){
     
        trgCand++;
   
        //the following line does a copy of muon with the name trgMatchedMuon (same properties different adress)
  	pat::Muon trgMatchedMuon(muon);
  
        // only take muons with valid tracks
  
        //store
  	trgMatchedMuon.addUserInt("muIdx",     muIdx    );
   	trgMatchedMuon.addUserInt("trgObjIdx", trgObjIdx);

        trgMatchedMuon.addUserInt("mu7_ip4",           pass_7_4);
        trgMatchedMuon.addUserInt("prescale_mu7_ip4",  prescale_7_4);

        trgMatchedMuon.addUserInt("mu8_ip3",           pass_8_3);
        trgMatchedMuon.addUserInt("prescale_mu8_ip3",  prescale_8_3);

        trgMatchedMuon.addUserInt("mu8_ip5",           pass_8_5);
        trgMatchedMuon.addUserInt("prescale_mu8_ip5",  prescale_8_5);

        trgMatchedMuon.addUserInt("mu8p5_ip3p5",           pass_8p5_3p5);
        trgMatchedMuon.addUserInt("prescale_mu8p5_ip3p5",  prescale_8p5_3p5);


        trgMatchedMuon.addUserInt("mu8_ip6",           pass_8_6);
        trgMatchedMuon.addUserInt("prescale_mu8_ip6",  prescale_8_6);

        trgMatchedMuon.addUserInt("mu9_ip4",           pass_9_4);
        trgMatchedMuon.addUserInt("prescale_mu9_ip4",  prescale_9_4);

        trgMatchedMuon.addUserInt("mu9_ip5",           pass_9_5);
        trgMatchedMuon.addUserInt("prescale_mu9_ip5",  prescale_9_5);

        trgMatchedMuon.addUserInt("mu9_ip6",           pass_9_6);
        trgMatchedMuon.addUserInt("prescale_mu9_ip6",  prescale_9_6);

        trgMatchedMuon.addUserInt("mu10p5_ip3p5",           pass_10p5_3p5);
        trgMatchedMuon.addUserInt("prescale_mu10p5_ip3p5",  prescale_10p5_3p5);

        trgMatchedMuon.addUserInt("mu12_ip6",              pass_12_6);
        trgMatchedMuon.addUserInt("prescale_mu12_ip6",  prescale_12_6);

        trgMatchedMuon.addUserInt("filter_mu7_4"         ,  filterLabelMu7_4_Vec      .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu8_3"         ,  filterLabelMu8_3_Vec      .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu8_5"         ,  filterLabelMu8_5_Vec      .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu8p5_3p5"     ,  filterLabelMu8p5_3p5_Vec  .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu8_6"         ,  filterLabelMu8_6_Vec      .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu9_4"         ,  filterLabelMu9_4_Vec      .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu9_5"         ,  filterLabelMu9_5_Vec      .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu9_6"         ,  filterLabelMu9_6_Vec      .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu10p5_3p5"    ,  filterLabelMu10p5_3p5_Vec .at(trgObjIdx));
        trgMatchedMuon.addUserInt("filter_mu12_6"        ,  filterLabelMu12_6_Vec     .at(trgObjIdx));

        //std::cout << "found trg obj!" << std::endl; 
        //if (!beamSpotHandle.isValid()) continue;
        //const reco::BeamSpot& beamSpot = *beamSpotHandle;
        //double dxyErr = muon.bestTrack()->dxyError();
        //double ipSig;
        //if (dxyErr > 0) {
        //  ipSig = std::abs(muon.bestTrack()->dxy(beamSpot.position()) / dxyErr);
        //  std::cout << "ip sig is: " << ipSig << std::endl;
        //} else {
        //  ipSig = 0;
        //  std::cout << "dxyError zero or negative" << std::endl;
        //}

        //std::cout << "it has impact parameter: " << ipSig << std::endl;
        counter++;
        trgMuons->emplace_back(trgMatchedMuon); 
  
      }
    } //close the loop over pat muons
  } //close if condition of HLT

  //store the new collections in the event
  iEvent.put(std::move(trgMuons),    "trgMuons");

} //close produce function

void Trigger::endJob(){
std::cout << "\n--------- TRIGGER MODULE ----------\n" << std::endl;
std::cout << "#Events                                                   : " << nEvents  << std::endl;
std::cout << "#Events where HLT fired                                   : " << passesHLT  << std::endl;
std::cout << "#Events where HLT fired and HLT Filter passed             : " << passesHLTFilter << std::endl;
std::cout << "#Muon candidates which possibly fire HLT (unmatched)      : " << patTriggerCand  << std::endl;
std::cout << "#Muon candidates which fired HLT (matched to trg object!) : " << trgCand         << std::endl;

}

DEFINE_FWK_MODULE(Trigger);

