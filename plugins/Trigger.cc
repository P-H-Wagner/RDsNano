//#include "RDs-tools-private/RDsNano/src/helper.h" 
//#include "RDs-tools-private/RDsNano/src/classes.h" 

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
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

class Trigger : public edm::EDProducer {

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
    const string trgFilterLabelMu8_5_ ;
    const string trgFilterLabelMu8_6_ ;
    const string trgFilterLabelMu9_4_ ;
    const string trgFilterLabelMu9_5_ ;
    const string trgFilterLabelMu9_6_ ;
    const string trgFilterLabelMu12_6_;

    const string hlt_7_4_p0_;
    const string hlt_7_4_p1_;
    const string hlt_7_4_p2_;
    const string hlt_7_4_p3_;
    const string hlt_7_4_p4_;

    const string hlt_8_3_p0_;
    const string hlt_8_3_p1_;
    const string hlt_8_3_p2_;
    const string hlt_8_3_p3_;
    const string hlt_8_3_p4_;

    const string hlt_8_5_p0_;
    const string hlt_8_5_p1_;
    const string hlt_8_5_p2_;
    const string hlt_8_5_p3_;
    const string hlt_8_5_p4_;

    const string hlt_8_6_p0_;
    const string hlt_8_6_p1_;
    const string hlt_8_6_p2_;
    const string hlt_8_6_p3_;
    const string hlt_8_6_p4_;

    const string hlt_9_4_p0_;
    const string hlt_9_4_p1_;
    const string hlt_9_4_p2_;
    const string hlt_9_4_p3_;
    const string hlt_9_4_p4_;

    const string hlt_9_5_p0_;
    const string hlt_9_5_p1_;
    const string hlt_9_5_p2_;
    const string hlt_9_5_p3_;
    const string hlt_9_5_p4_;

    const string hlt_9_6_p0_;
    const string hlt_9_6_p1_;
    const string hlt_9_6_p2_;
    const string hlt_9_6_p3_;
    const string hlt_9_6_p4_;

    const string hlt_12_6_p0_;
    const string hlt_12_6_p1_;
    const string hlt_12_6_p2_;
    const string hlt_12_6_p3_;
    const string hlt_12_6_p4_;


 
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
  trgFilterLabelMu8_5_ (iConfig.getParameter<string>                    ("trgFilterLabelMu8_5" )),
  trgFilterLabelMu8_6_ (iConfig.getParameter<string>                    ("trgFilterLabelMu8_6" )),
  trgFilterLabelMu9_4_ (iConfig.getParameter<string>                    ("trgFilterLabelMu9_4" )),
  trgFilterLabelMu9_5_ (iConfig.getParameter<string>                    ("trgFilterLabelMu9_5" )),
  trgFilterLabelMu9_6_ (iConfig.getParameter<string>                    ("trgFilterLabelMu9_6" )),
  trgFilterLabelMu12_6_(iConfig.getParameter<string>                    ("trgFilterLabelMu12_6")),

  hlt_7_4_p0_(iConfig.getParameter<string>                           ("hlt_7_4_p0")),
  hlt_7_4_p1_(iConfig.getParameter<string>                           ("hlt_7_4_p1")),
  hlt_7_4_p2_(iConfig.getParameter<string>                           ("hlt_7_4_p2")),
  hlt_7_4_p3_(iConfig.getParameter<string>                           ("hlt_7_4_p3")),
  hlt_7_4_p4_(iConfig.getParameter<string>                           ("hlt_7_4_p4")),

  hlt_8_3_p0_(iConfig.getParameter<string>                           ("hlt_8_3_p0")),
  hlt_8_3_p1_(iConfig.getParameter<string>                           ("hlt_8_3_p1")),
  hlt_8_3_p2_(iConfig.getParameter<string>                           ("hlt_8_3_p2")),
  hlt_8_3_p3_(iConfig.getParameter<string>                           ("hlt_8_3_p3")),
  hlt_8_3_p4_(iConfig.getParameter<string>                           ("hlt_8_3_p4")),

  hlt_8_5_p0_(iConfig.getParameter<string>                           ("hlt_8_5_p0")),
  hlt_8_5_p1_(iConfig.getParameter<string>                           ("hlt_8_5_p1")),
  hlt_8_5_p2_(iConfig.getParameter<string>                           ("hlt_8_5_p2")),
  hlt_8_5_p3_(iConfig.getParameter<string>                           ("hlt_8_5_p3")),
  hlt_8_5_p4_(iConfig.getParameter<string>                           ("hlt_8_5_p4")),

  hlt_8_6_p0_(iConfig.getParameter<string>                           ("hlt_8_6_p0")),
  hlt_8_6_p1_(iConfig.getParameter<string>                           ("hlt_8_6_p1")),
  hlt_8_6_p2_(iConfig.getParameter<string>                           ("hlt_8_6_p2")),
  hlt_8_6_p3_(iConfig.getParameter<string>                           ("hlt_8_6_p3")),
  hlt_8_6_p4_(iConfig.getParameter<string>                           ("hlt_8_6_p4")),

  hlt_9_4_p0_(iConfig.getParameter<string>                           ("hlt_9_4_p0")),
  hlt_9_4_p1_(iConfig.getParameter<string>                           ("hlt_9_4_p1")),
  hlt_9_4_p2_(iConfig.getParameter<string>                           ("hlt_9_4_p2")),
  hlt_9_4_p3_(iConfig.getParameter<string>                           ("hlt_9_4_p3")),
  hlt_9_4_p4_(iConfig.getParameter<string>                           ("hlt_9_4_p4")),

  hlt_9_5_p0_(iConfig.getParameter<string>                           ("hlt_9_5_p0")),
  hlt_9_5_p1_(iConfig.getParameter<string>                           ("hlt_9_5_p1")),
  hlt_9_5_p2_(iConfig.getParameter<string>                           ("hlt_9_5_p2")),
  hlt_9_5_p3_(iConfig.getParameter<string>                           ("hlt_9_5_p3")),
  hlt_9_5_p4_(iConfig.getParameter<string>                           ("hlt_9_5_p4")),

  hlt_9_6_p0_(iConfig.getParameter<string>                           ("hlt_9_6_p0")),
  hlt_9_6_p1_(iConfig.getParameter<string>                           ("hlt_9_6_p1")),
  hlt_9_6_p2_(iConfig.getParameter<string>                           ("hlt_9_6_p2")),
  hlt_9_6_p3_(iConfig.getParameter<string>                           ("hlt_9_6_p3")),
  hlt_9_6_p4_(iConfig.getParameter<string>                           ("hlt_9_6_p4")),

  hlt_12_6_p0_(iConfig.getParameter<string>                          ("hlt_12_6_p0")),
  hlt_12_6_p1_(iConfig.getParameter<string>                          ("hlt_12_6_p1")),
  hlt_12_6_p2_(iConfig.getParameter<string>                          ("hlt_12_6_p2")),
  hlt_12_6_p3_(iConfig.getParameter<string>                          ("hlt_12_6_p3")),
  hlt_12_6_p4_(iConfig.getParameter<string>                          ("hlt_12_6_p4")),

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

  ///////////// Debugging:
  std::cout << "Available trigger names:" << std::endl;
  for (unsigned int i = 0; i < names.size(); ++i) {
    std::cout << names.triggerName(i) << std::endl;
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
 
  // Getting the indices of the HLT paths
  unsigned int index_7_4_p0      = names.triggerIndex(hlt_7_4_p0_); 
  unsigned int index_7_4_p1      = names.triggerIndex(hlt_7_4_p1_); 
  unsigned int index_7_4_p2      = names.triggerIndex(hlt_7_4_p2_); 
  unsigned int index_7_4_p3      = names.triggerIndex(hlt_7_4_p3_); 
  unsigned int index_7_4_p4      = names.triggerIndex(hlt_7_4_p4_); 

  unsigned int index_8_3_p0      = names.triggerIndex(hlt_8_3_p0_); 
  unsigned int index_8_3_p1      = names.triggerIndex(hlt_8_3_p1_); 
  unsigned int index_8_3_p2      = names.triggerIndex(hlt_8_3_p2_); 
  unsigned int index_8_3_p3      = names.triggerIndex(hlt_8_3_p3_); 
  unsigned int index_8_3_p4      = names.triggerIndex(hlt_8_3_p4_); 

  unsigned int index_8_5_p0      = names.triggerIndex(hlt_8_5_p0_); 
  unsigned int index_8_5_p1      = names.triggerIndex(hlt_8_5_p1_); 
  unsigned int index_8_5_p2      = names.triggerIndex(hlt_8_5_p2_); 
  unsigned int index_8_5_p3      = names.triggerIndex(hlt_8_5_p3_); 
  unsigned int index_8_5_p4      = names.triggerIndex(hlt_8_5_p4_); 

  unsigned int index_8_6_p0      = names.triggerIndex(hlt_8_6_p0_); 
  unsigned int index_8_6_p1      = names.triggerIndex(hlt_8_6_p1_); 
  unsigned int index_8_6_p2      = names.triggerIndex(hlt_8_6_p2_); 
  unsigned int index_8_6_p3      = names.triggerIndex(hlt_8_6_p3_); 
  unsigned int index_8_6_p4      = names.triggerIndex(hlt_8_6_p4_); 

  unsigned int index_9_4_p0      = names.triggerIndex(hlt_9_4_p0_); 
  unsigned int index_9_4_p1      = names.triggerIndex(hlt_9_4_p1_); 
  unsigned int index_9_4_p2      = names.triggerIndex(hlt_9_4_p2_); 
  unsigned int index_9_4_p3      = names.triggerIndex(hlt_9_4_p3_); 
  unsigned int index_9_4_p4      = names.triggerIndex(hlt_9_4_p4_); 

  unsigned int index_9_5_p0      = names.triggerIndex(hlt_9_5_p0_); 
  unsigned int index_9_5_p1      = names.triggerIndex(hlt_9_5_p1_); 
  unsigned int index_9_5_p2      = names.triggerIndex(hlt_9_5_p2_); 
  unsigned int index_9_5_p3      = names.triggerIndex(hlt_9_5_p3_); 
  unsigned int index_9_5_p4      = names.triggerIndex(hlt_9_5_p4_); 

  unsigned int index_9_6_p0      = names.triggerIndex(hlt_9_6_p0_); 
  unsigned int index_9_6_p1      = names.triggerIndex(hlt_9_6_p1_); 
  unsigned int index_9_6_p2      = names.triggerIndex(hlt_9_6_p2_); 
  unsigned int index_9_6_p3      = names.triggerIndex(hlt_9_6_p3_); 
  unsigned int index_9_6_p4      = names.triggerIndex(hlt_9_6_p4_); 

  unsigned int index_12_6_p0      = names.triggerIndex(hlt_12_6_p0_); 
  unsigned int index_12_6_p1      = names.triggerIndex(hlt_12_6_p1_); 
  unsigned int index_12_6_p2      = names.triggerIndex(hlt_12_6_p2_); 
  unsigned int index_12_6_p3      = names.triggerIndex(hlt_12_6_p3_); 
  unsigned int index_12_6_p4      = names.triggerIndex(hlt_12_6_p4_); 

  
  std::cout << "trigger bits size: " << triggerBits->size() << std::endl;
  std::cout << "Mu 7: " << std::endl;
  std::cout <<"0 index  "<< index_7_4_p0 << std::endl; // << " and accepted? " <<      (triggerBits->accept(index_7_4_p0)) << std::endl;
  std::cout <<"1 index  "<< index_7_4_p1 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_7_4_p1)) << std::endl;
  std::cout <<"2 index  "<< index_7_4_p2 << std::endl; // << " and accepted? " <<    (triggerBits->accept(index_7_4_p2)) << std::endl;
  std::cout <<"3 index  "<< index_7_4_p3 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_7_4_p3)) << std::endl;
  std::cout <<"4 index  "<< index_7_4_p4 << std::endl; // << " and accepted? " <<    (triggerBits->accept(index_7_4_p4)) << std::endl;
  std::cout << "Mu 9 IP4 : " << std::endl;
  std::cout <<"0 index  "<< index_9_4_p0 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_4_p0)) << std::endl;
  std::cout <<"1 index  "<< index_9_4_p1 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_4_p1)) << std::endl;
  std::cout <<"2 index  "<< index_9_4_p2 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_4_p2)) << std::endl;
  std::cout <<"3 index  "<< index_9_4_p3 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_4_p3)) << std::endl;
  std::cout <<"4 index  "<< index_9_4_p4 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_4_p4)) << std::endl;
  std::cout << "Mu 9 IP5 : " << std::endl;
  std::cout <<"0 index  "<< index_9_5_p0 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_5_p0)) << std::endl;
  std::cout <<"1 index  "<< index_9_5_p1 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_5_p1)) << std::endl;
  std::cout <<"2 index  "<< index_9_5_p2 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_5_p2)) << std::endl;
  std::cout <<"3 index  "<< index_9_5_p3 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_5_p3)) << std::endl;
  std::cout <<"4 index  "<< index_9_5_p4 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_5_p4)) << std::endl;
  std::cout << "Mu 9 IP6 : " << std::endl;
  std::cout <<"0 index  "<< index_9_6_p0 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_6_p0)) << std::endl;
  std::cout <<"1 index  "<< index_9_6_p1 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_6_p1)) << std::endl;
  std::cout <<"2 index  "<< index_9_6_p2 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_6_p2)) << std::endl;
  std::cout <<"3 index  "<< index_9_6_p3 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_6_p3)) << std::endl;
  std::cout <<"4 index  "<< index_9_6_p4 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_9_6_p4)) << std::endl;
  std::cout << "Mu 8 IP3: " << std::endl;
  std::cout <<"0 index  "<< index_8_3_p0 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_8_3_p0)) << std::endl;
  std::cout <<"1 index  "<< index_8_3_p1 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_8_3_p1)) << std::endl;
  std::cout <<"2 index  "<< index_8_3_p2 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_8_3_p2)) << std::endl;
  std::cout <<"3 index  "<< index_8_3_p3 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_8_3_p3)) << std::endl;
  std::cout <<"4 index  "<< index_8_3_p4 << std::endl; // << " and accepted? " <<     (triggerBits->accept(index_8_3_p4)) << std::endl;



  //std::cout << "Mu 7: " << std::endl;
  //std::cout <<"0 index is valid? "<< (index_7_4_p0 < triggerBits->size())  << " and accepted? " <<      (triggerBits->accept(index_7_4_p0)) << std::endl;
  //std::cout <<"1 index is valid? "<< (index_7_4_p1 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_7_4_p1)) << std::endl;
  //std::cout <<"2 index is valid? "<< (index_7_4_p2 < triggerBits->size())  << " and accepted? " <<    (triggerBits->accept(index_7_4_p2)) << std::endl;
  //std::cout <<"3 index is valid? "<< (index_7_4_p3 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_7_4_p3)) << std::endl;
  //std::cout <<"4 index is valid? "<< (index_7_4_p4 < triggerBits->size())  << " and accepted? " <<    (triggerBits->accept(index_7_4_p4)) << std::endl;
  //std::cout << "Mu 9 IP4 : " << std::endl;
  //std::cout <<"0 index is valid? "<< (index_9_4_p0 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_4_p0)) << std::endl;
  //std::cout <<"1 index is valid? "<< (index_9_4_p1 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_4_p1)) << std::endl;
  //std::cout <<"2 index is valid? "<< (index_9_4_p2 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_4_p2)) << std::endl;
  //std::cout <<"3 index is valid? "<< (index_9_4_p3 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_4_p3)) << std::endl;
  //std::cout <<"4 index is valid? "<< (index_9_4_p4 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_4_p4)) << std::endl;
  //std::cout << "Mu 9 IP5 : " << std::endl;
  //std::cout <<"0 index is valid? "<< (index_9_5_p0 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_5_p0)) << std::endl;
  //std::cout <<"1 index is valid? "<< (index_9_5_p1 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_5_p1)) << std::endl;
  //std::cout <<"2 index is valid? "<< (index_9_5_p2 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_5_p2)) << std::endl;
  //std::cout <<"3 index is valid? "<< (index_9_5_p3 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_5_p3)) << std::endl;
  //std::cout <<"4 index is valid? "<< (index_9_5_p4 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_5_p4)) << std::endl;
  std::cout << "Mu 9 IP6 : " << std::endl;
  std::cout <<"0 index is valid? "<< (index_9_6_p0 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_6_p0)) << std::endl;
  std::cout <<"1 index is valid? "<< (index_9_6_p1 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_6_p1)) << std::endl;
  std::cout <<"2 index is valid? "<< (index_9_6_p2 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_6_p2)) << std::endl;
  std::cout <<"3 index is valid? "<< (index_9_6_p3 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_6_p3)) << std::endl;
  std::cout <<"4 index is valid? "<< (index_9_6_p4 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_9_6_p4)) << std::endl;
  std::cout << "Mu 8 IP3: " << std::endl;
  std::cout <<"0 index is valid? "<< (index_8_3_p0 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_8_3_p0)) << std::endl;
  std::cout <<"1 index is valid? "<< (index_8_3_p1 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_8_3_p1)) << std::endl;
  std::cout <<"2 index is valid? "<< (index_8_3_p2 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_8_3_p2)) << std::endl;
  std::cout <<"3 index is valid? "<< (index_8_3_p3 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_8_3_p3)) << std::endl;
  std::cout <<"4 index is valid? "<< (index_8_3_p4 < triggerBits->size())  << " and accepted? " <<     (triggerBits->accept(index_8_3_p4)) << std::endl;

  //default is false  
  bool pass_7_4_p0  = false, pass_7_4_p1  = false, pass_7_4_p2  = false, pass_7_4_p3  = false, pass_7_4_p4  = false;
  bool pass_8_3_p0  = false, pass_8_3_p1  = false, pass_8_3_p2  = false, pass_8_3_p3  = false, pass_8_3_p4  = false;
  bool pass_8_5_p0  = false, pass_8_5_p1  = false, pass_8_5_p2  = false, pass_8_5_p3  = false, pass_8_5_p4  = false;
  bool pass_8_6_p0  = false, pass_8_6_p1  = false, pass_8_6_p2  = false, pass_8_6_p3  = false, pass_8_6_p4  = false;
  bool pass_9_4_p0  = false, pass_9_4_p1  = false, pass_9_4_p2  = false, pass_9_4_p3  = false, pass_9_4_p4  = false;
  bool pass_9_5_p0  = false, pass_9_5_p1  = false, pass_9_5_p2  = false, pass_9_5_p3  = false, pass_9_5_p4  = false;
  bool pass_9_6_p0  = false, pass_9_6_p1  = false, pass_9_6_p2  = false, pass_9_6_p3  = false, pass_9_6_p4  = false;
  bool pass_12_6_p0 = false, pass_12_6_p1 = false, pass_12_6_p2 = false, pass_12_6_p3 = false, pass_12_6_p4 = false;

  // check first if the index is valid, i.e. if it is not out of range (maximum is givrn by triggerBits->size())
  // and if so, check if the trigger has fired with accept() 

  pass_7_4_p0      = ((index_7_4_p0 < triggerBits->size())   && (triggerBits->accept(index_7_4_p0)));
  pass_7_4_p1      = ((index_7_4_p1 < triggerBits->size())   && (triggerBits->accept(index_7_4_p1)));
  pass_7_4_p2      = ((index_7_4_p2 < triggerBits->size())   && (triggerBits->accept(index_7_4_p2)));
  pass_7_4_p3      = ((index_7_4_p3 < triggerBits->size())   && (triggerBits->accept(index_7_4_p3)));
  pass_7_4_p4      = ((index_7_4_p4 < triggerBits->size())   && (triggerBits->accept(index_7_4_p4)));

  pass_8_3_p0      = ((index_8_3_p0 < triggerBits->size())   && (triggerBits->accept(index_8_3_p0)));
  pass_8_3_p1      = ((index_8_3_p1 < triggerBits->size())   && (triggerBits->accept(index_8_3_p1)));
  pass_8_3_p2      = ((index_8_3_p2 < triggerBits->size())   && (triggerBits->accept(index_8_3_p2)));
  pass_8_3_p3      = ((index_8_3_p3 < triggerBits->size())   && (triggerBits->accept(index_8_3_p3)));
  pass_8_3_p4      = ((index_8_3_p4 < triggerBits->size())   && (triggerBits->accept(index_8_3_p4)));

  pass_8_5_p0      = ((index_8_5_p0 < triggerBits->size())   && (triggerBits->accept(index_8_5_p0)));
  pass_8_5_p1      = ((index_8_5_p1 < triggerBits->size())   && (triggerBits->accept(index_8_5_p1)));
  pass_8_5_p2      = ((index_8_5_p2 < triggerBits->size())   && (triggerBits->accept(index_8_5_p2)));
  pass_8_5_p3      = ((index_8_5_p3 < triggerBits->size())   && (triggerBits->accept(index_8_5_p3)));
  pass_8_5_p4      = ((index_8_5_p4 < triggerBits->size())   && (triggerBits->accept(index_8_5_p4)));

  pass_8_6_p0      = ((index_8_6_p0 < triggerBits->size())   && (triggerBits->accept(index_8_6_p0)));
  pass_8_6_p1      = ((index_8_6_p1 < triggerBits->size())   && (triggerBits->accept(index_8_6_p1)));
  pass_8_6_p2      = ((index_8_6_p2 < triggerBits->size())   && (triggerBits->accept(index_8_6_p2)));
  pass_8_6_p3      = ((index_8_6_p3 < triggerBits->size())   && (triggerBits->accept(index_8_6_p3)));
  pass_8_6_p4      = ((index_8_6_p4 < triggerBits->size())   && (triggerBits->accept(index_8_6_p4)));

  pass_9_4_p0      = ((index_9_4_p0 < triggerBits->size())   && (triggerBits->accept(index_9_4_p0)));
  pass_9_4_p1      = ((index_9_4_p1 < triggerBits->size())   && (triggerBits->accept(index_9_4_p1)));
  pass_9_4_p2      = ((index_9_4_p2 < triggerBits->size())   && (triggerBits->accept(index_9_4_p2)));
  pass_9_4_p3      = ((index_9_4_p3 < triggerBits->size())   && (triggerBits->accept(index_9_4_p3)));
  pass_9_4_p4      = ((index_9_4_p4 < triggerBits->size())   && (triggerBits->accept(index_9_4_p4)));

  pass_9_5_p0      = ((index_9_5_p0 < triggerBits->size())   && (triggerBits->accept(index_9_5_p0)));
  pass_9_5_p1      = ((index_9_5_p1 < triggerBits->size())   && (triggerBits->accept(index_9_5_p1)));
  pass_9_5_p2      = ((index_9_5_p2 < triggerBits->size())   && (triggerBits->accept(index_9_5_p2)));
  pass_9_5_p3      = ((index_9_5_p3 < triggerBits->size())   && (triggerBits->accept(index_9_5_p3)));
  pass_9_5_p4      = ((index_9_5_p4 < triggerBits->size())   && (triggerBits->accept(index_9_5_p4)));

  pass_9_6_p0      = ((index_9_6_p0 < triggerBits->size())   && (triggerBits->accept(index_9_6_p0)));
  pass_9_6_p1      = ((index_9_6_p1 < triggerBits->size())   && (triggerBits->accept(index_9_6_p1)));
  pass_9_6_p2      = ((index_9_6_p2 < triggerBits->size())   && (triggerBits->accept(index_9_6_p2)));
  pass_9_6_p3      = ((index_9_6_p3 < triggerBits->size())   && (triggerBits->accept(index_9_6_p3)));
  pass_9_6_p4      = ((index_9_6_p4 < triggerBits->size())   && (triggerBits->accept(index_9_6_p4)));

  pass_12_6_p0      = ((index_12_6_p0 < triggerBits->size())   && (triggerBits->accept(index_12_6_p0)));
  pass_12_6_p1      = ((index_12_6_p1 < triggerBits->size())   && (triggerBits->accept(index_12_6_p1)));
  pass_12_6_p2      = ((index_12_6_p2 < triggerBits->size())   && (triggerBits->accept(index_12_6_p2)));
  pass_12_6_p3      = ((index_12_6_p3 < triggerBits->size())   && (triggerBits->accept(index_12_6_p3)));
  pass_12_6_p4      = ((index_12_6_p4 < triggerBits->size())   && (triggerBits->accept(index_12_6_p4)));


  bool pass_7_4    = (pass_7_4_p0 || pass_7_4_p1 || pass_7_4_p2 || pass_7_4_p3 || pass_7_4_p4) ;
  bool pass_8_3    = (pass_8_3_p0 || pass_8_3_p1 || pass_8_3_p2 || pass_8_3_p3 || pass_8_3_p4) ;
  bool pass_8_5    = (pass_8_5_p0 || pass_8_5_p1 || pass_8_5_p2 || pass_8_5_p3 || pass_8_5_p4) ;
  bool pass_8_6    = (pass_8_6_p0 || pass_8_6_p1 || pass_8_6_p2 || pass_8_6_p3 || pass_8_6_p4) ;
  bool pass_9_4    = (pass_9_4_p0 || pass_9_4_p1 || pass_9_4_p2 || pass_9_4_p3 || pass_9_4_p4) ;
  bool pass_9_5    = (pass_9_5_p0 || pass_9_5_p1 || pass_9_5_p2 || pass_9_5_p3 || pass_9_5_p4) ;
  bool pass_9_6    = (pass_9_6_p0 || pass_9_6_p1 || pass_9_6_p2 || pass_9_6_p3 || pass_9_6_p4) ;
  bool pass_12_6   = (pass_12_6_p0 || pass_12_6_p1 || pass_12_6_p2 || pass_12_6_p3 || pass_12_6_p4) ;

 
  //std::cout << "before passing " << muons->size() << std::endl;
  //only continue when we the event passes one of the triggers
  if (pass_7_4 || pass_8_3 || pass_8_5 || pass_8_6 || pass_9_4 || pass_9_5 || pass_9_6 || pass_12_6) {

    //std::cout << "passed :) " << muons->size() << std::endl;
    
    // get the prescale (if one part is on, all should be, but leet's be on the safe side)
    int prescale_7_4_p0 = triggerPrescales->getPrescaleForIndex(index_7_4_p0);
    int prescale_7_4_p1 = triggerPrescales->getPrescaleForIndex(index_7_4_p1);
    int prescale_7_4_p2 = triggerPrescales->getPrescaleForIndex(index_7_4_p2);
    int prescale_7_4_p3 = triggerPrescales->getPrescaleForIndex(index_7_4_p3);
    int prescale_7_4_p4 = triggerPrescales->getPrescaleForIndex(index_7_4_p4);

    int prescale_8_3_p0 = triggerPrescales->getPrescaleForIndex(index_8_3_p0);
    int prescale_8_3_p1 = triggerPrescales->getPrescaleForIndex(index_8_3_p1);
    int prescale_8_3_p2 = triggerPrescales->getPrescaleForIndex(index_8_3_p2);
    int prescale_8_3_p3 = triggerPrescales->getPrescaleForIndex(index_8_3_p3);
    int prescale_8_3_p4 = triggerPrescales->getPrescaleForIndex(index_8_3_p4);

    int prescale_8_5_p0 = triggerPrescales->getPrescaleForIndex(index_8_5_p0);
    int prescale_8_5_p1 = triggerPrescales->getPrescaleForIndex(index_8_5_p1);
    int prescale_8_5_p2 = triggerPrescales->getPrescaleForIndex(index_8_5_p2);
    int prescale_8_5_p3 = triggerPrescales->getPrescaleForIndex(index_8_5_p3);
    int prescale_8_5_p4 = triggerPrescales->getPrescaleForIndex(index_8_5_p4);

    int prescale_8_6_p0 = triggerPrescales->getPrescaleForIndex(index_8_6_p0);
    int prescale_8_6_p1 = triggerPrescales->getPrescaleForIndex(index_8_6_p1);
    int prescale_8_6_p2 = triggerPrescales->getPrescaleForIndex(index_8_6_p2);
    int prescale_8_6_p3 = triggerPrescales->getPrescaleForIndex(index_8_6_p3);
    int prescale_8_6_p4 = triggerPrescales->getPrescaleForIndex(index_8_6_p4);

    int prescale_9_4_p0 = triggerPrescales->getPrescaleForIndex(index_9_4_p0);
    int prescale_9_4_p1 = triggerPrescales->getPrescaleForIndex(index_9_4_p1);
    int prescale_9_4_p2 = triggerPrescales->getPrescaleForIndex(index_9_4_p2);
    int prescale_9_4_p3 = triggerPrescales->getPrescaleForIndex(index_9_4_p3);
    int prescale_9_4_p4 = triggerPrescales->getPrescaleForIndex(index_9_4_p4);

    int prescale_9_5_p0 = triggerPrescales->getPrescaleForIndex(index_9_5_p0);
    int prescale_9_5_p1 = triggerPrescales->getPrescaleForIndex(index_9_5_p1);
    int prescale_9_5_p2 = triggerPrescales->getPrescaleForIndex(index_9_5_p2);
    int prescale_9_5_p3 = triggerPrescales->getPrescaleForIndex(index_9_5_p3);
    int prescale_9_5_p4 = triggerPrescales->getPrescaleForIndex(index_9_5_p4);

    int prescale_9_6_p0 = triggerPrescales->getPrescaleForIndex(index_9_6_p0);
    int prescale_9_6_p1 = triggerPrescales->getPrescaleForIndex(index_9_6_p1);
    int prescale_9_6_p2 = triggerPrescales->getPrescaleForIndex(index_9_6_p2);
    int prescale_9_6_p3 = triggerPrescales->getPrescaleForIndex(index_9_6_p3);
    int prescale_9_6_p4 = triggerPrescales->getPrescaleForIndex(index_9_6_p4);

    int prescale_12_6_p0 = triggerPrescales->getPrescaleForIndex(index_12_6_p0);
    int prescale_12_6_p1 = triggerPrescales->getPrescaleForIndex(index_12_6_p1);
    int prescale_12_6_p2 = triggerPrescales->getPrescaleForIndex(index_12_6_p2);
    int prescale_12_6_p3 = triggerPrescales->getPrescaleForIndex(index_12_6_p3);
    int prescale_12_6_p4 = triggerPrescales->getPrescaleForIndex(index_12_6_p4);


    int prescale_7_4    = 5*(prescale_7_4_p0  || prescale_7_4_p1  || prescale_7_4_p2  || prescale_7_4_p3  || prescale_7_4_p4 );
    int prescale_8_3    = 5*(prescale_8_3_p0  || prescale_8_3_p1  || prescale_8_3_p2  || prescale_8_3_p3  || prescale_8_3_p4 );
    int prescale_8_5    = 5*(prescale_8_5_p0  || prescale_8_5_p1  || prescale_8_5_p2  || prescale_8_5_p3  || prescale_8_5_p4 );
    int prescale_8_6    = 5*(prescale_8_6_p0  || prescale_8_6_p1  || prescale_8_6_p2  || prescale_8_6_p3  || prescale_8_6_p4 );
    int prescale_9_4    = 5*(prescale_9_4_p0  || prescale_9_4_p1  || prescale_9_4_p2  || prescale_9_4_p3  || prescale_9_4_p4 );
    int prescale_9_5    = 5*(prescale_9_5_p0  || prescale_9_5_p1  || prescale_9_5_p2  || prescale_9_5_p3  || prescale_9_5_p4 );
    int prescale_9_6    = 5*(prescale_9_6_p0  || prescale_9_6_p1  || prescale_9_6_p2  || prescale_9_6_p3  || prescale_9_6_p4 );
    int prescale_12_6   = 5*(prescale_12_6_p0 || prescale_12_6_p1 || prescale_12_6_p2 || prescale_12_6_p3 || prescale_12_6_p4 );

    passesHLT++;

    //////////////////////////////////////////////////////////////////////
    // Make sure that you can find a pat muon matching the HLT object   //
    // and take the best match in dR as candidate                       //
    //////////////////////////////////////////////////////////////////////
  
    // count the nr of pat::muons which fired the trigger -> reset to 0 in every event
    int counter = 0;
  
    //std::cout << "#muons in this event " << muons->size() << std::endl;

    for (unsigned int muIdx=0; muIdx<muons->size(); ++muIdx){
      //if (iEvent.id().event() != 31516) continue;  
      //if (iEvent.id().luminosityBlock() != 2007) continue;  
      //std::cout << iEvent.id().event() << std::endl;
      //std::cout << "event found! "<< std::endl;

      const pat::Muon& muon = (*muons)[muIdx];    
  
      //std::cout<<"found pat muon with pt:"<< muon.pt() << std::endl;
      // muon cuts

      if (!muSelection_(muon)) continue; 
      
      patTriggerCand++;
  
      // initialize start values
      float drMuonTrgObj   = 0.0;

      int   trgObjIdx      = -1;
      int   iTrg           = 0 ;
      int   iMatch         = 0 ;

      int   filterLabelMu7_4  = -1;
      int   filterLabelMu8_3  = -1;
      int   filterLabelMu8_5  = -1;
      int   filterLabelMu8_6  = -1;
      int   filterLabelMu9_4  = -1;
      int   filterLabelMu9_5  = -1;
      int   filterLabelMu9_6  = -1;
      int   filterLabelMu12_6 = -1;

      for(unsigned int objIdx=0; objIdx < triggerObjects->size(); ++objIdx){
  
        iTrg++;
  
        pat::TriggerObjectStandAlone trgObj = (*triggerObjects)[objIdx];      

        //unpack trigger labels
        trgObj.unpackFilterLabels(iEvent, *triggerBits);
  
        std::vector<std::string> filterLabels = trgObj.filterLabels();
    
        // check if the trigger object was actually firing our trigger

        //if (!((pass_7_4 && trgObj.hasFilterLabel(trgFilterLabelMu7_4_) ) || (pass_9_6 && trgObj.hasFilterLabel(trgFilterLabelMu9_6_)) || (pass_8_3 && trgObj.hasFilterLabel(trgFilterLabelMu8_3_)) )) continue; 
        
        const bool passTrigger =
            (pass_7_4  && trgObj.hasFilterLabel(trgFilterLabelMu7_4_))  ||
            (pass_8_3  && trgObj.hasFilterLabel(trgFilterLabelMu8_3_))  ||
            (pass_8_5  && trgObj.hasFilterLabel(trgFilterLabelMu8_5_))  ||
            (pass_8_6  && trgObj.hasFilterLabel(trgFilterLabelMu8_6_))  ||
            (pass_9_4  && trgObj.hasFilterLabel(trgFilterLabelMu9_4_))  ||
            (pass_9_5  && trgObj.hasFilterLabel(trgFilterLabelMu9_5_))  ||
            (pass_9_6  && trgObj.hasFilterLabel(trgFilterLabelMu9_6_))  ||
            (pass_12_6 && trgObj.hasFilterLabel(trgFilterLabelMu12_6_)) ;
        
        if (!passTrigger) continue;
 
        passesHLTFilter++;   
        iMatch++;

        float dr = reco::deltaR(trgObj, muon);
     
        if( ((dr < drMuonTrgObj) || (trgObjIdx < 0))  && (dr < maxdR_))
        {
              //modify dR and the reco and triggering muon index
  	      drMuonTrgObj = dr;
  	      trgObjIdx    = objIdx;
              filterLabelMu7_4    = trgObj.hasFilterLabel(trgFilterLabelMu7_4_) ;
              filterLabelMu8_3    = trgObj.hasFilterLabel(trgFilterLabelMu8_3_) ;
              filterLabelMu8_5    = trgObj.hasFilterLabel(trgFilterLabelMu8_5_) ;
              filterLabelMu8_6    = trgObj.hasFilterLabel(trgFilterLabelMu8_6_) ;
              filterLabelMu9_4    = trgObj.hasFilterLabel(trgFilterLabelMu9_4_) ;
              filterLabelMu9_5    = trgObj.hasFilterLabel(trgFilterLabelMu9_5_) ;
              filterLabelMu9_6    = trgObj.hasFilterLabel(trgFilterLabelMu9_6_) ;
              filterLabelMu12_6   = trgObj.hasFilterLabel(trgFilterLabelMu12_6_);

        }
  
      } // closing loop over trg muons
  
      //save pat muon if we found a matching candidate 
      if(trgObjIdx != -1){
     
        trgCand++;
   
        //the following line does a copy of muon with the name trgMatchedMuon (same properties different adress)
  	pat::Muon trgMatchedMuon(muon);
  
        // only take muons with valid tracks
        const reco::TransientTrack ttrackTrgMuon(*(muon.bestTrack()), paramField);  
        if (!ttrackTrgMuon.isValid()) continue;

  
        //store
  	trgMatchedMuon.addUserInt("muIdx",                muIdx    );
   	trgMatchedMuon.addUserInt("trgObjIdx",            trgObjIdx);

        trgMatchedMuon.addUserInt("mu7_ip4",              pass_7_4);

        trgMatchedMuon.addUserInt("mu7_ip4_p0",           pass_7_4_p0);
        trgMatchedMuon.addUserInt("mu7_ip4_p1",           pass_7_4_p1);
        trgMatchedMuon.addUserInt("mu7_ip4_p2",           pass_7_4_p2);
        trgMatchedMuon.addUserInt("mu7_ip4_p3",           pass_7_4_p3);
        trgMatchedMuon.addUserInt("mu7_ip4_p4",           pass_7_4_p4);

        trgMatchedMuon.addUserInt("prescale_mu7_ip4_p0",  prescale_7_4_p0);
        trgMatchedMuon.addUserInt("prescale_mu7_ip4_p1",  prescale_7_4_p1);
        trgMatchedMuon.addUserInt("prescale_mu7_ip4_p2",  prescale_7_4_p2);
        trgMatchedMuon.addUserInt("prescale_mu7_ip4_p3",  prescale_7_4_p3);
        trgMatchedMuon.addUserInt("prescale_mu7_ip4_p4",  prescale_7_4_p4);
        trgMatchedMuon.addUserInt("prescale_mu7_ip4",     prescale_7_4);

        trgMatchedMuon.addUserInt("mu8_ip3",              pass_8_3);

        trgMatchedMuon.addUserInt("mu8_ip3_p0",           pass_8_3_p0);
        trgMatchedMuon.addUserInt("mu8_ip3_p1",           pass_8_3_p1);
        trgMatchedMuon.addUserInt("mu8_ip3_p2",           pass_8_3_p2);
        trgMatchedMuon.addUserInt("mu8_ip3_p3",           pass_8_3_p3);
        trgMatchedMuon.addUserInt("mu8_ip3_p4",           pass_8_3_p4);

        trgMatchedMuon.addUserInt("prescale_mu8_ip3_p0",  prescale_8_3_p0);
        trgMatchedMuon.addUserInt("prescale_mu8_ip3_p1",  prescale_8_3_p1);
        trgMatchedMuon.addUserInt("prescale_mu8_ip3_p2",  prescale_8_3_p2);
        trgMatchedMuon.addUserInt("prescale_mu8_ip3_p3",  prescale_8_3_p3);
        trgMatchedMuon.addUserInt("prescale_mu8_ip3_p4",  prescale_8_3_p4);
        trgMatchedMuon.addUserInt("prescale_mu8_ip3",     prescale_8_3);

        trgMatchedMuon.addUserInt("mu8_ip5",              pass_8_5);

        trgMatchedMuon.addUserInt("mu8_ip5_p0",           pass_8_5_p0);
        trgMatchedMuon.addUserInt("mu8_ip5_p1",           pass_8_5_p1);
        trgMatchedMuon.addUserInt("mu8_ip5_p2",           pass_8_5_p2);
        trgMatchedMuon.addUserInt("mu8_ip5_p3",           pass_8_5_p3);
        trgMatchedMuon.addUserInt("mu8_ip5_p4",           pass_8_5_p4);

        trgMatchedMuon.addUserInt("prescale_mu8_ip5_p0",  prescale_8_5_p0);
        trgMatchedMuon.addUserInt("prescale_mu8_ip5_p1",  prescale_8_5_p1);
        trgMatchedMuon.addUserInt("prescale_mu8_ip5_p2",  prescale_8_5_p2);
        trgMatchedMuon.addUserInt("prescale_mu8_ip5_p3",  prescale_8_5_p3);
        trgMatchedMuon.addUserInt("prescale_mu8_ip5_p4",  prescale_8_5_p4);
        trgMatchedMuon.addUserInt("prescale_mu8_ip5",     prescale_8_5);

        trgMatchedMuon.addUserInt("mu8_ip6",              pass_8_6);

        trgMatchedMuon.addUserInt("mu8_ip6_p0",           pass_8_6_p0);
        trgMatchedMuon.addUserInt("mu8_ip6_p1",           pass_8_6_p1);
        trgMatchedMuon.addUserInt("mu8_ip6_p2",           pass_8_6_p2);
        trgMatchedMuon.addUserInt("mu8_ip6_p3",           pass_8_6_p3);
        trgMatchedMuon.addUserInt("mu8_ip6_p4",           pass_8_6_p4);

        trgMatchedMuon.addUserInt("prescale_mu8_ip6_p0",  prescale_8_6_p0);
        trgMatchedMuon.addUserInt("prescale_mu8_ip6_p1",  prescale_8_6_p1);
        trgMatchedMuon.addUserInt("prescale_mu8_ip6_p2",  prescale_8_6_p2);
        trgMatchedMuon.addUserInt("prescale_mu8_ip6_p3",  prescale_8_6_p3);
        trgMatchedMuon.addUserInt("prescale_mu8_ip6_p4",  prescale_8_6_p4);
        trgMatchedMuon.addUserInt("prescale_mu8_ip6",     prescale_8_6);


        trgMatchedMuon.addUserInt("mu9_ip4",              pass_9_4);
 
        trgMatchedMuon.addUserInt("mu9_ip4_p0",           pass_9_4_p0);
        trgMatchedMuon.addUserInt("mu9_ip4_p1",           pass_9_4_p1);
        trgMatchedMuon.addUserInt("mu9_ip4_p2",           pass_9_4_p2);
        trgMatchedMuon.addUserInt("mu9_ip4_p3",           pass_9_4_p3);
        trgMatchedMuon.addUserInt("mu9_ip4_p4",           pass_9_4_p4);

        trgMatchedMuon.addUserInt("prescale_mu9_ip4_p0",  prescale_9_4_p0);
        trgMatchedMuon.addUserInt("prescale_mu9_ip4_p1",  prescale_9_4_p1);
        trgMatchedMuon.addUserInt("prescale_mu9_ip4_p2",  prescale_9_4_p2);
        trgMatchedMuon.addUserInt("prescale_mu9_ip4_p3",  prescale_9_4_p3);
        trgMatchedMuon.addUserInt("prescale_mu9_ip4_p4",  prescale_9_4_p4);
        trgMatchedMuon.addUserInt("prescale_mu9_ip4"   ,  prescale_9_4);

        trgMatchedMuon.addUserInt("mu9_ip5",              pass_9_5);
 
        trgMatchedMuon.addUserInt("mu9_ip5_p0",           pass_9_5_p0);
        trgMatchedMuon.addUserInt("mu9_ip5_p1",           pass_9_5_p1);
        trgMatchedMuon.addUserInt("mu9_ip5_p2",           pass_9_5_p2);
        trgMatchedMuon.addUserInt("mu9_ip5_p3",           pass_9_5_p3);
        trgMatchedMuon.addUserInt("mu9_ip5_p4",           pass_9_5_p4);

        trgMatchedMuon.addUserInt("prescale_mu9_ip5_p0",  prescale_9_5_p0);
        trgMatchedMuon.addUserInt("prescale_mu9_ip5_p1",  prescale_9_5_p1);
        trgMatchedMuon.addUserInt("prescale_mu9_ip5_p2",  prescale_9_5_p2);
        trgMatchedMuon.addUserInt("prescale_mu9_ip5_p3",  prescale_9_5_p3);
        trgMatchedMuon.addUserInt("prescale_mu9_ip5_p4",  prescale_9_5_p4);
        trgMatchedMuon.addUserInt("prescale_mu9_ip5"   ,  prescale_9_5);

        trgMatchedMuon.addUserInt("mu9_ip6",              pass_9_6);
 
        trgMatchedMuon.addUserInt("mu9_ip6_p0",           pass_9_6_p0);
        trgMatchedMuon.addUserInt("mu9_ip6_p1",           pass_9_6_p1);
        trgMatchedMuon.addUserInt("mu9_ip6_p2",           pass_9_6_p2);
        trgMatchedMuon.addUserInt("mu9_ip6_p3",           pass_9_6_p3);
        trgMatchedMuon.addUserInt("mu9_ip6_p4",           pass_9_6_p4);

        trgMatchedMuon.addUserInt("prescale_mu9_ip6_p0",  prescale_9_6_p0);
        trgMatchedMuon.addUserInt("prescale_mu9_ip6_p1",  prescale_9_6_p1);
        trgMatchedMuon.addUserInt("prescale_mu9_ip6_p2",  prescale_9_6_p2);
        trgMatchedMuon.addUserInt("prescale_mu9_ip6_p3",  prescale_9_6_p3);
        trgMatchedMuon.addUserInt("prescale_mu9_ip6_p4",  prescale_9_6_p4);
        trgMatchedMuon.addUserInt("prescale_mu9_ip6"   ,  prescale_9_6);

        trgMatchedMuon.addUserInt("mu12_ip6",              pass_12_6);
 
        trgMatchedMuon.addUserInt("mu12_ip6_p0",           pass_12_6_p0);
        trgMatchedMuon.addUserInt("mu12_ip6_p1",           pass_12_6_p1);
        trgMatchedMuon.addUserInt("mu12_ip6_p2",           pass_12_6_p2);
        trgMatchedMuon.addUserInt("mu12_ip6_p3",           pass_12_6_p3);
        trgMatchedMuon.addUserInt("mu12_ip6_p4",           pass_12_6_p4);

        trgMatchedMuon.addUserInt("prescale_mu12_ip6_p0",  prescale_12_6_p0);
        trgMatchedMuon.addUserInt("prescale_mu12_ip6_p1",  prescale_12_6_p1);
        trgMatchedMuon.addUserInt("prescale_mu12_ip6_p2",  prescale_12_6_p2);
        trgMatchedMuon.addUserInt("prescale_mu12_ip6_p3",  prescale_12_6_p3);
        trgMatchedMuon.addUserInt("prescale_mu12_ip6_p4",  prescale_12_6_p4);
        trgMatchedMuon.addUserInt("prescale_mu12_ip6"   ,  prescale_12_6);


        trgMatchedMuon.addUserInt("filter_mu7_4"         ,  filterLabelMu7_4);
        trgMatchedMuon.addUserInt("filter_mu8_3"         ,  filterLabelMu8_3);
        trgMatchedMuon.addUserInt("filter_mu8_5"         ,  filterLabelMu8_5);
        trgMatchedMuon.addUserInt("filter_mu8_6"         ,  filterLabelMu8_6);
        trgMatchedMuon.addUserInt("filter_mu9_4"         ,  filterLabelMu9_4);
        trgMatchedMuon.addUserInt("filter_mu9_5"         ,  filterLabelMu9_5);
        trgMatchedMuon.addUserInt("filter_mu9_6"         ,  filterLabelMu9_6);
        trgMatchedMuon.addUserInt("filter_mu12_6"        ,  filterLabelMu12_6);

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
std::cout << "#Events where HLT fired and HLT Filter passed             : " << passesHLTFilter  << std::endl;
std::cout << "#Muon candidates which possibly fire HLT (unmatched)      : " << patTriggerCand  << std::endl;
std::cout << "#Muon candidates which fired HLT (matched to trg object!) : " << patTriggerCand  << std::endl;

}

DEFINE_FWK_MODULE(Trigger);

