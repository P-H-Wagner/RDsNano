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

class Trigger : public edm::EDProducer {

public:

    //define the transienttrackcollection type, does not exist!
    //typedef std::vector<reco::TransientTrack> TransientTrackCollection;
    
    //constructor, takes a reference to the ParameterSet 
    explicit Trigger(const edm::ParameterSet&);
    //destructor
    ~Trigger() override {};
   
    int counter = 0;

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
    const string trgFilterLabel_;
    const double maxdR_; 
    //for filter wrt trigger ????
    //const double dzTrg_cleaning_; 

    //muon selection
    //const double ptMin_;          // min pT in all muons for B candidates
    //const double absEtaMax_;      //max eta ""
    //const bool   softMuonsOnly_;    //cuts muons without soft ID
};

//constructor definition (outsde class)
Trigger::Trigger(const edm::ParameterSet& iConfig):
 
  //for the muons
  muSelection_(iConfig.getParameter<std::string>("muSelection")), 
  muonTag(iConfig.getParameter<edm::InputTag>("muonCollection")),
  muonSrc_(consumes<std::vector<pat::Muon>>(muonTag)),
  //for trigger info
  triggerBitTag(iConfig.getParameter<edm::InputTag>("trgResultsCollection")),
  triggerBits_(consumes<edm::TriggerResults>(triggerBitTag)),

  triggerObjectTag(iConfig.getParameter<edm::InputTag>("trgObjectsCollection")),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerObjectTag)),

  triggerPrescaleTag(iConfig.getParameter<edm::InputTag>("trgPrescaleCollection")),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(triggerPrescaleTag)),
  //vertex info
  vertexSrcTag(iConfig.getParameter<edm::InputTag>("vtxCollection")),
  vertexSrc_(consumes<reco::VertexCollection>(vertexSrcTag)), 

  //parameters

  trgFilterLabel_(iConfig.getParameter<string>("trgFilterLabel")),
  maxdR_(iConfig.getParameter<double>("maxdR_matching"))
  //dzTrg_cleaning_(iConfig.getParameter<double>("dzForCleaning_wrtTrgMuon")),
  //ptMin_(iConfig.getParameter<double>("ptMin")),
  //absEtaMax_(iConfig.getParameter<double>("absEtaMax")),
  // softMuonsOnly_(iConfig.getParameter<bool>("softMuonsOnly"))
{
  // produce 2 collections: 
  produces<pat::MuonCollection>("trgMuons");  
  produces<TransientTrackCollection>("ttracksTrgMuons"); // tracks of selected muons 
}



void Trigger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
 
  //Define handles
  //edm::ESHandle<MagneticField> bFieldHandle;
  //iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexSrc_, vertexHandle);

  //is used to store information about the results of trigger decisions in a given event
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  //is used to store information about the results of trigger decisions in a given event
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  // for every event, get the name list of the triggers
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  ///////////// Debugging:
  //std::cout << "Available trigger names:" << std::endl;
  //for (unsigned int i = 0; i < names.size(); ++i) {
  //  if("HLT_Mu7_IP4" in names.triggerName(i)){
  //  std::cout << names.triggerName(i) << std::endl;
  //}
  //}
  /////////////
  
  //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  //pat muons contain more info than reco muons, f.e. trigger info! 
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonSrc_, muons);

  // to save
  std::unique_ptr<pat::MuonCollection>      trgMuons   ( new pat::MuonCollection );
  std::unique_ptr<TransientTrackCollection> ttracksTrgMuons( new TransientTrackCollection);
 
  // Getting the indices of the HLT paths
  unsigned int index_7_4_p0      = names.triggerIndex("HLT_Mu7_IP4_part0_v2"); 
  unsigned int index_7_4_p1      = names.triggerIndex("HLT_Mu7_IP4_part1_v2"); 
  unsigned int index_7_4_p2      = names.triggerIndex("HLT_Mu7_IP4_part2_v2"); 
  unsigned int index_7_4_p3      = names.triggerIndex("HLT_Mu7_IP4_part3_v2"); 
  unsigned int index_7_4_p4      = names.triggerIndex("HLT_Mu7_IP4_part4_v2"); 
  //unsigned int index_8_3      = names.triggerIndex("HLT_Mu8_IP3");
  //unsigned int index_8_5      = names.triggerIndex("HLT_Mu8_IP5");
  //unsigned int index_8_6      = names.triggerIndex("HLT_Mu8_IP6");
  //unsigned int index_8p5_3p5  = names.triggerIndex("HLT_Mu8p5_IP3p5");
  //unsigned int index_9_4      = names.triggerIndex("HLT_Mu9_IP4");
  //unsigned int index_9_5      = names.triggerIndex("HLT_Mu9_IP5");
  //unsigned int index_10p5_3p5 = names.triggerIndex("HLT_Mu10p5_IP3p5");
  //unsigned int index_12_6     = names.triggerIndex("HLT_Mu12_IP6");

  //std::cout << (index_7_4_p0 < triggerBits->size())  << "and" <<      (triggerBits->accept(index_7_4_p0)) << std::endl;
  //std::cout << (index_7_4_p1 < triggerBits->size())  << "and" <<     (triggerBits->accept(index_7_4_p1)) << std::endl;
  //std::cout << (index_7_4_p2 < triggerBits->size())  << "and" <<    (triggerBits->accept(index_7_4_p2)) << std::endl;
  //std::cout << (index_7_4_p3 < triggerBits->size())  << "and" <<     (triggerBits->accept(index_7_4_p3)) << std::endl;
  //std::cout << (index_7_4_p4 < triggerBits->size())  << "and" <<    (triggerBits->accept(index_7_4_p4)) << std::endl;

  //default is false  
  bool pass_7_4_p0_path      = false;
  bool pass_7_4_p1_path      = false;
  bool pass_7_4_p2_path      = false;
  bool pass_7_4_p3_path      = false;
  bool pass_7_4_p4_path      = false;

  //bool pass_8_3_path      = false;
  //bool pass_8_5_path      = false;
  //bool pass_8_6_path      = false;
  //bool pass_8p5_3p5_path  = false;
  //bool pass_9_4_path      = false;
  //bool pass_9_5_path      = false;
  //bool pass_10p5_3p5_path = false;
  //bool pass_12_6_path     = false;

  // check first if the index is valid, i.e. if it is not out of range (maximum is givrn by triggerBits->size())
  // and if so, check if the trigger has fired with accept() 

  pass_7_4_p0_path      = ((index_7_4_p0 < triggerBits->size())      && (triggerBits->accept(index_7_4_p0)));
  pass_7_4_p1_path      = ((index_7_4_p1 < triggerBits->size())      && (triggerBits->accept(index_7_4_p1)));
  pass_7_4_p2_path      = ((index_7_4_p2 < triggerBits->size())      && (triggerBits->accept(index_7_4_p2)));
  pass_7_4_p3_path      = ((index_7_4_p3 < triggerBits->size())      && (triggerBits->accept(index_7_4_p3)));
  pass_7_4_p4_path      = ((index_7_4_p4 < triggerBits->size())      && (triggerBits->accept(index_7_4_p4)));
  //pass_8_3_path      = ((index_8_3 < triggerBits->size())      && (triggerBits->accept(index_8_3)));
  //pass_8_5_path      = ((index_8_5 < triggerBits->size())      && (triggerBits->accept(index_8_5)));
  //pass_8_6_path      = ((index_8_6 < triggerBits->size())      && (triggerBits->accept(index_8_6)));
  //pass_8p5_3p5_path  = ((index_8p5_3p5 < triggerBits->size())  && (triggerBits->accept(index_8p5_3p5)));
  //pass_9_4_path      = ((index_9_4 < triggerBits->size())      && (triggerBits->accept(index_9_4)));
  //pass_9_5_path      = ((index_9_5 < triggerBits->size())      && (triggerBits->accept(index_9_5)));
  //pass_10p5_3p5_path = ((index_10p5_3p5 < triggerBits->size()) && (triggerBits->accept(index_10p5_3p5)));
  //pass_12_6_path     = ((index_12_6 < triggerBits->size())     && (triggerBits->accept(index_12_6)));


  //define vector out of bools
  std::vector<bool> trgMuonFrom_7_4_flag;

  //only continue when we the event passes the HLT_Mu7_IP4
  if (pass_7_4_p0_path || pass_7_4_p1_path || pass_7_4_p2_path || pass_7_4_p3_path || pass_7_4_p4_path){
  //std::cout<<"found trigger!" << std::endl;
  // define vectors of ints of length muons->size(), all values set to 0
  std::vector<int> isTriggerMuon(muons->size(), 0);

  //////////////////////////////////////////////////////////////////////
  // Make sure that you can find a pat  muon matching the HLT object  //
  //////////////////////////////////////////////////////////////////////

  // count the nr of pat::muons which fired the trigger -> reset to 0 in every event
  int counter = 0;


  // now loop over all pat::muons
  for (unsigned int muIdx=0; muIdx<muons->size(); ++muIdx){
    if (iEvent.id().event() != 121643971) continue;  
    //if (iEvent.id().luminosityBlock() != 103) continue;  
    //std::cout << iEvent.id().event() << std::endl;
    //access the muon at the muIdx-position
    const pat::Muon& muon = (*muons)[muIdx];    

    std::cout<<"found pat muon with pt:"<< muon.pt() << std::endl;
    // muon cuts
    if (!muSelection_(muon)) continue; 

    //check if the pat muon is matched to some trigger object (by using the function triggerObjectMatchByPath()
    //bool isMatched = !(muon.triggerObjectMatchByPath("HLT_Mu7_IP4")==nullptr);
    std::cout<< "muon passed the muon selection!" << std::endl;

    // initialize start values
    float drMuonTrgObj = -1.;
    int muonIdx        = -1;
    int trgObjIdx      = -1;

    
    //////////////////////////////////////////////////////////////////////
    // Find the (pat muon, HLT)  pair with the smallest dR and save it  //
    //////////////////////////////////////////////////////////////////////
    int iTrg = 0;
    int iMatch = 0;
    for(unsigned int objIdx=0; objIdx < triggerObjects->size(); ++objIdx){

      iTrg++;

      pat::TriggerObjectStandAlone trgObj = (*triggerObjects)[objIdx];      
      //unpack trigger labels

      trgObj.unpackFilterLabels(iEvent, *triggerBits);

      std::vector<std::string> filterLabels = trgObj.filterLabels();
  
      //check if the triggermuon was actually firing the trigger
      if(!trgObj.hasFilterLabel(trgFilterLabel_)) continue;
    
      iMatch++;
      //std::cout<< "we have a trg object, is it matching?!!" << std::endl;
      //std::cout << "Filter labels for trigger object:" << std::endl;
      //for (const std::string& label : filterLabels) {
      //  std::cout << "- " << label << std::endl;
      //} 

      //save the dR between the triggering muon and the pat muon 
      float dr = reco::deltaR(trgObj, muon);
   
      //the first element of the loop can enter bc we allow drMuon.. = -1
      //the others can only enter if their dR is smaller (and they match the trigger)
      //---> we filter out the one (muon,trigger muon) pair, which has the smallest dR
      
      if( ( (dr < drMuonTrgObj) || (drMuonTrgObj == -1))  && (dr < maxdR_))
      {
              //modify dR and the reco and triggering muon index
	      drMuonTrgObj = dr;
	      muonIdx      = muIdx;
	      trgObjIdx    = objIdx;
      }

    } // closing loop over trg muons

    //save pat muon if we found a matching candidate 
    if(muonIdx != -1)
     {
            //the following line does a copy of muon with the name trgMatchedMuon (same properties different adress)
	    pat::Muon trgMatchedMuon(muon);

            std::cout<< "we found a trigger match!" << std::endl;
            //save also the tracks
            const reco::TransientTrack ttrackTrgMuon(*(muon.bestTrack()), paramField);  
            if (!ttrackTrgMuon.isValid()) continue;

            //add the index of the matching trigger muon in the event
	    trgMatchedMuon.addUserInt("muonIdx", muonIdx);
 
            //add the index of the original muon in the event
 	    trgMatchedMuon.addUserInt("trgObjIdx", trgObjIdx);
 	    trgMatchedMuon.addUserInt("trgPrescale", triggerPrescales->getPrescaleForIndex(trgObjIdx));

	    // append the pat::muon which could be matched to the output
            trgMuons->emplace_back(trgMatchedMuon); 
            ttracksTrgMuons->emplace_back(ttrackTrgMuon);
            //std::cout<<" Yes, it is matched and saved :)!" << std::endl;

            counter++;

      }
      else {std::cout << "did not find a trigger match" << std::endl;}
  } //close the loop over pat muons
  
  //std::cout <<  "wehave:" << counter << std::endl;   

  } //close if condition of HLT
  //store the new collections in the event
  iEvent.put(std::move(trgMuons),    "trgMuons");
  iEvent.put(std::move(ttracksTrgMuons), "ttracksTrgMuons");

//close the if condition for the HLT_Mu7_IP4
}//close produce function

DEFINE_FWK_MODULE(Trigger);

