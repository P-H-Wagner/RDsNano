
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class notEmpty : public edm::EDFilter {
public:
    explicit notEmpty(const edm::ParameterSet&);
    ~notEmpty() override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    bool filter(edm::Event&, const edm::EventSetup&) override;

    const edm::InputTag bsTag;
    const edm::EDGetTokenT<pat::CompositeCandidateCollection> bs_;
};

notEmpty::notEmpty(const edm::ParameterSet& iConfig): 
  bsTag(iConfig.getParameter<edm::InputTag>("bs")),
  bs_(consumes<pat::CompositeCandidateCollection>(bsTag)) {}

notEmpty::~notEmpty() {
}

bool notEmpty::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<pat::CompositeCandidateCollection> myCandidates;
    iEvent.getByToken(bs_, myCandidates);

    for (const auto& myCandidate : *myCandidates) {
        // Access the 'arrived' value for the current entry
        const int arrivedValue = myCandidate.userInt("arrived");

        // Check the condition: Only pass events where arrived == 1
        if (arrivedValue >0) {
            return true;  // Event passes the filter
        }
    }

    return false;  // Event does not pass the filter
}

void notEmpty::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // Define your parameters here
}

DEFINE_FWK_MODULE(notEmpty);
