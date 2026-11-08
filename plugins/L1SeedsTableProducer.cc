// Event table "L1" with the decision and the prescale of the L1 seeds of the 2018 B-parking paths:
//   L1_<seed>     final decision from gtStage2Digis (bit looked up in the L1 menu of the run);
//                 0 if the seed is not in that menu, -1 if gtStage2Digis is absent
//   L1_<seed>_ps  prescale of the lumi section from the brilcalc tables (MC: 1); the last lumi section of
//                 the run when beyond the table, 0 if the run is not in the table, NaN if the lumi section
//                 is missing inside the run
// Same logic and inputs as Bmmm/Analysis/python/L1Seeds2018.py; the inputs are converted once into
// data/l1_seeds_bparking2018.json by production/make_l1_payload.py.

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <nlohmann/json.hpp>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class L1SeedsTableProducer : public edm::global::EDProducer<> {
public:
  explicit L1SeedsTableProducer(const edm::ParameterSet& cfg)
      : algToken_(consumes<BXVector<GlobalAlgBlk>>(cfg.getParameter<edm::InputTag>("src"))),
        name_(cfg.getParameter<std::string>("name")) {
    std::ifstream fin(edm::FileInPath(cfg.getParameter<std::string>("payload")).fullPath());
    const nlohmann::json j = nlohmann::json::parse(fin);
    seeds_ = j.at("seeds").get<std::vector<std::string>>();
    for (const auto& s : seeds_)
      for (const auto& [run, segs] : j.at("prescales").at(s).items()) {
        auto& v = ps_[s][std::stoi(run)];
        for (const auto& g : segs) v.push_back({g[0].get<int>(), g[1].get<int>(), g[2].get<double>()});
      }
    for (const auto& [run, menu] : j.at("run_menu").items()) runMenu_[std::stoi(run)] = menu.get<std::string>();
    for (const auto& [menu, bits] : j.at("menus").items())
      for (const auto& [s, bit] : bits.items()) menus_[menu][s] = bit.get<int>();
    mcMenu_ = j.at("mc_menu").get<std::string>();
    runMin_ = runMenu_.begin()->first;
    runMax_ = runMenu_.rbegin()->first;
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    const bool mc = !iEvent.isRealData();
    const int run = iEvent.id().run();
    const int ls = iEvent.id().luminosityBlock();

    // menu of the run (MC: the MC menu); outside the 2018 map: none (decisions 0), as Bmmm
    const std::map<std::string, int>* menu = nullptr;
    std::string menuName;
    if (mc) menuName = mcMenu_;
    else if (run >= runMin_ && run <= runMax_) {
      auto it = runMenu_.find(run);
      if (it == runMenu_.end()) {   // closest run of the map
        auto hi = runMenu_.lower_bound(run);
        auto lo = hi == runMenu_.begin() ? hi : std::prev(hi);
        it = (hi == runMenu_.end() || std::abs(lo->first - run) <= std::abs(hi->first - run)) ? lo : hi;
      }
      menuName = it->second;
    }
    auto im = menus_.find(menuName);
    if (im != menus_.end()) menu = &im->second;

    edm::Handle<BXVector<GlobalAlgBlk>> alg;
    iEvent.getByToken(algToken_, alg);
    const GlobalAlgBlk* blk = (alg.isValid() && alg->size(0) > 0) ? &alg->at(0, 0) : nullptr;

    auto tab = std::make_unique<nanoaod::FlatTable>(1, name_, true);
    tab->setDoc("L1 seeds of the 2018 B-parking paths: final decision and brilcalc prescale (MC: 1); as Bmmm");
    for (const auto& s : seeds_) {
      int dec = -1;
      if (blk) {
        dec = 0;
        if (menu) {
          auto ib = menu->find(s);
          if (ib != menu->end()) dec = blk->getAlgoDecisionFinal(ib->second) ? 1 : 0;
        }
      }
      float ps = 1.f;
      if (!mc) {
        ps = 0.f;
        auto is = ps_.find(s);
        if (is != ps_.end()) {
          auto ir = is->second.find(run);
          if (ir != is->second.end()) {
            ps = NAN;
            const auto& segs = ir->second;
            for (const auto& g : segs)
              if (ls >= g.first && ls <= g.last) { ps = g.ps; break; }
            if (std::isnan(ps) && !segs.empty() && ls > segs.back().last) ps = segs.back().ps;
          }
        }
      }
      const std::string col = s.rfind("L1_", 0) == 0 ? s.substr(3) : s;
      tab->addColumnValue<int>(col, dec, "final decision of " + s + " (1/0; -1: no gtStage2Digis)");
      tab->addColumnValue<float>(col + "_ps", ps, "prescale of " + s + " in this lumi section (MC: 1)");
    }
    iEvent.put(std::move(tab));
  }

private:
  struct Seg { int first, last; double ps; };
  const edm::EDGetTokenT<BXVector<GlobalAlgBlk>> algToken_;
  const std::string name_;
  std::vector<std::string> seeds_;
  std::map<std::string, std::map<int, std::vector<Seg>>> ps_;
  std::map<int, std::string> runMenu_;
  std::map<std::string, std::map<std::string, int>> menus_;
  std::string mcMenu_;
  int runMin_ = 0, runMax_ = 0;
};

DEFINE_FWK_MODULE(L1SeedsTableProducer);
