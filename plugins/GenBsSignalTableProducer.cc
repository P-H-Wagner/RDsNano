// Event table "GenBs" (MC): the signal Bs decay of the event, as the Bmmm R(Ds) ntuplizer (tree_gen):
//   GenBs_code     0 Ds mu nu, 1 Ds* mu nu, 2 Ds tau nu, 3 Ds* tau nu, -1 none;  GenBs_n_bs
//   GenBs_pt, eta, phi, y, mass, pdgid (first copy: flavour at production), beta, gamma,
//   GenBs_lxyz, lxy, ct (production point genParticles:xyz0, else the first-copy vertex; decay point the
//   vertex of the daughters of the last copy), sv_x/y/z                     <- Bmmm bs_gen_info (gen_bs_*)
//   GenBs_ham_<p>_{pt,eta,phi,mass,pdgid}, p = dsx ds dsg lep nu taumu taunut taunum: the pre-FSR (first
//   copy) four-momenta of the semileptonic decay for Hammer                 <- Bmmm gen_hammer_leaves
// NaN where not defined. Definitions in genBmmm.h.

#include <memory>
#include <string>
#include <vector>
#include <cmath>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "genBmmm.h"

class GenBsSignalTableProducer : public edm::global::EDProducer<> {
public:
  explicit GenBsSignalTableProducer(const edm::ParameterSet& cfg)
      : genToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("src"))),
        pvToken_(consumes<math::XYZPointF>(cfg.getParameter<edm::InputTag>("genPV"))),
        name_(cfg.getParameter<std::string>("name")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    using namespace bmmmgen;
    edm::Handle<reco::GenParticleCollection> gen;
    iEvent.getByToken(genToken_, gen);
    std::vector<P> g;
    g.reserve(gen->size());
    for (const auto& p : *gen) g.push_back(&p);
    const BsSignalDecay d = findSignalDecay(g);

    auto tab = std::make_unique<nanoaod::FlatTable>(1, name_, true);
    tab->setDoc("signal Bs decay of the event and its pre-FSR decay products (Hammer), as Bmmm");
    tab->addColumnValue<int>("code", d.code, "0 Ds mu nu, 1 Ds* mu nu, 2 Ds tau nu, 3 Ds* tau nu, -1 none");
    tab->addColumnValue<int>("n_bs", d.nBs, "number of decaying Bs (last copies)");

    float pt = NAN, eta = NAN, phi = NAN, y = NAN, m = NAN, pdg = NAN, beta = NAN, gamma = NAN;
    float lxyz = NAN, lxy = NAN, ct = NAN, svx = NAN, svy = NAN, svz = NAN;
    if (d.code >= 0 && d.bs) {
      P bs = d.bs;
      pt = bs->pt(); eta = bs->eta(); phi = bs->phi();
      y = bs->energy() > std::abs(bs->pz()) ? 0.5 * std::log((bs->energy() + bs->pz()) / (bs->energy() - bs->pz())) : NAN;
      m = mass(bs); pdg = d.bsFirst->pdgId();
      beta = bs->energy() > 0 ? bs->p() / bs->energy() : NAN;
      gamma = m > 0 ? bs->energy() / m : NAN;
      if (bs->numberOfDaughters() > 0) {
        P dau = bs->daughter(0);
        svx = dau->vx(); svy = dau->vy(); svz = dau->vz();
        double pv[3] = {d.bsFirst->vx(), d.bsFirst->vy(), d.bsFirst->vz()};
        edm::Handle<math::XYZPointF> gpv;
        iEvent.getByToken(pvToken_, gpv);
        if (gpv.isValid() && std::isfinite(gpv->x()) && std::isfinite(gpv->y()) && std::isfinite(gpv->z())) {
          pv[0] = gpv->x(); pv[1] = gpv->y(); pv[2] = gpv->z();
        }
        const double dx = svx - pv[0], dy = svy - pv[1], dz = svz - pv[2];
        lxyz = std::sqrt(dx * dx + dy * dy + dz * dz);
        lxy = std::hypot(dx, dy);
        ct = bs->p() > 0 ? lxyz * m / bs->p() : NAN;
      }
    }
    tab->addColumnValue<float>("pt", pt, "Bs pT (last copy)");
    tab->addColumnValue<float>("eta", eta, "Bs eta");
    tab->addColumnValue<float>("phi", phi, "Bs phi");
    tab->addColumnValue<float>("y", y, "Bs rapidity");
    tab->addColumnValue<float>("mass", m, "Bs mass");
    tab->addColumnValue<float>("pdgid", pdg, "pdgId of the first copy (flavour at production)");
    tab->addColumnValue<float>("beta", beta, "Bs beta");
    tab->addColumnValue<float>("gamma", gamma, "Bs gamma");
    tab->addColumnValue<float>("lxyz", lxyz, "3D flight distance [cm]");
    tab->addColumnValue<float>("lxy", lxy, "transverse flight distance [cm]");
    tab->addColumnValue<float>("ct", ct, "proper decay length c*t [cm] (lifetime reweighting)");
    tab->addColumnValue<float>("sv_x", svx, "decay vertex x [cm]");
    tab->addColumnValue<float>("sv_y", svy, "decay vertex y [cm]");
    tab->addColumnValue<float>("sv_z", svz, "decay vertex z [cm]");

    // Hammer leaves (first copies = momenta at the decay vertex, before FSR)
    std::vector<std::pair<std::string, P>> ham = {{"dsx", nullptr}, {"ds", nullptr}, {"dsg", nullptr}, {"lep", nullptr},
                                                  {"nu", nullptr}, {"taumu", nullptr}, {"taunut", nullptr}, {"taunum", nullptr}};
    if (d.code >= 0) {
      ham[0].second = d.dsx;
      if (d.dsg) { ham[1].second = d.ds; ham[2].second = d.dsg; }
      ham[3].second = d.lep;
      ham[4].second = d.nu;
      if (d.tauNuT) { ham[5].second = d.mu; ham[6].second = d.tauNuT; ham[7].second = d.tauNuM; }
    }
    for (const auto& [n, p] : ham) {
      tab->addColumnValue<float>("ham_" + n + "_pt",    p ? p->pt() : NAN,     "Hammer " + n + " pT (pre-FSR)");
      tab->addColumnValue<float>("ham_" + n + "_eta",   p ? p->eta() : NAN,    "Hammer " + n + " eta");
      tab->addColumnValue<float>("ham_" + n + "_phi",   p ? p->phi() : NAN,    "Hammer " + n + " phi");
      tab->addColumnValue<float>("ham_" + n + "_mass",  p ? mass(p) : NAN,     "Hammer " + n + " mass");
      tab->addColumnValue<float>("ham_" + n + "_pdgid", p ? p->pdgId() : NAN,  "Hammer " + n + " pdgId");
    }
    iEvent.put(std::move(tab));
  }

private:
  const edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  const edm::EDGetTokenT<math::XYZPointF> pvToken_;
  const std::string name_;
};

DEFINE_FWK_MODULE(GenBsSignalTableProducer);
