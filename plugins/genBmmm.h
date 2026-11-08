#ifndef PhysicsTools_RDsNano_genBmmm_h
#define PhysicsTools_RDsNano_genBmmm_h
// Generator truth of the Bmmm R(Ds) ntuplizer (Bmmm/Analysis/python/RDsGenHistory.py), ported with the same
// definitions so that the two ntuples carry identical truth:
//   classifyCandidate  <- classify_candidate / _bs_code   (per candidate: gen_cand, gen_ds_b_ct, ...)
//   BsSignalDecay      <- BsSignalDecay.from_genparticles (per event: the signal Bs decay)
//   hammer leaves      <- gen_hammer_leaves                (pre-FSR four-momenta of the semileptonic decay)
//   bsGenInfo          <- bs_gen_info                      (Bs kinematics and flight, lifetime reweighting)

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace bmmmgen {

constexpr int BS = 531, B0 = 511, BP = 521, LB = 5122, DS = 431, DSST = 433, PHI = 333;
constexpr int KAON = 321, PION = 211, MU = 13, TAU = 15, GAMMA = 22, PI0 = 111, NU_MU = 14, NU_TAU = 16;

typedef const reco::Candidate* P;

inline bool isB(int pdg) {
  const int a = std::abs(pdg) % 10000;
  if (a >= 1000) return (a / 1000) % 10 == 5;
  return (a / 100) % 10 == 5 && (a / 10) % 10 != 5;
}
inline bool isCharm(int pdg) {
  const int a = std::abs(pdg) % 10000;
  if (a >= 1000) return (a / 1000) % 10 == 4;
  return (a / 100) % 10 == 4 && (a / 10) % 10 != 4;
}
// identity of a gen particle (pruned and packed copies of the same particle compare equal), as Bmmm key()
inline bool same(P a, P b) {
  if (!a || !b) return false;
  if (a == b) return true;
  return a->pdgId() == b->pdgId() && a->px() == b->px() && a->py() == b->py() && a->pz() == b->pz() &&
         a->vx() == b->vx() && a->vy() == b->vy() && a->vz() == b->vz();
}
inline P firstCopy(P p) {
  for (int i = 0; i < 200 && p; ++i) {
    if (p->numberOfMothers() == 0) return p;
    P m = p->mother(0);
    if (!m || std::abs(m->pdgId()) != std::abs(p->pdgId())) return p;
    p = m;
  }
  return p;
}
inline P lastCopy(P p) {
  for (int i = 0; i < 200 && p; ++i) {
    P nxt = nullptr;
    for (size_t d = 0; d < p->numberOfDaughters(); ++d)
      if (std::abs(p->daughter(d)->pdgId()) == std::abs(p->pdgId())) { nxt = p->daughter(d); break; }
    if (!nxt) return p;
    p = nxt;
  }
  return p;
}
inline P motherOf(P p) {
  P fc = firstCopy(p);
  return (fc && fc->numberOfMothers() > 0) ? fc->mother(0) : nullptr;
}
inline P bAncestor(P p) {
  P node = p;
  for (int i = 0; i < 200 && node; ++i) {
    if (node->numberOfMothers() == 0) return nullptr;
    node = node->mother(0);
    if (node && isB(node->pdgId())) return lastCopy(node);
  }
  return nullptr;
}
inline std::vector<P> directDaughters(P p, bool keepPhotons = false) {
  std::vector<P> out;
  P lc = lastCopy(p);
  for (size_t i = 0; i < lc->numberOfDaughters(); ++i) {
    P d = lc->daughter(i);
    if (!keepPhotons && std::abs(d->pdgId()) == GAMMA) continue;
    out.push_back(d);
  }
  return out;
}
inline std::vector<int> signature(P p) {
  std::vector<int> s;
  for (P d : directDaughters(p)) s.push_back(std::abs(d->pdgId()));
  std::sort(s.begin(), s.end());
  return s;
}
// 0 Ds mu nu, 1 Ds* mu nu, 2 Ds tau nu, 3 Ds* tau nu, -1 otherwise
inline int signalCode(const std::vector<int>& s) {
  if (s == std::vector<int>{MU, NU_MU, DS})   return 0;
  if (s == std::vector<int>{MU, NU_MU, DSST}) return 1;
  if (s == std::vector<int>{TAU, NU_TAU, DS}) return 2;
  if (s == std::vector<int>{TAU, NU_TAU, DSST}) return 3;
  return -1;
}
inline double mass(P p) { return std::sqrt(std::max(p->energy() * p->energy() - p->p() * p->p(), 0.)); }
// c*t [cm] of a b hadron: L m / |p|, L from the production vertex of its first copy to the vertex of the
// daughters of its last copy (B mixing included)
inline float bCt(P b) {
  if (!b) return NAN;
  P fc = firstCopy(b), lc = lastCopy(b);
  if (lc->numberOfDaughters() == 0) return NAN;
  P d = lc->daughter(0);
  const double L = std::sqrt(std::pow(d->vx() - fc->vx(), 2) + std::pow(d->vy() - fc->vy(), 2) + std::pow(d->vz() - fc->vz(), 2));
  const double p = lc->p();
  return p > 0 ? L * mass(lc) / p : NAN;
}
// closest charm hadron (at its last copy) that is an ancestor of all of gens
inline P commonCharmAncestor(const std::vector<P>& gens) {
  std::vector<std::vector<P>> chains;
  for (P g : gens) {
    std::vector<P> chain;
    P node = g;
    for (int i = 0; i < 200 && node && node->numberOfMothers() > 0; ++i) {
      node = node->mother(0);
      if (node && isCharm(node->pdgId())) chain.push_back(lastCopy(node));
    }
    if (chain.empty()) return nullptr;
    chains.push_back(chain);
  }
  P node = gens[0];
  for (int i = 0; i < 200 && node && node->numberOfMothers() > 0; ++i) {
    node = node->mother(0);
    if (!node || !isCharm(node->pdgId())) continue;
    P k = lastCopy(node);
    bool all = true;
    for (size_t c = 1; c < chains.size() && all; ++c)
      all = std::any_of(chains[c].begin(), chains[c].end(), [&](P x) { return same(x, k); });
    if (all) return k;
  }
  return nullptr;
}
// closest same-charge gen particle within drMax
inline P matchTrack(double eta, double phi, int q, const std::vector<P>& gen, double drMax) {
  P best = nullptr;
  double best2 = drMax * drMax;
  for (P g : gen) {
    if (g->charge() != q) continue;
    const double d2 = reco::deltaR2(eta, phi, g->eta(), g->phi());
    if (d2 < best2) { best = g; best2 = d2; }
  }
  return best;
}

inline int bsCode(P bs, P ds, P gmu, P muMother) {
  bs = lastCopy(bs);
  P lep = (muMother && same(lastCopy(muMother), bs)) ? gmu : nullptr;
  if (!lep && muMother && std::abs(muMother->pdgId()) == TAU) {
    P tm = motherOf(muMother);
    if (tm && same(lastCopy(tm), bs)) lep = muMother;
  }
  if (!lep) return 11;
  const int code = signalCode(signature(bs));
  if (code < 0) return 10;
  P dsm = motherOf(ds);
  if (dsm && same(lastCopy(dsm), bs))
    return (std::abs(firstCopy(ds)->pdgId()) == DS && (code == 0 || code == 2)) ? code : 10;
  if (dsm && std::abs(dsm->pdgId()) == DSST) {
    P dsstm = motherOf(dsm);
    if (dsstm && same(lastCopy(dsstm), bs) && (code == 1 || code == 3)) return code;
  }
  return 10;
}

struct CandTruth {
  float cand = NAN, ds_true = NAN, charm_pdgid = NAN, ds_mother_pdgid = NAN, ds_b_pdgid = NAN, mu_b_pdgid = NAN,
        mu_mother_pdgid = NAN, mu_from_tau = NAN, same_b = NAN, n_matched = NAN, ds_b_ct = NAN, mu_b_ct = NAN;
};
// gk1, gk2, gpi, gmu: matched gen particles (nullptr when unmatched)
inline CandTruth classifyCandidate(P gk1, P gk2, P gpi, P gmu) {
  CandTruth o;
  int nm = 0;
  for (P g : {gk1, gk2, gpi, gmu}) nm += g != nullptr;
  o.n_matched = nm;
  P charm = (gk1 && gk2 && gpi) ? commonCharmAncestor({gk1, gk2, gpi}) : nullptr;
  P ds = (charm && std::abs(charm->pdgId()) == DS) ? charm : nullptr;
  if (charm) o.charm_pdgid = charm->pdgId();
  o.ds_true = ds != nullptr;
  if (ds) { P m = motherOf(ds); o.ds_mother_pdgid = m ? m->pdgId() : 0; }
  P bDs = ds ? bAncestor(ds) : nullptr;
  if (bDs) { o.ds_b_pdgid = firstCopy(bDs)->pdgId(); o.ds_b_ct = bCt(bDs); }
  P bMu = nullptr, muMother = nullptr;
  if (gmu) {
    muMother = motherOf(gmu);
    o.mu_mother_pdgid = muMother ? muMother->pdgId() : 0;
    o.mu_from_tau = muMother && std::abs(muMother->pdgId()) == TAU;
    bMu = bAncestor(gmu);
    if (bMu) { o.mu_b_pdgid = firstCopy(bMu)->pdgId(); o.mu_b_ct = bCt(bMu); }
  }
  const bool sameB = bDs && bMu && same(bDs, bMu);
  if (bDs && gmu) o.same_b = sameB;
  int code;
  if (nm < 4 || !gmu || std::abs(gmu->pdgId()) != MU) code = -1;
  else if (!ds) code = 30;
  else if (!sameB) code = bMu ? 20 : 21;
  else {
    const int b = std::abs(bDs->pdgId());
    code = b == BS ? bsCode(bDs, ds, gmu, muMother) : b == B0 ? 12 : b == BP ? 13 : b == LB ? 14 : 15;
  }
  o.cand = code;
  return o;
}

// ---------------------------------------------------------------- event level: the signal Bs decay
struct BsSignalDecay {
  int code = -1, nBs = 0;
  P bs = nullptr, bsFirst = nullptr, dsx = nullptr, ds = nullptr, dsg = nullptr, lep = nullptr, nu = nullptr,
    mu = nullptr, tauNuT = nullptr, tauNuM = nullptr, phi = nullptr;
};
inline BsSignalDecay findSignalDecay(const std::vector<P>& gen) {
  BsSignalDecay out;
  std::vector<P> cands;
  for (P p : gen) {
    if (std::abs(p->pdgId()) != BS || p->numberOfDaughters() == 0) continue;
    P lc = lastCopy(p);
    if (std::none_of(cands.begin(), cands.end(), [&](P x) { return same(x, lc); })) cands.push_back(lc);
  }
  out.nBs = cands.size();
  int bestRank = -1;
  for (P bs : cands) {
    const int code = signalCode(signature(bs));
    if (code < 0) continue;
    P dsx = nullptr, lep = nullptr, nu = nullptr;
    for (P d : directDaughters(bs)) {
      const int a = std::abs(d->pdgId());
      if (!dsx && (a == DS || a == DSST)) dsx = d;
      if (!lep && (a == MU || a == TAU)) lep = d;
      if (!nu && (a == 12 || a == 14 || a == 16)) nu = d;
    }
    P ds = dsx, dsg = nullptr;
    if (std::abs(dsx->pdgId()) == DSST) {
      ds = nullptr;
      for (P d : directDaughters(dsx, true)) {
        const int a = std::abs(d->pdgId());
        if (!ds && a == DS) ds = d;
        if (!dsg && (a == GAMMA || a == PI0)) dsg = d;
      }
      if (!ds) continue;
    }
    // Ds -> phi(K K) pi
    P phi = nullptr;
    {
      auto dd = directDaughters(ds);
      P ph = nullptr, pi = nullptr; int nph = 0, npi = 0;
      for (P d : dd) { if (std::abs(d->pdgId()) == PHI) { ph = d; ++nph; } if (std::abs(d->pdgId()) == PION) { pi = d; ++npi; } }
      if (dd.size() == 2 && nph == 1 && npi == 1 && pi) {
        auto kd = directDaughters(ph);
        int nk = 0; for (P k : kd) nk += std::abs(k->pdgId()) == KAON;
        if (nk == 2 && kd.size() == 2) phi = ph;
      }
    }
    P mu = nullptr, nut = nullptr, num = nullptr;
    if (std::abs(lep->pdgId()) == TAU) {
      auto td = directDaughters(lep);
      std::vector<int> ids; for (P d : td) ids.push_back(std::abs(d->pdgId()));
      std::sort(ids.begin(), ids.end());
      if (ids == std::vector<int>{MU, NU_MU, NU_TAU})
        for (P d : td) { const int a = std::abs(d->pdgId()); if (a == MU) mu = d; if (a == NU_TAU) nut = d; if (a == NU_MU) num = d; }
    } else mu = lep;
    const int rank = 2 * (phi != nullptr) + (mu != nullptr);
    if (rank > bestRank) {
      bestRank = rank;
      out.bs = bs; out.code = code; out.dsx = dsx; out.ds = ds; out.dsg = dsg; out.lep = lep; out.nu = nu;
      out.mu = mu; out.tauNuT = nut; out.tauNuM = num; out.phi = phi;
    }
  }
  if (out.bs) out.bsFirst = firstCopy(out.bs);
  return out;
}

}  // namespace bmmmgen
#endif
