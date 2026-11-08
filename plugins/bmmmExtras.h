#ifndef PhysicsTools_RDsNano_bmmmExtras_h
#define PhysicsTools_RDsNano_bmmmExtras_h
// Quantities ported from the Bmmm R(Ds) ntuplizer (Bmmm/Analysis/python), with the same
// definitions, so that the two ntuples can be compared branch by branch:
//   HitContent      <- TrackHitContent.compute_hit_content
//   massError       <- RDsCandidate.mass_error
//   IsoInputs/iso   <- RDsCandidate._isolation (+ JpsiChargedCandidate constants)

#include <array>
#include <cmath>
#include <vector>
#include <algorithm>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace bmmm {

// ---------------------------------------------------------------- hit content
struct HitContent {
  int n_pix_hit = -1, n_pix_b_hit = -1, n_pix_e_hit = -1, n_pix_layer = -1, n_trk_layer = -1;
  int pix_first_b_layer = 0, pix_first_e_disk = 0, pix_first_layer = 0;
  int n_pix_miss_inner = -1, n_pix_inact_inner = 0;
  int pix_valid_mask = 0, pix_miss_mask = 0, pix_inact_mask = 0, pix_hit_count = 0;
};

// Phase-1 pixel: 4 BPix layers (bits 0-3), 3 FPix disks (bits 4-6); 2 bits per layer in pix_hit_count
inline HitContent hitContent(const reco::Track& trk) {
  HitContent h;
  const reco::HitPattern& hp = trk.hitPattern();
  int counts[7] = {0, 0, 0, 0, 0, 0, 0};
  const reco::HitPattern::HitCategory cats[3] = {reco::HitPattern::TRACK_HITS, reco::HitPattern::MISSING_INNER_HITS,
                                                 reco::HitPattern::MISSING_OUTER_HITS};
  for (auto cat : cats) {
    for (int i = 0; i < hp.numberOfAllHits(cat); ++i) {
      const uint16_t pattern = hp.getHitPattern(cat, i);
      if (!reco::HitPattern::pixelHitFilter(pattern)) continue;
      const bool barrel = reco::HitPattern::pixelBarrelHitFilter(pattern);
      const int layer = static_cast<int>(reco::HitPattern::getLayer(pattern));
      if (barrel ? (layer < 1 || layer > 4) : (layer < 1 || layer > 3)) continue;   // outside Phase-1: skip
      const int bit = barrel ? layer - 1 : 4 + layer - 1;
      if (reco::HitPattern::validHitFilter(pattern)) {
        if (cat != reco::HitPattern::TRACK_HITS) continue;
        if (barrel) h.pix_first_b_layer = h.pix_first_b_layer == 0 ? layer : std::min(h.pix_first_b_layer, layer);
        else        h.pix_first_e_disk  = h.pix_first_e_disk  == 0 ? layer : std::min(h.pix_first_e_disk,  layer);
        h.pix_valid_mask |= 1 << bit;
        counts[bit] = std::min(counts[bit] + 1, 3);
      } else if (reco::HitPattern::inactiveHitFilter(pattern) || reco::HitPattern::badHitFilter(pattern)) {
        h.pix_inact_mask |= 1 << bit;
      } else if (reco::HitPattern::missingHitFilter(pattern)) {
        h.pix_miss_mask |= 1 << bit;
      }
    }
  }
  h.pix_miss_mask  &= ~h.pix_valid_mask;
  h.pix_inact_mask &= ~h.pix_valid_mask;
  for (int b = 0; b < 7; ++b) h.pix_hit_count |= counts[b] << (2 * b);
  for (int i = 0; i < hp.numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS); ++i) {
    const uint16_t pattern = hp.getHitPattern(reco::HitPattern::MISSING_INNER_HITS, i);
    if (reco::HitPattern::pixelHitFilter(pattern) && reco::HitPattern::inactiveHitFilter(pattern)) h.n_pix_inact_inner++;
  }
  h.n_pix_hit   = hp.numberOfValidPixelHits();
  h.n_pix_b_hit = hp.numberOfValidPixelBarrelHits();
  h.n_pix_e_hit = hp.numberOfValidPixelEndcapHits();
  h.n_pix_layer = hp.pixelLayersWithMeasurement();
  h.n_trk_layer = hp.trackerLayersWithMeasurement();
  h.pix_first_layer = h.pix_first_b_layer ? h.pix_first_b_layer : (h.pix_first_e_disk ? 4 + h.pix_first_e_disk : 0);
  h.n_pix_miss_inner = hp.numberOfLostPixelHits(reco::HitPattern::MISSING_INNER_HITS);
  return h;
}

// ---------------------------------------------------------------- mass error
// sigma(m) of a system of tracks by linear propagation of the (q/p, lambda, phi) block of each track's
// covariance (tracks independent). p: calibrated 3-momentum and mass of each track; cov: the covariance of
// the track as in the fits; f: calibration factor (q/p row/column scaled by 1/f as Bmmm's fit_track; 1 in
// data); scal: relative q/p resolution of the calibration without the random number (NaN if none).
// Returns {sigma, sigma_cal}: sigma_cal adds (J_qp * scal * q/p)^2 per track (NaN if any scal is NaN).
struct MassErrTrack {
  double px, py, pz, m;
  int q;
  reco::TrackBase::CovarianceMatrix cov;
  double f;
  double scal;
};
inline std::pair<float, float> massError(const std::vector<MassErrTrack>& trks) {
  double E = 0, Px = 0, Py = 0, Pz = 0;
  for (const auto& t : trks) {
    E += std::sqrt(t.px * t.px + t.py * t.py + t.pz * t.pz + t.m * t.m);
    Px += t.px; Py += t.py; Pz += t.pz;
  }
  const double m2 = E * E - Px * Px - Py * Py - Pz * Pz;
  if (!(m2 > 0)) return {NAN, NAN};
  const double M = std::sqrt(m2);
  double var = 0, varCal = 0;
  bool calOk = true;
  for (const auto& t : trks) {
    const double e = std::sqrt(t.px * t.px + t.py * t.py + t.pz * t.pz + t.m * t.m);
    const double g[3] = {(E * t.px / e - Px) / M, (E * t.py / e - Py) / M, (E * t.pz / e - Pz) / M};
    const double pmag = std::sqrt(t.px * t.px + t.py * t.py + t.pz * t.pz);
    const double ptr = std::hypot(t.px, t.py);
    const double lam = std::atan2(t.pz, ptr), ph = std::atan2(t.py, t.px);
    const double eLam[3] = {-std::sin(lam) * std::cos(ph), -std::sin(lam) * std::sin(ph), std::cos(lam)};
    const double ePhi[3] = {-std::sin(ph), std::cos(ph), 0.};
    const double gp = g[0] * t.px + g[1] * t.py + g[2] * t.pz;
    const double J[3] = {-t.q * pmag * gp, pmag * (g[0] * eLam[0] + g[1] * eLam[1] + g[2] * eLam[2]),
                         ptr * (g[0] * ePhi[0] + g[1] * ePhi[1] + g[2] * ePhi[2])};
    const double s[3] = {1. / t.f, 1., 1.};   // D = diag(1/f, 1, 1, ...) on the q/p row and column
    double v = 0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) v += J[i] * s[i] * t.cov(i, j) * s[j] * J[j];
    var += v;
    if (std::isfinite(t.scal)) {
      const double a = J[0] * t.scal * t.q / pmag;
      varCal += v + a * a;
    } else calOk = false;
  }
  return {var >= 0 ? static_cast<float>(std::sqrt(var)) : NAN,
          (calOk && varCal >= 0) ? static_cast<float>(std::sqrt(varCal)) : NAN};
}

// ---------------------------------------------------------------- isolation
// Muon-POG style PF isolation of an arbitrary direction, against the candidate PV, as Bmmm:
// charged PF candidates associated to the vertex closest in z (the candidate PV z replacing the offline
// vertex it corresponds to); charged hadrons from the PV in the sum, charged candidates from other vertices
// (pt > 0.5) as the PU term; neutrals with pt > 0.5; "clean" removes the four signal tracks.
struct IsoInputs {
  std::vector<float> chEta, chPhi, chPt, chVz; std::vector<int> chQ; std::vector<char> chHad;
  std::vector<float> ntEta, ntPhi, ntPt; std::vector<char> ntNh;
  void fill(const std::vector<pat::PackedCandidate>& pf) {
    for (const auto& c : pf) {
      const int pdg = std::abs(c.pdgId());
      if (pdg == 211 || pdg == 321 || pdg == 11 || pdg == 13 || pdg == 2212 || pdg == 999211) {
        chEta.push_back(c.eta()); chPhi.push_back(c.phi()); chPt.push_back(c.pt()); chVz.push_back(c.vz());
        chQ.push_back(c.charge()); chHad.push_back(pdg == 211 || pdg == 321 || pdg == 2212 || pdg == 999211);
      } else if (pdg == 111 || pdg == 130 || pdg == 310 || pdg == 2112 || pdg == 22) {
        if (c.pt() > 0.5) {
          ntEta.push_back(c.eta()); ntPhi.push_back(c.phi()); ntPt.push_back(c.pt());
          ntNh.push_back(pdg != 22);
        }
      }
    }
  }
};
struct IsoResult { float ch, ch_clean, pu, nh, ph, iso, iso_clean, reliso, reliso_clean; };
struct SignalTrk { float eta, phi, pt; int q; };

// fromPV[i]: charged candidate i associated to the candidate PV; isSig[i]: it is one of the signal tracks
inline void isoAssociation(const IsoInputs& in, const std::vector<float>& vtxZ, int pvIdx, float pvZ,
                           const std::vector<SignalTrk>& sig, std::vector<char>& fromPV, std::vector<char>& isSig) {
  const size_t n = in.chPt.size();
  fromPV.assign(n, 0); isSig.assign(n, 0);
  const bool hasPV = pvIdx >= 0 && pvIdx < static_cast<int>(vtxZ.size());
  for (size_t i = 0; i < n; ++i) {
    if (hasPV) {
      int best = -1; float bestDz = 0;
      for (size_t v = 0; v < vtxZ.size(); ++v) {
        const float z = (static_cast<int>(v) == pvIdx) ? pvZ : vtxZ[v];
        const float dz = std::abs(z - in.chVz[i]);
        if (best < 0 || dz < bestDz) { best = v; bestDz = dz; }
      }
      fromPV[i] = best == pvIdx;
    }
    for (const auto& s : sig)
      if (in.chQ[i] == s.q && reco::deltaR(s.eta, s.phi, in.chEta[i], in.chPhi[i]) < 0.01 &&
          std::abs(in.chPt[i] - s.pt) < 0.05 * s.pt) { isSig[i] = 1; break; }
  }
}

inline IsoResult isolation(const IsoInputs& in, const std::vector<char>& fromPV, const std::vector<char>& isSig,
                           float eta, float phi, float pt, float cone) {
  double ch = 0, chClean = 0, pu = 0, puClean = 0, nh = 0, ph = 0;
  for (size_t i = 0; i < in.chPt.size(); ++i) {
    const float dr = reco::deltaR(eta, phi, in.chEta[i], in.chPhi[i]);
    if (!(dr < cone && dr > 1e-4)) continue;
    if (fromPV[i] && in.chHad[i]) { ch += in.chPt[i]; if (!isSig[i]) chClean += in.chPt[i]; }
    if (!fromPV[i] && in.chPt[i] > 0.5 && dr > 0.01) { pu += in.chPt[i]; if (!isSig[i]) puClean += in.chPt[i]; }
  }
  for (size_t i = 0; i < in.ntPt.size(); ++i) {
    const float dr = reco::deltaR(eta, phi, in.ntEta[i], in.ntPhi[i]);
    if (!(dr < cone && dr > 1e-2)) continue;
    if (in.ntNh[i]) nh += in.ntPt[i]; else ph += in.ntPt[i];
  }
  IsoResult r;
  r.ch = ch; r.ch_clean = chClean; r.pu = pu; r.nh = nh; r.ph = ph;
  r.iso       = ch      + std::max(0., nh + ph - 0.5 * pu);
  r.iso_clean = chClean + std::max(0., nh + ph - 0.5 * puClean);
  r.reliso       = pt > 0 ? r.iso / pt : NAN;
  r.reliso_clean = pt > 0 ? r.iso_clean / pt : NAN;
  return r;
}

}  // namespace bmmm
#endif
