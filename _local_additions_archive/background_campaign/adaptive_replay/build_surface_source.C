#include "TFile.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace {
constexpr int kMaxHit = 300000;

struct State {
  Long64_t entry = -1;
  int pid = 0;
  int trid = 0;
  double rate = 0.0;
  double x = 0.0, y = 0.0, z = 0.0;
  double px = 0.0, py = 0.0, pz = 0.0;
  double p = 0.0, k = 0.0;
  double strength = 1.0;
};

struct BankState : State {
  double mass = 0.0;
  double affinity = 0.0;
  double physical_pdf = 0.0;
  double proposal_pdf = 0.0;
};

double LogK(const State& h) { return std::log(std::max(h.k, 1.0e-9)); }
double Tx(const State& h) { return h.pz != 0.0 ? h.px / h.pz : 0.0; }
double Ty(const State& h) { return h.pz != 0.0 ? h.py / h.pz : 0.0; }

State MakeState(Long64_t entry, int j, double rate, const int* pid,
                const int* trid, const double* x, const double* y,
                const double* z, const double* px, const double* py,
                const double* pz, const double* p, const double* k) {
  State h;
  h.entry = entry; h.pid = pid[j]; h.trid = trid[j]; h.rate = rate;
  h.x = x[j]; h.y = y[j]; h.z = z[j];
  h.px = px[j]; h.py = py[j]; h.pz = pz[j]; h.p = p[j]; h.k = k[j];
  return h;
}

template <class Feature>
double WeightedScale(const std::vector<BankState>& bank, Feature feature,
                     double floor) {
  long double sw = 0.0L, s1 = 0.0L, s2 = 0.0L;
  for (const auto& b : bank) {
    const double value = feature(b);
    sw += b.mass; s1 += b.mass * value; s2 += b.mass * value * value;
  }
  const long double variance = sw > 0.0L
      ? s2 / sw - (s1 / sw) * (s1 / sw) : 0.0L;
  return std::max(floor, static_cast<double>(
      0.35L * std::sqrt(std::max(0.0L, variance))));
}

template <class Feature>
double LearnedScale(const std::vector<State>& centers, Feature feature,
                    double floor) {
  long double sw = 0.0L, s1 = 0.0L, s2 = 0.0L;
  for (const auto& center : centers) {
    const double weight = std::max(center.strength, 1.0e-30);
    const double value = feature(center);
    sw += weight; s1 += weight * value; s2 += weight * value * value;
  }
  const long double variance = sw > 0.0L
      ? s2 / sw - (s1 / sw) * (s1 / sw) : 0.0L;
  return std::max(floor, static_cast<double>(
      0.20L * std::sqrt(std::max(0.0L, variance))));
}
}

// Construct one full-support replay proposal from an untouched remoll file.
// The physical bank contains every unique forward track crossing at one common
// surface, including every particle species. Only the targeted component is
// concentrated near parents/ancestors of the requested ShowerMax response.
// A secondary's own crossing is deliberately ineligible: replaying it cannot
// improve the probability of the interaction that created it.
void build_surface_source(const char* input, int target_detector,
                          const char* observable, int surface_detector,
                          const char* output, int replay_events = 1000,
                          double physical_fraction = 0.20,
                          unsigned seed = 2202701, int top_histories = 400,
                          const char* learned_file = "",
                          const char* learned_tree = "") {
  if (replay_events <= 0 || physical_fraction <= 0.0 || physical_fraction > 1.0) {
    std::cerr << "SURFACE_ERROR invalid replay_events or physical_fraction\n";
    return;
  }
  const bool energy_observable = std::string(observable) == "rate_energy";
  TFile source(input, "READ");
  auto* tree = static_cast<TTree*>(source.Get("T"));
  if (source.IsZombie() || !tree) {
    std::cerr << "SURFACE_ERROR cannot open input tree\n";
    return;
  }
  tree->SetMakeClass(1);
  tree->SetBranchStatus("*", 0);
  int nhit = 0;
  double rate = 0.0;
  auto det = std::make_unique<int[]>(kMaxHit);
  auto pid = std::make_unique<int[]>(kMaxHit);
  auto trid = std::make_unique<int[]>(kMaxHit);
  auto mtrid = std::make_unique<int[]>(kMaxHit);
  auto x = std::make_unique<double[]>(kMaxHit);
  auto y = std::make_unique<double[]>(kMaxHit);
  auto z = std::make_unique<double[]>(kMaxHit);
  auto px = std::make_unique<double[]>(kMaxHit);
  auto py = std::make_unique<double[]>(kMaxHit);
  auto pz = std::make_unique<double[]>(kMaxHit);
  auto p = std::make_unique<double[]>(kMaxHit);
  auto k = std::make_unique<double[]>(kMaxHit);
  for (const char* branch : {"rate", "hit", "hit.det", "hit.pid", "hit.trid",
                             "hit.mtrid", "hit.x", "hit.y", "hit.z", "hit.px",
                             "hit.py", "hit.pz", "hit.p", "hit.k"}) {
    tree->SetBranchStatus(branch, 1);
  }
  tree->SetBranchAddress("rate", &rate);
  tree->SetBranchAddress("hit", &nhit);
  tree->SetBranchAddress("hit.det", det.get());
  tree->SetBranchAddress("hit.pid", pid.get());
  tree->SetBranchAddress("hit.trid", trid.get());
  tree->SetBranchAddress("hit.mtrid", mtrid.get());
  tree->SetBranchAddress("hit.x", x.get());
  tree->SetBranchAddress("hit.y", y.get());
  tree->SetBranchAddress("hit.z", z.get());
  tree->SetBranchAddress("hit.px", px.get());
  tree->SetBranchAddress("hit.py", py.get());
  tree->SetBranchAddress("hit.pz", pz.get());
  tree->SetBranchAddress("hit.p", p.get());
  tree->SetBranchAddress("hit.k", k.get());

  const Long64_t entries = tree->GetEntries();
  std::vector<std::pair<double, Long64_t>> ranked;
  long double target_importance = 0.0L;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    if (nhit > kMaxHit) nhit = kMaxHit;
    double score = 0.0;
    for (int j = 0; j < nhit; ++j) {
      if (det[j] != target_detector) continue;
      score += std::abs(rate) * (energy_observable ? std::max(k[j], 0.0) : 1.0);
    }
    target_importance += score;
    if (score > 0.0) ranked.emplace_back(score, entry);
  }
  std::sort(ranked.begin(), ranked.end(), std::greater<>());
  if (static_cast<int>(ranked.size()) > top_histories) ranked.resize(top_histories);

  std::vector<State> centers;
  long double ranked_importance = 0.0L, covered_importance = 0.0L;
  long long target_hits = 0, covered_hits = 0;
  for (const auto& ranked_entry : ranked) {
    const Long64_t entry = ranked_entry.second;
    tree->GetEntry(entry);
    if (nhit > kMaxHit) nhit = kMaxHit;
    std::map<int, int> mothers;
    std::map<int, int> surface_hit;
    for (int j = 0; j < nhit; ++j) {
      if (trid[j] > 0 && !mothers.count(trid[j])) mothers[trid[j]] = mtrid[j];
      if (det[j] == surface_detector && pz[j] > 0.0 && trid[j] > 0)
        surface_hit[trid[j]] = j;
    }
    for (int j = 0; j < nhit; ++j) {
      if (det[j] != target_detector) continue;
      const long double importance = std::abs(rate) *
          (energy_observable ? std::max(k[j], 0.0) : 1.0);
      ranked_importance += importance;
      ++target_hits;
      // A primary response may be replayed from its own upstream state.  A
      // secondary response must begin at its mother, never at the already-made
      // secondary.  Continue up the recorded ancestry only when necessary.
      int track = mtrid[j] > 0 ? mtrid[j] : trid[j], found = -1;
      std::set<int> seen;
      for (int depth = 0; track > 0 && depth < 1000; ++depth) {
        if (!seen.insert(track).second) break;
        const auto surface = surface_hit.find(track);
        if (surface != surface_hit.end() && z[surface->second] < z[j]) {
          found = surface->second;
          break;
        }
        const auto mother = mothers.find(track);
        track = mother == mothers.end() ? 0 : mother->second;
      }
      if (found < 0) continue;
      covered_importance += importance;
      ++covered_hits;
      centers.push_back(MakeState(entry, found, rate, pid.get(), trid.get(),
                                  x.get(), y.get(), z.get(), px.get(), py.get(),
                                  pz.get(), p.get(), k.get()));
    }
  }
  if (centers.empty()) {
    std::cerr << "SURFACE_ERROR no target parent/ancestor crosses requested "
                 "surface; child-only crossings are not valid production "
                 "biasing states\n";
    return;
  }

  bool learned_used = false;
  if (learned_file && learned_file[0] && learned_tree && learned_tree[0]) {
    TFile learning(learned_file, "READ");
    auto* learned = static_cast<TTree*>(learning.Get(learned_tree));
    if (!learned) {
      std::cerr << "SURFACE_ERROR requested learned-response tree is missing\n";
      return;
    }
    State center;
    long long trials = 0, successes = 0;
    double response_sum = 0.0;
    learned->SetBranchAddress("entry", &center.entry);
    learned->SetBranchAddress("pid", &center.pid);
    learned->SetBranchAddress("trid", &center.trid);
    learned->SetBranchAddress("physical_rate", &center.rate);
    learned->SetBranchAddress("x", &center.x); learned->SetBranchAddress("y", &center.y);
    learned->SetBranchAddress("z", &center.z); learned->SetBranchAddress("px", &center.px);
    learned->SetBranchAddress("py", &center.py); learned->SetBranchAddress("pz", &center.pz);
    learned->SetBranchAddress("p", &center.p); learned->SetBranchAddress("k", &center.k);
    learned->SetBranchAddress("trials", &trials);
    learned->SetBranchAddress("successes", &successes);
    learned->SetBranchAddress("response_sum", &response_sum);
    std::vector<State> successful;
    for (Long64_t i = 0; i < learned->GetEntries(); ++i) {
      learned->GetEntry(i);
      if (successes <= 0 || trials <= 0 || response_sum <= 0.0) continue;
      center.strength = response_sum / trials;
      successful.push_back(center);
    }
    if (successful.empty()) {
      std::cerr << "SURFACE_ERROR learned-response tree has no successful states\n";
      return;
    }
    centers.swap(successful);
    learned_used = true;
  }

  std::vector<BankState> bank;
  std::map<int, long long> bank_species;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    if (nhit > kMaxHit) nhit = kMaxHit;
    std::set<int> seen;
    for (int j = 0; j < nhit; ++j) {
      if (det[j] != surface_detector || pz[j] <= 0.0 || trid[j] <= 0) continue;
      if (!seen.insert(trid[j]).second) continue;
      BankState b;
      static_cast<State&>(b) = MakeState(entry, j, rate, pid.get(), trid.get(),
                                         x.get(), y.get(), z.get(), px.get(),
                                         py.get(), pz.get(), p.get(), k.get());
      b.mass = std::abs(rate);
      if (b.mass > 0.0 && std::isfinite(b.mass)) {
        bank.push_back(b);
        ++bank_species[b.pid];
      }
    }
  }
  if (bank.empty()) {
    std::cerr << "SURFACE_ERROR physical surface bank is empty\n";
    return;
  }
  const double sx = learned_used ? LearnedScale(centers, [](const auto& h) { return h.x; }, 2.0)
                                 : WeightedScale(bank, [](const auto& h) { return h.x; }, 5.0);
  const double sy = learned_used ? LearnedScale(centers, [](const auto& h) { return h.y; }, 2.0)
                                 : WeightedScale(bank, [](const auto& h) { return h.y; }, 5.0);
  const double sk = learned_used ? LearnedScale(centers, [](const auto& h) { return LogK(h); }, 0.05)
                                 : WeightedScale(bank, [](const auto& h) { return LogK(h); }, 0.15);
  const double stx = learned_used ? LearnedScale(centers, [](const auto& h) { return Tx(h); }, 0.0005)
                                  : WeightedScale(bank, [](const auto& h) { return Tx(h); }, 0.002);
  const double sty = learned_used ? LearnedScale(centers, [](const auto& h) { return Ty(h); }, 0.0005)
                                  : WeightedScale(bank, [](const auto& h) { return Ty(h); }, 0.002);

  std::map<int, std::vector<const State*>> centers_by_pid;
  for (const auto& center : centers) centers_by_pid[center.pid].push_back(&center);
  double maximum_strength = 0.0;
  for (const auto& center : centers) maximum_strength = std::max(maximum_strength, center.strength);
  long double physical_norm = 0.0L, targeted_norm = 0.0L;
  for (auto& b : bank) {
    double affinity = 0.0;
    const auto found = centers_by_pid.find(b.pid);
    if (found != centers_by_pid.end()) {
      for (const State* center : found->second) {
        const double d2 = std::pow((b.x - center->x) / sx, 2) +
                          std::pow((b.y - center->y) / sy, 2) +
                          std::pow((LogK(b) - LogK(*center)) / sk, 2) +
                          std::pow((Tx(b) - Tx(*center)) / stx, 2) +
                          std::pow((Ty(b) - Ty(*center)) / sty, 2);
        const double relative_strength = maximum_strength > 0.0
            ? center->strength / maximum_strength : 1.0;
        affinity = std::max(affinity, relative_strength * std::exp(-0.5 * d2));
      }
    }
    b.affinity = affinity;
    physical_norm += b.mass;
    targeted_norm += b.mass * b.affinity;
  }
  if (targeted_norm <= 0.0L) {
    std::cerr << "SURFACE_ERROR targeted density has no bank overlap\n";
    return;
  }
  std::vector<double> cumulative;
  double cumulative_q = 0.0;
  for (auto& b : bank) {
    b.physical_pdf = static_cast<double>(b.mass / physical_norm);
    const double targeted_pdf = static_cast<double>(b.mass * b.affinity / targeted_norm);
    b.proposal_pdf = physical_fraction * b.physical_pdf +
                     (1.0 - physical_fraction) * targeted_pdf;
    cumulative_q += b.proposal_pdf;
    cumulative.push_back(cumulative_q);
  }
  for (double& value : cumulative) value /= cumulative_q;
  for (auto& b : bank) b.proposal_pdf /= cumulative_q;

  TFile replay_file(output, "RECREATE");
  TTree density("density", "Complete observed surface bank and current proposal");
  BankState density_row;
  density.Branch("entry", &density_row.entry);
  density.Branch("pid", &density_row.pid);
  density.Branch("trid", &density_row.trid);
  density.Branch("rate", &density_row.rate);
  density.Branch("x", &density_row.x); density.Branch("y", &density_row.y);
  density.Branch("z", &density_row.z);
  density.Branch("px", &density_row.px); density.Branch("py", &density_row.py);
  density.Branch("pz", &density_row.pz); density.Branch("k", &density_row.k);
  density.Branch("affinity", &density_row.affinity);
  density.Branch("physical_pdf", &density_row.physical_pdf);
  density.Branch("proposal_pdf", &density_row.proposal_pdf);
  for (const auto& b : bank) { density_row = b; density.Fill(); }
  density.Write();

  TTree replay("T", "Adaptive full-support surface replay source");
  remollEvent_t replay_event;
  std::vector<remollGenericDetectorHit_t> replay_hits;
  double replay_rate = 0.0, proposal_pdf = 0.0, physical_pdf = 0.0, affinity = 0.0;
  Long64_t source_entry = -1;
  int source_track = 0, source_pid = 0, source_det = surface_detector;
  int proposal_target = target_detector;
  replay.Branch("rate", &replay_rate, "rate/D");
  replay.Branch("ev", &replay_event);
  replay.Branch("hit", &replay_hits);
  replay.Branch("source_entry", &source_entry, "source_entry/L");
  replay.Branch("source_track", &source_track, "source_track/I");
  replay.Branch("source_pid", &source_pid, "source_pid/I");
  replay.Branch("source_det", &source_det, "source_det/I");
  replay.Branch("proposal_target", &proposal_target, "proposal_target/I");
  replay.Branch("proposal_pdf", &proposal_pdf, "proposal_pdf/D");
  replay.Branch("physical_pdf", &physical_pdf, "physical_pdf/D");
  replay.Branch("affinity", &affinity, "affinity/D");
  TRandom3 rng(seed);
  constexpr int replay_tag = 990035;
  for (int draw = 0; draw < replay_events; ++draw) {
    const double uniform = rng.Uniform();
    size_t index = std::lower_bound(cumulative.begin(), cumulative.end(), uniform)
        - cumulative.begin();
    if (index >= bank.size()) index = bank.size() - 1;
    const auto& b = bank[index];
    std::memset(&replay_event, 0, sizeof(replay_event));
    remollGenericDetectorHit_t hit;
    std::memset(&hit, 0, sizeof(hit));
    hit.det = replay_tag; hit.pid = b.pid; hit.trid = 1;
    hit.x = b.x; hit.y = b.y; hit.z = b.z;
    hit.px = b.px; hit.py = b.py; hit.pz = b.pz; hit.p = b.p; hit.k = b.k;
    const double pmag = std::sqrt(hit.px*hit.px + hit.py*hit.py + hit.pz*hit.pz);
    if (pmag > 0.0) {
      hit.x += 0.01 * hit.px / pmag;
      hit.y += 0.01 * hit.py / pmag;
      hit.z += 0.01 * hit.pz / pmag;
    }
    replay_hits.assign(1, hit);
    replay_rate = b.rate / (replay_events * b.proposal_pdf);
    source_entry = b.entry; source_track = b.trid; source_pid = b.pid;
    proposal_pdf = b.proposal_pdf; physical_pdf = b.physical_pdf;
    affinity = b.affinity;
    replay.Fill();
  }
  replay.Write();
  replay_file.Close();

  std::cout << "SURFACE_BUILD target=" << target_detector
            << " observable=" << observable
            << " surface=" << surface_detector
            << " input_events=" << entries
            << " bank=" << bank.size()
            << " species=" << bank_species.size()
            << " centers=" << centers.size()
            << " learned=" << learned_used
            << " ranked_hits=" << target_hits
            << " covered_hits=" << covered_hits
            << " ranked_coverage="
            << (ranked_importance > 0.0L ? static_cast<double>(covered_importance / ranked_importance) : 0.0)
            << " physical_fraction=" << physical_fraction
            << " replay_events=" << replay_events << "\n";
}
