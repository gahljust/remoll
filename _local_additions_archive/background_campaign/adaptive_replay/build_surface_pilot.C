#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace {
constexpr int kMaxHit = 300000;
struct PilotState {
  Long64_t entry = -1;
  int pid = 0, trid = 0;
  double physical_rate = 0.0, importance = 0.0;
  double x = 0.0, y = 0.0, z = 0.0;
  double px = 0.0, py = 0.0, pz = 0.0, p = 0.0, k = 0.0;
};
}

// Produce a training source made only from exact recorded parents/ancestors of
// the requested detector response. A secondary's own surface crossing is not
// eligible. The pilot is for learning and carries rate=1; it is never included
// in the physical estimator.
void build_surface_pilot(const char* input, int target_detector,
                         const char* observable, int surface_detector,
                         const char* output, int pilot_events = 1000,
                         int unique_states = 200, unsigned seed = 2202701,
                         int top_histories = 400) {
  const bool energy = std::string(observable) == "rate_energy";
  TFile source(input, "READ");
  auto* tree = static_cast<TTree*>(source.Get("T"));
  if (source.IsZombie() || !tree || pilot_events <= 0 || unique_states <= 0) {
    std::cerr << "PILOT_ERROR invalid input or controls\n";
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
                             "hit.py", "hit.pz", "hit.p", "hit.k"})
    tree->SetBranchStatus(branch, 1);
  tree->SetBranchAddress("rate", &rate); tree->SetBranchAddress("hit", &nhit);
  tree->SetBranchAddress("hit.det", det.get()); tree->SetBranchAddress("hit.pid", pid.get());
  tree->SetBranchAddress("hit.trid", trid.get()); tree->SetBranchAddress("hit.mtrid", mtrid.get());
  tree->SetBranchAddress("hit.x", x.get()); tree->SetBranchAddress("hit.y", y.get());
  tree->SetBranchAddress("hit.z", z.get()); tree->SetBranchAddress("hit.px", px.get());
  tree->SetBranchAddress("hit.py", py.get()); tree->SetBranchAddress("hit.pz", pz.get());
  tree->SetBranchAddress("hit.p", p.get()); tree->SetBranchAddress("hit.k", k.get());

  const Long64_t entries = tree->GetEntries();
  std::vector<std::pair<double, Long64_t>> ranked;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    double score = 0.0;
    for (int j = 0; j < std::min(nhit, kMaxHit); ++j)
      if (det[j] == target_detector)
        score += std::abs(rate) * (energy ? std::max(k[j], 0.0) : 1.0);
    if (score > 0.0) ranked.emplace_back(score, entry);
  }
  std::sort(ranked.begin(), ranked.end(), std::greater<>());
  if (static_cast<int>(ranked.size()) > top_histories) ranked.resize(top_histories);

  std::map<std::pair<Long64_t, int>, PilotState> candidates;
  long long target_hits = 0, covered_hits = 0;
  long double total_importance = 0.0L, covered_importance = 0.0L;
  for (const auto& ranked_entry : ranked) {
    const Long64_t entry = ranked_entry.second;
    tree->GetEntry(entry);
    const int count = std::min(nhit, kMaxHit);
    std::map<int, int> mothers, surface_hits;
    for (int j = 0; j < count; ++j) {
      if (trid[j] > 0 && !mothers.count(trid[j])) mothers[trid[j]] = mtrid[j];
      if (det[j] == surface_detector && trid[j] > 0 && pz[j] > 0.0)
        surface_hits[trid[j]] = j;
    }
    for (int j = 0; j < count; ++j) {
      if (det[j] != target_detector) continue;
      const double importance = std::abs(rate) * (energy ? std::max(k[j], 0.0) : 1.0);
      ++target_hits; total_importance += importance;
      int track = mtrid[j] > 0 ? mtrid[j] : trid[j], found = -1;
      std::set<int> seen;
      for (int depth = 0; track > 0 && depth < 1000; ++depth) {
        if (!seen.insert(track).second) break;
        const auto at_surface = surface_hits.find(track);
        if (at_surface != surface_hits.end() && z[at_surface->second] < z[j]) {
          found = at_surface->second; break;
        }
        const auto mother = mothers.find(track);
        track = mother == mothers.end() ? 0 : mother->second;
      }
      if (found < 0) continue;
      ++covered_hits; covered_importance += importance;
      const auto key = std::make_pair(entry, trid[found]);
      auto& state = candidates[key];
      state.entry = entry; state.pid = pid[found]; state.trid = trid[found];
      state.physical_rate = rate; state.importance += importance;
      state.x = x[found]; state.y = y[found]; state.z = z[found];
      state.px = px[found]; state.py = py[found]; state.pz = pz[found];
      state.p = p[found]; state.k = k[found];
    }
  }
  std::vector<PilotState> selected;
  for (const auto& item : candidates) selected.push_back(item.second);
  std::sort(selected.begin(), selected.end(), [](const auto& a, const auto& b) {
    return a.importance > b.importance;
  });
  if (static_cast<int>(selected.size()) > unique_states) selected.resize(unique_states);
  if (selected.empty()) {
    std::cerr << "PILOT_ERROR no exact parent/ancestor surface states "
                 "recovered; child-only crossings are ineligible\n";
    return;
  }

  TFile output_file(output, "RECREATE");
  TTree replay("T", "Exact surface-state response training source");
  remollEvent_t event;
  std::vector<remollGenericDetectorHit_t> hits;
  double replay_rate = 1.0, physical_rate = 0.0, center_importance = 0.0;
  Long64_t source_entry = -1;
  int source_track = 0, source_pid = 0, pilot_state = 0;
  replay.Branch("rate", &replay_rate, "rate/D"); replay.Branch("ev", &event);
  replay.Branch("hit", &hits); replay.Branch("source_entry", &source_entry, "source_entry/L");
  replay.Branch("source_track", &source_track, "source_track/I");
  replay.Branch("source_pid", &source_pid, "source_pid/I");
  replay.Branch("pilot_state", &pilot_state, "pilot_state/I");
  replay.Branch("physical_rate", &physical_rate, "physical_rate/D");
  replay.Branch("center_importance", &center_importance, "center_importance/D");
  std::mt19937 rng(seed);
  std::vector<int> order(selected.size());
  for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
  constexpr int replay_tag = 990035;
  for (int draw = 0; draw < pilot_events; ++draw) {
    if (draw % static_cast<int>(selected.size()) == 0)
      std::shuffle(order.begin(), order.end(), rng);
    const int selected_index = order[draw % selected.size()];
    const auto& state = selected[selected_index];
    std::memset(&event, 0, sizeof(event));
    remollGenericDetectorHit_t hit; std::memset(&hit, 0, sizeof(hit));
    hit.det = replay_tag; hit.pid = state.pid; hit.trid = 1;
    hit.x = state.x; hit.y = state.y; hit.z = state.z;
    hit.px = state.px; hit.py = state.py; hit.pz = state.pz;
    hit.p = state.p; hit.k = state.k;
    const double pmag = std::sqrt(hit.px*hit.px + hit.py*hit.py + hit.pz*hit.pz);
    if (pmag > 0.0) {
      hit.x += 0.01*hit.px/pmag; hit.y += 0.01*hit.py/pmag; hit.z += 0.01*hit.pz/pmag;
    }
    hits.assign(1, hit);
    source_entry = state.entry; source_track = state.trid; source_pid = state.pid;
    pilot_state = selected_index; physical_rate = state.physical_rate;
    center_importance = state.importance;
    replay.Fill();
  }
  replay.Write(); output_file.Close();
  std::cout << "PILOT_BUILD target=" << target_detector << " observable=" << observable
            << " surface=" << surface_detector << " exact_states=" << selected.size()
            << " target_hits=" << target_hits << " covered_hits=" << covered_hits
            << " ranked_coverage=" << (total_importance > 0.0L
                 ? static_cast<double>(covered_importance/total_importance) : 0.0)
            << " pilot_events=" << pilot_events << "\n";
}
