#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace {
constexpr int kMaxHit = 300000;
struct Learned {
  Long64_t entry = -1;
  int pid = 0, trid = 0;
  double physical_rate = 0.0, importance = 0.0;
  double x = 0.0, y = 0.0, z = 0.0;
  double px = 0.0, py = 0.0, pz = 0.0, p = 0.0, k = 0.0;
  long long trials = 0, successes = 0;
  double response_sum = 0.0, response_sum2 = 0.0;
};
}

// Match a single-thread exact-state pilot result entry-by-entry to its source,
// then update one learned-response tree inside a shared learning ROOT file.
void analyze_surface_pilot(const char* pilot_source, const char* pilot_result,
                           int target_detector, const char* observable,
                           const char* learning_file, const char* tree_name) {
  TFile source_file(pilot_source, "READ"), result_file(pilot_result, "READ");
  auto* source = static_cast<TTree*>(source_file.Get("T"));
  auto* result = static_cast<TTree*>(result_file.Get("T"));
  if (!source || !result || source->GetEntries() != result->GetEntries()) {
    std::cerr << "PILOT_ERROR source/result trees are missing or misaligned\n";
    return;
  }
  const bool energy = std::string(observable) == "rate_energy";
  source->SetMakeClass(1); result->SetMakeClass(1);
  source->SetBranchStatus("*", 0); result->SetBranchStatus("*", 0);
  int source_nhit = 0, result_nhit = 0;
  auto source_pid_hit = std::make_unique<int[]>(kMaxHit);
  auto source_x = std::make_unique<double[]>(kMaxHit);
  auto source_y = std::make_unique<double[]>(kMaxHit);
  auto source_z = std::make_unique<double[]>(kMaxHit);
  auto source_px = std::make_unique<double[]>(kMaxHit);
  auto source_py = std::make_unique<double[]>(kMaxHit);
  auto source_pz = std::make_unique<double[]>(kMaxHit);
  auto source_p = std::make_unique<double[]>(kMaxHit);
  auto source_k = std::make_unique<double[]>(kMaxHit);
  auto result_det = std::make_unique<int[]>(kMaxHit);
  auto result_k = std::make_unique<double[]>(kMaxHit);
  Long64_t source_entry = -1;
  int source_track = 0, source_pid = 0;
  double physical_rate = 0.0, importance = 0.0;
  for (const char* branch : {"hit", "hit.pid", "hit.x", "hit.y", "hit.z",
                             "hit.px", "hit.py", "hit.pz", "hit.p", "hit.k",
                             "source_entry", "source_track", "source_pid",
                             "physical_rate", "center_importance"})
    source->SetBranchStatus(branch, 1);
  source->SetBranchAddress("hit", &source_nhit);
  source->SetBranchAddress("hit.pid", source_pid_hit.get());
  source->SetBranchAddress("hit.x", source_x.get()); source->SetBranchAddress("hit.y", source_y.get());
  source->SetBranchAddress("hit.z", source_z.get()); source->SetBranchAddress("hit.px", source_px.get());
  source->SetBranchAddress("hit.py", source_py.get()); source->SetBranchAddress("hit.pz", source_pz.get());
  source->SetBranchAddress("hit.p", source_p.get()); source->SetBranchAddress("hit.k", source_k.get());
  source->SetBranchAddress("source_entry", &source_entry);
  source->SetBranchAddress("source_track", &source_track);
  source->SetBranchAddress("source_pid", &source_pid);
  source->SetBranchAddress("physical_rate", &physical_rate);
  source->SetBranchAddress("center_importance", &importance);
  result->SetBranchStatus("hit", 1); result->SetBranchStatus("hit.det", 1);
  result->SetBranchStatus("hit.k", 1);
  result->SetBranchAddress("hit", &result_nhit);
  result->SetBranchAddress("hit.det", result_det.get());
  result->SetBranchAddress("hit.k", result_k.get());

  std::map<std::pair<Long64_t, int>, Learned> learned;
  // Merge any earlier pilots for this target.
  {
    TFile previous(learning_file, "READ");
    auto* old = static_cast<TTree*>(previous.Get(tree_name));
    if (old) {
      Learned row;
      old->SetBranchAddress("entry", &row.entry); old->SetBranchAddress("pid", &row.pid);
      old->SetBranchAddress("trid", &row.trid); old->SetBranchAddress("physical_rate", &row.physical_rate);
      old->SetBranchAddress("importance", &row.importance); old->SetBranchAddress("x", &row.x);
      old->SetBranchAddress("y", &row.y); old->SetBranchAddress("z", &row.z);
      old->SetBranchAddress("px", &row.px); old->SetBranchAddress("py", &row.py);
      old->SetBranchAddress("pz", &row.pz); old->SetBranchAddress("p", &row.p);
      old->SetBranchAddress("k", &row.k); old->SetBranchAddress("trials", &row.trials);
      old->SetBranchAddress("successes", &row.successes);
      old->SetBranchAddress("response_sum", &row.response_sum);
      old->SetBranchAddress("response_sum2", &row.response_sum2);
      for (Long64_t i = 0; i < old->GetEntries(); ++i) {
        old->GetEntry(i); learned[{row.entry, row.trid}] = row;
      }
    }
  }

  long long pilot_successes = 0;
  for (Long64_t entry = 0; entry < source->GetEntries(); ++entry) {
    source->GetEntry(entry); result->GetEntry(entry);
    if (source_nhit <= 0) continue;
    const auto key = std::make_pair(source_entry, source_track);
    auto& row = learned[key];
    row.entry = source_entry; row.pid = source_pid; row.trid = source_track;
    row.physical_rate = physical_rate; row.importance = std::max(row.importance, importance);
    row.x = source_x[0]; row.y = source_y[0]; row.z = source_z[0];
    row.px = source_px[0]; row.py = source_py[0]; row.pz = source_pz[0];
    row.p = source_p[0]; row.k = source_k[0];
    double response = 0.0;
    for (int j = 0; j < std::min(result_nhit, kMaxHit); ++j)
      if (result_det[j] == target_detector)
        response += energy ? std::max(result_k[j], 0.0) : 1.0;
    ++row.trials; row.response_sum += response; row.response_sum2 += response*response;
    if (response > 0.0) { ++row.successes; ++pilot_successes; }
  }

  TFile output(learning_file, "UPDATE");
  output.Delete((std::string(tree_name) + ";*").c_str());
  TTree tree(tree_name, "Measured conditional response of exact surface states");
  Learned row;
  tree.Branch("entry", &row.entry); tree.Branch("pid", &row.pid); tree.Branch("trid", &row.trid);
  tree.Branch("physical_rate", &row.physical_rate); tree.Branch("importance", &row.importance);
  tree.Branch("x", &row.x); tree.Branch("y", &row.y); tree.Branch("z", &row.z);
  tree.Branch("px", &row.px); tree.Branch("py", &row.py); tree.Branch("pz", &row.pz);
  tree.Branch("p", &row.p); tree.Branch("k", &row.k); tree.Branch("trials", &row.trials);
  tree.Branch("successes", &row.successes); tree.Branch("response_sum", &row.response_sum);
  tree.Branch("response_sum2", &row.response_sum2);
  long long total_trials = 0, total_successes = 0, successful_states = 0;
  for (const auto& item : learned) {
    row = item.second; tree.Fill(); total_trials += row.trials; total_successes += row.successes;
    if (row.successes > 0) ++successful_states;
  }
  tree.Write(); output.Close();
  std::cout << "PILOT_ANALYSIS target=" << target_detector << " observable=" << observable
            << " pilot_events=" << source->GetEntries() << " pilot_successes=" << pilot_successes
            << " total_trials=" << total_trials << " total_successes=" << total_successes
            << " successful_states=" << successful_states << " tree=" << tree_name << "\n";
}
