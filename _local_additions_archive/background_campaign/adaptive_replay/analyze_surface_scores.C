#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {
constexpr int kMaxHit = 300000;
constexpr std::size_t kTail = 201;

struct Moment {
  long long signal = 0;
  long double sum = 0.0L, sum2 = 0.0L, sum3 = 0.0L, sum4 = 0.0L;
  long double maximum = 0.0L;
  std::vector<long double> largest;
  void Add(long double value) {
    if (value != 0.0L) ++signal;
    sum += value; sum2 += value*value; sum3 += value*value*value;
    sum4 += value*value*value*value;
    const long double magnitude = std::abs(value);
    maximum = std::max(maximum, magnitude);
    if (magnitude == 0.0L) return;
    auto position = std::lower_bound(largest.begin(), largest.end(), magnitude,
                                     std::greater<long double>());
    largest.insert(position, magnitude);
    if (largest.size() > kTail) largest.pop_back();
  }
};

std::string Group(int detector) {
  if (detector < 70030 || detector > 72730 || (detector - 70030) % 100 != 0)
    return "";
  const int stack = (detector - 70030) / 100;
  if (stack % 4 == 0) return "showermax_open";
  if (stack % 4 == 2) return "showermax_closed";
  return "showermax_transition";
}
}

// mode=original: rate is already divided by the original N, so multiply by
// tree entries to recover one independent primary-history score.
// mode=replay: rate contains 1/replay_events from source generation, so use
// score_scale=replay_events to recover one independent importance score.
void analyze_surface_scores(const char* input, const char* mode,
                            int score_scale, const char* output) {
  TFile file(input, "READ");
  auto* tree = static_cast<TTree*>(file.Get("T"));
  if (file.IsZombie() || !tree) {
    std::cerr << "SCORE_ERROR cannot open input tree\n";
    return;
  }
  tree->SetMakeClass(1);
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("rate", 1);
  tree->SetBranchStatus("hit", 1);
  tree->SetBranchStatus("hit.det", 1);
  tree->SetBranchStatus("hit.k", 1);
  double rate = 0.0;
  int nhit = 0;
  auto det = std::make_unique<int[]>(kMaxHit);
  auto k = std::make_unique<double[]>(kMaxHit);
  tree->SetBranchAddress("rate", &rate);
  tree->SetBranchAddress("hit", &nhit);
  tree->SetBranchAddress("hit.det", det.get());
  tree->SetBranchAddress("hit.k", k.get());

  const Long64_t entries = tree->GetEntries();
  const long double scale = std::string(mode) == "original"
      ? static_cast<long double>(entries) : static_cast<long double>(score_scale);
  std::map<std::string, Moment> moments;
  for (int detector = 70030; detector <= 72730; detector += 100) {
    moments["tile_" + std::to_string(detector) + "/rate"];
    moments["tile_" + std::to_string(detector) + "/rate_energy"];
  }
  for (const char* group : {"showermax_open", "showermax_closed",
                            "showermax_transition"}) {
    moments[std::string(group) + "/rate"];
    moments[std::string(group) + "/rate_energy"];
  }

  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    if (nhit > kMaxHit) nhit = kMaxHit;
    std::map<std::string, long double> event_scores;
    for (int j = 0; j < nhit; ++j) {
      const std::string group = Group(det[j]);
      if (group.empty()) continue;
      const long double weighted_rate = rate * scale;
      const std::string tile = "tile_" + std::to_string(det[j]);
      event_scores[tile + "/rate"] += weighted_rate;
      event_scores[tile + "/rate_energy"] += weighted_rate * std::max(k[j], 0.0);
      event_scores[group + "/rate"] += weighted_rate;
      event_scores[group + "/rate_energy"] += weighted_rate * std::max(k[j], 0.0);
    }
    for (auto& item : moments) item.second.Add(event_scores[item.first]);
  }

  std::ofstream out(output);
  out << "key\tevents\tsignal\tsum\tsum2\tsum3\tsum4\tmaximum\tlargest\n";
  out << std::setprecision(20);
  for (const auto& item : moments) {
    const Moment& value = item.second;
    out << item.first << '\t' << entries << '\t' << value.signal << '\t'
        << value.sum << '\t' << value.sum2 << '\t' << value.sum3 << '\t'
        << value.sum4 << '\t' << value.maximum << '\t';
    if (value.largest.empty()) out << '-';
    else {
      for (std::size_t index = 0; index < value.largest.size(); ++index) {
        if (index) out << ',';
        out << value.largest[index];
      }
    }
    out << '\n';
  }
  out.close();
  std::cout << "SCORE_ANALYSIS mode=" << mode << " entries=" << entries
            << " statistics=" << moments.size() << " output=" << output << "\n";
}
