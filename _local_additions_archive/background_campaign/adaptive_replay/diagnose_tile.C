#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr int kMaxHit = 300000;
constexpr int kString = 255;

struct History {
  Long64_t entry{-1};
  long double score{0.0L};
  double rate{0.0};
};

struct Breakdown {
  long double score{0.0L};
  long long hits{0};
  long long histories{0};
  Long64_t last_entry{-1};
};

double RelativeError(long long n, long double sum, long double sum2) {
  if (n < 2 || sum == 0.0L) return 0.0;
  const long double variance = std::max(
      0.0L, (sum2 - sum * sum / static_cast<long double>(n)) /
                static_cast<long double>(n - 1));
  return static_cast<double>(std::sqrt(variance / n) /
                             std::abs(sum / n));
}

std::string Clean(const char* text) {
  std::string value = text == nullptr || text[0] == '\0' ? "-" : text;
  std::replace(value.begin(), value.end(), '\t', ' ');
  std::replace(value.begin(), value.end(), '\n', ' ');
  return value;
}

}  // namespace

// Diagnose one selected tile observable without modifying the input file.
// The output prefix receives .evolution.tsv, .histories.tsv, .hits.tsv, and
// .sources.tsv files. Scores use one generated remoll event as the independent
// statistical history.
void diagnose_tile(const char* input, int detector = 71630,
                   const char* observable = "rate",
                   const char* output_prefix = "tile_diagnostic",
                   int evolution_stride = 1000, int top_histories = 100) {
  const bool energy = std::string(observable) == "rate_energy";
  TFile file(input, "READ");
  auto* tree = static_cast<TTree*>(file.Get("T"));
  if (file.IsZombie() || tree == nullptr || evolution_stride <= 0 ||
      top_histories <= 0) {
    std::cerr << "TILE_DIAGNOSTIC_ERROR invalid input or controls\n";
    return;
  }

  tree->SetMakeClass(1);
  tree->SetBranchStatus("*", 0);
  for (const char* branch : {
           "rate", "hit", "hit.det", "hit.pid", "hit.trid", "hit.mtrid",
           "hit.k", "hit.x", "hit.y", "hit.z", "hit.vx", "hit.vy",
           "hit.vz", "hit.px", "hit.py", "hit.pz",
           "hit.creator_physvol_name[255]", "hit.creator_process_name[255]",
           "hit.creator_material_name[255]"}) {
    tree->SetBranchStatus(branch, 1);
  }
  double rate = 0.0;
  int nhit = 0;
  auto det = std::make_unique<int[]>(kMaxHit);
  auto pid = std::make_unique<int[]>(kMaxHit);
  auto trid = std::make_unique<int[]>(kMaxHit);
  auto mtrid = std::make_unique<int[]>(kMaxHit);
  auto k = std::make_unique<double[]>(kMaxHit);
  auto x = std::make_unique<double[]>(kMaxHit);
  auto y = std::make_unique<double[]>(kMaxHit);
  auto z = std::make_unique<double[]>(kMaxHit);
  auto vx = std::make_unique<double[]>(kMaxHit);
  auto vy = std::make_unique<double[]>(kMaxHit);
  auto vz = std::make_unique<double[]>(kMaxHit);
  auto px = std::make_unique<double[]>(kMaxHit);
  auto py = std::make_unique<double[]>(kMaxHit);
  auto pz = std::make_unique<double[]>(kMaxHit);
  auto creator_volume = std::make_unique<char[]>(kMaxHit * kString);
  auto creator_process = std::make_unique<char[]>(kMaxHit * kString);
  auto creator_material = std::make_unique<char[]>(kMaxHit * kString);
  tree->SetBranchAddress("rate", &rate);
  tree->SetBranchAddress("hit", &nhit);
  tree->SetBranchAddress("hit.det", det.get());
  tree->SetBranchAddress("hit.pid", pid.get());
  tree->SetBranchAddress("hit.trid", trid.get());
  tree->SetBranchAddress("hit.mtrid", mtrid.get());
  tree->SetBranchAddress("hit.k", k.get());
  tree->SetBranchAddress("hit.x", x.get());
  tree->SetBranchAddress("hit.y", y.get());
  tree->SetBranchAddress("hit.z", z.get());
  tree->SetBranchAddress("hit.vx", vx.get());
  tree->SetBranchAddress("hit.vy", vy.get());
  tree->SetBranchAddress("hit.vz", vz.get());
  tree->SetBranchAddress("hit.px", px.get());
  tree->SetBranchAddress("hit.py", py.get());
  tree->SetBranchAddress("hit.pz", pz.get());
  tree->SetBranchAddress("hit.creator_physvol_name[255]", creator_volume.get());
  tree->SetBranchAddress("hit.creator_process_name[255]", creator_process.get());
  tree->SetBranchAddress("hit.creator_material_name[255]", creator_material.get());
  const Long64_t entries = tree->GetEntries();
  const long double scale = static_cast<long double>(entries);

  std::ofstream evolution(std::string(output_prefix) + ".evolution.tsv");
  evolution << "events\tsignal_histories\tmean\trse\tneff\tmaximum_fraction\n";
  evolution << std::setprecision(20);
  std::vector<History> histories;
  histories.reserve(entries);
  long double sum = 0.0L, sum2 = 0.0L, maximum = 0.0L;
  long long signal = 0;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    long double score = 0.0L;
    const int count = std::min(nhit, kMaxHit);
    for (int index = 0; index < count; ++index) {
      if (det[index] != detector) continue;
      score += rate * scale * (energy ? std::max(k[index], 0.0) : 1.0);
    }
    if (score != 0.0L) {
      ++signal;
      histories.push_back({entry, score, rate});
    }
    sum += score;
    sum2 += score * score;
    maximum = std::max(maximum, std::abs(score));
    const long long n = entry + 1;
    if (n % evolution_stride == 0 || n == entries) {
      evolution << n << '\t' << signal << '\t' << sum / n << '\t'
                << RelativeError(n, sum, sum2) << '\t'
                << (sum2 > 0.0L ? sum * sum / sum2 : 0.0L) << '\t'
                << (sum != 0.0L ? maximum / std::abs(sum) : 0.0L) << '\n';
    }
  }
  evolution.close();

  std::sort(histories.begin(), histories.end(), [](const auto& a, const auto& b) {
    return std::abs(a.score) > std::abs(b.score);
  });
  if (static_cast<int>(histories.size()) > top_histories)
    histories.resize(top_histories);

  std::ofstream history_out(std::string(output_prefix) + ".histories.tsv");
  std::ofstream hit_out(std::string(output_prefix) + ".hits.tsv");
  history_out << "rank\tentry\tscore\tfraction_of_total\trate\n";
  hit_out << "rank\tentry\tevent_score\thit_index\tpid\ttrid\tmtrid\tk_mev"
             "\tx_mm\ty_mm\tz_mm\tvx_mm\tvy_mm\tvz_mm\tpx_mev\tpy_mev"
             "\tpz_mev\tcreator_volume\tcreator_process\tcreator_material\n";
  history_out << std::setprecision(20);
  hit_out << std::setprecision(20);
  std::map<std::tuple<std::string, std::string, std::string>, Breakdown> sources;
  for (std::size_t rank = 0; rank < histories.size(); ++rank) {
    const auto& history = histories[rank];
    history_out << rank + 1 << '\t' << history.entry << '\t' << history.score
                << '\t' << (sum != 0.0L ? history.score / sum : 0.0L) << '\t'
                << history.rate << '\n';
    tree->GetEntry(history.entry);
    const int count = std::min(nhit, kMaxHit);
    for (int index = 0; index < count; ++index) {
      if (det[index] != detector) continue;
      const long double hit_score = history.rate * scale *
          (energy ? std::max(k[index], 0.0) : 1.0);
      const bool primary = mtrid[index] == 0;
      const std::string kind = primary ? "primary" : "secondary";
      const std::string volume = primary ? "-" :
          Clean(creator_volume.get() + index * kString);
      const std::string process = primary ? "-" :
          Clean(creator_process.get() + index * kString);
      auto& source = sources[{kind, volume, process}];
      source.score += hit_score;
      ++source.hits;
      if (source.last_entry != history.entry) {
        ++source.histories;
        source.last_entry = history.entry;
      }
      hit_out << rank + 1 << '\t' << history.entry << '\t' << history.score
              << '\t' << index << '\t' << pid[index] << '\t' << trid[index] << '\t'
              << mtrid[index] << '\t' << k[index] << '\t' << x[index] << '\t' << y[index]
              << '\t' << z[index] << '\t' << vx[index] << '\t' << vy[index] << '\t'
              << vz[index] << '\t' << px[index] << '\t' << py[index] << '\t' << pz[index]
              << '\t' << volume << '\t' << process << '\t'
              << (primary ? "-" : Clean(
                  creator_material.get() + index * kString)) << '\n';
    }
  }
  history_out.close();
  hit_out.close();

  std::vector<std::pair<decltype(sources)::key_type, Breakdown>> ranked_sources(
      sources.begin(), sources.end());
  std::sort(ranked_sources.begin(), ranked_sources.end(), [](const auto& a, const auto& b) {
    return std::abs(a.second.score) > std::abs(b.second.score);
  });
  std::ofstream source_out(std::string(output_prefix) + ".sources.tsv");
  source_out << "kind\tcreator_volume\tcreator_process\tscore"
                "\tfraction_of_total\thits\thistories\n";
  source_out << std::setprecision(20);
  for (const auto& item : ranked_sources) {
    source_out << std::get<0>(item.first) << '\t' << std::get<1>(item.first)
               << '\t' << std::get<2>(item.first) << '\t' << item.second.score
               << '\t' << (sum != 0.0L ? item.second.score / sum : 0.0L)
               << '\t' << item.second.hits << '\t' << item.second.histories << '\n';
  }
  source_out.close();

  std::cout << "TILE_DIAGNOSTIC detector=" << detector
            << " observable=" << observable << " events=" << entries
            << " signal=" << signal << " rse=" << RelativeError(entries, sum, sum2)
            << " neff=" << (sum2 > 0.0L ? sum * sum / sum2 : 0.0L)
            << " top_histories=" << histories.size()
            << " output_prefix=" << output_prefix << '\n';
}
