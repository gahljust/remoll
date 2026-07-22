// ancestry_attribution.C
//
// Answers "where do the secondaries that reach this detector come from?"
// for any remoll output file — the measured importance map that biasing,
// surface placement, and capture-threshold choices should be driven by.
//
// For every hit in the selected detector range, its rate weight is
// attributed to (creation volume, creator process, species class), plus
// a direct-primary category. Outputs, ranked by rate:
//   1. per creation volume (the "places that need secondary showers")
//   2. per (creation volume x process x species)
//   3. the kinetic-energy spectrum of contributing hits per volume
//      decade, to inform capture/recording thresholds
//
// Usage (any remoll file, e.g. the June sieve-out runs or
// showermax_acceptance.root):
//
//   root -l -b -q '_local_additions_archive/background_campaign/
//       adaptive_replay/ancestry_attribution.C(
//       "file.root", 30, 30, "attribution.tsv")'
//   root -l -b -q '_local_additions_archive/background_campaign/
//       adaptive_replay/ancestry_attribution.C(
//       "file.root", 70030, 72730, "attribution_stack.tsv")'
//
// Rate convention: one event-rate weight per selected hit (the standard
// remoll virtual-plane practice). Energy-weighted rate uses the hit
// kinetic energy in MeV.

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "include/remolltypes.hh"

namespace {

const char* SpeciesClass(int pid)
{
  switch (pid) {
    case 11: return "e-";
    case -11: return "e+";
    case 22: return "gamma";
    case 13: case -13: return "mu";
    case 211: case -211: return "pi";
    case 2112: return "neutron";
    case 2212: return "proton";
    default: return "other";
  }
}

struct Accumulator {
  double rate{0.0};
  double rate_energy_mev{0.0};
  unsigned long long hits{0ULL};
  // log10(KE/MeV) decade counters from <0 (index 0) to >=4 (index 5)
  double rate_by_decade[6]{};
};

void Fill(Accumulator& accumulator, double rate, double ke_mev)
{
  accumulator.rate += rate;
  accumulator.rate_energy_mev += rate * ke_mev;
  ++accumulator.hits;
  int decade = ke_mev <= 0.0 ? 0
      : static_cast<int>(std::floor(std::log10(ke_mev))) + 1;
  if (decade < 0) decade = 0;
  if (decade > 5) decade = 5;
  accumulator.rate_by_decade[decade] += rate;
}

} // namespace

void ancestry_attribution(
    const char* filename,
    int det_lo = 30,
    int det_hi = 30,
    const char* output = "attribution.tsv")
{
  TFile* file = TFile::Open(filename, "READ");
  if (file == nullptr || file->IsZombie()) {
    printf("ancestry_attribution: cannot open %s\n", filename);
    return;
  }
  TTree* tree = nullptr;
  file->GetObject("T", tree);
  if (tree == nullptr) {
    printf("ancestry_attribution: no tree T in %s\n", filename);
    return;
  }

  remollEvent_t* event = nullptr;
  double rate = 0.0;
  std::vector<remollGenericDetectorHit_t>* hits = nullptr;
  tree->SetBranchAddress("hit", &hits);
  if (tree->GetBranch("rate") != nullptr) {
    tree->SetBranchAddress("rate", &rate);
  }
  if (tree->GetBranch("ev") != nullptr) {
    tree->SetBranchAddress("ev", &event);
  }

  // volume -> accumulator; and (volume, process, species) -> accumulator
  std::map<std::string, Accumulator> by_volume;
  std::map<std::tuple<std::string, std::string, std::string>,
           Accumulator> detailed;
  Accumulator total;

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    if (hits == nullptr) continue;
    const double weight =
        (std::isfinite(rate) && rate != 0.0) ? rate : 1.0;
    for (const remollGenericDetectorHit_t& hit : *hits) {
      if (hit.det < det_lo || hit.det > det_hi) continue;
      // Primaries have mtrid == 0 and no creator process.
      const bool primary = hit.mtrid == 0;
      std::string volume = primary
          ? std::string("(direct primary)")
          : std::string(hit.creator_physvol_name);
      if (volume.empty()) volume = "(unknown volume)";
      std::string process = primary
          ? std::string("-")
          : std::string(hit.creator_process_name);
      if (process.empty()) process = "(unknown process)";
      const std::string species = SpeciesClass(hit.pid);

      Fill(total, weight, hit.k);
      Fill(by_volume[volume], weight, hit.k);
      Fill(detailed[std::make_tuple(volume, process, species)],
           weight, hit.k);
    }
  }

  if (total.hits == 0ULL) {
    printf("ancestry_attribution: no hits in det range %d..%d\n",
           det_lo, det_hi);
    return;
  }

  // Rank volumes by rate.
  std::vector<std::pair<std::string, const Accumulator*>> ranked;
  ranked.reserve(by_volume.size());
  for (const auto& entry : by_volume) {
    ranked.emplace_back(entry.first, &entry.second);
  }
  std::sort(ranked.begin(), ranked.end(),
            [](const auto& a, const auto& b) {
              return a.second->rate > b.second->rate;
            });

  std::ofstream out(output);
  out.precision(10);
  out << "# format=remoll-ancestry-attribution-v1\n"
      << "# file=" << filename
      << " det_lo=" << det_lo << " det_hi=" << det_hi << "\n"
      << "# total_rate_hz=" << total.rate
      << " total_rate_energy_mev_hz=" << total.rate_energy_mev
      << " total_hits=" << total.hits << "\n";

  out << "table\tvolume\tprocess\tspecies\trate_hz\trate_fraction"
      << "\tcumulative_fraction\trate_energy_mev_hz\thits"
      << "\trate_ke_lt1mev\trate_ke_1_10mev\trate_ke_10_100mev"
      << "\trate_ke_100mev_1gev\trate_ke_1_10gev\trate_ke_ge10gev\n";

  printf("\n=== Rate attribution by creation volume "
         "(det %d..%d, total %.4g Hz) ===\n", det_lo, det_hi, total.rate);
  double cumulative = 0.0;
  for (const auto& entry : ranked) {
    const Accumulator& accumulator = *entry.second;
    const double fraction = accumulator.rate / total.rate;
    cumulative += fraction;
    printf("%-45s %12.4g Hz  %6.2f%%  (cum %6.2f%%)\n",
           entry.first.c_str(), accumulator.rate,
           100.0 * fraction, 100.0 * cumulative);
    out << "volume\t" << entry.first << "\t-\t-\t"
        << accumulator.rate << '\t' << fraction << '\t' << cumulative
        << '\t' << accumulator.rate_energy_mev << '\t'
        << accumulator.hits;
    for (double decade : accumulator.rate_by_decade) {
      out << '\t' << decade;
    }
    out << '\n';
  }

  // Detailed rows, ranked.
  std::vector<std::pair<std::tuple<std::string, std::string, std::string>,
                        const Accumulator*>> detail_ranked;
  detail_ranked.reserve(detailed.size());
  for (const auto& entry : detailed) {
    detail_ranked.emplace_back(entry.first, &entry.second);
  }
  std::sort(detail_ranked.begin(), detail_ranked.end(),
            [](const auto& a, const auto& b) {
              return a.second->rate > b.second->rate;
            });
  cumulative = 0.0;
  printf("\n=== Top 25 (volume x process x species) contributions ===\n");
  int printed = 0;
  for (const auto& entry : detail_ranked) {
    const Accumulator& accumulator = *entry.second;
    const double fraction = accumulator.rate / total.rate;
    cumulative += fraction;
    if (printed++ < 25) {
      printf("%-38s %-16s %-8s %11.4g Hz  %6.2f%%\n",
             std::get<0>(entry.first).c_str(),
             std::get<1>(entry.first).c_str(),
             std::get<2>(entry.first).c_str(),
             accumulator.rate, 100.0 * fraction);
    }
    out << "detail\t" << std::get<0>(entry.first) << '\t'
        << std::get<1>(entry.first) << '\t'
        << std::get<2>(entry.first) << '\t'
        << accumulator.rate << '\t' << fraction << '\t' << cumulative
        << '\t' << accumulator.rate_energy_mev << '\t'
        << accumulator.hits;
    for (double decade : accumulator.rate_by_decade) {
      out << '\t' << decade;
    }
    out << '\n';
  }
  out.close();
  printf("\nancestry_attribution: wrote %s\n", output);
  printf("Use the per-volume ranking to (a) place surface-source "
         "recording planes\njust downstream of the top volumes, "
         "(b) assign importance/weight windows,\nand (c) read the KE "
         "decade columns to justify capture and recording\nthresholds "
         "(a threshold is harmless if the rate below it is "
         "negligible).\n");
}
