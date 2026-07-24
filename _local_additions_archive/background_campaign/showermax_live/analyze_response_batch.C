#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

namespace {
struct Moment {
  long long nonzero = 0;
  long double sum = 0, sum2 = 0, sum3 = 0, sum4 = 0, maximum = 0;
  std::vector<long double> largest;
  void add(long double x) {
    if (x != 0) ++nonzero;
    sum += x; sum2 += x*x; sum3 += x*x*x; sum4 += x*x*x*x;
    x = std::abs(x);
    maximum = std::max(maximum, x);
    if (x == 0) return;
    auto at = std::lower_bound(largest.begin(), largest.end(), x,
                               std::greater<long double>());
    largest.insert(at, x);
    if (largest.size() > 201) largest.pop_back();
  }
};

struct Spectrum {
  static constexpr int bins = 60;
  long long hits = 0;
  long long ignored = 0;
  long double rate = 0;
  long double ignored_rate = 0;
  std::array<long long,bins> bin_hits{};
  std::array<long long,bins> bin_ignored{};
  std::array<long double,bins> bin_rate{};
  std::array<long double,bins> bin_ignored_rate{};
  void add(double energy, double weight, bool is_ignored) {
    constexpr double low = -3.0;       // 0.001 MeV
    constexpr double high = 4.30103;   // 20000 MeV
    const double log_energy = std::log10(std::max(1.0e-3, energy));
    const int bin = std::clamp(
        static_cast<int>((log_energy-low)/(high-low)*bins), 0, bins-1);
    ++hits;
    rate += weight;
    ++bin_hits[bin];
    bin_rate[bin] += weight;
    if (is_ignored) {
      ++ignored;
      ignored_rate += weight;
      ++bin_ignored[bin];
      bin_ignored_rate[bin] += weight;
    }
  }
};

const std::array<int, 8> pids = {11, -11, 22, 13, -13, 211, -211, 2112};
const std::array<const char*, 8> names =
    {"e-", "e+", "gamma", "mu-", "mu+", "pi-", "pi+", "neutron"};
using Coeff = std::array<double, 5>;
std::map<int, std::map<double, Coeff>> fits;

bool load_fits(const std::string& data_dir) {
  for (std::size_t particle = 0; particle < pids.size(); ++particle) {
    std::ifstream input(data_dir + "/fit_param_xy_" + names[particle] + "_ifarm.csv");
    if (!input) return false;
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream row(line);
      double energy;
      Coeff c{};
      if (row >> energy >> c[0] >> c[1] >> c[2] >> c[3] >> c[4])
        fits[pids[particle]][energy] = c;
    }
  }
  return true;
}

bool supported(int pid) { return fits.count(pid) != 0; }

bool response_energy_available(int pid, double energy) {
  if (!supported(pid)) return false;
  const auto& table = fits.at(pid);
  return energy >= table.begin()->first && energy <= table.rbegin()->first;
}

std::pair<double,double> local_xy(double x, double y) {
  constexpr double radius = 1100.0;
  constexpr int planes = 28;
  const double twopi = 2.0 * std::acos(-1.0);
  double phi = std::atan2(y, x);
  if (phi < 0) phi += twopi;
  const int plane = static_cast<int>(std::floor(phi / (twopi / planes) + 0.5)) % planes;
  const double angle = plane * twopi / planes;
  return {(x - radius*std::cos(angle))*std::cos(angle)
             + (y - radius*std::sin(angle))*std::sin(angle),
          -(x - radius*std::cos(angle))*std::sin(angle)
             + (y - radius*std::sin(angle))*std::cos(angle)};
}

double response(int det, int pid, double energy, double x, double y) {
  if (!response_energy_available(pid, energy)) return 0.0;
  const auto& table = fits.at(pid);
  auto upper = table.lower_bound(energy);
  if (upper == table.end()) upper = std::prev(table.end());
  auto lower = upper;
  if (upper != table.begin() && upper->first != energy) lower = std::prev(upper);
  const double low = lower->first;
  const double high = upper->first;
  const auto xy = local_xy(x, y);
  const double lx = xy.first;
  const double ly = xy.second;
  auto eval = [&](double e) {
    const auto& c = table.at(e);
    return 0.5*(c[0] + lx*c[1]) + 0.5*(c[2] + ly*c[3] + ly*ly*c[4]);
  };
  double value = eval(low);
  if (high != low) value += (eval(high) - value) * (energy - low) / (high - low);
  const double long_pass = ((det - 73000) % 4 == 3) ? 0.22 : 0.33;
  return std::max(0.0, value * long_pass);
}

std::string tail(const Moment& value) {
  std::ostringstream out;
  out << std::setprecision(20);
  for (std::size_t i = 0; i < value.largest.size(); ++i) {
    if (i) out << ',';
    out << value.largest[i];
  }
  return out.str().empty() ? "-" : out.str();
}

int ring_index(int detector) {
  if (!((detector >= 110000 && detector < 112000) ||
        (detector >= 150000 && detector < 152000))) return -1;
  const int suffix = detector % 100;
  if (suffix == 10) return 0;
  if (suffix == 20) return 1;
  if (suffix == 30) return 2;
  if (suffix == 40) return 3;
  if (suffix == 51 || suffix == 52 || suffix == 53) return 4;
  if (suffix == 60) return 5;
  return -1;
}

// Match the established remoll main-detector analysis partition exactly:
// within each 2*pi/7 period the regions are
// closed, transition, open, transition, closed.
int azimuth_region(double x, double y) {
  const double pi = std::acos(-1.0);
  const double period = 2.0*pi/7.0;
  double phi = std::atan2(y, x);
  if (phi < 0.0) phi += 2.0*pi;
  double local = std::fmod(phi, period);
  if (local < 0.0) local += period;
  if (local < pi/28.0 || local >= 7.0*pi/28.0) return 0; // closed
  if (local < 3.0*pi/28.0 || local >= 5.0*pi/28.0) return 1; // transition
  return 2; // open
}
}

void analyze_response_batch(const char* filename, int septant,
                            const char* response_data_dir,
                            long long start_entry = 0,
                            long long requested_events = -1) {
  if (septant < 0 || septant > 6 || !load_fits(response_data_dir)) {
    std::cout << "LIVE_ERROR response_configuration\n";
    return;
  }
  TFile input(filename, "READ");
  auto* tree = dynamic_cast<TTree*>(input.Get("T"));
  if (!tree) { std::cout << "LIVE_ERROR missing_tree\n"; return; }
  double rate = 0;
  std::vector<remollGenericDetectorHit_t>* hits = nullptr;
  tree->SetBranchAddress("rate", &rate);
  tree->SetBranchAddress("hit", &hits);
  const long long normalization_events = tree->GetEntries();
  start_entry = std::clamp(start_entry, 0LL, normalization_events);
  const long long events = requested_events < 0
      ? normalization_events - start_entry
      : std::min(requested_events, normalization_events - start_entry);
  if (events <= 0) { std::cout << "LIVE_ERROR empty_range\n"; return; }
  const long double display_scale =
      static_cast<long double>(normalization_events) / events;
  std::array<Moment,28> shower_tiles;
  std::array<Moment,4> group; // closed, transition, open, total
  std::array<long double,28> shower_display_weight{};
  std::array<long long,28> shower_hit_count{};
  std::map<int,Spectrum> spectra;
  static constexpr int global_bins = 120;
  std::map<int,long double> circle_heat;
  std::map<int,long double> full_plane_heat;
  std::array<Moment,7> rate_moments;
  std::array<std::array<Moment,3>,7> rate_region_moments;
  std::array<std::array<Moment,14>,6> ring_tile_moments;
  std::array<std::map<int,long double>,7> rate_heat;
  std::array<long long,7> rate_hits{};
  long long circle_hits = 0;
  long long full_plane_hits = 0;
  long long unsupported = 0;

  for (long long offset = 0; offset < events; ++offset) {
    tree->GetEntry(start_entry + offset);
    std::array<long double,28> event_pe{};
    std::array<long double,7> event_rate{};
    std::array<std::array<long double,3>,7> event_rate_region{};
    std::array<std::array<long double,14>,6> event_ring_tile{};
    for (const auto& hit : *hits) {
      const int ring = ring_index(hit.det);
      if (ring >= 0) {
        const int rx = std::clamp(
            static_cast<int>((hit.x + 2000.0)/4000.0*global_bins),
            0, global_bins-1);
        const int ry = std::clamp(
            static_cast<int>((hit.y + 2000.0)/4000.0*global_bins),
            0, global_bins-1);
        const int bin = ry*global_bins + rx;
        const long double weighted_rate = rate * display_scale;
        event_rate[ring] += rate * normalization_events;
        event_rate[6] += rate * normalization_events;
        const int region = azimuth_region(hit.x, hit.y);
        event_rate_region[ring][region] += rate * normalization_events;
        event_rate_region[6][region] += rate * normalization_events;
        const int base = hit.det >= 150000 ? 150000 : 110000;
        const int segment = (hit.det - base)/100 - 1;
        if (segment >= 0 && segment < 14)
          event_ring_tile[ring][segment] += rate * normalization_events;
        rate_heat[ring][bin] += weighted_rate;
        rate_heat[6][bin] += weighted_rate;
        ++rate_hits[ring];
        ++rate_hits[6];
      }
      if (hit.pz <= 0) continue;
      const int gx = std::clamp(
          static_cast<int>((hit.x + 2000.0)/4000.0*global_bins),
          0, global_bins-1);
      const int gy = std::clamp(
          static_cast<int>((hit.y + 2000.0)/4000.0*global_bins),
          0, global_bins-1);
      const int global_bin = gy*global_bins + gx;
      if (hit.det == 30) {
        full_plane_heat[global_bin] += rate * display_scale;
        ++full_plane_hits;
        continue;
      }
      if (hit.det < 73001 || hit.det > 73028) continue;
      const bool ignored_response = !response_energy_available(hit.pid, hit.e);
      spectra[hit.pid].add(hit.e, rate * display_scale, ignored_response);
      if (ignored_response) {
        ++unsupported;
        continue;
      }
      const double pe = response(hit.det, hit.pid, hit.e, hit.x, hit.y);
      circle_heat[global_bin] += rate * display_scale * pe;
      ++circle_hits;
      const int index = hit.det - 73001;
      event_pe[index] += pe;
      shower_display_weight[index] += rate * display_scale * pe;
      ++shower_hit_count[index];
    }
    std::array<long double,4> group_score{};
    for (int i = 0; i < 28; ++i) {
      const long double score = rate * normalization_events * event_pe[i];
      shower_tiles[i].add(score);
      const int sector = i + 1;
      if (sector % 4 == 1) group_score[0] += score;
      else if (sector % 4 == 3) group_score[2] += score;
      else group_score[1] += score;
      group_score[3] += score;
    }
    for (int index = 0; index < 4; ++index)
      group[index].add(group_score[index]);
    for (int index = 0; index < 7; ++index)
      rate_moments[index].add(event_rate[index]);
    for (int index = 0; index < 7; ++index)
      for (int region = 0; region < 3; ++region)
        rate_region_moments[index][region].add(
            event_rate_region[index][region]);
    for (int ring = 0; ring < 6; ++ring)
      for (int segment = 0; segment < 14; ++segment)
        ring_tile_moments[ring][segment].add(event_ring_tile[ring][segment]);
  }

  std::cout << std::setprecision(20) << "LIVE_BATCH " << events << ' '
            << unsupported << ' ' << normalization_events << ' '
            << start_entry << "\n";
  const std::array<const char*,4> group_names =
      {"closed", "transition", "open", "showermax_total"};
  auto print = [&](const char* kind, const char* name, int det, const Moment& m) {
    std::cout << "LIVE_STAT " << kind << ' ' << name << ' ' << det << ' '
              << m.nonzero << ' ' << m.sum << ' ' << m.sum2 << ' ' << m.sum3
              << ' ' << m.sum4 << ' ' << m.maximum << ' ' << tail(m) << "\n";
  };
  for (int i = 0; i < 28; ++i) {
    const std::string name = "shmax_" + std::to_string(73001+i);
    print("tile", name.c_str(), 73001+i, shower_tiles[i]);
    std::cout << "LIVE_TILE " << name << ' ' << 73001+i << ' '
              << shower_hit_count[i] << ' ' << shower_display_weight[i] << "\n";
  }
  for (int i = 0; i < 4; ++i)
    print("group", group_names[i], i, group[i]);
  for (const auto& [pid, spectrum] : spectra) {
    std::cout << "LIVE_SPECIES " << pid << ' ' << spectrum.hits << ' '
              << spectrum.ignored << ' ' << spectrum.rate << ' '
              << spectrum.ignored_rate << "\n";
    for (int bin = 0; bin < Spectrum::bins; ++bin)
      if (spectrum.bin_hits[bin] != 0)
        std::cout << "LIVE_ENERGY " << pid << ' ' << bin << ' '
                  << spectrum.bin_hits[bin] << ' '
                  << spectrum.bin_ignored[bin] << ' '
                  << spectrum.bin_rate[bin] << ' '
                  << spectrum.bin_ignored_rate[bin] << "\n";
  }
  std::cout << "LIVE_GLOBAL circle " << circle_hits << "\n";
  for (const auto& [bin, weight] : circle_heat)
    std::cout << "LIVE_GLOBAL_BIN circle " << bin%global_bins << ' '
              << bin/global_bins << ' ' << weight << "\n";
  std::cout << "LIVE_GLOBAL full_plane " << full_plane_hits << "\n";
  for (const auto& [bin, weight] : full_plane_heat)
    std::cout << "LIVE_GLOBAL_BIN full_plane " << bin%global_bins << ' '
              << bin/global_bins << ' ' << weight << "\n";
  const std::array<const char*,7> rate_names =
      {"ring1","ring2","ring3","ring4","ring5","ring6","main_detector"};
  const std::array<const char*,3> region_names =
      {"closed","transition","open"};
  for (int index = 0; index < 7; ++index) {
    print("group", rate_names[index], index, rate_moments[index]);
    for (int region = 0; region < 3; ++region) {
      const std::string name = std::string(rate_names[index]) + "_"
          + region_names[region];
      print("group", name.c_str(), index, rate_region_moments[index][region]);
    }
    std::cout << "LIVE_GLOBAL " << rate_names[index] << ' '
              << rate_hits[index] << "\n";
    for (const auto& [bin, weight] : rate_heat[index])
      std::cout << "LIVE_GLOBAL_BIN " << rate_names[index] << ' '
                << bin%global_bins << ' ' << bin/global_bins << ' '
                << weight << "\n";
  }
  for (int ring = 0; ring < 6; ++ring)
    for (int segment = 0; segment < 14; ++segment) {
      const std::string name = "ring" + std::to_string(ring+1)
          + "_segment" + (segment < 9 ? "0" : "") + std::to_string(segment+1);
      print("tile", name.c_str(), segment+1,
            ring_tile_moments[ring][segment]);
    }
}
