#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "TFile.h"
#include "TTree.h"

namespace {
struct Moment {
  static constexpr std::size_t tail_size = 201;
  long long signal_events = 0;
  long double sum = 0.0L;
  long double sum2 = 0.0L;
  long double sum3 = 0.0L;
  long double sum4 = 0.0L;
  long double maximum = 0.0L;
  std::vector<long double> largest;

  void add(long double value) {
    if (value != 0.0L) ++signal_events;
    sum += value;
    sum2 += value * value;
    sum3 += value * value * value;
    sum4 += value * value * value * value;
    const long double magnitude = std::abs(value);
    maximum = std::max(maximum, magnitude);
    if (magnitude == 0.0L) return;
    const auto position = std::lower_bound(largest.begin(), largest.end(), magnitude,
                                           std::greater<long double>());
    largest.insert(position, magnitude);
    if (largest.size() > tail_size) largest.pop_back();
  }
};

struct CellMoment {
  long long sample_events = 0;
  Moment response;
};

struct Cell {
  int index;
  double low;
  double high;
};

std::string detector_group(int detector) {
  if (detector >= 110000 && detector < 152000) {
    const int suffix = detector % 100;
    if (suffix == 10) return "ring1";
    if (suffix == 20) return "ring2";
    if (suffix == 30) return "ring3";
    if (suffix == 40) return "ring4";
    if (suffix == 51 || suffix == 52 || suffix == 53) return "ring5";
    if (suffix == 60) return "ring6";
  }
  if (detector >= 70030 && detector <= 72730 && detector % 100 == 30) {
    const int stack = (detector - 70030) / 100;
    if (stack % 4 == 0) return "showermax_open";
    if (stack % 4 == 2) return "showermax_closed";
    return "showermax_transition";
  }
  if (detector >= 32 && detector <= 35) return "gem_diagnostic";
  if (detector >= 120 && detector <= 131) return "ring_virtual_diagnostic";
  return "";
}

bool precision_detector(int detector) {
  const std::string group = detector_group(detector);
  return group.rfind("ring", 0) == 0 || group.rfind("showermax_", 0) == 0;
}

std::vector<int> precision_detectors() {
  std::vector<int> detectors;
  static constexpr std::array<int, 8> suffixes = {10, 20, 30, 40, 51, 52, 53, 60};
  for (const int base : {110000, 150000}) {
    for (int segment = 1; segment <= 14; ++segment) {
      for (const int suffix : suffixes) detectors.push_back(base + 100 * segment + suffix);
    }
  }
  for (int detector = 70030; detector <= 72730; detector += 100) {
    detectors.push_back(detector);
  }
  return detectors;
}

std::string tail_values(const Moment& moment) {
  if (moment.largest.empty()) return "-";
  std::ostringstream output;
  output << std::setprecision(20);
  for (std::size_t index = 0; index < moment.largest.size(); ++index) {
    if (index != 0) output << ',';
    output << moment.largest[index];
  }
  return output.str();
}

Cell equal_cell(double value, double minimum, double maximum, int count) {
  const double clipped = std::clamp(value, minimum, std::nextafter(maximum, minimum));
  const int index = std::clamp(
      static_cast<int>((clipped - minimum) / (maximum - minimum) * count), 0, count - 1);
  const double width = (maximum - minimum) / count;
  return {index, minimum + index * width, minimum + (index + 1) * width};
}

Cell beam_cell(double value, double beam_gev) {
  static constexpr std::array<double, 7> edges = {0.0, 0.05, 0.20, 0.50, 0.80, 0.95, 1.0};
  const double fraction = std::clamp(value / beam_gev, 0.0, std::nextafter(1.0, 0.0));
  int index = 0;
  while (index + 1 < static_cast<int>(edges.size()) - 1 && fraction >= edges[index + 1]) {
    ++index;
  }
  return {index, edges[index] * beam_gev, edges[index + 1] * beam_gev};
}
}

void analyze_batch(const char* filename, double beam_gev, double theta_min_deg,
                   double theta_max_deg, double z_min_mm, double z_max_mm) {
  TFile input(filename, "READ");
  auto* tree = dynamic_cast<TTree*>(input.Get("T"));
  if (!tree) {
    std::cout << "CAMPAIGN_ERROR missing_tree\n";
    return;
  }

  double rate = 0.0;
  remollEvent_t* event = nullptr;
  remollBeamTarget_t* beam_target = nullptr;
  std::vector<remollEventParticle_t>* particles = nullptr;
  std::vector<remollGenericDetectorHit_t>* hits = nullptr;
  tree->SetBranchAddress("rate", &rate);
  tree->SetBranchAddress("ev", &event);
  tree->SetBranchAddress("bm", &beam_target);
  tree->SetBranchAddress("part", &particles);
  tree->SetBranchAddress("hit", &hits);

  using Key = std::tuple<std::string, std::string, std::string>;
  using CellKey = std::tuple<std::string, int, std::string, std::string>;
  std::map<Key, Moment> moments;
  std::map<CellKey, CellMoment> cells;
  std::map<std::pair<std::string, int>, Cell> cell_bounds;
  const long long events = tree->GetEntries();
  const double pi = std::acos(-1.0);
  const std::vector<std::string> groups = {
    "ring1", "ring2", "ring3", "ring4", "ring5", "ring6",
    "showermax_open", "showermax_closed", "showermax_transition",
    "gem_diagnostic", "ring_virtual_diagnostic"
  };
  const std::vector<int> detector_ids = precision_detectors();
  for (const int detector : detector_ids) {
    const std::string tile = "tile_" + std::to_string(detector);
    moments[{tile, "all", "rate"}];
    if (detector_group(detector).rfind("showermax_", 0) == 0) {
      moments[{tile, "all", "rate_energy"}];
    }
  }
  std::vector<std::pair<std::string, Cell>> defined_cells;
  for (int index = 0; index < 8; ++index) {
    const double width = (theta_max_deg - theta_min_deg) / 8.0;
    defined_cells.push_back({"theta", {index, theta_min_deg + index * width,
                                        theta_min_deg + (index + 1) * width}});
  }
  for (int index = 0; index < 14; ++index) {
    defined_cells.push_back({"phi", {index, index * 360.0 / 14.0,
                                      (index + 1) * 360.0 / 14.0}});
  }
  static constexpr std::array<double, 7> beam_edges = {0.0, 0.05, 0.20, 0.50, 0.80, 0.95, 1.0};
  for (int index = 0; index < 6; ++index) {
    defined_cells.push_back({"beamp", {index, beam_edges[index] * beam_gev,
                                        beam_edges[index + 1] * beam_gev}});
  }
  for (int index = 0; index < 4; ++index) {
    defined_cells.push_back({"vertexz", {index, index / 4.0, (index + 1) / 4.0}});
  }
  for (int index = 0; index < 6; ++index) {
    defined_cells.push_back({"outgoinge", {index, beam_edges[index], beam_edges[index + 1]}});
  }
  for (const auto& [axis, cell] : defined_cells) {
    cell_bounds[{axis, cell.index}] = cell;
    for (const auto& group : groups) {
      for (const auto& observable : {"rate", "rate_energy"}) {
        cells[{axis, cell.index, group, observable}];
      }
    }
  }

  for (long long entry = 0; entry < events; ++entry) {
    tree->GetEntry(entry);
    std::map<Key, long double> response;
    for (const auto& hit : *hits) {
      const std::string group = detector_group(hit.det);
      if (group.empty()) continue;
      const std::string ancestry = hit.mtrid == 0 ? "primary" : "secondary";
      for (const std::string& category : {std::string("all"), ancestry}) {
        response[{group, category, "rate"}] += rate * events;
        response[{group, category, "rate_energy"}] += rate * events * hit.k;
      }
      if (precision_detector(hit.det)) {
        const std::string tile = "tile_" + std::to_string(hit.det);
        response[{tile, "all", "rate"}] += rate * events;
        if (group.rfind("showermax_", 0) == 0) {
          response[{tile, "all", "rate_energy"}] += rate * events * hit.k;
        }
      }
    }

    const double theta_deg = event->thcom * 180.0 / pi;
    double phi_deg = particles->empty() ? 0.0 : particles->front().ph * 180.0 / pi;
    phi_deg = std::fmod(phi_deg, 360.0);
    if (phi_deg < 0.0) phi_deg += 360.0;
    const double vertex_fraction = z_max_mm > z_min_mm
        ? (beam_target->z - z_min_mm) / (z_max_mm - z_min_mm) : 0.5;
    const double proton_mass_mev = 938.27208816;
    const double outgoing_maximum = proton_mass_mev * event->beamp
        / (proton_mass_mev + event->beamp * (1.0 - std::cos(event->thcom)));
    const double outgoing_fraction = particles->empty() || outgoing_maximum <= 0.0
        ? 0.0 : particles->front().p / outgoing_maximum;
    const std::array<std::pair<std::string, Cell>, 5> event_cells = {{
      {"theta", equal_cell(theta_deg, theta_min_deg, theta_max_deg, 8)},
      {"phi", equal_cell(phi_deg, 0.0, 360.0, 14)},
      {"beamp", beam_cell(event->beamp / 1000.0, beam_gev)},
      {"vertexz", equal_cell(vertex_fraction, 0.0, 1.0, 4)},
      {"outgoinge", beam_cell(outgoing_fraction, 1.0)},
    }};

    for (const auto& group : groups) {
      for (const auto& ancestry : {"all", "primary", "secondary"}) {
        for (const auto& observable : {"rate", "rate_energy"}) {
          const Key key{group, ancestry, observable};
          moments[key].add(response[key]);
        }
      }
      for (const auto& observable : {"rate", "rate_energy"}) {
        const long double value = response[{group, "all", observable}];
        for (const auto& [axis, cell] : event_cells) {
          auto& target = cells[{axis, cell.index, group, observable}];
          ++target.sample_events;
          target.response.add(value);
          cell_bounds[{axis, cell.index}] = cell;
        }
      }
    }
    for (const int detector : detector_ids) {
      const std::string tile = "tile_" + std::to_string(detector);
      moments[{tile, "all", "rate"}].add(response[{tile, "all", "rate"}]);
      if (detector_group(detector).rfind("showermax_", 0) == 0) {
        moments[{tile, "all", "rate_energy"}].add(
            response[{tile, "all", "rate_energy"}]);
      }
    }
  }

  std::cout << std::setprecision(20);
  std::cout << "CAMPAIGN_BATCH " << events << "\n";
  for (const auto& item : moments) {
    const auto& [group, ancestry, observable] = item.first;
    const auto& value = item.second;
    std::cout << "CAMPAIGN_STAT " << group << ' ' << ancestry << ' '
              << observable << ' ' << value.signal_events << ' '
              << value.sum << ' ' << value.sum2 << ' ' << value.sum3 << ' '
              << value.sum4 << ' ' << value.maximum << ' ' << tail_values(value) << "\n";
  }
  for (const auto& item : cells) {
    const auto& [axis, index, group, observable] = item.first;
    const auto& value = item.second;
    const auto& bounds = cell_bounds.at({axis, index});
    std::cout << "CAMPAIGN_CELL " << axis << ' ' << index << ' '
              << bounds.low << ' ' << bounds.high << ' ' << group << ' '
              << observable << ' ' << value.sample_events << ' '
              << value.response.signal_events << ' ' << value.response.sum << ' '
              << value.response.sum2 << ' ' << value.response.maximum << "\n";
  }
}
