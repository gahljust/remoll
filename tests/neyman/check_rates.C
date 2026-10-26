// Use a Geant4 MT build with per-event reseeding and identical seed/setup/beamOn.
// Sequential builds require seed replay; quota rejection changes their RNG stream.
// Load the built libremoll before executing this macro.
#include <TFile.h>
#include <TTree.h>
#include <TObjString.h>
#include <TSystem.h>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>

void check(bool condition, const char* message) {
    if (!condition) {
        std::cerr << message << std::endl;
        gSystem->Exit(1);
    }
}
struct Event { double rate, xs, p, theta; };
std::map<int, Event> events(TFile& file) {
    auto tree = file.Get<TTree>("T");
    check(tree != nullptr, "missing event tree");
    tree->SetEstimate(tree->GetEntries()*4+1);
    check(tree->Draw("rate:ev.xs:seed.fEvtNo", "", "goff") == tree->GetEntries(), "failed to read events");
    std::map<int, Event> result;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i)
        check(result.emplace(int(tree->GetV3()[i]), Event{tree->GetV1()[i], tree->GetV2()[i], NAN, NAN}).second,
              "duplicate event ID");
    const auto particles = tree->Draw("part.p:atan2(sqrt(part.px*part.px+part.py*part.py),part.pz)*1000:seed.fEvtNo:part.pz",
                                      "part.pid==11", "goff");
    check(particles >= 0 && particles < tree->GetEstimate(), "failed to read particle coordinates");
    std::map<int, bool> forward;
    for (Long64_t i = 0; i < particles; ++i) {
        const int id = tree->GetV3()[i];
        auto& e = result.at(id);
        const bool f = tree->GetV4()[i] > 0;
        if (std::isnan(e.p) || (f && !forward[id]) || (f == forward[id] && tree->GetV1()[i] > e.p)) {
            e.p = tree->GetV1()[i];
            e.theta = tree->GetV2()[i];
            forward[id] = f;
        }
    }
    return result;
}
void check_rates(const char* ordinary, const char* sampled) {
    TFile baseline(ordinary), selected(sampled);
    check(!baseline.IsZombie() && !selected.IsZombie(), "cannot open inputs");
    auto a = baseline.Get<remollRunData>("run_data");
    auto b = selected.Get<remollRunData>("run_data");
    check(a && b && a->GetNthrown() == b->GetNthrown(), "candidate counts differ");
    const double n = a->GetNthrown();
    auto plan = selected.Get<TObjString>("neyman_plan");
    auto pa = baseline.Get<TObjString>("generator_proposal");
    auto pb = selected.Get<TObjString>("generator_proposal");
    check(plan && pa && pb && pa->GetString() == pb->GetString(), "missing plan or changed proposal");
    struct Cell { int group; double plo, phi, tlo, thi; };
    std::vector<Cell> cells;
    std::map<int, std::pair<double, long long>> groups;
    std::map<int, long long> counts;
    std::istringstream input(plan->GetString().Data());
    std::string line;
    for (int i = 0; i < 4; ++i) std::getline(input, line);
    int outside = -1;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        int ip, it, g;
        std::string plo, phi, tlo, thi;
        double probability;
        long long quota, pilot;
        std::istringstream row(line);
        check(bool(row >> ip >> it >> plo >> phi >> tlo >> thi >> g >> probability >> quota >> pilot), "bad saved plan");
        groups[g] = {probability, quota};
        if (ip < 0) outside = g;
        else cells.push_back({g, std::stod(plo), std::stod(phi), std::stod(tlo), std::stod(thi)});
    }
    const auto reference = events(baseline), kept = events(selected);
    check(reference.size() == n && !kept.empty(), "ordinary histories missing or empty selection");
    for (const auto& event : kept) {
        auto original = reference.find(event.first);
        check(original != reference.end(), "unknown candidate ID");
        const auto& e = event.second;
        int g = outside;
        for (const auto& c : cells)
            if (e.p >= c.plo && e.p < c.phi && e.theta >= c.tlo && e.theta < c.thi) { g = c.group; break; }
        check(groups.count(g), "unmapped event");
        ++counts[g];
        const double expected = original->second.rate * n * groups[g].first/groups[g].second;
        check(e.xs == original->second.xs, "ev.xs changed");
        check(std::isfinite(e.rate) && std::abs(e.rate-expected) <= 1e-12*std::abs(expected), "incorrect saved rate");
    }
    for (const auto& g : groups) check(counts[g.first] == g.second.second, "incorrect group quota");
    std::cout << "PASS " << kept.size() << " actual ROOT histories: cross sections, quota rates and group counts" << std::endl;
}
