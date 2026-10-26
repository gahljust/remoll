#include "remollNeymanPlan.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace {
void require(bool test, const std::string& message) {
    if (!test) throw std::runtime_error("Neyman: " + message);
}
std::string line(std::istream& input) {
    std::string value;
    require(bool(std::getline(input, value)), "incomplete plan");
    if (!value.empty() && value.back() == '\r') value.pop_back();
    return value;
}
struct Row {
    int ip, it, group;
    double plo, phi, tlo, thi, probability;
    long long quota, pilotCount;
};
}

remollNeymanPlan& remollNeymanPlan::Instance() {
    static remollNeymanPlan plan;
    return plan;
}

void remollNeymanPlan::Load(const std::string& path) {
    std::lock_guard<std::mutex> lock(fMutex);
    require(!fRunning, "cannot load a plan during a run");
    std::ifstream input(path);
    require(bool(input), "cannot open " + path);
    std::ostringstream contents;
    contents << input.rdbuf();
    std::istringstream parsed(contents.str());
    require(line(parsed) == "remoll-neyman-quota-v1", "expected remoll-neyman-quota-v1");
    // Read how the probabilities were obtained and which generator setup they belong to.
    std::string key, expected, source, provenance, extra;
    std::istringstream context(line(parsed));
    context >> key >> std::ws;
    require(key == "proposal" && context.peek() == '"' && bool(context >> std::quoted(expected)) &&
            !expected.empty() && !(context >> extra), "invalid proposal identity");
    long long pilotTotal;
    std::istringstream origin(line(parsed));
    require(bool(origin >> key >> source >> pilotTotal) && key == "probabilities", "invalid probability provenance");
    origin >> std::ws;
    require(origin.peek() == '"' && bool(origin >> std::quoted(provenance)) && !provenance.empty() &&
            !(origin >> extra), "missing quoted probability source");
    require((source == "known" && pilotTotal == 0) || (source == "pilot" && pilotTotal >= 2),
            "probabilities must be known or pilot-derived with a history count");
    require(line(parsed) == "ip\titheta\tp_lo_mev\tp_hi_mev\ttheta_lo_mrad\ttheta_hi_mrad\tgroup\tstratum_probability\ttransport_quota\tpilot_count",
            "unexpected TSV header");
    // Keep the old plan until the new one passes all checks.
    std::vector<Row> rows;
    int np = 0, nt = 0, ng = 0;
    std::string text;
    while (std::getline(parsed, text)) {
        if (text.empty() || text == "\r") continue;
        std::istringstream stream(text);
        std::vector<std::string> fields;
        std::string token;
        while (stream >> token) fields.push_back(token);
        require(fields.size() == 10, "expected ten fields per row");
        auto real = [&](int i) {
            size_t used;
            const double value = std::stod(fields[i], &used);
            require(used == fields[i].size(), "invalid number");
            return value;
        };
        auto integer = [&](int i) {
            size_t used;
            const long long value = std::stoll(fields[i], &used);
            require(used == fields[i].size(), "invalid integer");
            return value;
        };
        const auto ip = integer(0), it = integer(1), group = integer(6);
        require(ip >= -1 && ip < 1000000 && it >= -1 && it < 1000000 && group >= 0 && group <= 1000000,
                "invalid index");
        Row row{int(ip), int(it), int(group), real(2), real(3), real(4), real(5), real(7), integer(8), integer(9)};
        require((ip == -1) == (it == -1), "outside row must have both indices -1");
        // Each group needs at least two events to estimate how much its results vary.
        require(std::isfinite(row.probability) && row.probability > 0 && row.probability <= 1 && row.quota >= 2,
                "invalid probability or quota (<2)");
        // Check that the fractions match the pilot counts. This does not include uncertainty in those estimates.
        if (source == "pilot") {
            require(row.pilotCount > 0 && row.pilotCount <= pilotTotal &&
                    std::abs(row.probability-double(row.pilotCount)/pilotTotal) <= 1e-12*row.probability,
                    "probability does not match pilot count/total");
        } else {
            require(row.pilotCount == 0, "known probabilities require pilot_count=0");
        }
        if (ip >= 0) {
            require(std::isfinite(row.plo) && std::isfinite(row.phi) && row.plo >= 0 && row.plo < row.phi &&
                    std::isfinite(row.tlo) && std::isfinite(row.thi) && row.tlo >= 0 && row.tlo < row.thi, "invalid cell bounds");
            np = std::max(np, row.ip+1);
            nt = std::max(nt, row.it+1);
        }
        ng = std::max(ng, row.group+1);
        rows.push_back(row);
        require(rows.size() <= 1000001, "grid exceeds one million cells");
    }
    require(np > 0 && nt > 0 && static_cast<long long>(np)*nt <= 1000000, "empty or oversized grid");
    std::vector<double> p(np+1, NAN), theta(nt+1, NAN), probability(ng, -1);
    std::vector<long long> quota(ng, -1), pilotCounts(ng, -1);
    std::vector<int> cell(np*nt, -1);
    int outside = -1;
    auto edge = [](double& old, double value) {
        require(std::isnan(old) || old == value, "inconsistent grid edges");
        old = value;
    };
    // Cells may share a group: its probability, quota and pilot count are repeated, not added.
    for (const auto& row : rows) {
        const int g = row.group;
        require(probability[g] < 0 || (probability[g] == row.probability && quota[g] == row.quota &&
                pilotCounts[g] == row.pilotCount), "inconsistent group");
        probability[g] = row.probability;
        quota[g] = row.quota;
        pilotCounts[g] = row.pilotCount;
        if (row.ip == -1) {
            require(outside < 0, "duplicate outside row");
            outside = g;
            continue;
        }
        const int index = row.ip*nt+row.it;
        require(cell[index] < 0, "duplicate cell");
        cell[index] = g;
        edge(p[row.ip], row.plo);
        edge(p[row.ip+1], row.phi);
        edge(theta[row.it], row.tlo);
        edge(theta[row.it+1], row.thi);
    }
    // Every candidate must map to a group, including points outside the rectangular grid.
    require(outside >= 0 && std::find(cell.begin(), cell.end(), -1) == cell.end(), "missing cells or outside group");
    require(std::all_of(probability.begin(), probability.end(), [](double x){return x > 0;}), "missing group");
    require(std::abs(std::accumulate(probability.begin(), probability.end(), 0.0L)-1) < 1e-12,
            "probabilities do not sum to one");
    long long totalQuota = 0, totalPilot = 0;
    for (size_t g = 0; g < quota.size(); ++g) {
        require(quota[g] <= std::numeric_limits<long long>::max()-totalQuota, "quota overflow");
        require(pilotCounts[g] <= std::numeric_limits<long long>::max()-totalPilot, "pilot count overflow");
        totalQuota += quota[g];
        totalPilot += pilotCounts[g];
    }
    require(totalPilot == pilotTotal, "pilot counts do not sum to pilot total");
    // Use the checked plan and keep a copy to save in the ROOT file.
    fP = p;
    fTheta = theta;
    fProbability = probability;
    fQuota = quota;
    fCell = cell;
    fOutside = outside;
    fExpectedProposal = expected;
    fText = contents.str();
    fLoaded = true;
}

std::string remollNeymanPlan::Text() const {
    std::lock_guard<std::mutex> lock(fMutex);
    return fText;
}
std::string remollNeymanPlan::Proposal() const {
    std::lock_guard<std::mutex> lock(fMutex);
    return fCurrentProposal;
}
void remollNeymanPlan::Enable(bool enabled) {
    std::lock_guard<std::mutex> lock(fMutex);
    require(!fRunning, "cannot toggle sampling during a run");
    require(!enabled || fLoaded, "load a plan before enabling sampling");
    fEnabled = enabled;
}
void remollNeymanPlan::BeginRun(long long candidates) {
    std::lock_guard<std::mutex> lock(fMutex);
    require(!fRunning, "run already active");
    // Start the counts over for each beamOn. Even a large run can leave rare groups short of their quotas.
    if (fEnabled) {
        require(fLoaded && candidates >= std::accumulate(fQuota.begin(), fQuota.end(), 0LL),
                "candidate budget smaller than total quota");
        fAccepted.assign(fQuota.size(), 0);
    }
    fCandidates = candidates;
    fProcessed = 0;
    fCurrentProposal.clear();
    fRunning = true;
}
void remollNeymanPlan::RegisterProposal(const std::string& proposal) {
    std::lock_guard<std::mutex> lock(fMutex);
    // All workers must use the same generator settings, and those settings must match the plan.
    require(fRunning && !proposal.empty(), "missing active proposal");
    require(fCurrentProposal.empty() || fCurrentProposal == proposal, "workers have different proposals");
    require(!fEnabled || fExpectedProposal == proposal, "plan does not match the generator/beam/target proposal");
    fCurrentProposal = proposal;
}
bool remollNeymanPlan::Accept(double p, double theta, double& correction) {
    correction = 1;
    if (!fEnabled) return true;
    std::lock_guard<std::mutex> lock(fMutex);
    require(fRunning && !fCurrentProposal.empty() && fProcessed < fCandidates, "selection outside validated run");
    ++fProcessed;
    // A bin includes its lower edge but not its upper edge. Missing or out-of-range values use the outside group.
    int group = fOutside;
    if (std::isfinite(p) && std::isfinite(theta) && p >= fP.front() && p < fP.back() && theta >= fTheta.front() && theta < fTheta.back()) {
        const auto ip = std::upper_bound(fP.begin(), fP.end(), p)-fP.begin()-1;
        const auto it = std::upper_bound(fTheta.begin(), fTheta.end(), theta)-fTheta.begin()-1;
        group = fCell[ip*(fTheta.size()-1)+it];
    }
    // Keep events in this group until its quota is full, whether or not they hit a detector.
    if (fAccepted[group] == fQuota[group]) return false;
    ++fAccepted[group];
    // Estimate how many events ordinary sampling would put in this group: total events times
    // the group's expected fraction. Divide by the quota to make up for sampling it more or less often.
    // GeneratePrimaries() in remollPrimaryGeneratorAction.cc applies this correction to the rate.
    correction = double(fCandidates)*fProbability[group]/fQuota[group];
    return true;
}
void remollNeymanPlan::EndRun(bool complete) {
    std::lock_guard<std::mutex> lock(fMutex);
    // The rate correction divides by the planned quota, so every quota must be filled.
    // Reject incomplete runs before finalizing their ROOT output.
    const bool valid = fRunning && complete && fProcessed == fCandidates && fAccepted == fQuota;
    fRunning = false;
    require(!fEnabled || valid, "incomplete candidates or underfilled quotas; output is invalid");
}
