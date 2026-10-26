// Run with build/reroot from the repository root; see example.C.
// Use complete, independent ordinary runs from this branch, with identical settings.
// Each response expression gives one number per event; the event rate is applied below.
// Keep misses (response zero). No event cuts or already-sampled inputs are allowed.
// Pool sparse bins, reserve two events per group,
// then quotas chosen to balance the relative variance of the requested measurements.
// Pilot fractions are estimates. Uncertainty in those fractions is not included here.
#ifndef REMOLL_MAKE_NEYMAN_PLAN_C
#define REMOLL_MAKE_NEYMAN_PLAN_C
#ifndef __CLING__
#include "../../include/remollRunData.hh"
#endif
#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TUUID.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace neyman_pilot {
// Stop with an explanation instead of building a plan from invalid inputs.
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error("Neyman pilot: " + message);
}
// Choose the momentum/angle grid and how much pilot data a bin needs to stand on its own.
// Bins with too little data share a quota. includeAsymmetry adds asymmetry to the rate goals.
struct Settings {
    int momentumBins = 24, angleBins = 24;
    double maxMomentumMeV = 11100, minAngleMrad = 1e-4, maxAngleMrad = 2000;
    double minimumFraction = .001;
    long long minimumPilotCount = 100;
    bool includeAsymmetry = true;
};
// Keep means and centered sums for each response, so large rates do not hide small variations.
// d and n are the rate contribution and its asymmetry numerator. Their centered sums
// measure how each varies, and how they vary together; add() combines bins without rereading events.
struct Moment {
    long long count = 0;
    std::vector<double> d, n, d2, n2, dn;
    explicit Moment(size_t size) : d(size), n(size), d2(size), n2(size), dn(size) {}
    void fill(const std::vector<double>& contribution, double asymmetry) {
        ++count;
        for (size_t j = 0; j < d.size(); ++j) {
            const double x = contribution[j], y = x*asymmetry;
            require(std::isfinite(x) && std::isfinite(y), "nonfinite event contribution");
            const double dx = x-d[j], dy = y-n[j];
            d[j] += dx/count; n[j] += dy/count;
            d2[j] += dx*(x-d[j]); n2[j] += dy*(y-n[j]); dn[j] += dx*(y-n[j]);
        }
    }
    // Merge centered sums, including the difference between the two bins' means.
    void add(const Moment& other) {
        if (other.count == 0) return;
        if (count == 0) { *this = other; return; }
        const double fraction = double(other.count)/(count+other.count);
        const double cross = count*fraction;
        for (size_t j = 0; j < d.size(); ++j) {
            const double dx = other.d[j]-d[j], dy = other.n[j]-n[j];
            d2[j] += other.d2[j]+dx*dx*cross;
            n2[j] += other.n2[j]+dy*dy*cross;
            dn[j] += other.dn[j]+dx*dy*cross;
            d[j] += dx*fraction; n[j] += dy*fraction;
        }
        count += other.count;
    }
};
// probability is each group's fraction of pilot events. variance measures spread across
// the whole pilot; within measures spread inside each group, for each requested measurement.
struct Objectives {
    std::vector<double> probability, variance;
    std::vector<std::vector<double>> within;
};
// Turn the accumulated totals into the fractions and spreads used by allocate() below.
Objectives objectives(const std::vector<Moment>& groups, bool asymmetry) {
    const size_t responses = groups.front().d.size();
    const size_t size = responses * (asymmetry ? 2 : 1);
    Moment total(responses);
    for (const auto& group : groups) total.add(group);
    Objectives result;
    result.within.resize(groups.size(), std::vector<double>(size));
    result.variance.resize(size);
    // For rate, measure the spread of event contributions. For asymmetry, subtract the
    // overall numerator/rate ratio times each rate contribution from its numerator.
    // This accounts for the numerator and denominator changing together.
    auto variance = [](const Moment& m, size_t j, double a, bool ratio) {
        // These sums are already centered on the mean; dividing by count gives the variance.
        const double square = ratio ? m.n2[j] - 2*a*m.dn[j] + a*a*m.d2[j] : m.d2[j];
        const double v = square/m.count;
        require(std::isfinite(v), "nonfinite variance; check the responses and rates");
        return std::max(0.0, v);
    };
    // Use the same overall ratio for every group: they contribute to one asymmetry measurement.
    for (size_t j = 0; j < size; ++j) {
        const size_t response = j % responses;
        const double a = total.d[response] ? total.n[response]/total.d[response] : 0;
        result.variance[j] = variance(total, response, a, j >= responses);
        for (size_t h = 0; h < groups.size(); ++h)
            result.within[h][j] = variance(groups[h], response, a, j >= responses);
    }
    // Count all pilot events, including misses. These fractions also enter remoll's rate correction.
    for (const auto& group : groups) result.probability.push_back(double(group.count)/total.count);
    return result;
}
// Split the tracking budget between groups. weight sets the emphasis on each measurement.
std::vector<int> allocate(const Objectives& data, int budget, const std::vector<double>& weight) {
    std::vector<double> score(data.probability.size());
    for (size_t h = 0; h < score.size(); ++h) {
        double combined = 0;
        // Compare each measurement with its own overall spread, so its units do not dominate.
        for (size_t j = 0; j < weight.size(); ++j)
            if (data.variance[j] > 0 && weight[j] > 0)
                combined += weight[j]*data.within[h][j]/data.variance[j];
        // Multiply the group's pilot fraction by its combined spread to set its share of events.
        // For one measurement, the shares are proportional to pilot fraction times standard deviation.
        score[h] = data.probability[h]*std::sqrt(std::max(0.0, combined));
    }
    double sum = std::accumulate(score.begin(), score.end(), 0.0);
    // If the pilot sees no spread within any group, use its fractions for the spare slots.
    if (sum == 0) { score = data.probability; sum = 1; }
    // Reserve two events per group so production can estimate its variation.
    // Divide the remaining budget according to the scores above.
    std::vector<int> quota(score.size(), 2);
    const int remaining = budget - 2*int(score.size());
    require(remaining >= 0, "budget cannot give every group two events");
    std::vector<std::pair<double,int>> remainder;
    int assigned = 0;
    for (size_t h = 0; h < score.size(); ++h) {
        const double exact = remaining*score[h]/sum;
        const int whole = int(std::floor(exact));
        quota[h] += whole; assigned += whole;
        remainder.emplace_back(exact-whole, h);
    }
    // Quotas must be whole numbers. Give leftover slots to the largest rounding remainders
    // so the quotas add up to exactly the requested budget.
    std::sort(remainder.begin(), remainder.end(), std::greater<std::pair<double,int>>());
    for (int i = 0; i < remaining-assigned; ++i) ++quota[remainder[i].second];
    return quota;
}
// Start with equal emphasis on measurements with observed variation, then balance their expected precision.
std::vector<int> balanced_allocation(const Objectives& data, int budget) {
    std::vector<double> weight(data.variance.size());
    for (size_t j = 0; j < weight.size(); ++j) weight[j] = data.variance[j] > 0 ? 1 : 0;
    double sum = std::accumulate(weight.begin(), weight.end(), 0.0);
    require(sum > 0, "pilot has no measured variation; check responses or collect more events");
    for (double& w : weight) w /= sum;
    // Shift emphasis toward the least-improved measurement.
    // With one measurement this reduces to ordinary Neyman allocation plus minimum quotas.
    for (int iteration = 0; iteration < 60; ++iteration) {
        const auto quota = allocate(data, budget, weight);
        // Predict each measurement's variance with these quotas, relative to ordinary sampling
        // with the same tracking budget. A larger ratio means less improvement.
        std::vector<double> ratio(weight.size());
        for (size_t j = 0; j < weight.size(); ++j) {
            if (weight[j] == 0) continue;
            for (size_t h = 0; h < quota.size(); ++h)
                ratio[j] += data.probability[h]*data.probability[h]*data.within[h][j]/quota[h];
            ratio[j] *= budget/data.variance[j];
        }
        // Give the least-improved measurements more emphasis in the next pass.
        // Use 60 passes, with 0.20 setting the strength of each adjustment; this is an iterative balance.
        const double worst = *std::max_element(ratio.begin(), ratio.end());
        if (worst <= 0) break;
        sum = 0;
        for (size_t j = 0; j < weight.size(); ++j) {
            weight[j] *= std::exp(.20*ratio[j]/worst);
            sum += weight[j];
        }
        for (double& w : weight) w /= sum;
    }
    return allocate(data, budget, weight);
}
}

// Read ordinary pilots, group their events, calculate quotas, and write the plan remoll loads.
// budget is the number of events to track in production, not the number of generator attempts.
void make_plan(const std::vector<std::string>& pilots, const char* output, int budget,
               const std::vector<std::string>& responses,
               neyman_pilot::Settings settings = neyman_pilot::Settings()) {
    using namespace neyman_pilot;
    // Check the requested grid and budget before reading data; do not overwrite an existing plan.
    require(!pilots.empty() && !responses.empty() && budget >= 4, "supply pilots, responses and a budget >= 4");
    require(settings.momentumBins > 0 && settings.angleBins > 0 &&
            1LL*settings.momentumBins*settings.angleBins <= 1000000 &&
            std::isfinite(settings.maxMomentumMeV) && settings.maxMomentumMeV > 0 &&
            settings.minAngleMrad > 0 && settings.minAngleMrad < settings.maxAngleMrad &&
            std::isfinite(settings.maxAngleMrad) && settings.minimumPilotCount >= 2 &&
            settings.minimumFraction > 0 && settings.minimumFraction <= 1, "invalid binning or grouping settings");
    require(output && *output && !std::ifstream(output).good(), "output already exists or has no name");
    // Use evenly spaced momentum bins and logarithmically spaced angle bins, which give
    // finer detail at small angles. The extra outside bin keeps missing/out-of-range electrons.
    std::vector<double> p(settings.momentumBins+1), theta(settings.angleBins+1);
    for (size_t i = 0; i < p.size(); ++i) p[i] = settings.maxMomentumMeV*(double(i)/settings.momentumBins);
    for (size_t i = 0; i < theta.size(); ++i)
        theta[i] = settings.minAngleMrad*std::pow(settings.maxAngleMrad/settings.minAngleMrad,
                                               double(i)/settings.angleBins);
    // Reject collapsed or unrepresentable edges before writing a plan the reader cannot load.
    for (const auto* edges : {&p, &theta})
        for (size_t i = 0; i < edges->size(); ++i)
            require(((*edges)[i] == 0 || std::isnormal((*edges)[i])) &&
                    (i == 0 || (*edges)[i] > (*edges)[i-1]), "bin edges must be distinct and representable");
    const int outside = settings.momentumBins*settings.angleBins;
    std::vector<Moment> cells(outside+1, Moment(responses.size()));
    std::string proposal;
    std::set<std::string> files;
    long long count = 0;
    double referenceAsymmetry = 0;
    // Save the response definitions and input file identities so the plan can be traced back.
    std::ostringstream provenance;
    provenance << "responses=";
    for (const auto& expression : responses) provenance << '[' << expression << ']';
    provenance << "; asymmetry=" << settings.includeAsymmetry << "; pilots=";
    // Only combine complete ordinary pilots with matching generator settings. Counting a file
    // twice or using already-sampled events would give the wrong fractions and uncertainty estimates.
    for (const auto& path : pilots) {
        TFile file(path.c_str(), "READ");
        require(!file.IsZombie() && files.insert(file.GetUUID().AsString()).second, "missing or duplicate pilot: " + path);
        require(!file.Get("neyman_plan"), "pilot is already sampled: " + path);
        auto* tree = file.Get<TTree>("T");
        auto* run = file.Get<remollRunData>("run_data");
        auto* identity = file.Get<TObjString>("generator_proposal");
        require(tree && run && identity && tree->GetEntries() >= 2 &&
                static_cast<unsigned long long>(tree->GetEntries()) == run->GetNthrown(),
                "need a complete ordinary pilot with generator_proposal: " + path);
        const std::string current = identity->GetString().Data();
        require(!current.empty() && (proposal.empty() || proposal == current), "pilot generator settings differ");
        proposal = current;
        provenance << path << " (" << file.GetUUID().AsString() << ");";
        double rate = 0;
        remollEvent_t* event = nullptr;
        std::vector<remollEventParticle_t>* particles = nullptr;
        require(tree->SetBranchAddress("rate", &rate) >= 0 && tree->SetBranchAddress("ev", &event) >= 0 &&
                tree->SetBranchAddress("part", &particles) >= 0, "missing event branches");
        // Prepare each response once, then evaluate it for every event below.
        // A response can be a hit count, PE yield, or a 0/1 flag for reaching a detector.
        std::vector<std::unique_ptr<TTreeFormula>> formulas;
        for (size_t j = 0; j < responses.size(); ++j) {
            formulas.emplace_back(new TTreeFormula(("response"+std::to_string(j)).c_str(), responses[j].c_str(), tree));
            require(formulas.back()->GetNdim() > 0, "invalid response expression: " + responses[j]);
        }
        std::vector<double> contribution(responses.size());
        // Read every generated event. A detector miss still contributes to the pilot counts.
        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            require(tree->GetEntry(i) > 0 && event && particles, "unreadable pilot event");
            require(std::isfinite(rate) && rate >= 0 && (!settings.includeAsymmetry || std::isfinite(event->A)),
                    "invalid event rate or asymmetry");
            // Match GeneratePrimaries() in remollPrimaryGeneratorAction.cc: use the highest-momentum
            // forward electron, or the highest-momentum electron if none is forward. This only locates the bin.
            const remollEventParticle_t* forward = nullptr;
            const remollEventParticle_t* any = nullptr;
            for (const auto& particle : *particles) if (particle.pid == 11) {
                if (!any || particle.p > any->p) any = &particle;
                if (particle.pz > 0 && (!forward || particle.p > forward->p)) forward = &particle;
            }
            const auto* electron = forward ? forward : any;
            const double momentum = electron ? electron->p : NAN;
            const double angle = electron ? std::atan2(std::hypot(electron->px, electron->py), electron->pz)*1000 : NAN;
            // Match Accept() in remollNeymanPlan.cc: include the lower bin edge, exclude the upper.
            // Leave missing or out-of-range electrons in the outside bin.
            int cell = outside;
            if (std::isfinite(momentum) && std::isfinite(angle) && momentum >= p.front() && momentum < p.back() &&
                    angle >= theta.front() && angle < theta.back()) {
                const int ip = std::upper_bound(p.begin(), p.end(), momentum)-p.begin()-1;
                const int it = std::upper_bound(theta.begin(), theta.end(), angle)-theta.begin()-1;
                cell = ip*settings.angleBins+it;
            }
            auto& m = cells[cell];
            // Subtract a common asymmetry reference to avoid large offsets. objectives() accounts
            // for that shift, so the variation used for the asymmetry allocation stays the same.
            if (count == 0 && settings.includeAsymmetry) referenceAsymmetry = event->A;
            ++count;
            for (size_t j = 0; j < formulas.size(); ++j) {
                require(formulas[j]->GetNdata() == 1, "response must give one value per event; use Sum$ for hits");
                const double response = formulas[j]->EvalInstance();
                require(std::isfinite(response) && response >= 0, "response must be finite and nonnegative");
                // Undo the division by total events applied in GeneratePrimaries() before saving rate.
                // Otherwise a larger pilot file would appear to have smaller event contributions.
                // Multiply by the detector response here; multiply by asymmetry as well for its numerator.
                contribution[j] = rate*run->GetNthrown()*response;
            }
            // fill() keeps the centered sums used by objectives() above.
            m.fill(contribution, settings.includeAsymmetry ? event->A-referenceAsymmetry : 0);
        }
    }
    // Give well-populated bins their own groups; pool the rest so a few pilot events
    // do not have to determine a separate quota for every rare bin.
    std::vector<int> separate;
    for (int cell = 0; cell <= outside; ++cell)
        if (cells[cell].count >= settings.minimumPilotCount && double(cells[cell].count)/count >= settings.minimumFraction)
            separate.push_back(cell);
    std::sort(separate.begin(), separate.end(), [&](int a, int b) { return cells[a].count > cells[b].count; });
    // Limit separate groups to leave room for the pooled group and for slots beyond the minimums.
    const size_t limit = std::max(1, budget/4-1);
    if (separate.size() > limit) separate.resize(limit);
    long long pooled = count;
    for (int cell : separate) pooled -= cells[cell].count;
    // Empty bins still need support. If the pool cannot estimate its own variation,
    // merge the smallest observed group into it rather than giving unseen regions zero probability.
    if (pooled < 2 && !separate.empty()) separate.pop_back();
    // Map every bin to a group, then add its counts and response totals to that group.
    std::vector<int> groupOfCell(cells.size(), separate.size());
    for (size_t h = 0; h < separate.size(); ++h) groupOfCell[separate[h]] = h;
    std::vector<Moment> groups(separate.size()+1, Moment(responses.size()));
    for (size_t cell = 0; cell < cells.size(); ++cell) groups[groupOfCell[cell]].add(cells[cell]);
    // Calculate the pilot fractions and variation, then turn those into integer transport quotas.
    const auto data = objectives(groups, settings.includeAsymmetry);
    const auto quota = balanced_allocation(data, budget);
    // Write the format read by Load() in remollNeymanPlan.cc, including generator settings
    // and pilot counts. Production uses these to check the plan and correct event rates.
    require(proposal.find_first_of("\r\n") == std::string::npos &&
            provenance.str().find_first_of("\r\n") == std::string::npos, "plan metadata must fit on one line");
    std::ofstream plan(output);
    require(bool(plan), "cannot write plan");
    plan << std::setprecision(17) << "remoll-neyman-quota-v1\nproposal " << std::quoted(proposal)
         << "\nprobabilities pilot " << count << ' ' << std::quoted(provenance.str())
         << "\nip\titheta\tp_lo_mev\tp_hi_mev\ttheta_lo_mrad\ttheta_hi_mrad\tgroup\tstratum_probability\ttransport_quota\tpilot_count\n";
    // Write a row for every bin, including outside. Bins sharing a group repeat its fraction
    // and quota; remoll counts that quota once for the group, not once per row.
    for (int cell = 0; cell <= outside; ++cell) {
        const int group = groupOfCell[cell];
        if (cell == outside) plan << "-1\t-1\tnan\tnan\tnan\tnan";
        else {
            const int ip = cell/settings.angleBins, it = cell%settings.angleBins;
            plan << ip << '\t' << it << '\t' << p[ip] << '\t' << p[ip+1] << '\t' << theta[it] << '\t' << theta[it+1];
        }
        plan << '\t' << group << '\t' << data.probability[group] << '\t' << quota[group] << '\t' << groups[group].count << '\n';
    }
    plan.close();
    require(bool(plan), "failed to finish writing plan");
    // Divide each quota by its pilot fraction to estimate the generator attempts it needs.
    // The largest sets the scale for beamOn; allow extra attempts for random fluctuations.
    double expected = 0;
    for (size_t h = 0; h < groups.size(); ++h) expected = std::max(expected, quota[h]/data.probability[h]);
    std::cout << "Wrote " << output << ": " << count << " pilot events, " << groups.size() << " groups, "
              << budget << " events to track.\nLargest quota / pilot fraction: " << std::ceil(expected)
              << " candidates. Use a larger beamOn count; this is not a guarantee of filling quotas.\n";
}
#endif
