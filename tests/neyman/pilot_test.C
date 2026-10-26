// build/reroot -l -b -q 'tests/neyman/pilot_test.C("/tmp/neyman-pilot-test")'
// Use a fresh directory. Creates small ROOT fixtures and tests the real C++ plan reader.
#include "../../analysis/neyman/make_plan.C"
#include "../../include/remollNeymanPlan.hh"
#include <TSystem.h>
#include <functional>

namespace pilot_test_detail {
void check(bool ok) { if (!ok) throw std::runtime_error("pilot test failed"); }
void rejects(const std::function<void()>& action) {
    try { action(); } catch (const std::exception&) { return; }
    throw std::runtime_error("expected rejection");
}
void fixture(const std::string& path, int repeats = 1, bool sampled = false,
             bool incomplete = false, const char* proposal = "test proposal") {
    TFile file(path.c_str(), "NEW");
    check(!file.IsZombie());
    TTree tree("T", "ordinary pilot");
    remollRunData run;
    const int total = 200*repeats;
    run.SetNthrown(total + (incomplete ? 1 : 0));
    double rate = 0;
    remollEvent_t event{};
    std::vector<remollEventParticle_t> particles;
    std::vector<remollGenericDetectorHit_t> hits;
    tree.Branch("rate", &rate);
    tree.Branch("ev", &event);
    tree.Branch("part", &particles);
    tree.Branch("hit", &hits);
    for (int i = 0; i < total; ++i) {
        const bool high = i%200 >= 100;
        const double p = high ? 7500 : 2500;
        remollEventParticle_t electron{};
        electron.pid = 11; electron.p = p;
        electron.px = p*std::sin(.01); electron.pz = p*std::cos(.01);
        particles = {electron};
        hits.clear();
        if (i%2) {
            remollGenericDetectorHit_t hit{};
            hit.det = 28; hit.pid = 11; hit.pz = 1;
            hits.push_back(hit);
        }
        rate = (high ? 6.0 : 2.0)/total;
        event.A = i%3;
        tree.Fill();
    }
    tree.Write(); run.Write("run_data");
    TObjString(proposal).Write("generator_proposal");
    if (sampled) TObjString("already sampled").Write("neyman_plan");
}
}
void pilot_test(const char* directory) {
    using namespace pilot_test_detail;
    try {
        // Constant rates within each group, but asymmetry residual spreads of 2 and 4.
        // Only the ratio objective should drive the spare quotas, in a 1:2 split.
        neyman_pilot::Moment lowMoment(1), highMoment(1);
        lowMoment.fill({2}, 1); lowMoment.fill({2}, 3);
        highMoment.fill({4}, 1); highMoment.fill({4}, 3);
        // A small spread around a large rate must survive averaging and merging bins.
        neyman_pilot::Moment large(1), other(1);
        large.fill({1e12-1}, 0); other.fill({1e12+1}, 0); large.add(other);
        check(neyman_pilot::objectives({large}, false).variance[0] == 1);
        const auto data = neyman_pilot::objectives({lowMoment, highMoment}, true);
        check(data.variance == std::vector<double>({1, 10}));
        check(data.within[0] == std::vector<double>({0, 4}) && data.within[1] == std::vector<double>({0, 16}));
        check(neyman_pilot::balanced_allocation(data, 100) == std::vector<int>({34, 66}));
        check(gSystem->mkdir(directory, true) == 0);
        const std::string base = std::string(directory)+"/";
        fixture(base+"a.root"); fixture(base+"b.root", 2);
        fixture(base+"sampled.root", 1, true);
        fixture(base+"incomplete.root", 1, false, true);
        fixture(base+"changed.root", 1, false, false, "different generator");
        neyman_pilot::Settings settings;
        settings.momentumBins = 2; settings.angleBins = 1;
        settings.maxMomentumMeV = 10000; settings.minAngleMrad = 1; settings.maxAngleMrad = 20;
        settings.includeAsymmetry = false;
        const std::vector<std::string> response = {"Sum$(hit.det==28)>0"};
        make_plan({base+"a.root", base+"b.root"}, (base+"plan.tsv").c_str(), 100, response, settings);
        auto& plan = remollNeymanPlan::Instance();
        plan.Load(base+"plan.tsv"); plan.Enable(true); plan.BeginRun(200); plan.RegisterProposal("test proposal");
        int low = 0, high = 0;
        // Equal region fractions, with three times the spread in the high region:
        // after reserving two events each, split the other 96 as 24 / 72.
        for (int i = 0; i < 200; ++i) {
            double correction = 0;
            const bool upper = i >= 100;
            if (!plan.Accept(upper ? 7500 : 2500, 10, correction)) continue;
            (upper ? high : low)++;
            check(std::abs(correction-100.0/(upper ? 74 : 26)) < 1e-12);
        }
        check(low == 26 && high == 74); plan.EndRun(true); plan.Enable(false);
        for (const auto& inputs : std::vector<std::vector<std::string>>{
                {base+"a.root", base+"a.root"}, {base+"sampled.root"}, {base+"incomplete.root"},
                {base+"a.root", base+"changed.root"}, {base+"missing.root"}})
            rejects([&] { make_plan(inputs, (base+"bad.tsv").c_str(), 100, response, settings); });
        rejects([&] { make_plan({base+"a.root"}, (base+"plan.tsv").c_str(), 100, response, settings); });
        rejects([&] { make_plan({base+"a.root"}, (base+"bad.tsv").c_str(), 100, {"hit.det"}, settings); });
        rejects([&] { make_plan({base+"a.root"}, (base+"bad.tsv").c_str(), 100, {"-1"}, settings); });
        rejects([&] { make_plan({base+"a.root"}, (base+"bad.tsv").c_str(), 100, {"0"}, settings); });
        settings.minimumPilotCount = 1000; // All bins pool together, including unseen/outside bins.
        make_plan({base+"a.root"}, (base+"pooled.tsv").c_str(), 20, response, settings);
        plan.Load(base+"pooled.tsv"); plan.Enable(true); plan.BeginRun(20); plan.RegisterProposal("test proposal");
        for (int i = 0; i < 20; ++i) { double correction; check(plan.Accept(NAN, NAN, correction) && correction == 1); }
        plan.EndRun(true); plan.Enable(false);
        // Extreme finite inputs must either make a readable plan or fail before writing one.
        settings.maxMomentumMeV = 1e-320;
        rejects([&] { make_plan({base+"a.root"}, (base+"bad.tsv").c_str(), 20, response, settings); });
        settings.maxMomentumMeV = 1e308;
        make_plan({base+"a.root"}, (base+"large-grid.tsv").c_str(), 20, response, settings);
        plan.Load(base+"large-grid.tsv");
        settings.minAngleMrad = 1e-320;
        rejects([&] { make_plan({base+"a.root"}, (base+"bad.tsv").c_str(), 20, response, settings); });
        std::cout << "PASS rate/asymmetry allocation, pilot moments, unequal file sizes, exact quotas/corrections, outside support and invalid inputs\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; gSystem->Exit(1);
    }
}
