#include "remollNeymanPlan.hh"
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <thread>

void check(bool condition) {
    if (!condition) throw std::runtime_error("test failed");
}
void fails(const std::function<void()>& action) {
    try {
        action();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error("expected failure");
}
int main(int argc, char** argv) {
    try {
        check(argc == 2);
        const std::string path = argv[1];
        const std::string columns = "ip\titheta\tp_lo_mev\tp_hi_mev\ttheta_lo_mrad\ttheta_hi_mrad\tgroup\tstratum_probability\ttransport_quota\tpilot_count\n";
        const std::string header = "remoll-neyman-quota-v1\nproposal \"test\"\nprobabilities pilot 100 \"independent pilot fixture\"\n";
        const std::string rows = "0 0 0 10 0 20 0 0.9 50 90\n-1 -1 nan nan nan nan 1 0.1 50 10\n";
        auto write = [&](const std::string& text) {
            std::ofstream file(path);
            file << text;
        };
        auto& plan = remollNeymanPlan::Instance();
        double correction;
        check(plan.Accept(NAN, NAN, correction) && correction == 1);
        fails([&] { plan.Enable(true); });
        write(header + columns + rows);
        plan.Load(path);
        plan.Enable(true);
        fails([&] { plan.BeginRun(99); });
        for (int run = 0; run < 2; ++run) {
            plan.BeginRun(400);
            fails([&] { plan.Accept(1, 1, correction); });
            plan.RegisterProposal("test");
            fails([&] { plan.Load(path); });
            fails([&] { plan.Enable(false); });
            std::atomic<int> inside{0}, outside{0};
            std::vector<std::thread> workers;
            for (int i = 0; i < 4; ++i) {
                workers.emplace_back([&] {
                    for (int j = 0; j < 100; ++j) {
                        double factor;
                        const bool inner = j%2 == 0;
                        if (!plan.Accept(inner ? 5 : 10, 1, factor)) continue;
                        check(std::abs(factor-(inner ? 7.2 : 0.8)) < 1e-12);
                        if (inner) ++inside; else ++outside;
                    }
                });
            }
            for (auto& worker : workers) worker.join();
            check(inside == 50 && outside == 50);
            plan.EndRun(true);
        }
        plan.BeginRun(400);
        fails([&] { plan.RegisterProposal("changed target"); });
        fails([&] { plan.EndRun(false); });
        plan.BeginRun(400);
        plan.RegisterProposal("test");
        for (int i = 0; i < 400; ++i) plan.Accept(1, 1, correction);
        fails([&] { plan.EndRun(true); });
        plan.Enable(false);
        for (const auto& text : {"remoll-bernoulli-v1\n", "ip\titheta\n"}) {
            write(text);
            fails([&] { plan.Load(path); });
        }
        for (const auto& badRows : {
                "0 0 0 10 0 20 0 .9 50 90\n",
                "0 0 0 10 0 20 0 .9 50 89\n-1 -1 nan nan nan nan 1 .1 50 10\n",
                "0 0 0 10 0 20 0 nan 50 90\n-1 -1 nan nan nan nan 1 .1 50 10\n",
                "0 0 0 10 0 20 0 .9 1 90\n-1 -1 nan nan nan nan 1 .1 50 10\n"}) {
            write(header + columns + badRows);
            fails([&] { plan.Load(path); });
        }
        // One documented cell-count limit, no separate 10,000-axis restriction.
        std::string large = "remoll-neyman-quota-v1\nproposal \"test\"\nprobabilities known 0 \"single complete stratum\"\n" + columns;
        for (int i = 0; i < 10001; ++i)
            large += std::to_string(i)+" 0 "+std::to_string(i)+" "+std::to_string(i+1)+" 0 1 0 1 2 0\n";
        write(large + "-1 -1 nan nan nan nan 0 1 2 0\n");
        plan.Load(path);
        plan.Enable(true);
        plan.BeginRun(2);
        plan.RegisterProposal("test");
        check(plan.Accept(10000.5, .5, correction) && correction == 1);
        check(plan.Accept(NAN, NAN, correction) && correction == 1);
        plan.EndRun(true);
        std::cout << "PASS exact quotas/weights, proposal checks, concurrency, reset, underfill, provenance/parser checks\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
