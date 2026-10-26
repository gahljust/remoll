#ifndef REMOLL_NEYMAN_PLAN_HH
#define REMOLL_NEYMAN_PLAN_HH

#include <atomic>
#include <mutex>
#include <string>
#include <vector>

// Use a supplied plan; analyzing pilot runs and choosing quotas happens outside remoll.
// Workers share the same counts, so each quota applies to the whole run.
class remollNeymanPlan {
public:
    static remollNeymanPlan& Instance();
    void Load(const std::string& path);
    void Enable(bool enabled);
    bool Enabled() const { return fEnabled.load(); }
    void BeginRun(long long candidates);
    void RegisterProposal(const std::string& proposal);
    // Decide whether to track this event and return its rate correction.
    bool Accept(double momentumMeV, double thetaMrad, double& correction);
    void EndRun(bool complete);
    std::string Text() const;
    std::string Proposal() const;
private:
    mutable std::mutex fMutex;
    std::atomic<bool> fEnabled{false};
    bool fLoaded = false, fRunning = false;
    // fP/fTheta set the momentum/angle bins (MeV/mrad); fProbability gives each group's expected fraction.
    std::vector<double> fP, fTheta, fProbability;
    std::vector<int> fCell; // Maps each momentum/angle bin to its quota group; angle bins vary fastest.
    std::vector<long long> fQuota, fAccepted; // Events requested and accepted in each group.
    std::string fText, fExpectedProposal, fCurrentProposal;
    int fOutside = -1; // Group for missing electrons, nonfinite coordinates or points outside the grid.
    long long fCandidates = 0, fProcessed = 0; // Requested and processed candidates, including rejections.
};
#endif
