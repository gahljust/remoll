// From the repository root, after building remoll:
//   build/reroot -l -b -q 'analysis/neyman/example.C("pilot.root","plan.tsv",1000,28)'
// First make an ordinary pilot with your usual geometry, fields and generator settings:
//   /remoll/evgen/transportNeymanEnable false
//   /remoll/seed 12345
//   /remoll/filename pilot.root
//   /run/beamOn 100000
// Keep every event, including detector misses. Each input must be a complete run
// from this branch. Multiple files need independent seeds and the same settings.
// For production, keep those settings and the detector response unchanged, then use:
//   /remoll/seed 67890
//   /remoll/filename sampled.root
//   /remoll/evgen/transportNeymanFile plan.tsv
//   /remoll/evgen/transportNeymanEnable true
//   /run/beamOn <candidate count larger than the estimate printed by make_plan>
// The quota budget is events to track, not beamOn's number of generator candidates.
// Underfilled runs fail. Increase beamOn and rerun; do not use incomplete output.
// Saved rates already include the correction. Pilot-probability uncertainty and
// an independent check of precision are still needed for a final physics result.
#include "make_plan.C"
void example(const char* pilot = "pilot.root", const char* plan = "plan.tsv",
             int budget = 1000, int detector = 28) {
    neyman_pilot::Settings settings;
    // This example optimizes the rate of events with a forward electron hit, plus
    // their rate-weighted asymmetry. It does not include our study's PE response model.
    // Use your detector ID and response: a hit count, a PE sum, or a 0/1 event flag.
    // Do not multiply the expression by rate; make_plan does that for every event.
    const std::string response = "Sum$(hit.det=="+std::to_string(detector)+" && hit.pid==11 && hit.pz>0)>0";
    // Set false for rate alone; pass several expressions to balance several detectors.
    settings.includeAsymmetry = true;
    make_plan({pilot}, plan, budget, {response}, settings);
}
