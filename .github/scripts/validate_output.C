// Check the saved file itself: remoll can exit successfully after a macro error.
// Compile with reroot validate_output.C+ and the image's own remoll headers.
#include "remollRunData.hh"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
std::string required(const char* key) {
  const char* value = std::getenv(key);
  if (!value || !*value) throw std::runtime_error(std::string("Missing ") + key);
  return value;
}
unsigned long long positive(const char* key) {
  const auto text = required(key);
  unsigned long long value = 0;
  for (unsigned char c : text) {
    if (c < '0' || c > '9' ||
        value > (std::numeric_limits<unsigned long long>::max() - (c-'0')) / 10)
      throw std::runtime_error(std::string("Invalid positive integer: ") + key);
    value = value * 10 + (c-'0');
  }
  if (!value) throw std::runtime_error(std::string("Zero is invalid: ") + key);
  return value;
}
}

void validate_output() {
  try {
    const auto path = required("EXPECTED_ROOT");
    const auto seed = positive("EXPECTED_SEED");
    const auto events = positive("EXPECTED_EVENTS");
    const auto filtered = required("EXPECTED_FILTERED");
    if (filtered != "0" && filtered != "1")
      throw std::runtime_error("EXPECTED_FILTERED must be 0 or 1");
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
      throw std::runtime_error("ROOT file is missing, invalid, zombie, or recovered");
    auto* tree = file->Get<TTree>("T");
    auto* run = file->Get<remollRunData>("run_data");
    if (!tree || !run) throw std::runtime_error("Missing T or remollRunData run_data");
    if (run->GetSeed() != seed) throw std::runtime_error("run_data seed mismatch");
    if (run->GetNthrown() != events) throw std::runtime_error("run_data Nthrown mismatch");
    const auto entries = tree->GetEntries();
    if (entries <= 0 || static_cast<unsigned long long>(entries) > events)
      throw std::runtime_error("T must be nonempty and cannot exceed Nthrown");
    // Rejection sampling writes fewer entries; old images may not save a plan.
    // Check ordinary counts exactly, without claiming to verify old quota plans.
    if (filtered == "0" && static_cast<unsigned long long>(entries) != events)
      throw std::runtime_error("Unfiltered T entries do not match Nthrown");
    tree->SetBranchStatus("*", 1);
    for (Long64_t i = 0; i < entries; ++i)
      if (tree->GetEntry(i) <= 0) throw std::runtime_error("Unreadable T entry " + std::to_string(i));
    std::cout << "PASS: seed=" << seed << " Nthrown=" << events << " entries=" << entries;
    if (filtered == "1") std::cout << "; filtered output: quotas NOT verified";
    std::cout << std::endl;
    file->Close();
    gSystem->Exit(0);
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
}
