#ifndef __REMOLLINTERACTIONRECORDER_HH
#define __REMOLLINTERACTIONRECORDER_HH

#include "G4GenericMessenger.hh"
#include "G4String.hh"

#include <cstdint>
#include <fstream>
#include <mutex>

class G4Step;

// Optional step-level training-data recorder for material response maps.
// Each interaction is written as one continuation row plus one row for every
// secondary created in that step.  It deliberately records continuous states,
// not a pre-binned transfer matrix or a database of whole showers.
class remollInteractionRecorder {
  public:
    static remollInteractionRecorder& GetInstance();

    bool IsEnabled() const { return fEnabled; }
    void BeginRun();
    void RecordStep(const G4Step* step);
    void EndRun();

  private:
    remollInteractionRecorder();
    remollInteractionRecorder(const remollInteractionRecorder&) = delete;
    remollInteractionRecorder& operator=(
        const remollInteractionRecorder&) = delete;

    bool fEnabled{false};
    bool fRunActive{false};
    bool fIncludeTransportation{false};
    G4String fOutputPath{"interactions.tsv"};
    G4String fMaterial{""};
    G4String fPhysicalVolume{""};
    G4double fMinKineticEnergyMeV{0.0};
    std::ofstream fOutput;
    std::uint64_t fNextInteraction{0U};
    std::uint64_t fRecordedInteractions{0U};
    std::uint64_t fRecordedSecondaries{0U};
    mutable std::mutex fMutex;

    G4GenericMessenger fMessenger{
      this,
      "/remoll/interaction/",
      "Material interaction training-data recording properties"};
};

#endif // __REMOLLINTERACTIONRECORDER_HH
