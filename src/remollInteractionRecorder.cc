#include "remollInteractionRecorder.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

#include "remollSystemOfUnits.hh"

#include <iomanip>
#include <string>
#include <vector>

namespace {

std::string Clean(const G4String& value)
{
  std::string result = value;
  for (char& character : result) {
    if (character == '\t' || character == '\n' || character == '\r') {
      character = ' ';
    }
  }
  return result;
}

} // namespace

remollInteractionRecorder& remollInteractionRecorder::GetInstance()
{
  static remollInteractionRecorder instance;
  return instance;
}

remollInteractionRecorder::remollInteractionRecorder()
{
  fMessenger.DeclareProperty(
      "enable", fEnabled, "Enable material interaction recording");
  fMessenger.DeclareProperty(
      "output", fOutputPath, "Interaction table output file (.tsv)");
  fMessenger.DeclareProperty(
      "material", fMaterial,
      "Only record this Geant4 material name; empty records all materials");
  fMessenger.DeclareProperty(
      "includeTransportation", fIncludeTransportation,
      "Also record Transportation and step-limiter steps");
  fMessenger.DeclareProperty(
      "minKineticEnergyMeV", fMinKineticEnergyMeV,
      "Do not record parent steps below this pre-step kinetic energy [MeV]");
}

void remollInteractionRecorder::BeginRun()
{
  std::lock_guard<std::mutex> lock(fMutex);
  if (!fEnabled || fRunActive) return;
  fOutput.open(fOutputPath.data(), std::ios::out | std::ios::trunc);
  if (!fOutput) {
    G4cerr << "remoll/interaction: cannot open " << fOutputPath
           << "; disabling" << G4endl;
    fEnabled = false;
    return;
  }
  fOutput << "event\tinteraction\ttrack\tparent\trole\tsecondary_index"
             "\tchild_track\tpdg"
             "\tmaterial\tprocess\tprocess_type\tprocess_subtype\tstatus"
             "\tpre_x_mm\tpre_y_mm\tpre_z_mm\tpre_px_mev\tpre_py_mev"
             "\tpre_pz_mev\tpre_ke_mev\tpost_x_mm\tpost_y_mm\tpost_z_mm"
             "\tpost_px_mev\tpost_py_mev\tpost_pz_mev\tpost_ke_mev"
             "\ttime_ns\tstep_mm\tedep_mev\tweight\n";
  fOutput << std::setprecision(17);
  fNextInteraction = 0U;
  fRecordedInteractions = 0U;
  fRecordedSecondaries = 0U;
  fRunActive = true;
  G4cout << "remoll/interaction: recording to " << fOutputPath;
  if (!fMaterial.empty()) G4cout << " for material " << fMaterial;
  G4cout << G4endl;
}

void remollInteractionRecorder::RecordStep(const G4Step* step)
{
  if (!fEnabled || !fRunActive || step == nullptr) return;
  const G4StepPoint* pre = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  const G4Track* parent = step->GetTrack();
  if (pre == nullptr || post == nullptr || parent == nullptr) return;
  const G4Material* material = pre->GetMaterial();
  if (material == nullptr) return;
  if (!fMaterial.empty() && material->GetName() != fMaterial) return;
  if (pre->GetKineticEnergy() / MeV < fMinKineticEnergyMeV) return;

  const G4VProcess* process = post->GetProcessDefinedStep();
  const G4String process_name = process == nullptr
      ? G4String("none") : process->GetProcessName();
  if (!fIncludeTransportation
      && (process_name == "Transportation"
          || process_name == "CoupledTransportation"
          || process_name == "StepLimiter")) return;

  const G4Event* event =
      G4EventManager::GetEventManager()->GetConstCurrentEvent();
  const G4int event_id = event == nullptr ? -1 : event->GetEventID();
  const G4ThreeVector pre_x = pre->GetPosition();
  const G4ThreeVector pre_p = pre->GetMomentum();
  const G4ThreeVector post_x = post->GetPosition();
  const G4ThreeVector post_p = post->GetMomentum();
  const G4int process_type = process == nullptr ? -1 : process->GetProcessType();
  const G4int process_subtype =
      process == nullptr ? -1 : process->GetProcessSubType();
  const std::vector<const G4Track*>* secondaries =
      step->GetSecondaryInCurrentStep();

  std::lock_guard<std::mutex> lock(fMutex);
  if (!fRunActive) return;
  const std::uint64_t interaction = fNextInteraction++;
  auto write_prefix = [&](const char* role, std::size_t secondary_index,
                          G4int child_track, G4int pdg,
                          const G4ThreeVector& out_x,
                          const G4ThreeVector& out_p, G4double out_ke,
                          G4double out_time, G4double weight) {
    fOutput << event_id << '\t' << interaction << '\t'
            << parent->GetTrackID() << '\t' << parent->GetParentID() << '\t'
            << role << '\t' << secondary_index << '\t' << child_track
            << '\t' << pdg << '\t'
            << Clean(material->GetName()) << '\t' << Clean(process_name)
            << '\t' << process_type << '\t' << process_subtype << '\t'
            << parent->GetTrackStatus() << '\t'
            << pre_x.x() / mm << '\t' << pre_x.y() / mm << '\t'
            << pre_x.z() / mm << '\t' << pre_p.x() / MeV << '\t'
            << pre_p.y() / MeV << '\t' << pre_p.z() / MeV << '\t'
            << pre->GetKineticEnergy() / MeV << '\t' << out_x.x() / mm
            << '\t' << out_x.y() / mm << '\t' << out_x.z() / mm << '\t'
            << out_p.x() / MeV << '\t' << out_p.y() / MeV << '\t'
            << out_p.z() / MeV << '\t' << out_ke / MeV << '\t'
            << out_time / ns << '\t' << step->GetStepLength() / mm << '\t'
            << step->GetTotalEnergyDeposit() / MeV << '\t' << weight << '\n';
  };

  write_prefix("continuation", 0U, parent->GetTrackID(),
               parent->GetParticleDefinition()->GetPDGEncoding(), post_x,
               post_p, post->GetKineticEnergy(), post->GetGlobalTime(),
               parent->GetWeight());
  if (secondaries != nullptr) {
    std::size_t secondary_index = 0U;
    for (const G4Track* child : *secondaries) {
      if (child == nullptr) continue;
      write_prefix("secondary", ++secondary_index, child->GetTrackID(),
                   child->GetParticleDefinition()->GetPDGEncoding(),
                   child->GetPosition(), child->GetMomentum(),
                   child->GetKineticEnergy(), child->GetGlobalTime(),
                   child->GetWeight());
      ++fRecordedSecondaries;
    }
  }
  ++fRecordedInteractions;
}

void remollInteractionRecorder::EndRun()
{
  std::lock_guard<std::mutex> lock(fMutex);
  if (!fRunActive) return;
  fOutput.close();
  fRunActive = false;
  G4cout << "remoll/interaction: wrote " << fRecordedInteractions
         << " interactions and " << fRecordedSecondaries
         << " secondary rows to " << fOutputPath << G4endl;
}
