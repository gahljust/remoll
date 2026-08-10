#include "remollGenFixedInelasticC12.hh"

#include "remollEvent.hh"
#include "remollVertex.hh"

#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4ElectroVDNuclearModel.hh"
#include "G4Element.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4HadSecondary.hh"
#include "G4HadronicProcess.hh"
#include "G4Material.hh"
#include "G4Nucleus.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysListUtil.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace {

constexpr G4int kCarbonZ = 6;
constexpr G4int kCarbonA = 12;

const G4Element* GetCarbonElement(const G4Material* material) {
  if (material == nullptr || material->GetNumberOfElements() != 1) return nullptr;

  const G4Element* element = material->GetElement(0);
  if (element == nullptr || element->GetZasInt() != kCarbonZ) return nullptr;
  return element;
}

G4double MomentumFromKineticEnergy(const G4double kineticEnergy,
                                   const G4double mass) {
  return std::sqrt(std::max(0.0, kineticEnergy * (kineticEnergy + 2.0 * mass)));
}

void AddPrimary(remollEvent* event, const G4ParticleDefinition* definition,
                const G4ThreeVector& momentum) {
  G4PrimaryParticle particle(definition, momentum.x(), momentum.y(), momentum.z());
  event->ProduceNewParticle(G4ThreeVector(), &particle);
}

} // namespace

remollGenFixedInelasticC12::remollGenFixedInelasticC12()
    : remollVEventGen("fixedinelasticC12"),
      fElectronNuclearProcess(nullptr),
      fElectronNuclearModel(nullptr) {
  // Geant4 returns all final-state directions in the incident-electron frame.
  // PolishEvent applies the sampled beam/multiple-scattering rotation afterward.
  fApplyMultScatt = true;
}

remollGenFixedInelasticC12::~remollGenFixedInelasticC12() = default;

void remollGenFixedInelasticC12::InitializeG4Model() {
  if (fElectronNuclearProcess != nullptr && fElectronNuclearModel != nullptr) return;

  // This method is deliberately lazy.  The generator objects are constructed
  // before /run/initialize, while the reference physics list creates and
  // registers electronNuclear during /run/initialize.
  fElectronNuclearProcess =
      G4PhysListUtil::FindInelasticProcess(G4Electron::Electron());
  if (fElectronNuclearProcess == nullptr
      || fElectronNuclearProcess->GetProcessName() != "electronNuclear") {
    G4cerr << "fixedinelasticC12 requires the Geant4 electronNuclear process. "
           << "Select a reference physics list that registers G4EmExtraPhysics "
           << "before /run/initialize." << G4endl;
    std::exit(1);
  }

  fElectronNuclearModel = dynamic_cast<G4ElectroVDNuclearModel*>(
      fElectronNuclearProcess->GetHadronicModel("G4ElectroVDNuclearModel"));
  if (fElectronNuclearModel == nullptr) {
    G4cerr << "fixedinelasticC12 found electronNuclear, but it does not use "
           << "G4ElectroVDNuclearModel." << G4endl;
    std::exit(1);
  }

  G4cout << "fixedinelasticC12: borrowing "
         << fElectronNuclearProcess->GetProcessName() << " / "
         << fElectronNuclearModel->GetModelName()
         << " from the active Geant4 physics list; target nucleus is C12"
         << G4endl;
}

void remollGenFixedInelasticC12::SamplePhysics(remollVertex* vertex,
                                                remollEvent* event) {
  InitializeG4Model();

  const G4Material* material = vertex->GetMaterial();
  const G4Element* carbon = GetCarbonElement(material);
  if (carbon == nullptr) {
    G4cerr << "fixedinelasticC12 requires an elemental carbon target (Z=6); "
           << "active material is "
           << (material == nullptr ? "<null>" : material->GetName())
           << G4endl;
    std::exit(1);
  }

  const G4double beamKineticEnergy = vertex->GetBeamEnergy();
  G4DynamicParticle incoming(G4Electron::Electron(), G4ThreeVector(0.0, 0.0, 1.0),
                             beamKineticEnergy);

  // Microscopic cross section per carbon atom from the exact data store used
  // by the transported electronNuclear process.
  const G4double totalCrossSection =
      fElectronNuclearProcess->GetElementCrossSection(&incoming, carbon, material);
  if (!(totalCrossSection >= 0.0) || !std::isfinite(totalCrossSection)) {
    G4cerr << "fixedinelasticC12 received an invalid Geant4 cross section at "
           << beamKineticEnergy / GeV << " GeV" << G4endl;
    std::exit(1);
  }

  G4HadProjectile projectile(incoming);
  G4Nucleus target(kCarbonA, kCarbonZ);
  G4HadFinalState* finalState =
      fElectronNuclearModel->ApplyYourself(projectile, target);
  if (finalState == nullptr) {
    G4cerr << "fixedinelasticC12: G4ElectroVDNuclearModel returned no final state"
           << G4endl;
    std::exit(1);
  }

  const G4double finalElectronKineticEnergy = finalState->GetEnergyChange();
  const G4ThreeVector finalElectronDirection =
      finalState->GetMomentumChange().unit();
  const G4bool acceptedInteraction =
      finalElectronKineticEnergy >= 0.0
      && finalElectronKineticEnergy < beamKineticEnergy;

  // CalculateEMVertex contains a rejection correction between the electro- and
  // photo-nuclear parameterizations.  Preserve that rejection as a zero-weight
  // trial instead of resampling; conditioning on successes would bias the rate.
  event->SetEffCrossSection(acceptedInteraction ? totalCrossSection : 0.0);
  event->SetAsymmetry(0.0);

  if (acceptedInteraction) {
    const G4double electronMass = G4Electron::Electron()->GetPDGMass();
    const G4double initialTotalEnergy = beamKineticEnergy + electronMass;
    const G4double finalTotalEnergy =
        finalElectronKineticEnergy + electronMass;
    const G4double initialMomentum =
        MomentumFromKineticEnergy(beamKineticEnergy, electronMass);
    const G4double finalMomentum =
        MomentumFromKineticEnergy(finalElectronKineticEnergy, electronMass);
    const G4double cosTheta = std::clamp(finalElectronDirection.z(), -1.0, 1.0);
    // Evaluate Q2 without subtracting two O(E^2) terms.  The electronuclear
    // spectrum is sharply forward and reaches Q2 values for which the direct
    // expression 2(E E' - p p' cos(theta) - m^2) loses most of its digits.
    // The two terms below are algebraically identical but remain positive and
    // well conditioned in the collinear and small-energy-loss limits.
    const G4double energyDifference =
        initialTotalEnergy - finalTotalEnergy;
    const G4double massTerm =
        electronMass * electronMass * energyDifference * energyDifference
        / (initialTotalEnergy * finalTotalEnergy
           + initialMomentum * finalMomentum
           - electronMass * electronMass);
    const G4double sinThetaSquared =
        finalElectronDirection.x() * finalElectronDirection.x()
        + finalElectronDirection.y() * finalElectronDirection.y();
    const G4double angularTerm = cosTheta > -1.0
        ? initialMomentum * finalMomentum * sinThetaSquared / (1.0 + cosTheta)
        : 2.0 * initialMomentum * finalMomentum;
    const G4double q2 = std::max(0.0, 2.0 * (massTerm + angularTerm));
    const G4double scatteringAngle =
        std::atan2(std::sqrt(sinThetaSquared), cosTheta);
    const G4double energyTransfer =
        beamKineticEnergy - finalElectronKineticEnergy;
    const G4double w2 = proton_mass_c2 * proton_mass_c2
        + 2.0 * proton_mass_c2 * energyTransfer - q2;

    event->SetQ2(q2);
    event->SetW2(w2);
    event->SetThCoM(scatteringAngle);
    if (energyTransfer > 0.0) {
      event->SetXbj(q2 / (2.0 * proton_mass_c2 * energyTransfer));
    }

    // G4HadronicProcess::FillResult keeps the continuing projectile as the
    // first track.  Do the same so remoll's polarization handling remains
    // attached to the scattered electron.
    AddPrimary(event, G4Electron::Electron(),
               finalMomentum * finalElectronDirection);

  } else {
    // A zero-cross-section, unchanged electron keeps the Geant4 event valid
    // without adding entries to weighted detector distributions.
    AddPrimary(event, G4Electron::Electron(), incoming.GetMomentum());
    event->SetQ2(0.0);
    event->SetW2(0.0);
    event->SetThCoM(0.0);
  }

  for (std::size_t index = 0;
       index < finalState->GetNumberOfSecondaries(); ++index) {
    G4DynamicParticle* secondary =
        finalState->GetSecondary(index)->GetParticle();
    if (secondary == nullptr) continue;

    if (acceptedInteraction && secondary->GetDefinition() != nullptr) {
      // Match G4HadronicProcess::FillResult's correction for off-shell model
      // products before turning them into new primary tracks.
      const G4double nominalMass = secondary->GetDefinition()->GetPDGMass();
      const G4double dynamicMass = secondary->GetMass();
      if (std::abs(dynamicMass - nominalMass) > 1.0 * keV) {
        secondary->SetKineticEnergy(std::max(
            secondary->GetKineticEnergy() + dynamicMass - nominalMass,
            0.001 * eV));
        secondary->SetMass(nominalMass);
      }

      AddPrimary(event, secondary->GetDefinition(), secondary->GetMomentum());
    }

    // The model transfers ownership of dynamic secondaries to its caller.
    delete secondary;
  }
  finalState->ClearSecondaries();
}
