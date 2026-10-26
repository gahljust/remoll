// Build with the same ROOT/Geant4 flags as remoll and link libremoll.
// Material names matter: remollGenpElastic selects G4_lH2 by name for radiation loss.
#include "remollBeamTarget.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
int main() {
    auto a = new G4Material("G4_lH2", 1., 1.008*g/mole, .07*g/cm3);
    auto b = new G4Material("other_hydrogen", 1., 1.008*g/mole, .14*g/cm3);
    auto world = new G4LogicalVolume(new G4Box("world", 100*mm, 100*mm, 100*mm), a, "world");
    auto mother = new G4PVPlacement(nullptr, {}, world, "mother", nullptr, false, 0);
    remollBeamTarget beam;
    beam.AddMotherVolume(mother, "LH2");
    for (auto material : {a, b}) {
        auto volume = new G4LogicalVolume(new G4Box("target", 10*mm, 10*mm, 10*mm), material, "target");
        auto placement = new G4PVPlacement(nullptr, {}, volume, "target", world, false, 0);
        beam.AddTargetVolume(placement, material == a ? "upstream" : "LH2");
    }
    const auto before = beam.SamplingProposal();
    a->SetName("other_hydrogen"); b->SetName("G4_lH2");
    if (before == beam.SamplingProposal()) {
        std::cerr << "Changed elastic radiation-loss material was not detected\n";
        return 1;
    }
    std::cout << "PASS proposal includes material names\n";
}
