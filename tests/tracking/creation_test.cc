// Build with Geant4 flags and src/remollTrackingAction.cc.
// Resuming after optical secondaries must not change the stored birthplace.
#include "remollTrackingAction.hh"
#include "remollUserTrackInformation.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
int main() {
 auto particle=G4Geantino::Definition();
 G4ParticleTable::GetParticleTable()->SetReadiness(true);
 auto material=new G4Material("birthMaterial",1.,1.*g/mole,1.*g/cm3);
 auto other=new G4Material("laterMaterial",1.,1.*g/mole,2.*g/cm3);
 auto make=[&](const char* name,G4Material* material) {
   auto lv=new G4LogicalVolume(new G4Box(name,1*mm,1*mm,1*mm),material,name);
   return new G4PVPlacement(nullptr,{},lv,name,nullptr,false,0);
 };
 auto born=make("birthVolume",material), later=make("laterVolume",other);
 G4Track track(new G4DynamicParticle(particle,{0,0,1},1*MeV),0,{});
 auto locate=[&](G4VPhysicalVolume* pv){G4NavigationHistory h;h.SetFirstEntry(pv);track.SetTouchableHandle(G4TouchableHandle(new G4TouchableHistory(h)));};
 remollTrackingAction action;
 locate(born); action.PreUserTrackingAction(&track);
 auto info=dynamic_cast<remollUserTrackInformation*>(track.GetUserInformation());
 std::cout << "First: " << info->GetCreatorPhysicalVolumeName() << std::endl;
 track.IncrementCurrentStepNumber(); track.SetTrackStatus(fSuspend); locate(later);
 action.PreUserTrackingAction(&track);
 std::cout << "Resumed: " << info->GetCreatorPhysicalVolumeName() << " " << info->GetCreatorMaterialName() << std::endl;
 return info->GetCreatorPhysicalVolumeName()=="birthVolume" && info->GetCreatorMaterialName()=="birthMaterial" ? 0:1;
}
