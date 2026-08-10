#ifndef __REMOLLGENFIXEDINELASTICC12_HH
#define __REMOLLGENFIXEDINELASTICC12_HH

#include "remollVEventGen.hh"

class G4ElectroVDNuclearModel;
class G4HadronicProcess;

// Force one Geant4 electron-nuclear interaction on C12 at the remoll target
// vertex.  The process cross section and final-state model are borrowed from
// the active reference physics list after /run/initialize, so this generator
// follows the same model stack as normal Geant4 transport.
class remollGenFixedInelasticC12 : public remollVEventGen {
public:
  remollGenFixedInelasticC12();
  virtual ~remollGenFixedInelasticC12();

private:
  void InitializeG4Model();
  void SamplePhysics(remollVertex*, remollEvent*) override;

  G4HadronicProcess* fElectronNuclearProcess;
  G4ElectroVDNuclearModel* fElectronNuclearModel;
};

#endif // __REMOLLGENFIXEDINELASTICC12_HH
