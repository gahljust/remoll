#ifndef remollPrimaryGeneratorAction_h
#define remollPrimaryGeneratorAction_h 1

#include "remollBeamTarget.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4GenericMessenger.hh"
#include "G4String.hh"

#include <map>
#include <memory>
#include <vector>

class G4ParticleGun;
class G4Event;
class remollVEventGen;
class remollEvent;

class remollPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    remollPrimaryGeneratorAction();
    virtual ~remollPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

    const remollEvent* GetEvent() const { return fEvent; }
    // Fixed-quota sampling generates many cheap candidates after a stratum is
    // already full.  Those rejected candidates must count in Nthrown but need
    // no full ROOT record.  Accepted zero-weight trials remain writable.
    G4bool ShouldWriteEvent() const { return fWriteEvent; }

    void SetGenerator(G4String&);
    void LoadTransportGate(G4String);
    void LoadNeymanTransportPlan(G4String);

  private:
    std::map<G4String,std::shared_ptr<remollVEventGen>> fEvGenMap;
    std::shared_ptr<remollVEventGen> fEventGen;
    G4String fEventGenName;

    std::map<G4String,std::shared_ptr<G4VPrimaryGenerator>> fPriGenMap;
    std::shared_ptr<G4VPrimaryGenerator> fPriGen;
    G4String fPriGenName;

    G4ParticleGun* fParticleGun;

    remollBeamTarget fBeamTarg;


    remollEvent *fEvent;

    G4int fRateCopy;
    G4bool fGeneratorOnly;
    G4bool fTransportGateEnabled;
    G4bool fNeymanTransportEnabled;
    G4bool fWriteEvent;
    G4String fTransportGateFile;
    G4int fTransportGatePBins;
    G4int fTransportGateThetaBins;
    std::vector<G4double> fTransportGatePEdges;
    std::vector<G4double> fTransportGateThetaEdges;
    std::vector<G4double> fTransportGateKeep;
    G4GenericMessenger fEvGenMessenger{this,"/remoll/evgen/","Remoll event generator properties"};

    G4double TransportGateProbability() const;
    G4bool ApplyNeymanTransportPlan(G4double& correction) const;

    G4double fEffCrossSection;
};

#endif
