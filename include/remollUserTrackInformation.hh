#ifndef remollUserTrackInformation_h_
#define remollUserTrackInformation_h_

#include "G4StepStatus.hh"
#include "G4VUserTrackInformation.hh"

class remollUserTrackInformation : public G4VUserTrackInformation
{
  public:
    remollUserTrackInformation() { fStepStatus = fUndefined; };
    virtual ~remollUserTrackInformation() { };
    G4StepStatus GetStepStatus() const { return fStepStatus; };
    void SetStepStatus(G4StepStatus stepstatus) { fStepStatus = stepstatus; };

    G4int GetCreatorPhysicalVolumeID() const { return fCreatorPhysicalVolumeID; };
    void SetCreatorPhysicalVolumeID(G4int id) { fCreatorPhysicalVolumeID = id; };
    G4int GetCreatorProcessID() const { return fCreatorProcessID; };
    void SetCreatorProcessID(G4int id) { fCreatorProcessID = id; };
    G4int GetCreatorMaterialID() const { return fCreatorMaterialID; };
    void SetCreatorMaterialID(G4int id) { fCreatorMaterialID = id; };
  private:
    G4StepStatus fStepStatus;
    G4int fCreatorPhysicalVolumeID{-1};
    G4int fCreatorProcessID{-1};
    G4int fCreatorMaterialID{-1};
};

#endif /* remollUserTrackInformation_h_ */
