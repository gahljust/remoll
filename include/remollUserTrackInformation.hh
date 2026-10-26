#ifndef remollUserTrackInformation_h_
#define remollUserTrackInformation_h_

#include "G4StepStatus.hh"
#include "G4VUserTrackInformation.hh"

#include <string>

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
    const std::string& GetCreatorPhysicalVolumeName() const { return fCreatorPhysicalVolumeName; };
    void SetCreatorPhysicalVolumeName(const std::string& name) { fCreatorPhysicalVolumeName = name; };
    const std::string& GetCreatorProcessName() const { return fCreatorProcessName; };
    void SetCreatorProcessName(const std::string& name) { fCreatorProcessName = name; };
    const std::string& GetCreatorMaterialName() const { return fCreatorMaterialName; };
    void SetCreatorMaterialName(const std::string& name) { fCreatorMaterialName = name; };
  private:
    G4StepStatus fStepStatus;
    G4int fCreatorPhysicalVolumeID{-1};
    G4int fCreatorProcessID{-1};
    G4int fCreatorMaterialID{-1};
    std::string fCreatorPhysicalVolumeName;
    std::string fCreatorProcessName;
    std::string fCreatorMaterialName;
};

#endif /* remollUserTrackInformation_h_ */
