#include "remollVEventGen.hh"

#include <cassert>
#include <algorithm>
#include <cmath>

#include "G4ParticleGun.hh"
#include "G4RotationMatrix.hh"
#include "Randomize.hh"

#include "remollBeamTarget.hh"
#include "remollVertex.hh"
#include "remollEvent.hh"
#include "remollRun.hh"
#include "remollRunData.hh"

G4double remollVEventGen::fTh_min = 0.0*deg;
G4double remollVEventGen::fTh_max = 180.0*deg;
G4double remollVEventGen::fThCoM_min = 0.0*deg;
G4double remollVEventGen::fThCoM_max = 180.0*deg;
G4double remollVEventGen::fPh_min = 0.0*deg;
G4double remollVEventGen::fPh_max = 360.0*deg;
G4double remollVEventGen::fE_min = 0.0*deg;
G4double remollVEventGen::fE_max = 11.0*GeV;

remollVEventGen::remollVEventGen(const G4String name)
: fName(name),
  fThCoMBiasMode("physical"),
  fThCoMBiasMin(-1.0),
  fThCoMBiasMax(-1.0),
  fBeamPol("0"),
  fNumberOfParticles(1),fParticleGun(0),
  fBeamTarg(0),
  fThisGenMessenger(this,"/remoll/evgen/" + name + "/","Remoll " + name + " generator properties")
{
    // Set initial number of particles and create particle gun
    SetNumberOfParticles(fNumberOfParticles);

    // Create event generator messenger
    fEvGenMessenger.DeclarePropertyWithUnit("emax","GeV",fE_max,"Maximum generation energy");
    fEvGenMessenger.DeclarePropertyWithUnit("emin","GeV",fE_min,"Minimum generation energy");
    fEvGenMessenger.DeclarePropertyWithUnit("thmax","deg",fTh_max,"Maximum generation theta angle");
    fEvGenMessenger.DeclarePropertyWithUnit("thmin","deg",fTh_min,"Minimum generation theta angle");
    fEvGenMessenger.DeclarePropertyWithUnit("phmax","deg",fPh_max,"Maximum generation phi angle");
    fEvGenMessenger.DeclarePropertyWithUnit("phmin","deg",fPh_min,"Minimum generation phi angle");
    fEvGenMessenger.DeclarePropertyWithUnit("thcommax","deg",fThCoM_max,"Maximum CoM generation theta angle");
    fEvGenMessenger.DeclarePropertyWithUnit("thcommin","deg",fThCoM_min,"Minimum CoM generation theta angle");
    fEvGenMessenger.DeclareProperty("beamPolarization",fBeamPol,"Polarization direction: +L, +H, +V, -L, -H, -V, 0");
    fEvGenMessenger.DeclareMethod(
        "printlimits",
        &remollVEventGen::PrintEventGen,
        "Print the event generator limits");

    fThCoMBiasMessenger.DeclareProperty("mode",fThCoMBiasMode,"Thrown-angle bias mode: physical or uniform");
    fThCoMBiasMessenger.DeclarePropertyWithUnit("min","deg",fThCoMBiasMin,"Minimum biased thrown angle");
    fThCoMBiasMessenger.DeclarePropertyWithUnit("max","deg",fThCoMBiasMax,"Maximum biased thrown angle");

    fSamplingType = kActiveTargetVolume;
    fApplyMultScatt = false;
}

remollVEventGen::~remollVEventGen()
{
}

G4double remollVEventGen::SampleThCoMWithBias(G4double& bias_weight)
{
    return SampleThetaWithBias(fThCoM_min, fThCoM_max, bias_weight);
}

G4bool remollVEventGen::HasThetaBias() const
{
    return !(fThCoMBiasMode == "physical" || fThCoMBiasMode == "none");
}

G4double remollVEventGen::SampleThetaWithBias(G4double min_limit, G4double max_limit, G4double& bias_weight)
{
    bias_weight = 1.0;

    if (fThCoMBiasMode == "physical" || fThCoMBiasMode == "none") {
        return acos(G4RandFlat::shoot(cos(max_limit), cos(min_limit)));
    }

    if (fThCoMBiasMode == "uniform") {
        G4double min_th = fThCoMBiasMin >= 0.0 ? std::max(fThCoMBiasMin, min_limit) : min_limit;
        G4double max_th = fThCoMBiasMax >= 0.0 ? std::min(fThCoMBiasMax, max_limit) : max_limit;

        if (max_th <= min_th) {
            G4cerr << "ERROR: invalid /remoll/bias/thcom uniform range: "
                   << min_th/deg << " to " << max_th/deg << " deg" << G4endl;
            exit(1);
        }

        G4double th = G4RandFlat::shoot(min_th, max_th);
        G4double physical_pdf = sin(th)/(cos(min_limit) - cos(max_limit));
        G4double sample_pdf = 1.0/(max_th - min_th);
        bias_weight = physical_pdf/sample_pdf;
        return th;
    }

    G4cerr << "ERROR: unknown /remoll/bias/thcom/mode " << fThCoMBiasMode << G4endl;
    exit(1);
}

void remollVEventGen::PrintEventGen()
{
  G4cout << "Event generator: " << fName << G4endl;
  G4cout << "E =     [" << fE_min/GeV  << "," << fE_max/GeV  << "] GeV" << G4endl;
  G4cout << "phi =   [" << fPh_min/deg << "," << fPh_max/deg << "] deg" << G4endl;
  G4cout << "theta = [" << fTh_min/deg << "," << fTh_max/deg << "] deg" << G4endl;
  G4cout << "theta (COM) = [" << fThCoM_min/deg << "," << fThCoM_max/deg << "] deg" << G4endl;
}

void remollVEventGen::SetNumberOfParticles(G4int n)
{
  // Store new number of particles
  fNumberOfParticles = n;

  // Delete old particle gun
  if (fParticleGun != nullptr) {
    delete fParticleGun;
    fParticleGun = 0;
  }
  // Create new particle gun
  fParticleGun = new G4ParticleGun(fNumberOfParticles);
}

remollEvent* remollVEventGen::GenerateEvent()
{
    // Set up beam/target vertex
    remollVertex vert   = fBeamTarg->SampleVertex(fSamplingType);

    /////////////////////////////////////////////////////////////////////
    // Create and initialize values for event
    remollEvent *thisev = new remollEvent();
    thisev->SetBeamTarget(fBeamTarg);
    thisev->fBeamBiasWeight = fBeamTarg->fBeamBiasWeight;
    thisev->fBeamBiasPhysicalPdf = fBeamTarg->fBeamBiasPhysicalPdf;
    thisev->fBeamBiasSamplePdf = fBeamTarg->fBeamBiasSamplePdf;
    thisev->fBiasWeight = thisev->fBeamBiasWeight;

    thisev->fVertexPos    = fBeamTarg->fVer;
    if( fApplyMultScatt ) {
        thisev->fBeamMomentum = fBeamTarg->fSampledEnergy*(fBeamTarg->fDir.unit());
    } else {
        thisev->fBeamMomentum = fBeamTarg->fSampledEnergy*G4ThreeVector(0.0, 0.0, 1.0);
    }
    /////////////////////////////////////////////////////////////////////

    SamplePhysics(&vert, thisev);

    PolishEvent(thisev);

    return thisev;
}


void remollVEventGen::PolishEvent(remollEvent *ev) {
    /*!
       Here it's our job to:
          Make sure the event is sane
          Apply multiple scattering effects to the final
        products if applicable
      Calculate rates from our given luminosity
      Calculate measured asymmetry from polarization
      Calculate vertex offsets
     */

    if( !ev->EventIsSane() ) {
        G4cerr << __FILE__ << " line " << __LINE__ << ":  Event check failed for generator " << fName << ".  Aborting" << G4endl;
        ev->Print();
        exit(1);
    }

    G4ThreeVector rotax      = (-1)*(fBeamTarg->fDir.cross(G4ThreeVector(0.0, 0.0, 1.0))).unit();
    G4RotationMatrix msrot;
    msrot.rotate(fBeamTarg->fDir.theta(), rotax);

    std::vector<G4ThreeVector>::iterator iter;

    if( fApplyMultScatt ) {
        for( iter = ev->fPartRealMom.begin(); iter != ev->fPartRealMom.end(); iter++ ) {
            //  rotate direction vectors based on multiple scattering
            (*iter) *= msrot;
        }

        // Rotate position offsets due to multiple scattering
        for( iter = ev->fPartPos.begin(); iter != ev->fPartPos.end(); iter++ ) {
            //  rotate direction vectors based on multiple scattering
            (*iter) *= msrot;
        }
    }

    // Add base vertex
    for( iter = ev->fPartPos.begin(); iter != ev->fPartPos.end(); iter++ ) {
        (*iter) += ev->fVertexPos;
    }

    ev->fmAsym = ev->fAsym*fBeamTarg->fBeamPolarization;
}
