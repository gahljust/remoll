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
G4String remollVEventGen::fThCoMBiasMode = "physical";
G4double remollVEventGen::fThCoMBiasMin = -1.0;
G4double remollVEventGen::fThCoMBiasMax = -1.0;
G4double remollVEventGen::fThCoMBiasPhysicalFraction = 1.0;
G4String remollVEventGen::fPhiBiasMode = "physical";
G4double remollVEventGen::fPhiBiasMin = -1.0;
G4double remollVEventGen::fPhiBiasMax = -1.0;
G4double remollVEventGen::fPhiBiasPhysicalFraction = 1.0;
G4String remollVEventGen::fOutgoingEnergyBiasMode = "physical";
G4double remollVEventGen::fOutgoingEnergyBiasMinFraction = 0.0;
G4double remollVEventGen::fOutgoingEnergyBiasMaxFraction = 1.0;
G4double remollVEventGen::fOutgoingEnergyBiasPhysicalFraction = 1.0;

remollVEventGen::remollVEventGen(const G4String name)
: fName(name),
  fLastThetaPhysicalPdf(1.0),
  fLastThetaSamplePdf(1.0),
  fLastPhiPhysicalPdf(1.0),
  fLastPhiSamplePdf(1.0),
  fLastOutgoingEnergyPhysicalPdf(1.0),
  fLastOutgoingEnergySamplePdf(1.0),
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

    fThCoMBiasMessenger.DeclareProperty("mode",fThCoMBiasMode,"Thrown-angle bias mode: physical, uniform, or mixture");
    fThCoMBiasMessenger.DeclarePropertyWithUnit("min","deg",fThCoMBiasMin,"Minimum biased thrown angle");
    fThCoMBiasMessenger.DeclarePropertyWithUnit("max","deg",fThCoMBiasMax,"Maximum biased thrown angle");
    fThCoMBiasMessenger.DeclareProperty("physicalFraction",fThCoMBiasPhysicalFraction,"Physical component probability in mixture mode");
    fPhiBiasMessenger.DeclareProperty("mode",fPhiBiasMode,"Azimuthal bias mode: physical, uniform, or mixture");
    fPhiBiasMessenger.DeclarePropertyWithUnit("min","deg",fPhiBiasMin,"Minimum biased azimuthal angle");
    fPhiBiasMessenger.DeclarePropertyWithUnit("max","deg",fPhiBiasMax,"Maximum biased azimuthal angle");
    fPhiBiasMessenger.DeclareProperty("physicalFraction",fPhiBiasPhysicalFraction,"Physical component probability in mixture mode");
    fOutgoingEnergyBiasMessenger.DeclareProperty("mode",fOutgoingEnergyBiasMode,"Outgoing-energy bias mode: physical or mixture");
    fOutgoingEnergyBiasMessenger.DeclareProperty("minFraction",fOutgoingEnergyBiasMinFraction,"Minimum fraction of the kinematic maximum");
    fOutgoingEnergyBiasMessenger.DeclareProperty("maxFraction",fOutgoingEnergyBiasMaxFraction,"Maximum fraction of the kinematic maximum");
    fOutgoingEnergyBiasMessenger.DeclareProperty("physicalFraction",fOutgoingEnergyBiasPhysicalFraction,"Physical component probability in mixture mode");

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

G4double remollVEventGen::ThrowThetaProposal(G4double min_limit, G4double max_limit,
                                             ThetaProposal_t proposal) const
{
    if (proposal == kInverseOneMinusCosSquared) {
        const G4double icth_a = 1.0/(1.0 - cos(min_limit));
        const G4double icth_b = 1.0/(1.0 - cos(max_limit));
        const G4double sampv = 1.0/G4RandFlat::shoot(icth_b, icth_a);
        return acos(1.0 - sampv);
    }
    return acos(G4RandFlat::shoot(cos(max_limit), cos(min_limit)));
}

G4double remollVEventGen::ThetaProposalPdf(G4double th, G4double min_limit, G4double max_limit,
                                           ThetaProposal_t proposal) const
{
    const G4double cthmin = cos(min_limit);
    const G4double cthmax = cos(max_limit);
    if (proposal == kInverseOneMinusCosSquared) {
        const G4double icth_a = 1.0/(1.0 - cthmin);
        const G4double icth_b = 1.0/(1.0 - cthmax);
        const G4double u = 1.0 - cos(th);
        if (u <= 0.0) return 0.0;
        // u = 1/X with X uniform on [icth_b,icth_a] gives q(u) = 1/((icth_a-icth_b) u^2);
        // multiplying by |du/dth| = sin(th) puts it in the same measure as the target.
        return sin(th)/((icth_a - icth_b)*u*u);
    }
    return sin(th)/(cthmin - cthmax);
}

G4double remollVEventGen::SampleThetaWithBias(G4double min_limit, G4double max_limit,
                                              G4double& bias_weight, ThetaProposal_t proposal)
{
    bias_weight = 1.0;
    const G4double normalization = cos(min_limit) - cos(max_limit);

    if (fThCoMBiasMode == "physical" || fThCoMBiasMode == "none") {
        G4double th = ThrowThetaProposal(min_limit, max_limit, proposal);
        fLastThetaPhysicalPdf = sin(th)/normalization;
        fLastThetaSamplePdf = ThetaProposalPdf(th, min_limit, max_limit, proposal);
        // Unity for the flat-in-cos proposal; the 1/(1-cos)^2 sampling factor
        // otherwise, which is what the generators used to fold in by hand.
        bias_weight = fLastThetaSamplePdf > 0.0
            ? fLastThetaPhysicalPdf/fLastThetaSamplePdf : 0.0;
        return th;
    }

    if (fThCoMBiasMode == "uniform" || fThCoMBiasMode == "mixture") {
        G4double min_th = fThCoMBiasMin >= 0.0 ? std::max(fThCoMBiasMin, min_limit) : min_limit;
        G4double max_th = fThCoMBiasMax >= 0.0 ? std::min(fThCoMBiasMax, max_limit) : max_limit;

        if (max_th <= min_th) {
            G4cerr << "ERROR: invalid /remoll/bias/thcom uniform range: "
                   << min_th/deg << " to " << max_th/deg << " deg" << G4endl;
            exit(1);
        }

        G4double th;
        if (fThCoMBiasMode == "mixture") {
            if (fThCoMBiasPhysicalFraction <= 0.0 || fThCoMBiasPhysicalFraction > 1.0) {
                G4cerr << "ERROR: /remoll/bias/thcom/physicalFraction must be in (0,1]" << G4endl;
                exit(1);
            }
            if (G4UniformRand() < fThCoMBiasPhysicalFraction)
                th = ThrowThetaProposal(min_limit, max_limit, proposal);
            else
                th = G4RandFlat::shoot(min_th, max_th);
        } else {
            th = G4RandFlat::shoot(min_th, max_th);
        }
        fLastThetaPhysicalPdf = sin(th)/normalization;
        G4double uniform_pdf = (th >= min_th && th <= max_th) ? 1.0/(max_th - min_th) : 0.0;
        // The physical component of the mixture is the *proposal* density, not
        // the target measure.  Using the target measure here would silently
        // discard the 1/(1-cos)^2 importance sampling.
        G4double physical_pdf = ThetaProposalPdf(th, min_limit, max_limit, proposal);
        fLastThetaSamplePdf = fThCoMBiasMode == "mixture"
            ? fThCoMBiasPhysicalFraction*physical_pdf
              + (1.0-fThCoMBiasPhysicalFraction)*uniform_pdf
            : uniform_pdf;
        bias_weight = fLastThetaSamplePdf > 0.0
            ? fLastThetaPhysicalPdf/fLastThetaSamplePdf : 0.0;
        return th;
    }

    G4cerr << "ERROR: unknown /remoll/bias/thcom/mode " << fThCoMBiasMode << G4endl;
    exit(1);
}

G4double remollVEventGen::SamplePhiWithBias(G4double& bias_weight)
{
    bias_weight = 1.0;
    const G4double physical_pdf = 1.0/(fPh_max - fPh_min);
    if (fPhiBiasMode == "physical" || fPhiBiasMode == "none") {
        fLastPhiPhysicalPdf = physical_pdf;
        fLastPhiSamplePdf = physical_pdf;
        return G4RandFlat::shoot(fPh_min, fPh_max);
    }
    G4double min_ph = fPhiBiasMin >= 0.0 ? std::max(fPhiBiasMin, fPh_min) : fPh_min;
    G4double max_ph = fPhiBiasMax >= 0.0 ? std::min(fPhiBiasMax, fPh_max) : fPh_max;
    if (max_ph <= min_ph) {
        G4cerr << "ERROR: invalid /remoll/bias/phi range" << G4endl;
        exit(1);
    }
    G4double ph;
    if (fPhiBiasMode == "mixture") {
        if (fPhiBiasPhysicalFraction <= 0.0 || fPhiBiasPhysicalFraction > 1.0) {
            G4cerr << "ERROR: /remoll/bias/phi/physicalFraction must be in (0,1]" << G4endl;
            exit(1);
        }
        ph = G4UniformRand() < fPhiBiasPhysicalFraction
            ? G4RandFlat::shoot(fPh_min, fPh_max)
            : G4RandFlat::shoot(min_ph, max_ph);
    } else if (fPhiBiasMode == "uniform") {
        ph = G4RandFlat::shoot(min_ph, max_ph);
    } else {
        G4cerr << "ERROR: unknown /remoll/bias/phi/mode " << fPhiBiasMode << G4endl;
        exit(1);
    }
    G4double uniform_pdf = (ph >= min_ph && ph <= max_ph) ? 1.0/(max_ph-min_ph) : 0.0;
    fLastPhiPhysicalPdf = physical_pdf;
    fLastPhiSamplePdf = fPhiBiasMode == "mixture"
        ? fPhiBiasPhysicalFraction*physical_pdf + (1.0-fPhiBiasPhysicalFraction)*uniform_pdf
        : uniform_pdf;
    bias_weight = fLastPhiPhysicalPdf/fLastPhiSamplePdf;
    return ph;
}

G4double remollVEventGen::SampleOutgoingEnergyWithBias(G4double maximum,
                                                       G4double& bias_weight)
{
    bias_weight = 1.0;
    const G4double physical_pdf = 1.0/maximum;
    if (fOutgoingEnergyBiasMode == "physical" || fOutgoingEnergyBiasMode == "none") {
        fLastOutgoingEnergyPhysicalPdf = physical_pdf;
        fLastOutgoingEnergySamplePdf = physical_pdf;
        return G4RandFlat::shoot(0.0, maximum);
    }
    if (fOutgoingEnergyBiasMode != "mixture") {
        G4cerr << "ERROR: unknown /remoll/bias/outgoinge/mode "
               << fOutgoingEnergyBiasMode << G4endl;
        exit(1);
    }
    const G4double min_fraction = std::max(0.0, fOutgoingEnergyBiasMinFraction);
    const G4double max_fraction = std::min(1.0, fOutgoingEnergyBiasMaxFraction);
    if (max_fraction <= min_fraction || fOutgoingEnergyBiasPhysicalFraction <= 0.0
            || fOutgoingEnergyBiasPhysicalFraction > 1.0) {
        G4cerr << "ERROR: invalid /remoll/bias/outgoinge mixture settings" << G4endl;
        exit(1);
    }
    const G4double minimum = min_fraction*maximum;
    const G4double target_maximum = max_fraction*maximum;
    const G4double energy = G4UniformRand() < fOutgoingEnergyBiasPhysicalFraction
        ? G4RandFlat::shoot(0.0, maximum)
        : G4RandFlat::shoot(minimum, target_maximum);
    const G4double target_pdf = energy >= minimum && energy <= target_maximum
        ? 1.0/(target_maximum-minimum) : 0.0;
    fLastOutgoingEnergyPhysicalPdf = physical_pdf;
    fLastOutgoingEnergySamplePdf = fOutgoingEnergyBiasPhysicalFraction*physical_pdf
        + (1.0-fOutgoingEnergyBiasPhysicalFraction)*target_pdf;
    bias_weight = fLastOutgoingEnergyPhysicalPdf/fLastOutgoingEnergySamplePdf;
    return energy;
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
    thisev->fVertexBiasWeight = fBeamTarg->fVertexBiasWeight;
    thisev->fVertexBiasPhysicalPdf = fBeamTarg->fVertexBiasPhysicalPdf;
    thisev->fVertexBiasSamplePdf = fBeamTarg->fVertexBiasSamplePdf;
    thisev->fBiasWeight = thisev->fBeamBiasWeight*thisev->fVertexBiasWeight;

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
