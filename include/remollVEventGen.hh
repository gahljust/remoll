#ifndef __REMOLLVEVENTGEN_HH
#define __REMOLLVEVENTGEN_HH

#include "remolltypes.hh"
#include "remollglobs.hh"
#include "remollVertex.hh"

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4GenericMessenger.hh"

/*!
   Generic base class for event generators
   This provides an interface for everyone to
   derive from.

   Ultimately this returns a remollEvent which is
   what the PrimaryGeneratorAction is going to use and
   contains information that will go in the output.

   It needs to be aware of remollBeamTarget and remollRunData,
   take a generically generated event assuming ideal beam
   and transform it into what is going to be simulated.
*/

class G4ParticleGun;

class remollEvent;
class remollBeamTarget;
class remollRunData;

class remollVEventGen {
    public:
	remollVEventGen(G4String name);
	virtual ~remollVEventGen();

	virtual void PrintEventGen();

	remollEvent *GenerateEvent();

	G4String GetName() { return fName; }

	void SetBeamTarget(remollBeamTarget* bt) {
	  fBeamTarg = bt;
	}

	void SetSamplingType(SamplingType_t type) { fSamplingType = type; }
	SamplingType_t GetSamplingType() const { return fSamplingType; }

	void SetDoMultScatt( G4bool multscatt ){ fApplyMultScatt = multscatt; }

        void SetEmin(double emin) { fE_min = emin; }
        void SetEmax(double emax) { fE_max = emax; }
        void SetPhmin(double phmin) { fPh_min = phmin; }
        void SetPhmax(double phmax) { fPh_max = phmax; }
        void SetThmin(double thmin) { fTh_min = thmin; }
        void SetThmax(double thmax) { fTh_max = thmax; }

    protected:
        G4double SampleThCoMWithBias(G4double& bias_weight);
        G4double SampleThetaWithBias(G4double min_th, G4double max_th, G4double& bias_weight);
        G4double SamplePhiWithBias(G4double& bias_weight);
        G4double SampleOutgoingEnergyWithBias(G4double maximum, G4double& bias_weight);
        G4bool HasThetaBias() const;

	// Generator name
        G4String fName;

    protected:
	// Generation limits
	static G4double fThCoM_min, fThCoM_max;
	static G4double fTh_min, fTh_max;
	static G4double fPh_min, fPh_max;
	static G4double fE_min, fE_max;
        static G4String fThCoMBiasMode;
        static G4double fThCoMBiasMin;
        static G4double fThCoMBiasMax;
        static G4double fThCoMBiasPhysicalFraction;
        G4double fLastThetaPhysicalPdf;
        G4double fLastThetaSamplePdf;
        static G4String fPhiBiasMode;
        static G4double fPhiBiasMin;
        static G4double fPhiBiasMax;
        static G4double fPhiBiasPhysicalFraction;
        G4double fLastPhiPhysicalPdf;
        G4double fLastPhiSamplePdf;
        static G4String fOutgoingEnergyBiasMode;
        static G4double fOutgoingEnergyBiasMinFraction;
        static G4double fOutgoingEnergyBiasMaxFraction;
        static G4double fOutgoingEnergyBiasPhysicalFraction;
        G4double fLastOutgoingEnergyPhysicalPdf;
        G4double fLastOutgoingEnergySamplePdf;

  G4String fBeamPol;
public:
  const G4String GetBeamPolarization(){return fBeamPol;}

    protected:
	// Number of particles
	G4int fNumberOfParticles;
	// Particle gun
        G4ParticleGun* fParticleGun;
    public:
	// Set the number of particles
	void SetNumberOfParticles(G4int n);
	// Get a new particle gun for this generator
        G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    protected:
	remollBeamTarget* fBeamTarg;

	void PolishEvent(remollEvent *);

	// Pure virtual function that needs to be filled out
	virtual void SamplePhysics(remollVertex *, remollEvent *) = 0;

    protected:

	SamplingType_t fSamplingType;
	G4bool     fApplyMultScatt;

    private:
	// Generic messenger as protected to be used in derived classes
	G4GenericMessenger fEvGenMessenger{
            this,
            "/remoll/evgen/",
            "Remoll event generator properties"};
	G4GenericMessenger fThCoMBiasMessenger{
            this,
            "/remoll/bias/thcom/",
            "Thrown-angle biasing"};
	G4GenericMessenger fPhiBiasMessenger{
            this,
            "/remoll/bias/phi/",
            "Azimuthal-angle biasing"};
	G4GenericMessenger fOutgoingEnergyBiasMessenger{
            this,
            "/remoll/bias/outgoinge/",
            "Outgoing-energy biasing"};
    protected:
	G4GenericMessenger fThisGenMessenger;
};


#endif//__REMOLLVEVENTGEN_HH
