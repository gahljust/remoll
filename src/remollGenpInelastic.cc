#include "remollGenpInelastic.hh"

#include "christy_bosted_inelastic.h"

#include "Randomize.hh"

#include "remollEvent.hh"
#include "remollVertex.hh"
#include "G4Material.hh"

#include "remolltypes.hh"

remollGenpInelastic::remollGenpInelastic()
: remollVEventGen("inelastic") {
    fTh_min =     0.1*deg;
    fTh_max =     5.0*deg;

    fApplyMultScatt = true;
}

remollGenpInelastic::~remollGenpInelastic(){
}

void remollGenpInelastic::SamplePhysics(remollVertex *vert, remollEvent *evt){
    // Generate inelastic event

    double beamE = vert->GetBeamEnergy();
    double mp    = proton_mass_c2;

    // 1/(1-cos)^2 proposal: the inelastic cross section still carries the
    // forward Mott-like rise, so flat-in-cos throwing wastes most of the
    // sample on unreachable small angles.
    G4double th_bias_weight = 1.0;
    double th = SampleThetaWithBias(fTh_min, fTh_max, th_bias_weight,
                                    kInverseOneMinusCosSquared);
    evt->fThCoMBiasWeight = th_bias_weight;
    evt->fThCoMBiasPhysicalPdf = fLastThetaPhysicalPdf;
    evt->fThCoMBiasSamplePdf = fLastThetaSamplePdf;
    evt->fBiasWeight *= th_bias_weight;
    G4double ph_bias_weight = 1.0;
    double ph = SamplePhiWithBias(ph_bias_weight);
    evt->fPhiBiasWeight = ph_bias_weight;
    evt->fPhiBiasPhysicalPdf = fLastPhiPhysicalPdf;
    evt->fPhiBiasSamplePdf = fLastPhiSamplePdf;
    evt->fBiasWeight *= ph_bias_weight;
    double efmax = mp*beamE/(mp + beamE*(1.0-cos(th)));;
    G4double energy_bias_weight = 1.0;
    double ef = SampleOutgoingEnergyWithBias(efmax, energy_bias_weight);
    evt->fOutgoingEnergyBiasWeight = energy_bias_weight;
    evt->fOutgoingEnergyBiasPhysicalPdf = fLastOutgoingEnergyPhysicalPdf;
    evt->fOutgoingEnergyBiasSamplePdf = fLastOutgoingEnergySamplePdf;
    evt->fBiasWeight *= energy_bias_weight;

    double thissigma_p = sigma_p( beamE/GeV, th, ef/GeV )*nanobarn/GeV;
    double thissigma_n = sigma_n( beamE/GeV, th, ef/GeV )*nanobarn/GeV;

    double sigmatot = thissigma_p*vert->GetMaterial()->GetZ() +
	//  Effective neutron number...  I don't like it either  SPR 2/14/2013
	thissigma_n*(vert->GetMaterial()->GetA()*mole/g - vert->GetMaterial()->GetZ());

    double V = (fPh_max - fPh_min) * (cos(fTh_min) - cos(fTh_max)) * efmax;

    evt->SetEffCrossSection(sigmatot*V);

    if( vert->GetMaterial()->GetNumberOfElements() != 1 ){
	G4cerr << __FILE__ << " line " << __LINE__ << 
	    ": Error!  Some lazy programmer didn't account for complex materials in the moller process!" << G4endl;
	exit(1);
    }


    double Q2 = 2.0*beamE*ef*(1.0-cos(th));
    evt->SetQ2( Q2 );
    evt->SetThCoM(th);

    G4double APV = -1.0*Q2*0.8e-4/GeV/GeV; // Empirical APV value, 
                                      // stolen from mollerClass.C in mollersim
	// R-L asymmetry for ep inelastic should be negative

    evt->SetAsymmetry(APV);


    evt->SetW2( mp*mp + 2.0*mp*(beamE-ef) - Q2 );

    evt->ProduceNewParticle( G4ThreeVector(0.0, 0.0, 0.0), 
	                     G4ThreeVector(ef*sin(th)*cos(ph), ef*sin(th)*sin(ph), ef*cos(th) ), 
			     "e-" );

    return;
}
