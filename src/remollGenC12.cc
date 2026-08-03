#include "remollGenC12.hh"

#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "remollEvent.hh"
#include "remollVertex.hh"
#include "remolltypes.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

extern G4int F1F2IN09(G4int Z, G4int IA, G4double qsq,
		      G4double wsq, G4double &F1, G4double &F2);

extern G4int F1F2QE09(G4int Z, G4int IA, G4double qsq,
		      G4double wsq, G4double &F1, G4double &F2);

std::map<G4int,G4String> remollGenC12::fNamesMap = remollGenC12::CreateNamesMap();

std::map<G4int,G4String> remollGenC12::CreateNamesMap() {
  std::map<G4int,G4String> names;
  names[0] = "elasticC12";
  names[1] = "quasielasticC12";
  names[2] = "inelasticC12";
  return names;
}

remollGenC12::remollGenC12(G4int physicsType)
: remollVEventGen(fNamesMap[physicsType]) {
    type = physicsType;

    fTh_min =     0.1*deg;
    fTh_max =     5.0*deg;
    fApplyMultScatt = true;
}

remollGenC12::~remollGenC12() {
}

namespace {
struct InelasticEnergyDomain {
    G4double low;
    G4double high;
    G4double fullWidth;
};

InelasticEnergyDomain GetInelasticEnergyDomain(G4double beamE, G4double theta)
{
    // Resolve the interval that the old rejection loop sampled conditionally.
    // Sampling this interval directly is equivalent in physical mode and lets
    // the shared outgoing-energy mixture evaluate its proposal density exactly.
    const G4double q2_per_energy = 4.0*beamE*std::pow(std::sin(theta/2.0), 2);
    const G4double denominator = 2.0*proton_mass_c2 + q2_per_energy;
    const G4double low = std::max(
        electron_mass_c2,
        (proton_mass_c2*proton_mass_c2 + 2.0*proton_mass_c2*beamE
         - 10.0*GeV*GeV)/denominator);
    G4double high = std::min(
        beamE,
        (proton_mass_c2*proton_mass_c2 + 2.0*proton_mass_c2*beamE)
        /denominator);
    if (q2_per_energy > 0.0)
        high = std::min(high, 10.0*GeV*GeV/q2_per_energy);
    return {low, high, beamE-electron_mass_c2};
}
}

void remollGenC12::SamplePhysics(remollVertex *vert, remollEvent *evt) {

    G4double beamE = vert->GetBeamEnergy(); // in MeV (it can be modified by beam loss)
    // Sample with 1/(1-cos)^2, which very nearly cancels the Mott divergence
    // in the C12 elastic cross section.  Throwing flat in cos(th) here (as
    // this generator did) spreads the event weight over ~6 decades and is the
    // origin of the pathological rate tail seen in long carbon runs.
    G4double th_bias_weight = 1.0;
    G4double th = SampleThetaWithBias(fTh_min, fTh_max, th_bias_weight,
                                      kInverseOneMinusCosSquared); // radians
    evt->fThCoMBiasWeight = th_bias_weight;
    evt->fThCoMBiasPhysicalPdf = fLastThetaPhysicalPdf;
    evt->fThCoMBiasSamplePdf = fLastThetaSamplePdf;
    evt->fBiasWeight *= th_bias_weight;

    // Full solid angle.  The 1/(1-cos)^2 reweighting factor that used to be
    // folded in here by hand (the commented-out samp_fact) is now returned by
    // SampleThetaWithBias and applied to fEffXs through fBiasWeight, so it
    // must not be applied twice.
    G4double phaseSpaceFactor = (fPh_max - fPh_min) * (cos(fTh_min) - cos(fTh_max));
    G4double ph_bias_weight = 1.0;
    G4double ph = SamplePhiWithBias(ph_bias_weight);
    evt->fPhiBiasWeight = ph_bias_weight;
    evt->fPhiBiasPhysicalPdf = fLastPhiPhysicalPdf;
    evt->fPhiBiasSamplePdf = fLastPhiSamplePdf;
    evt->fBiasWeight *= ph_bias_weight;
    G4double eOut=0;
    G4double fWeight=0;
    G4double Q2=0;
    G4double effectiveXsection=0;
    G4double W2=0;
    G4double asym=0;

    switch (type) {
    case 0:
      GenElastic(beamE,th,Q2,W2,effectiveXsection,fWeight,eOut,asym);
      break;
    case 1:
      GenQuasiElastic(beamE,th,Q2,W2,effectiveXsection,fWeight,eOut,asym);
      phaseSpaceFactor=phaseSpaceFactor*(beamE/GeV - electron_mass_c2/GeV);
      break;
    case 2: {
      const InelasticEnergyDomain domain = GetInelasticEnergyDomain(beamE, th);
      const G4double accepted_width = domain.high-domain.low;
      if (domain.fullWidth <= 0.0 || accepted_width <= 0.0) {
        G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << G4endl;
        G4cerr << "  empty C12 inelastic outgoing-energy domain: beamE theta "
               << beamE/GeV << " " << th/deg << G4endl;
        exit(1);
      }

      // Sample the same accepted interval as the former rejection loop.  The
      // shared helper returns p_conditional/q_conditional; multiplying by the
      // accepted/full width ratio gives p_full/q_conditional and preserves the
      // full-range phase-space normalization below.
      G4double energy_bias_weight = 1.0;
      eOut = domain.low
          + SampleOutgoingEnergyWithBias(accepted_width, energy_bias_weight);
      energy_bias_weight *= accepted_width/domain.fullWidth;
      evt->fOutgoingEnergyBiasWeight = energy_bias_weight;
      evt->fOutgoingEnergyBiasPhysicalPdf = 1.0/domain.fullWidth;
      evt->fOutgoingEnergyBiasSamplePdf = fLastOutgoingEnergySamplePdf;
      evt->fBiasWeight *= energy_bias_weight;

      GenInelastic(beamE,th,eOut,Q2,W2,effectiveXsection,fWeight,asym);
      phaseSpaceFactor=phaseSpaceFactor*(beamE/GeV - electron_mass_c2/GeV);
      break;
    }
    default:
      G4cerr<<"Unknown C12 event type"<<G4endl;
      exit(1);
    }
    
    const G4Material* material = vert->GetMaterial();
    if (material == nullptr || material->GetNumberOfElements() != 1
        || std::abs(material->GetZ() - 6.0) > 1.0e-9) {
        G4cerr << __FILE__ << " line " << __LINE__ <<
               "  : Error!  The C12 generator requires an elemental carbon "
               "target (Z=6), but the active target material is "
               << (material == nullptr ? "<null>" : material->GetName())
               << ". Select a carbon target before /run/beamOn, for example "
               "/control/execute macros/target/Optics1.mac." << G4endl;
        exit(1);
    }

    evt->SetThCoM(th);
    evt->SetEffCrossSection(phaseSpaceFactor*effectiveXsection*microbarn);
    evt->SetQ2( Q2 ); //MeV^2
    evt->SetAsymmetry(asym);
    evt->SetW2( W2 );
    evt->ProduceNewParticle( G4ThreeVector(0.0, 0.0, 0.0),
                             G4ThreeVector(eOut*sin(th)*cos(ph), eOut*sin(th)*sin(ph), eOut*cos(th) ),
                             "e-" );

    return;
}

void remollGenC12::GenInelastic(G4double beamE,G4double theta,G4double eOut,
                               G4double &Q2,G4double &W2,G4double &effectiveXsection,
                               G4double &fWeight,G4double &asym) {

    G4double F1 = 0.0;
    G4double F2 = 0.0;
    G4double W1 = 0.0;
    G4double xsect = 0.0;

    const G4double CTH = cos(theta/2.0);
    const G4double STH = sin(theta/2.0);
    const G4double T2THE = STH*STH/CTH/CTH;
    const G4double Nu = beamE - eOut;
    Q2 = 4.0*beamE*eOut*STH*STH;
    W2 = pow(proton_mass_c2,2) + 2.0*proton_mass_c2*Nu - Q2;
    
    // Mott scattering
    G4double MOTT = pow((hbarc/GeV/m*fine_structure_const/(2.*beamE/GeV)*CTH/STH/STH),2)/microbarn; // units: ubarn

    G4int A=12;
    G4int Z=6;
    G4int bad=F1F2IN09(Z, A, Q2/GeV/GeV, W2/GeV/GeV, F1, F2);  
    if(bad != 0){
      G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << G4endl;	
      G4cerr << "  result was (-1 = Q2,W2 out of range in resmodd | -2 = A<3 | 1 inf/nan F1/F2) : " << bad <<G4endl;
      G4cerr << "     Q2 W2 A "<<Q2/GeV/GeV<<" "<<W2/GeV/GeV<<" "<<A<< G4endl;
      G4cerr << "     STH eOut beamE "<<STH<<" "<<eOut<<" "<<beamE<<G4endl;
      asym=0.99999;
      effectiveXsection=0;
      exit(1);
    }
    W1 = F1/proton_mass_c2;
 
    xsect = MOTT*(F2/(Nu/GeV) + 2.0*T2THE*(W1*GeV));

    if (xsect < 0 || std::isnan(xsect) ) {
        G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << G4endl;	
        G4cerr<<"Inelatic xsection problem: "<<xsect<<" "<<F1<<" "<<F2
              <<" "<<theta<<" "<<Q2/GeV/GeV<<" "<<W2/GeV/GeV<<" "<<G4endl;
        asym=0.99999;
        effectiveXsection=0;
        exit(1);
    }
//---Yuxiang Zhao, R-L asymmetry is negative
    asym=-1.0*Q2*0.8e-4/GeV/GeV;
    fWeight = xsect*sin(theta);
    effectiveXsection=xsect;
}




void remollGenC12::GenQuasiElastic(G4double beamE,G4double theta,
				  G4double &Q2,G4double &W2,G4double &effectiveXsection,
				  G4double &fWeight,G4double &eOut,G4double &asym) {  
  G4double F1 = 0.0;
  G4double F2 = 0.0;
  G4double w1 = 0.0;
  G4double w2 = 0.0;  
  G4double xsect = 0.0;
  
  // Generate flat energy distribution of outgoing electron
  eOut =  electron_mass_c2 + G4UniformRand()*(beamE - electron_mass_c2);
  
  G4double CTH = cos(theta/2.0);
  G4double STH = sin(theta/2.0);
  G4double T2THE = STH*STH/CTH/CTH;
  G4double Nu = beamE - eOut;//[MeV]
  Q2 = 4.0*beamE*eOut*STH*STH;
  W2 = proton_mass_c2*proton_mass_c2 + 2.0*proton_mass_c2*Nu - Q2;
  
  // Mott scattering
  G4double MOTT = pow((hbarc/GeV/m*fine_structure_const/(2.*beamE/GeV)*CTH/STH/STH),2)/microbarn; // units: ubarn

  G4int A=12;
  G4int Z=6;
  F1F2QE09(Z, A, Q2/GeV/GeV, W2/GeV/GeV, F1, F2);
  
  w1 = F1/(proton_mass_c2/GeV);
  w2 = F2/(Nu/GeV);
  
  xsect = MOTT*(w2 + 2.0*T2THE*w1);
  
  // In some cases a negative F2 is returned giving a negative cross section
  if (xsect <= 0) xsect = 0.0;

  fWeight = xsect*sin(theta);
  effectiveXsection = xsect;

  ///~~~ Aymmetry calculation
  asym= -GF/(4.*pi*fine_structure_const*sqrt(2.)) * Q2 * (QWp);
}

// copy from QweakSimEPEvent::Elastic_Cross_Section_Carbon
// Calculate the Elastic Carbon cross-section using PWBA Form Factors (FF) 
// deterimined using the method described in HANS F. EHRBNBERG et. all Phys. Rev. 113,2 (1959)
void remollGenC12::GenElastic(G4double beamE,G4double theta,
			     G4double &Q2,G4double &W2,G4double &effectiveXsection,
			     G4double &fWeight,G4double &eOut,G4double &asym) {
  // Physical Constants      
  G4double M_A = 931.494 * MeV;                 // a.m.u. in MeV
  G4double Z = 6.0;                             // # of protons
  G4double A = 12.0;                            // Atomic Weight
  G4double M = M_A*A;                           // Nuclear Mass
  G4double a = 1.65;                            // well parameter [fm]
  G4double alpha_FS = 1.0/137.035999074;        // Fine Structure Constant
  G4double m_e = 0.511 * MeV;                   // electron mass in MeV

  // Kinematic variables
  G4double CTH = cos(theta/2.0);
  G4double STH = sin(theta/2.0);
  G4double ETA = 1.0+2.0*beamE*STH*STH/M;
  eOut = beamE/ETA;
  Q2 = 4*beamE*eOut*STH*STH;                    // [MeV^2]
  W2 = M*M;
  G4double myhbarc = hbarc / MeV / fermi;       // 197.3269631 MeV fm
  G4double q2 = Q2/(myhbarc*myhbarc);           //convert MeV^2 into fm^(-2)

  // Dynamics
  // Mott Cross Section
  G4double SigmaMott = pow(alpha_FS*myhbarc*Z/(2.0*beamE),2);    // [fm^2]
  SigmaMott *= (CTH*CTH)/pow(STH,4)*(eOut/beamE)*10000;         // [ub] 10000 is the conversion factor          
  G4double FF = (1.0 - (q2*a*a)/9.0)*exp(-(q2*a*a)/4.0);        // Form Factor
  G4double sigma = SigmaMott*FF*FF;
  // Schwinger correction from Mo and Tsai Rev.Mod.Phys. 41 205 (1969) eq. (II.2)
  G4double delta_Schwinger = 0.0;
  G4double SchwingerDeltaE = 15*MeV;  // this has to be checked --- hanjie
  if(SchwingerDeltaE!=0) {
    G4double FunctionofTheta = log (STH*STH) * log (CTH*CTH);
    delta_Schwinger = (-2.0*alpha_FS/pi) * ((log(beamE/SchwingerDeltaE) - 13.0/12.0)
                             * (log(Q2/(m_e*m_e)) - 1.0) + 17.0/36.0 + FunctionofTheta/2.0);
  }
  sigma = sigma*(1.0 + delta_Schwinger);
  effectiveXsection = sigma; 
  fWeight = sigma*sin(theta);
  asym= -GF/(4.*pi*fine_structure_const*sqrt(2.)) * Q2 * (QWp + QWn*(A-Z)/Z);

  return;
}
