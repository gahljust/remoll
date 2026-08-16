#include "remollPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Version.hh"
#include "G4Exception.hh"
#include "Randomize.hh"

#include "remollHEPEvtInterface.hh"
#ifdef G4LIB_USE_HEPMC
#include "HepMCG4AsciiInterface.hh"
#ifdef G4LIB_USE_PYTHIA
#include "HepMCG4PythiaInterface.hh"
#endif
#endif

#include "remollBeamTarget.hh"
#include "remollVEventGen.hh"
#include "remollEvent.hh"
#include "remollRun.hh"
#include "remollRunData.hh"
#include "remolltypes.hh"
#include "globals.hh"

#include "remollGenMoller.hh"
#include "remollGenpElastic.hh"
#include "remollGenpInelastic.hh"
#include "remollGenPion.hh"
#include "remollGenBeam.hh"
#include "remollGenC12.hh"
#include "remollGenFixedInelasticC12.hh"
#include "remollGenFlat.hh"
#include "remollGenExternal.hh"
#include "remollGenAl.hh"
#include "remollGenLUND.hh"
#include "remollGenHyperon.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <mutex>
#include <sstream>

namespace {

struct NeymanTransportPlanState {
  std::mutex mutex;
  G4String file;
  G4int pBins = 0;
  G4int thetaBins = 0;
  G4int invalidGroup = -1;
  std::vector<G4double> pEdges;
  std::vector<G4double> thetaEdges;
  std::vector<G4int> cellGroup;
  std::vector<G4double> probability;
  std::vector<long long> quota;
  std::vector<long long> accepted;
  long long remaining = 0;
  long long processed = 0;
  G4bool loaded = false;
};

NeymanTransportPlanState& NeymanPlanState() {
  static NeymanTransportPlanState state;
  return state;
}

void SamplerFatal(const char* code, const G4String& message) {
  G4Exception("remollPrimaryGeneratorAction", code, FatalException,
              message.c_str());
}

}  // namespace

remollPrimaryGeneratorAction::remollPrimaryGeneratorAction()
  : fEventGen(0),fPriGen(0),fParticleGun(0),fEvent(0),fRateCopy(0),fGeneratorOnly(false),
    fTransportGateEnabled(false),fNeymanTransportEnabled(false),
    fWriteEvent(true),
    fTransportGatePBins(0),fTransportGateThetaBins(0),
    fEffCrossSection(0)
{
    static bool has_been_warned = false;
    if (! has_been_warned) {
      G4cout << "remoll: All possible event generators are instantiated every time." << G4endl;
      G4cout << "remoll: This means some will not find necessary input files or" << G4endl;
      G4cout << "remoll: print other information in the next few lines." << G4endl;
      has_been_warned = true;
    }

    // Populate map with all possible event generators
    fEvGenMap["moller"] = std::make_shared<remollGenMoller>();
    fEvGenMap["elastic"] = std::make_shared<remollGenpElastic>();
    fEvGenMap["inelastic"] = std::make_shared<remollGenpInelastic>();
    fEvGenMap["pion"] = std::make_shared<remollGenPion>();
    fEvGenMap["beam"] = std::make_shared<remollGenBeam>();
    fEvGenMap["flat"] = std::make_shared<remollGenFlat>();
    fEvGenMap["elasticAl"] = std::make_shared<remollGenAl>(0);
    fEvGenMap["quasielasticAl"] = std::make_shared<remollGenAl>(1);
    fEvGenMap["inelasticAl"] = std::make_shared<remollGenAl>(2);
    fEvGenMap["external"] = std::make_shared<remollGenExternal>();
    fEvGenMap["pion_LUND"] = std::make_shared<remollGenLUND>();
    fEvGenMap["elasticC12"] = std::make_shared<remollGenC12>(0);
    fEvGenMap["quasielasticC12"] = std::make_shared<remollGenC12>(1);
    fEvGenMap["inelasticC12"] = std::make_shared<remollGenC12>(2);
    fEvGenMap["fixedinelasticC12"] = std::make_shared<remollGenFixedInelasticC12>();
    fEvGenMap["hyperon"] = std::make_shared<remollGenHyperon>();

    // Populate map with all possible primary generators
    fPriGenMap["particlegun"] = std::make_shared<G4ParticleGun>();
    fPriGenMap["HEPEvt"] = std::make_shared<remollHEPEvtInterface>();
    #ifdef G4LIB_USE_HEPMC
    fPriGenMap["hepmcAscii"] = std::make_shared<HepMCG4AsciiInterface>();
    #ifdef G4LIB_USE_PYTHIA
    fPriGenMap["hepmcPythia"] = std::make_shared<HepMCG4PythiaInterface>();
    #endif
    #endif

    // Default generator
    G4String default_generator = "moller";
    SetGenerator(default_generator);

    // Create event generator messenger
    fEvGenMessenger.DeclareMethod("set",&remollPrimaryGeneratorAction::SetGenerator,"Select physics generator");
    fEvGenMessenger.DeclarePropertyWithUnit("sigma","picobarn",fEffCrossSection,"Set effective cross section");
    fEvGenMessenger.DeclareProperty("copyRate",fRateCopy,"ExtGen: copy rate from previous sim");
    fEvGenMessenger.DeclareProperty("generatorOnly",fGeneratorOnly,"Generate and write ev/bm/part/rate branches without Geant4 tracking");
    fEvGenMessenger.DeclareMethod("transportGateFile",&remollPrimaryGeneratorAction::LoadTransportGate,
        "Load a p',theta transport keep-probability TSV file");
    fEvGenMessenger.DeclareProperty("transportGateEnable",fTransportGateEnabled,
        "Apply the loaded unbiased p',theta gate before Geant4 transport");
    fEvGenMessenger.DeclareMethod("transportNeymanFile",
        &remollPrimaryGeneratorAction::LoadNeymanTransportPlan,
        "Load a fixed-quota p',theta Neyman transport plan");
    fEvGenMessenger.DeclareProperty("transportNeymanEnable",fNeymanTransportEnabled,
        "Apply the loaded fixed-quota Neyman plan before Geant4 transport");
}

void remollPrimaryGeneratorAction::LoadTransportGate(G4String path)
{
    std::ifstream input(path.c_str());
    if (!input) {
      fTransportGateEnabled = false;
      SamplerFatal("remoll-sampler-001",
                   "Unable to open transport gate " + path);
      return;
    }

    struct Row { G4int ip, it; G4double plo, phi, tlo, thi, keep; };
    std::vector<Row> rows;
    G4String line;
    std::getline(input, line);
    G4int pBins = 0, thetaBins = 0;
    while (std::getline(input, line)) {
      std::stringstream stream(line);
      std::vector<G4String> field;
      G4String value;
      while (std::getline(stream, value, '\t')) field.push_back(value);
      if (field.size() < 8) continue;
      Row row{std::stoi(field[0]), std::stoi(field[1]),
              std::stod(field[2]), std::stod(field[3]),
              std::stod(field[4]), std::stod(field[5]), std::stod(field[7])};
      if (!std::isfinite(row.keep) || row.keep <= 0.0 || row.keep > 1.0) continue;
      rows.push_back(row);
      pBins = std::max(pBins, row.ip + 1);
      thetaBins = std::max(thetaBins, row.it + 1);
    }
    if (rows.empty() || pBins <= 0 || thetaBins <= 0) {
      fTransportGateEnabled = false;
      SamplerFatal("remoll-sampler-002",
                   "No valid bins in transport gate " + path);
      return;
    }

    std::vector<G4double> pEdges(pBins + 1, 0.0), thetaEdges(thetaBins + 1, 0.0);
    std::vector<G4double> keep(pBins * thetaBins, 1.0);
    std::vector<G4bool> seen(pBins * thetaBins, false);
    for (const auto& row : rows) {
      if (row.ip < 0 || row.ip >= pBins || row.it < 0 || row.it >= thetaBins
          || seen[row.ip * thetaBins + row.it]) {
        fTransportGateEnabled = false;
        SamplerFatal("remoll-sampler-003",
                     "Duplicate or invalid cell in transport gate " + path);
        return;
      }
      pEdges[row.ip] = row.plo;
      pEdges[row.ip + 1] = row.phi;
      thetaEdges[row.it] = row.tlo;
      thetaEdges[row.it + 1] = row.thi;
      keep[row.ip * thetaBins + row.it] = row.keep;
      seen[row.ip * thetaBins + row.it] = true;
    }
    if (std::find(seen.begin(), seen.end(), false) != seen.end()
        || !std::is_sorted(pEdges.begin(), pEdges.end())
        || !std::is_sorted(thetaEdges.begin(), thetaEdges.end())) {
      fTransportGateEnabled = false;
      SamplerFatal("remoll-sampler-004",
                   "Incomplete or unordered transport gate " + path);
      return;
    }
    fTransportGateFile = path;
    fTransportGatePBins = pBins;
    fTransportGateThetaBins = thetaBins;
    fTransportGatePEdges.swap(pEdges);
    fTransportGateThetaEdges.swap(thetaEdges);
    fTransportGateKeep.swap(keep);
    G4cout << "Loaded transport gate " << path << " (" << pBins << " x "
           << thetaBins << " bins)" << G4endl;
}

G4double remollPrimaryGeneratorAction::TransportGateProbability() const
{
    if (!fTransportGateEnabled || fTransportGateKeep.empty() || fEvent == nullptr)
      return 1.0;

    G4int selected = -1;
    G4int selectedAny = -1;
    for (size_t index = 0; index < fEvent->fPartType.size(); ++index) {
      if (fEvent->fPartType[index] == nullptr
          || fEvent->fPartType[index]->GetPDGEncoding() != 11) continue;
      if (selectedAny < 0 || fEvent->fPartRealMom[index].mag()
          > fEvent->fPartRealMom[selectedAny].mag()) selectedAny = index;
      if (fEvent->fPartRealMom[index].z() > 0.0
          && (selected < 0 || fEvent->fPartRealMom[index].mag()
          > fEvent->fPartRealMom[selected].mag())) selected = index;
    }
    if (selected < 0) selected = selectedAny;
    if (selected < 0) return 1.0;

    const auto& momentum = fEvent->fPartRealMom[selected];
    const G4double p = momentum.mag() / MeV;
    const G4double theta = std::atan2(momentum.perp(), momentum.z()) * 1000.0;
    if (!std::isfinite(p) || !std::isfinite(theta)
        || p < fTransportGatePEdges.front() || p >= fTransportGatePEdges.back()
        || theta < fTransportGateThetaEdges.front()
        || theta >= fTransportGateThetaEdges.back()) return 1.0;
    const G4int ip = std::upper_bound(fTransportGatePEdges.begin(),
        fTransportGatePEdges.end(), p) - fTransportGatePEdges.begin() - 1;
    const G4int it = std::upper_bound(fTransportGateThetaEdges.begin(),
        fTransportGateThetaEdges.end(), theta) - fTransportGateThetaEdges.begin() - 1;
    if (ip < 0 || ip >= fTransportGatePBins || it < 0 || it >= fTransportGateThetaBins)
      return 1.0;
    return fTransportGateKeep[ip * fTransportGateThetaBins + it];
}

void remollPrimaryGeneratorAction::LoadNeymanTransportPlan(G4String path)
{
    auto& state = NeymanPlanState();
    std::lock_guard<std::mutex> lock(state.mutex);
    if (state.loaded && state.file == path) return;
    std::ifstream input(path.c_str());
    if (!input) {
      fNeymanTransportEnabled = false;
      SamplerFatal("remoll-sampler-005",
                   "Unable to open Neyman transport plan " + path);
      return;
    }
    struct Row {
      G4int ip, it, group;
      G4double plo, phi, tlo, thi, probability;
      long long quota;
    };
    std::vector<Row> rows;
    G4String line;
    std::getline(input, line);
    G4int pBins = 0, thetaBins = 0, groups = 0;
    while (std::getline(input, line)) {
      std::stringstream stream(line);
      std::vector<G4String> field;
      G4String value;
      while (std::getline(stream, value, '\t')) field.push_back(value);
      if (field.size() < 9) continue;
      Row row{std::stoi(field[0]), std::stoi(field[1]), std::stoi(field[6]),
              std::stod(field[2]), std::stod(field[3]),
              std::stod(field[4]), std::stod(field[5]),
              std::stod(field[7]), std::stoll(field[8])};
      if (row.group < 0 || !std::isfinite(row.probability)
          || row.probability < 0.0 || row.quota <= 0) continue;
      rows.push_back(row);
      if (row.ip >= 0 && row.it >= 0) {
        pBins = std::max(pBins, row.ip + 1);
        thetaBins = std::max(thetaBins, row.it + 1);
      }
      groups = std::max(groups, row.group + 1);
    }
    if (rows.empty() || pBins <= 0 || thetaBins <= 0 || groups <= 0) {
      fNeymanTransportEnabled = false;
      SamplerFatal("remoll-sampler-006",
                   "No valid bins in Neyman transport plan " + path);
      return;
    }
    state.file = path;
    state.pBins = pBins;
    state.thetaBins = thetaBins;
    state.invalidGroup = -1;
    state.pEdges.assign(pBins + 1, 0.0);
    state.thetaEdges.assign(thetaBins + 1, 0.0);
    state.cellGroup.assign(pBins * thetaBins, -1);
    state.probability.assign(groups, -1.0);
    state.quota.assign(groups, -1);
    for (const auto& row : rows) {
      if ((state.probability[row.group] >= 0.0
           && (std::abs(state.probability[row.group] - row.probability) > 1.0e-12
               || state.quota[row.group] != row.quota))) {
        state.loaded = false;
        fNeymanTransportEnabled = false;
        SamplerFatal("remoll-sampler-007",
                     "Inconsistent Neyman stratum in " + path);
        return;
      }
      state.probability[row.group] = row.probability;
      state.quota[row.group] = row.quota;
      if (row.ip < 0 || row.it < 0) {
        if (state.invalidGroup >= 0 && state.invalidGroup != row.group) {
          state.loaded = false;
          fNeymanTransportEnabled = false;
          SamplerFatal("remoll-sampler-008",
                       "Multiple invalid-event strata in " + path);
          return;
        }
        state.invalidGroup = row.group;
      } else {
        if (row.ip >= pBins || row.it >= thetaBins
            || state.cellGroup[row.ip * thetaBins + row.it] >= 0) {
          state.loaded = false;
          fNeymanTransportEnabled = false;
          SamplerFatal("remoll-sampler-009",
                       "Duplicate or invalid Neyman cell in " + path);
          return;
        }
        state.pEdges[row.ip] = row.plo; state.pEdges[row.ip + 1] = row.phi;
        state.thetaEdges[row.it] = row.tlo;
        state.thetaEdges[row.it + 1] = row.thi;
        state.cellGroup[row.ip * thetaBins + row.it] = row.group;
      }
    }
    state.remaining = 0;
    G4double probabilityTotal = 0.0;
    for (G4int group = 0; group < groups; ++group) {
      if (state.probability[group] < 0.0 || state.quota[group] <= 0) {
        state.loaded = false;
        fNeymanTransportEnabled = false;
        SamplerFatal("remoll-sampler-010",
                     "Incomplete Neyman transport group in " + path);
        return;
      }
      probabilityTotal += state.probability[group];
      state.remaining += state.quota[group];
    }
    if (state.invalidGroup < 0
        || std::find(state.cellGroup.begin(), state.cellGroup.end(), -1)
             != state.cellGroup.end()
        || !std::is_sorted(state.pEdges.begin(), state.pEdges.end())
        || !std::is_sorted(state.thetaEdges.begin(), state.thetaEdges.end())
        || std::abs(probabilityTotal - 1.0) > 5.0e-6) {
      state.loaded = false;
      fNeymanTransportEnabled = false;
      SamplerFatal("remoll-sampler-011",
                   "Incomplete or unnormalized Neyman plan " + path);
      return;
    }
    state.accepted.assign(groups, 0);
    state.processed = 0;
    state.loaded = true;
    G4cout << "Loaded Neyman transport plan " << path << " (" << groups
           << " strata, " << state.remaining << " transports)" << G4endl;
}

G4bool remollPrimaryGeneratorAction::ApplyNeymanTransportPlan(
    G4double& correction) const
{
    correction = 1.0;
    if (!fNeymanTransportEnabled || fEvent == nullptr) return true;
    auto& state = NeymanPlanState();
    G4int selected = -1, selectedAny = -1;
    for (size_t index = 0; index < fEvent->fPartType.size(); ++index) {
      if (fEvent->fPartType[index] == nullptr
          || fEvent->fPartType[index]->GetPDGEncoding() != 11) continue;
      if (selectedAny < 0 || fEvent->fPartRealMom[index].mag()
          > fEvent->fPartRealMom[selectedAny].mag()) selectedAny = index;
      if (fEvent->fPartRealMom[index].z() > 0.0
          && (selected < 0 || fEvent->fPartRealMom[index].mag()
          > fEvent->fPartRealMom[selected].mag())) selected = index;
    }
    if (selected < 0) selected = selectedAny;

    std::lock_guard<std::mutex> lock(state.mutex);
    if (!state.loaded) return true;
    const G4double nthrown = remollRun::GetRunData()->GetNthrown();
    ++state.processed;
    auto report_if_last = [&state, nthrown]() {
      if (state.processed != static_cast<long long>(nthrown)) return;
      if (state.remaining == 0) {
        G4cout << "Neyman transport validation: all quotas filled" << G4endl;
      } else {
        std::ostringstream message;
        message << "Neyman transport plan underfilled by " << state.remaining
                << " histories; output is invalid";
        SamplerFatal("remoll-sampler-012", message.str());
      }
    };
    G4int group = state.invalidGroup;
    if (selected >= 0) {
      const auto& momentum = fEvent->fPartRealMom[selected];
      const G4double p = momentum.mag() / MeV;
      const G4double theta = std::atan2(momentum.perp(), momentum.z()) * 1000.0;
      if (std::isfinite(p) && std::isfinite(theta)
          && p >= state.pEdges.front() && p < state.pEdges.back()
          && theta >= state.thetaEdges.front() && theta < state.thetaEdges.back()) {
        const G4int ip = std::upper_bound(state.pEdges.begin(), state.pEdges.end(), p)
            - state.pEdges.begin() - 1;
        const G4int it = std::upper_bound(state.thetaEdges.begin(), state.thetaEdges.end(), theta)
            - state.thetaEdges.begin() - 1;
        group = state.cellGroup[ip * state.thetaBins + it];
      }
    }
    if (group < 0 || group >= G4int(state.quota.size())) {
      SamplerFatal("remoll-sampler-013",
                   "Generated event does not map to a Neyman stratum");
      return false;
    }
    if (state.accepted[group] >= state.quota[group]) {
      report_if_last();
      return false;
    }
    ++state.accepted[group];
    --state.remaining;
    correction = nthrown * state.probability[group] / state.quota[group];
    if (state.remaining == 0)
      G4cout << "Neyman transport plan filled all quotas" << G4endl;
    report_if_last();
    return true;
}

remollPrimaryGeneratorAction::~remollPrimaryGeneratorAction()
{
}

void remollPrimaryGeneratorAction::SetGenerator(G4String& genname)
{
    // Set generator to null
    fEventGen = 0;
    fPriGen = 0;

    // Find event generator
    auto evgen = fEvGenMap.find(genname);
    if (evgen != fEvGenMap.end()) {
      G4cout << "Setting generator to " << genname << G4endl;
      fPriGen = 0;
      fPriGenName = "";
      fEventGen = evgen->second;
      fEventGenName = evgen->first;
      fParticleGun = fEventGen->GetParticleGun();
    }

    // Find primary generator
    auto prigen = fPriGenMap.find(genname);
    if (prigen != fPriGenMap.end()) {
      G4cout << "Setting generator to " << genname << G4endl;
      fPriGen = prigen->second;
      fPriGenName = prigen->first;
      fEventGen = 0;
      fEventGenName = "";
      fParticleGun = 0;
    }

    // No generator found
    if (!fEventGen && !fPriGen) {
      G4cerr << __FILE__ << " line " << __LINE__ << " - ERROR generator " << genname << " invalid" << G4endl;
      exit(1);
    }

    // Set the beam target
    if (fEventGen) {
      fEventGen->SetBeamTarget(&fBeamTarg);
    }
}

void remollPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    fWriteEvent = true;
    if (!fEventGen && !fPriGen) {
      G4cerr << __FILE__ << " line " << __LINE__ << " - No event generator found." << G4endl;
      exit(1);
    }

    // Delete old primary event
    if (fEvent != nullptr) {
      delete fEvent;
      fEvent = 0;
    }

    // 1. Using primary generator interface
    if (!fEventGen && fPriGen) {
      fPriGen->GeneratePrimaryVertex(anEvent);
      fEvent = new remollEvent(anEvent);
      fEvent->SetEffCrossSection(fEffCrossSection);
    }

    // 2. Using event generator interface
    if (fEventGen && !fPriGen) {

      // Helper function
      auto contains = [](const G4String& lhs, const G4String& rhs) {
        #if G4VERSION_NUMBER < 1100
          return lhs.contains(rhs);
        #else
          return G4StrUtil::contains(lhs, rhs);
        #endif
      };

      // Set beam polarization
      const G4String fBeamPol = fEventGen->GetBeamPolarization();
      G4ThreeVector cross(0,0,2);
      if (fBeamPol == "0") cross = G4ThreeVector(0,0,0);
      else {
        if     (contains(fBeamPol, "V")) cross = G4ThreeVector(1,0,0);
        else if(contains(fBeamPol, "H")) cross = G4ThreeVector(0,1,0);
        if (contains(fBeamPol, "-")) cross *= -1;
      }

      // Create new primary event
      fEvent = fEventGen->GenerateEvent();
      G4double transportCorrection = 1.0;
      G4bool transportEvent = true;
      if (fNeymanTransportEnabled) {
        transportEvent = ApplyNeymanTransportPlan(transportCorrection);
      } else {
        const G4double transportKeep = TransportGateProbability();
        // Do not perturb the legacy random stream when the gate is off or b == 1.
        transportEvent = transportKeep >= 1.0 || G4UniformRand() < transportKeep;
        transportCorrection = 1.0 / transportKeep;
      }
      if (transportEvent) {
        fEvent->fBiasWeight *= transportCorrection;
      } else {
        // Rejected histories remain part of Nthrown but carry zero contribution.
        fEvent->fEffXs = 0.0;
        fEvent->fRate = 0.0;
      }
      // A fixed-quota analysis knows Nthrown and n_h from its archived macro
      // and plan, so writing rejected candidates only duplicates bulky event,
      // geometry-map, and RNG-state branches.  Gate samples remain unchanged
      // because their IID analyzer currently uses the complete candidate tree.
      fWriteEvent = !fNeymanTransportEnabled || transportEvent;
      for (unsigned int pidx = 0; transportEvent && !fGeneratorOnly && pidx < fEvent->fPartType.size(); pidx++) {

        double p = fEvent->fPartRealMom[pidx].mag();
        double m = fEvent->fPartType[pidx]->GetPDGMass();
        double kinE = sqrt(p*p + m*m) - m;

        fParticleGun->SetParticleDefinition(fEvent->fPartType[pidx]);
        fParticleGun->SetParticleEnergy(kinE);
        fParticleGun->SetParticlePosition(fEvent->fPartPos[pidx]);
        fParticleGun->SetParticleMomentumDirection(fEvent->fPartRealMom[pidx].unit());

        G4ThreeVector pol(0,0,0);
        if (pidx == 0) {
          if (cross.mag() !=0) {
            if (cross.mag() == 1) //transverse polarization
              pol = G4ThreeVector( (fEvent->fPartRealMom[0].unit()).cross(cross));
            else if (contains(fBeamPol, "+") ) //positive helicity
              pol = fEvent->fPartRealMom[0].unit();
            else //negative helicity
              pol = - fEvent->fPartRealMom[0].unit();
          }
        }
        fParticleGun->SetParticlePolarization(pol);

        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }

    // Finally set the cross section and rate

    // Get number of thrown events
    G4double nthrown = remollRun::GetRunData()->GetNthrown();


    // Calculate rate
    SamplingType_t sampling_type = fEventGen->GetSamplingType();
    fEvent->fEffXs *= fEvent->fBiasWeight;
    if (fEvent->fRate != 0) fEvent->fRate *= fEvent->fBiasWeight;
    if (fEvent->fRate == 0) { // If the rate is set to 0 then calculate it using the cross section
        fEvent->fRate  = fEvent->fEffXs * fBeamTarg.GetEffLumin(sampling_type) / nthrown;

    } else if(!fRateCopy){ // For LUND - calculate rate and cross section
        fEvent->fEffXs = fEvent->fRate * nthrown / fBeamTarg.GetEffLumin(sampling_type);
        fEvent->fRate  = fEvent->fRate / nthrown;
    }

}
