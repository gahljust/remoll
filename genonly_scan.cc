/*!
  Fast remoll generator-only scanner.

  This executable reuses remollPrimaryGeneratorAction so the generator,
  radiative-loss, cross-section, asymmetry, and rate math remains remoll's
  own implementation, but it avoids the Geant4 run/tracking loop and writes a
  compact ROOT tree with generator-level scalars only.
*/

#include "remollDetectorConstruction.hh"
#include "remollBeamTarget.hh"
#include "remollEvent.hh"
#include "remollGenAl.hh"
#include "remollGenBeam.hh"
#include "remollGenC12.hh"
#include "remollGenExternal.hh"
#include "remollGenFlat.hh"
#include "remollGenHyperon.hh"
#include "remollGenLUND.hh"
#include "remollGenMoller.hh"
#include "remollGenPion.hh"
#include "remollGenpElastic.hh"
#include "remollGenpInelastic.hh"
#include "remollRun.hh"
#include "remollRunData.hh"
#include "remollVEventGen.hh"

#include "G4ParticleDefinition.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"

#ifdef __APPLE__
#include <unistd.h>
#endif

#include <algorithm>
#include <cctype>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Config {
  std::string macro;
  std::string geometry;
  std::string generator = "moller";
  std::string output = "genonly_scan.root";
  unsigned long long n_events = 0;
  unsigned long long progress_points = 100;
  G4long seed = 0;
  bool n_events_set = false;
  bool output_set = false;
};

struct MacroCommand {
  std::string command;
  std::string source;
  int line = 0;
};

std::string Trim(const std::string& text)
{
  const auto first = std::find_if_not(text.begin(), text.end(),
      [](unsigned char c) { return std::isspace(c); });
  const auto last = std::find_if_not(text.rbegin(), text.rend(),
      [](unsigned char c) { return std::isspace(c); }).base();
  return first < last ? std::string(first, last) : std::string();
}

std::string DirName(const std::string& path)
{
  const auto slash = path.find_last_of('/');
  return slash == std::string::npos ? "." : path.substr(0, slash);
}

std::string ResolvePath(const std::string& parent, const std::string& child)
{
  if (child.empty() || child.front() == '/') return child;

  const std::vector<std::string> candidates = {
      DirName(parent) + "/" + child,
      child};

  for (const auto& candidate : candidates) {
    std::ifstream test(candidate.c_str());
    if (test) return candidate;
  }

  return candidates.front();
}

bool StartsWith(const std::string& text, const std::string& prefix)
{
  return text.compare(0, prefix.size(), prefix) == 0;
}

void PrintUsage()
{
  G4cerr
      << "Usage: genonly_scan [-m macro] [-g geometry.gdml] [-n events] "
      << "[-o output.root] [-r seed] [--progress-points N] [macro]" << G4endl;
}

void ReadMacro(const std::string& path, std::vector<MacroCommand>& commands)
{
  std::ifstream input(path.c_str());
  if (!input) {
    G4cerr << "ERROR: could not open macro " << path << G4endl;
    exit(1);
  }

  std::string line;
  int line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    line = Trim(line);
    if (line.empty()) continue;

    if (StartsWith(line, "/control/execute ")) {
      const std::string nested = Trim(line.substr(std::string("/control/execute ").size()));
      ReadMacro(ResolvePath(path, nested), commands);
      continue;
    }

    commands.push_back({line, path, line_number});
  }
}

bool ParseBeamOn(const std::string& command, unsigned long long& n_events)
{
  if (!StartsWith(command, "/run/beamOn ")) return false;
  n_events = std::strtoull(Trim(command.substr(std::string("/run/beamOn ").size())).c_str(), nullptr, 10);
  return true;
}

bool ParseFilename(const std::string& command, std::string& output)
{
  if (!StartsWith(command, "/remoll/filename ")) return false;
  output = Trim(command.substr(std::string("/remoll/filename ").size()));
  return true;
}

bool ParseGeometry(const std::string& command, std::string& geometry)
{
  if (!StartsWith(command, "/remoll/setgeofile ")) return false;
  geometry = Trim(command.substr(std::string("/remoll/setgeofile ").size()));
  return true;
}

bool ParseGenerator(const std::string& command, std::string& generator)
{
  if (!StartsWith(command, "/remoll/evgen/set ")) return false;
  generator = Trim(command.substr(std::string("/remoll/evgen/set ").size()));
  return true;
}

bool ShouldSkipCommand(const std::string& command)
{
  return StartsWith(command, "/run/")
      || StartsWith(command, "/remoll/setgeofile ")
      || StartsWith(command, "/remoll/filename ")
      || StartsWith(command, "/remoll/evgen/set ")
      || StartsWith(command, "/remoll/evgen/generatorOnly ")
      || StartsWith(command, "/remoll/SD/")
      || StartsWith(command, "/remoll/physlist/")
      || StartsWith(command, "/remoll/parallel/")
      || StartsWith(command, "/remoll/kryptonite/")
      || StartsWith(command, "/vis/")
      || StartsWith(command, "/gui/")
      || StartsWith(command, "/tracking/")
      || StartsWith(command, "/process/")
      || StartsWith(command, "/random/");
}

void ApplyUsefulCommands(const std::vector<MacroCommand>& commands)
{
  G4UImanager* ui = G4UImanager::GetUIpointer();
  for (const auto& item : commands) {
    if (ShouldSkipCommand(item.command)) continue;
    const int status = ui->ApplyCommand(item.command);
    if (status != 0) {
      G4cerr << "ERROR: command failed in " << item.source << ":" << item.line
             << G4endl << "  " << item.command << G4endl;
      exit(1);
    }
  }
}

void InitializeParticles()
{
  G4BosonConstructor bosons;
  bosons.ConstructParticle();
  G4LeptonConstructor leptons;
  leptons.ConstructParticle();
  G4MesonConstructor mesons;
  mesons.ConstructParticle();
  G4BaryonConstructor baryons;
  baryons.ConstructParticle();
  G4IonConstructor ions;
  ions.ConstructParticle();
  G4ShortLivedConstructor short_lived;
  short_lived.ConstructParticle();
}

std::shared_ptr<remollVEventGen> MakeGenerator(const std::string& name)
{
  if (name == "moller") return std::make_shared<remollGenMoller>();
  if (name == "elastic") return std::make_shared<remollGenpElastic>();
  if (name == "inelastic") return std::make_shared<remollGenpInelastic>();
  if (name == "pion") return std::make_shared<remollGenPion>();
  if (name == "beam") return std::make_shared<remollGenBeam>();
  if (name == "flat") return std::make_shared<remollGenFlat>();
  if (name == "elasticAl") return std::make_shared<remollGenAl>(0);
  if (name == "quasielasticAl") return std::make_shared<remollGenAl>(1);
  if (name == "inelasticAl") return std::make_shared<remollGenAl>(2);
  if (name == "external") return std::make_shared<remollGenExternal>();
  if (name == "pion_LUND") return std::make_shared<remollGenLUND>();
  if (name == "elasticC12") return std::make_shared<remollGenC12>(0);
  if (name == "quasielasticC12") return std::make_shared<remollGenC12>(1);
  if (name == "inelasticC12") return std::make_shared<remollGenC12>(2);
  if (name == "hyperon") return std::make_shared<remollGenHyperon>();

  G4cerr << "ERROR: unsupported generator for genonly_scan: " << name << G4endl;
  exit(1);
}

void NormalizeRate(remollEvent* event, const remollBeamTarget& beam_target,
                   const remollVEventGen& generator, unsigned long long n_events)
{
  const G4double nthrown = static_cast<G4double>(n_events);
  const SamplingType_t sampling_type = generator.GetSamplingType();
  const G4double lumin = beam_target.GetEffLumin(sampling_type);

  if (event->fRate == 0) {
    event->fEffXs *= event->fBiasWeight;
    event->fRate = event->fEffXs * lumin / nthrown;
  } else {
    event->fRate *= event->fBiasWeight;
    event->fEffXs = event->fRate * nthrown / lumin;
    event->fRate = event->fRate / nthrown;
  }
}

G4ThreeVector ParticleMom(const remollEvent* event, size_t index)
{
  if (index >= event->fPartMom.size()) return G4ThreeVector();
  return event->fPartMom[index];
}

} // namespace

int main(int argc, char** argv)
{
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  gROOT->Reset();
#endif

  Config config;
  config.seed = time(0);
#ifdef __APPLE__
  config.seed += getpid();
#endif

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-m" && i + 1 < argc) config.macro = argv[++i];
    else if (arg == "-g" && i + 1 < argc) config.geometry = argv[++i];
    else if (arg == "-n" && i + 1 < argc) {
      config.n_events = std::strtoull(argv[++i], nullptr, 10);
      config.n_events_set = true;
    } else if (arg == "-o" && i + 1 < argc) {
      config.output = argv[++i];
      config.output_set = true;
    } else if (arg == "-r" && i + 1 < argc) {
      config.seed = std::atol(argv[++i]);
    } else if ((arg == "--progress-points" || arg == "-p") && i + 1 < argc) {
      config.progress_points = std::strtoull(argv[++i], nullptr, 10);
    } else if (!arg.empty() && arg.front() != '-') {
      config.macro = arg;
    } else {
      PrintUsage();
      return 1;
    }
  }

  std::vector<MacroCommand> commands;
  if (!config.macro.empty()) {
    ReadMacro(config.macro, commands);
    for (const auto& item : commands) {
      unsigned long long macro_events = 0;
      std::string macro_output;
      std::string macro_geometry;
      if (!config.n_events_set && ParseBeamOn(item.command, macro_events))
        config.n_events = macro_events;
      if (!config.output_set && ParseFilename(item.command, macro_output))
        config.output = macro_output;
      if (config.geometry.empty() && ParseGeometry(item.command, macro_geometry))
        config.geometry = macro_geometry;
      ParseGenerator(item.command, config.generator);
    }
  }

  if (config.n_events == 0) {
    G4cerr << "ERROR: set events with -n or /run/beamOn in the macro" << G4endl;
    return 1;
  }

  CLHEP::HepRandom::setTheSeed(config.seed);
  remollRun::GetRunData()->Init();
  remollRun::GetRunData()->SetSeed(config.seed);
  remollRun::GetRunData()->SetNthrown(config.n_events);
  if (!config.macro.empty())
    remollRun::GetRunData()->SetMacroFile(config.macro.c_str());

  remollDetectorConstruction detector("material", config.geometry.c_str());
  remollBeamTarget beam_target;
  InitializeParticles();
  std::shared_ptr<remollVEventGen> generator = MakeGenerator(config.generator);
  generator->SetBeamTarget(&beam_target);

  G4UImanager::GetUIpointer()->ApplyCommand("/remoll/geometry/overlapcheck false");
  detector.Construct();
  G4StateManager::GetStateManager()->SetNewState(G4State_Idle);

  ApplyUsefulCommands(commands);

  TFile output(config.output.c_str(), "RECREATE");
  if (output.IsZombie() || !output.IsOpen()) {
    G4cerr << "ERROR: could not open output file " << config.output << G4endl;
    return 1;
  }

  TTree tree("T", "remoll generator-only scan");

  ULong64_t entry = 0;
  Float_t xs = 0, rate = 0, beamp_GeV = 0, thcom_deg = 0;
  Float_t Q2 = 0, W2 = 0, A = 0;
  Float_t vx = 0, vy = 0, vz = 0;
  Float_t p0_GeV = 0, th0 = 0, ph0 = 0;
  Float_t radlen = 0, travelled = 0;
  Float_t bias_weight = 1, beam_bias_weight = 1, thcom_bias_weight = 1;

  tree.Branch("entry", &entry, "entry/l");
  tree.Branch("xs", &xs, "xs/F");
  tree.Branch("rate", &rate, "rate/F");
  tree.Branch("beamp_GeV", &beamp_GeV, "beamp_GeV/F");
  tree.Branch("thcom_deg", &thcom_deg, "thcom_deg/F");
  tree.Branch("Q2", &Q2, "Q2/F");
  tree.Branch("W2", &W2, "W2/F");
  tree.Branch("A", &A, "A/F");
  tree.Branch("vx", &vx, "vx/F");
  tree.Branch("vy", &vy, "vy/F");
  tree.Branch("vz", &vz, "vz/F");
  tree.Branch("p0_GeV", &p0_GeV, "p0_GeV/F");
  tree.Branch("th0", &th0, "th0/F");
  tree.Branch("ph0", &ph0, "ph0/F");
  tree.Branch("radlen", &radlen, "radlen/F");
  tree.Branch("travelled", &travelled, "travelled/F");
  tree.Branch("bias_weight", &bias_weight, "bias_weight/F");
  tree.Branch("beam_bias_weight", &beam_bias_weight, "beam_bias_weight/F");
  tree.Branch("thcom_bias_weight", &thcom_bias_weight, "thcom_bias_weight/F");

  const unsigned long long progress_step =
      config.progress_points == 0 ? 0 :
      std::max(1ULL, config.n_events / config.progress_points);

  for (entry = 0; entry < config.n_events; ++entry) {
    if (progress_step > 0 && entry % progress_step == 0) {
      G4cout << "Event " << entry << G4endl;
    }

    std::unique_ptr<remollEvent> remoll_event(generator->GenerateEvent());
    NormalizeRate(remoll_event.get(), beam_target, *generator, config.n_events);
    const remollBeamTarget* event_beam_target = remoll_event->GetBeamTarget();

    rate = static_cast<Float_t>(remoll_event->fRate * second);
    xs = static_cast<Float_t>(remoll_event->fEffXs / microbarn);
    beamp_GeV = static_cast<Float_t>(remoll_event->fBeamMomentum.mag() / GeV);
    thcom_deg = static_cast<Float_t>(remoll_event->fThCoM / deg);
    Q2 = static_cast<Float_t>(remoll_event->fQ2);
    W2 = static_cast<Float_t>(remoll_event->fW2);
    A = static_cast<Float_t>(remoll_event->fAsym / ppb);

    vx = static_cast<Float_t>(remoll_event->fVertexPos.x());
    vy = static_cast<Float_t>(remoll_event->fVertexPos.y());
    vz = static_cast<Float_t>(remoll_event->fVertexPos.z());

    radlen = static_cast<Float_t>(event_beam_target != nullptr ? event_beam_target->fRadiationLength : 0);
    travelled = static_cast<Float_t>(event_beam_target != nullptr ? event_beam_target->fTravelledLength : 0);
    bias_weight = static_cast<Float_t>(remoll_event->fBiasWeight);
    beam_bias_weight = static_cast<Float_t>(remoll_event->fBeamBiasWeight);
    thcom_bias_weight = static_cast<Float_t>(remoll_event->fThCoMBiasWeight);

    const auto mom0 = ParticleMom(remoll_event.get(), 0);
    p0_GeV = static_cast<Float_t>(mom0.mag() / GeV);
    th0 = static_cast<Float_t>(mom0.theta());
    ph0 = static_cast<Float_t>(mom0.phi());

    tree.Fill();
  }

  tree.Write();
  output.Close();

  if (progress_step > 0) {
    G4cout << "Event " << config.n_events << G4endl;
  }
  G4cout << "Wrote " << config.n_events << " generator-only events to "
         << config.output << G4endl;
  return 0;
}
