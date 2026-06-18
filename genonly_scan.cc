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
    event->fRate = event->fEffXs * lumin / nthrown;
  } else {
    event->fEffXs = event->fRate * nthrown / lumin;
    event->fRate = event->fRate / nthrown;
  }
}

int ParticlePid(const remollEvent* event, size_t index)
{
  if (index >= event->fPartType.size() || event->fPartType[index] == nullptr) return 0;
  return event->fPartType[index]->GetPDGEncoding();
}

G4ThreeVector ParticleMom(const remollEvent* event, size_t index)
{
  if (index >= event->fPartMom.size()) return G4ThreeVector();
  return event->fPartMom[index];
}

G4ThreeVector ParticleRealMom(const remollEvent* event, size_t index)
{
  if (index >= event->fPartRealMom.size()) return G4ThreeVector();
  return event->fPartRealMom[index];
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
  Int_t npart = 0;
  Int_t pid0 = 0, pid1 = 0;
  Double_t rate = 0, xs = 0, A = 0, Am = 0, Q2 = 0, W2 = 0, xbj = 0;
  Double_t thcom = 0, thcom_deg = 0;
  Double_t beamp = 0, beamp_GeV = 0, beamE = 0, beamE_GeV = 0, sampledE = 0, sampledE_GeV = 0;
  Double_t vx = 0, vy = 0, vz = 0;
  Double_t bmx = 0, bmy = 0, bmz = 0, bmth = 0, bmph = 0;
  Double_t radlen = 0, travelled = 0, eff_mat_len = 0;
  Double_t p0 = 0, p0_GeV = 0, th0 = 0, ph0 = 0, px0 = 0, py0 = 0, pz0 = 0;
  Double_t rp0 = 0, rth0 = 0, rph0 = 0, rpx0 = 0, rpy0 = 0, rpz0 = 0;
  Double_t p1 = 0, p1_GeV = 0, th1 = 0, ph1 = 0, px1 = 0, py1 = 0, pz1 = 0;
  Double_t rp1 = 0, rth1 = 0, rph1 = 0, rpx1 = 0, rpy1 = 0, rpz1 = 0;

  tree.Branch("entry", &entry, "entry/l");
  tree.Branch("npart", &npart, "npart/I");
  tree.Branch("pid0", &pid0, "pid0/I");
  tree.Branch("pid1", &pid1, "pid1/I");
  tree.Branch("rate", &rate, "rate/D");
  tree.Branch("xs", &xs, "xs/D");
  tree.Branch("A", &A, "A/D");
  tree.Branch("Am", &Am, "Am/D");
  tree.Branch("Q2", &Q2, "Q2/D");
  tree.Branch("W2", &W2, "W2/D");
  tree.Branch("xbj", &xbj, "xbj/D");
  tree.Branch("thcom", &thcom, "thcom/D");
  tree.Branch("thcom_deg", &thcom_deg, "thcom_deg/D");
  tree.Branch("beamp", &beamp, "beamp/D");
  tree.Branch("beamp_GeV", &beamp_GeV, "beamp_GeV/D");
  tree.Branch("beamE", &beamE, "beamE/D");
  tree.Branch("beamE_GeV", &beamE_GeV, "beamE_GeV/D");
  tree.Branch("sampledE", &sampledE, "sampledE/D");
  tree.Branch("sampledE_GeV", &sampledE_GeV, "sampledE_GeV/D");
  tree.Branch("vx", &vx, "vx/D");
  tree.Branch("vy", &vy, "vy/D");
  tree.Branch("vz", &vz, "vz/D");
  tree.Branch("bmx", &bmx, "bmx/D");
  tree.Branch("bmy", &bmy, "bmy/D");
  tree.Branch("bmz", &bmz, "bmz/D");
  tree.Branch("bmth", &bmth, "bmth/D");
  tree.Branch("bmph", &bmph, "bmph/D");
  tree.Branch("radlen", &radlen, "radlen/D");
  tree.Branch("travelled", &travelled, "travelled/D");
  tree.Branch("eff_mat_len", &eff_mat_len, "eff_mat_len/D");
  tree.Branch("p0", &p0, "p0/D");
  tree.Branch("p0_GeV", &p0_GeV, "p0_GeV/D");
  tree.Branch("th0", &th0, "th0/D");
  tree.Branch("ph0", &ph0, "ph0/D");
  tree.Branch("px0", &px0, "px0/D");
  tree.Branch("py0", &py0, "py0/D");
  tree.Branch("pz0", &pz0, "pz0/D");
  tree.Branch("rp0", &rp0, "rp0/D");
  tree.Branch("rth0", &rth0, "rth0/D");
  tree.Branch("rph0", &rph0, "rph0/D");
  tree.Branch("rpx0", &rpx0, "rpx0/D");
  tree.Branch("rpy0", &rpy0, "rpy0/D");
  tree.Branch("rpz0", &rpz0, "rpz0/D");
  tree.Branch("p1", &p1, "p1/D");
  tree.Branch("p1_GeV", &p1_GeV, "p1_GeV/D");
  tree.Branch("th1", &th1, "th1/D");
  tree.Branch("ph1", &ph1, "ph1/D");
  tree.Branch("px1", &px1, "px1/D");
  tree.Branch("py1", &py1, "py1/D");
  tree.Branch("pz1", &pz1, "pz1/D");
  tree.Branch("rp1", &rp1, "rp1/D");
  tree.Branch("rth1", &rth1, "rth1/D");
  tree.Branch("rph1", &rph1, "rph1/D");
  tree.Branch("rpx1", &rpx1, "rpx1/D");
  tree.Branch("rpy1", &rpy1, "rpy1/D");
  tree.Branch("rpz1", &rpz1, "rpz1/D");

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

    npart = static_cast<Int_t>(remoll_event->fPartType.size());
    pid0 = ParticlePid(remoll_event.get(), 0);
    pid1 = ParticlePid(remoll_event.get(), 1);

    rate = remoll_event->fRate * second;
    xs = remoll_event->fEffXs / microbarn;
    A = remoll_event->fAsym / ppb;
    Am = remoll_event->fmAsym / ppb;
    Q2 = remoll_event->fQ2;
    W2 = remoll_event->fW2;
    xbj = remoll_event->fXbj;
    thcom = remoll_event->fThCoM;
    thcom_deg = remoll_event->fThCoM / deg;

    beamp = remoll_event->fBeamMomentum.mag();
    beamp_GeV = beamp / GeV;
    beamE = remoll_event->fBeamE;
    beamE_GeV = beamE / GeV;
    sampledE = event_beam_target != nullptr ? event_beam_target->fSampledEnergy : 0;
    sampledE_GeV = sampledE / GeV;

    vx = remoll_event->fVertexPos.x();
    vy = remoll_event->fVertexPos.y();
    vz = remoll_event->fVertexPos.z();
    bmx = remoll_event->fBeamMomentum.x();
    bmy = remoll_event->fBeamMomentum.y();
    bmz = remoll_event->fBeamMomentum.z();
    bmth = remoll_event->fBeamMomentum.theta();
    bmph = remoll_event->fBeamMomentum.phi();

    radlen = event_beam_target != nullptr ? event_beam_target->fRadiationLength : 0;
    travelled = event_beam_target != nullptr ? event_beam_target->fTravelledLength : 0;
    eff_mat_len = event_beam_target != nullptr ? event_beam_target->fEffectiveMaterialLength : 0;

    const auto mom0 = ParticleMom(remoll_event.get(), 0);
    p0 = mom0.mag();
    p0_GeV = p0 / GeV;
    th0 = mom0.theta();
    ph0 = mom0.phi();
    px0 = mom0.x();
    py0 = mom0.y();
    pz0 = mom0.z();

    const auto real_mom0 = ParticleRealMom(remoll_event.get(), 0);
    rp0 = real_mom0.mag();
    rth0 = real_mom0.theta();
    rph0 = real_mom0.phi();
    rpx0 = real_mom0.x();
    rpy0 = real_mom0.y();
    rpz0 = real_mom0.z();

    const auto mom1 = ParticleMom(remoll_event.get(), 1);
    p1 = mom1.mag();
    p1_GeV = p1 / GeV;
    th1 = mom1.theta();
    ph1 = mom1.phi();
    px1 = mom1.x();
    py1 = mom1.y();
    pz1 = mom1.z();

    const auto real_mom1 = ParticleRealMom(remoll_event.get(), 1);
    rp1 = real_mom1.mag();
    rth1 = real_mom1.theta();
    rph1 = real_mom1.phi();
    rpx1 = real_mom1.x();
    rpy1 = real_mom1.y();
    rpz1 = real_mom1.z();

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
