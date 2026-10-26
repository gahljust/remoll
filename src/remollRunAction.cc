#include "remollRunAction.hh"
#include "remollNeymanPlan.hh"
#include "G4Exception.hh"

#include "G4RunManager.hh"
#include "G4Timer.hh"

#include "remollIO.hh"
#include "remollRun.hh"
#include "remollRunData.hh"
#include "remollBeamTarget.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"
namespace { G4Mutex remollRunActionMutex = G4MUTEX_INITIALIZER; }

remollRunAction::remollRunAction()
{
  fMessenger.DeclareMethod(
      "seed",
      &remollRunAction::UpdateSeed,
      "Set random engine seed")
      .SetParameterName("seed", false)
      .SetStates(G4State_PreInit,G4State_Idle);
  fMessenger.DeclareMethod(
      "interval",
      &remollRunAction::SetUpdateInterval,
      "Print this many progress points (i.e. 100 -> every 1%)")
      .SetParameterName("interval", false)
      .SetStates(G4State_PreInit,G4State_Idle);

  // Configure the process-wide plan on the master before BeginOfRunAction.
  // Worker-only generator commands are replayed too late for this shared state.
  if (!G4Threading::IsWorkerThread()) {
    fNeymanMessenger.reset(new G4GenericMessenger(this, "/remoll/evgen/", "Neyman sampling"));
    fNeymanMessenger->DeclareMethod("transportNeymanFile", &remollRunAction::LoadNeymanPlan,
        "Load a fixed-quota momentum-angle plan")
        .SetStates(G4State_PreInit, G4State_Idle).SetToBeBroadcasted(false);
    fNeymanMessenger->DeclareMethod("transportNeymanEnable", &remollRunAction::EnableNeyman,
        "Enable fixed-quota transport sampling (default false)")
        .SetStates(G4State_PreInit, G4State_Idle).SetToBeBroadcasted(false);
  }

  // Store random status before primary particle generation
  G4RunManager::GetRunManager()->StoreRandomNumberStatusToG4Event(1);
}

void remollRunAction::LoadNeymanPlan(G4String path) {
    try {
        remollNeymanPlan::Instance().Load(path);
    } catch (const std::exception& e) {
        G4Exception("LoadNeymanPlan", "remoll-neyman", FatalException, e.what());
    }
}
void remollRunAction::EnableNeyman(G4bool enabled) {
    try {
        remollNeymanPlan::Instance().Enable(enabled);
    } catch (const std::exception& e) {
        G4Exception("EnableNeyman", "remoll-neyman", FatalException, e.what());
    }
}

void remollRunAction::UpdateSeed(const G4long seed)
{
  G4Random::setTheSeed(seed);
  remollRun::GetRunData()->SetSeed(seed);
  G4cout << "Random seed set to " << seed << G4endl;
}

G4Run* remollRunAction::GenerateRun()
{
  return new remollRun();
}

void remollRunAction::BeginOfRunAction(const G4Run* run)
{
  // Cast into remollRun
  const remollRun* aRun = static_cast<const remollRun*>(run);

  // Print progress
  G4int evts_to_process = aRun->GetNumberOfEventToBeProcessed();
  G4RunManager::GetRunManager()->SetPrintProgress((evts_to_process > fInterval)
                                                  ? evts_to_process/ fInterval
                                                  : 1);

  if (IsMaster())
  {
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

    // Initialize shared quotas once on the master; the requested count includes rejected candidates.
    try {
      remollNeymanPlan::Instance().BeginRun(aRun->GetNumberOfEventToBeProcessed());
    } catch (const std::exception& e) {
      G4Exception("BeginOfRunAction", "remoll-neyman", FatalException, e.what());
      return;
    }
    fTimer.Start();

    G4AutoLock lock(&remollRunActionMutex);
    remollIO* io = remollIO::GetInstance();
    io->InitializeTree();

    remollRunData *rundata = remollRun::GetRunData();
    rundata->SetNthrown( aRun->GetNumberOfEventToBeProcessed() );
  }
}

void remollRunAction::EndOfRunAction(const G4Run* run)
{
  // Cast into remollRun
  const remollRun* aRun = static_cast<const remollRun*>(run);

  if (IsMaster())
  {
      // Require a complete run before WriteTree() writes the final ROOT output.
      try {
          remollNeymanPlan::Instance().EndRun(aRun->GetNumberOfEvent() == aRun->GetNumberOfEventToBeProcessed());
      } catch (const std::exception& e) {
          G4Exception("EndOfRunAction", "remoll-neyman", FatalException, e.what());
          return;
      }
      fTimer.Stop();

      G4cout << "### Run " << aRun->GetRunID() << " ended "
             << "(" << fTimer.GetUserElapsed() << "s)." << G4endl;

      G4AutoLock lock(&remollRunActionMutex);
      remollIO* io = remollIO::GetInstance();
      io->WriteTree();
  }
}

