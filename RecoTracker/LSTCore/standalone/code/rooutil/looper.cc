//  .
// ..: P. Chang, philip@physics.ucsd.edu

//#include "cpptqdm/tqdm.h"
#include "looper.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Event Looper (Looper) class template implementation
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////

//_________________________________________________________________________________________________
RooUtil::Looper::Looper()
    : tchain(0),
      listOfFiles(0),
      fileIter(0),
      tfile(0),
      ttree(0),
      ps(0),
      nEventsTotalInChain(0),
      nEventsTotalInTree(0),
      nEventsToProcess(-1),
      nEventsProcessed(0),
      indexOfEventInTTree(0),
      fastmode(true),
      treeclass(0),
      bar_id(0),
      print_rate(432),
      doskim(false),
      skimfilename(""),
      skimfile(0),
      skimtree(0),
      nEventsSkimmed(0),
      silent(false),
      isinit(false),
      use_treeclass_progress(false),
      //    use_tqdm_progress_bar( isatty(1) ),
      nskipped_batch(0),
      nskipped(0),
      nbatch_skip_threshold(500),
      nbatch_to_skip(5000),
      nskipped_threshold(100000),
      ncounter(0),
      teventlist(0) {
  bmark = new TBenchmark();
  //    bar.disable_colors();
}

//_________________________________________________________________________________________________
RooUtil::Looper::Looper(TChain* c, TreeUtil* t, int nevtToProc)
    : tchain(0),
      listOfFiles(0),
      fileIter(0),
      tfile(0),
      ttree(0),
      ps(0),
      nEventsTotalInChain(0),
      nEventsTotalInTree(0),
      nEventsToProcess(nevtToProc),
      nEventsProcessed(0),
      indexOfEventInTTree(0),
      fastmode(true),
      treeclass(0),
      bar_id(0),
      print_rate(432),
      doskim(false),
      skimfilename(""),
      skimfile(0),
      skimtree(0),
      nEventsSkimmed(0),
      silent(false),
      isinit(false),
      use_treeclass_progress(false),
      //    use_tqdm_progress_bar( isatty(1) ),
      nskipped_batch(0),
      nskipped(0),
      nbatch_skip_threshold(500),
      nbatch_to_skip(5000),
      nskipped_threshold(100000),
      ncounter(0),
      teventlist(0) {
  bmark = new TBenchmark();
  if (c && t)
    init(c, t, nevtToProc);
  //    bar.disable_colors();
}

//_________________________________________________________________________________________________
void RooUtil::Looper::init(TChain* c, TreeUtil* t, int nevtToProc) {
  listOfFiles = 0;
  if (fileIter)
    delete fileIter;
  fileIter = 0;
  tfile = 0;
  ttree = 0;
  ps = 0;
  nEventsTotalInChain = 0;
  nEventsTotalInTree = 0;
  nEventsToProcess = nevtToProc;
  nEventsProcessed = 0;
  indexOfEventInTTree = 0;
  fastmode = true;
  treeclass = 0;
  bar_id = 0;
  print_rate = 432;
  doskim = false;
  skimfilename = "";
  skimfile = 0;
  skimtree = 0;
  nEventsSkimmed = 0;
  silent = false;
  isinit = false;
  use_treeclass_progress = false;
  //    use_tqdm_progress_bar( isatty(1) ),
  nskipped_batch = 0;
  nskipped = 0;
  nbatch_skip_threshold = 500;
  nbatch_to_skip = 5000;
  nskipped_threshold = 100000;
  ncounter = 0;
  teventlist = 0;

  if (isinit)
    error(
        "The Looper is already initialized! Are you calling Looper::init(TChain* c, TreeUtil* t, int nevtToProcess) "
        "for the second time?",
        __FUNCTION__);

  initProgressBar();
  print("Start EventLooping");
  start();

  nEventsToProcess = nevtToProc;

  if (c)
    setTChain(c);

  if (t)
    setTreeClass(t);

  if (nEventsToProcess > 5000 || nEventsToProcess == -1)
    fastmode = true;

  c->GetEntry(0);

  // Check tree exists
  if (not c->GetTree()) {
    throw std::ios_base::failure("Failed to get TTree from input ntuple");
  }

  t->Init(c->GetTree());

  bmark->Start("benchmark");
  isinit = true;
}

//_________________________________________________________________________________________________
RooUtil::Looper::~Looper() {
  if (isinit) {
    end();

    // return
    using namespace std;
    bmark->Stop("benchmark");
    cout << endl;
    cout << "------------------------------" << endl;
    cout << "CPU  Time:	" << Form("%.01f", bmark->GetCpuTime("benchmark")) << endl;
    cout << "Real Time:	" << Form("%.01f", bmark->GetRealTime("benchmark")) << endl;
    cout << endl;
    // delete bmark;

    //        if ( fileIter )
    //            delete fileIter;
    //
    //        if ( tfile )
    //            delete tfile;
  }
}

//_________________________________________________________________________________________________
void RooUtil::Looper::setTChain(TChain* c) {
  if (c) {
    tchain = c;
    setNEventsToProcess();
    setFileList();
  } else
    error("You provided a null TChain pointer!", __FUNCTION__);
}

//_________________________________________________________________________________________________
void RooUtil::Looper::setTreeClass(TreeUtil* t) {
  if (t)
    treeclass = t;
  else
    error("You provided a null TreeClass pointer!", __FUNCTION__);
}

//_________________________________________________________________________________________________
void RooUtil::Looper::printCurrentEventIndex() {
  RooUtil::print(TString::Format("Current TFile = %s", tfile->GetName()));
  RooUtil::print(TString::Format("Current TTree = %s", ttree->GetName()));
  RooUtil::print(TString::Format("Current Entry # in TTree = %d", indexOfEventInTTree));
}

//_________________________________________________________________________________________________
bool RooUtil::Looper::nextTree() {
  if (!fileIter)
    error("fileIter not set but you are trying to access the next file", __FUNCTION__);

  // Get the TChainElement from TObjArrayIter.
  // If no more to run over, Next returns 0.
  TChainElement* chainelement = (TChainElement*)fileIter->Next();

  if (chainelement) {
    // If doskim is true and if this is the very first file being opened in the TChain,
    // flag it to create a tfile and ttree where the skimmed events will go to.
    bool createskimtree = false;

    if (!ttree && doskim)
      createskimtree = true;

    // // If there is already a TFile opened from previous iteration, close it.
    // if ( tfile )
    //     tfile->Close();

    // Open up a new file
    tfile = TFile::Open(chainelement->GetTitle());
    // Get the ttree
    ttree = (TTree*)tfile->Get(tchain->GetName());

    // If an eventindexmap has a key for this file then set the TEventList for this TTree
    // std::cout <<  " chainelement->GetTitle(): " << chainelement->GetTitle() <<  std::endl;
    // std::cout <<  " eventindexmap.hasEventList(chainelement->GetTitle()): " << eventindexmap.hasEventList(chainelement->GetTitle()) <<  std::endl;
    if (eventindexmap.hasEventList(chainelement->GetTitle())) {
      // std::cout << chainelement->GetTitle() << std::endl;
      teventlist = eventindexmap.getEventList(chainelement->GetTitle());
      ttree->SetEventList(teventlist);
    } else {
      teventlist = 0;
    }

    if (!ttree)
      error("TTree is null!??", __FUNCTION__);

    // Set some fast mode stuff
    if (fastmode) {
      TTreeCache::SetLearnEntries(1000);
      print("TTreeCache enabled");
    }

    if (fastmode)
      ttree->SetCacheSize(128 * 1024 * 1024);
    else
      ttree->SetCacheSize(-1);

    // Print some info to stdout
    print("Looping " + TString(tfile->GetName()) + "/TTree:" + TString(ttree->GetName()));
    //        printProgressBar(true);
    // Reset the nEventsTotalInTree in this tree
    nEventsTotalInTree = ttree->GetEntries();
    // Reset the event index as we got a new ttree
    indexOfEventInTTree = 0;
    // Set the ttree to the TreeUtil
    treeclass->Init(ttree);

    // If skimming create the skim tree after the treeclass inits it.
    // This is to make sure the branch addresses are correct.
    if (createskimtree)
      createSkimTree();
    else if (doskim)
      copyAddressesToSkimTree();

    //        // TTreePerfStats
    //        if ( ps )
    //            ps->SaveAs( "perf.root" );
    //
    //        ps = new TTreePerfStats( "ioperf", ttree );
    // Return that I got a good one
    return true;
  } else {
    // Announce that we are done with this chain
    //        print("");
    //        print("Done with all trees in this chain", __FUNCTION__);
    return false;
  }
}

//_________________________________________________________________________________________________
bool RooUtil::Looper::allEventsInTreeProcessed() {
  if (teventlist) {
    if (indexOfEventInTTree >= (unsigned int)teventlist->GetN()) {
      unsigned int curindex =
          teventlist->GetEntry(indexOfEventInTTree - 1);  // since I just increased by one a few lines before
      unsigned int previndex = indexOfEventInTTree >= 2 ? teventlist->GetEntry(indexOfEventInTTree - 2) : 0;
      nEventsToProcess += (curindex - previndex);
      return true;
    } else {
      return false;
    }
  } else {
    if (indexOfEventInTTree >= nEventsTotalInTree)
      return true;
    else
      return false;
  }
}

//_________________________________________________________________________________________________
bool RooUtil::Looper::allEventsInChainProcessed() {
  if (nEventsProcessed >= (unsigned int)nEventsToProcess)
    return true;
  else
    return false;
}

//_________________________________________________________________________________________________
bool RooUtil::Looper::nextEventInTree() {
  //    treeclass->progress(nEventsProcessed, nEventsToProcess);
  // Sanity check before loading the next event.
  if (!ttree)
    error("current ttree not set!", __FUNCTION__);

  if (!tfile)
    error("current tfile not set!", __FUNCTION__);

  if (!fileIter)
    error("fileIter not set!", __FUNCTION__);

  // Check whether I processed everything
  if (allEventsInTreeProcessed())
    return false;

  if (allEventsInChainProcessed())
    return false;

  // if fast mode do some extra
  if (fastmode)
    ttree->LoadTree(teventlist ? teventlist->GetEntry(indexOfEventInTTree) : indexOfEventInTTree);

  // Set the event index in TreeUtil
  treeclass->GetEntry(teventlist ? teventlist->GetEntry(indexOfEventInTTree) : indexOfEventInTTree);
  // Increment the counter for this ttree
  ++indexOfEventInTTree;
  // Increment the counter for the entire tchain
  // If there is teventlist set then the skipping is a bit more complex...
  if (teventlist) {
    unsigned int curindex =
        teventlist->GetEntry(indexOfEventInTTree - 1);  // since I just increased by one a few lines before
    unsigned int previndex = indexOfEventInTTree >= 2 ? teventlist->GetEntry(indexOfEventInTTree - 2) : 0;
    nEventsToProcess += (curindex - previndex);
  } else {
    ++nEventsProcessed;
  }
  // Print progress
  printProgressBar();
  // If all fine return true
  return true;
}

//_________________________________________________________________________________________________
bool RooUtil::Looper::nextEvent() {
  if (!isinit)
    error(
        "The Looper is not initialized! please call properly Looper::init(TChain* c, TreeUtil* t, int nevtToProcess) "
        "first!",
        __FUNCTION__);

  // If no tree it means this is the beginning of the loop.
  if (!ttree) {
    //        std::cout << " I think this is the first tree " << std::endl;
    // Load the next tree if it returns true, then proceed to next event in tree.
    while (nextTree()) {
      // If the next event in tree was successfully loaded return true, that it's good.
      if (nextEventInTree()) {
        //                std::cout << " I think this is the first event in first tree" << std::endl;
        // Set the boolean that a new file opened for this event
        isnewfileopened = true;
        return true;
      }
      // If the first event in this tree was not good, continue to the next tree
      else
        continue;
    }

    // If looping over all trees, we fail to find first event that's good,
    // return false and call it quits.
    // At this point it will exit the loop without processing any events.
    //        printProgressBar();
    // Set the boolean that a new file has not opened for this event
    isnewfileopened = false;
    return false;
  }
  // If tree exists, it means that we're in the middle of a loop
  else {
    // If next event is successfully loaded proceed.
    if (nextEventInTree()) {
      // Set the boolean that a new file has not opened for this event
      isnewfileopened = false;
      return true;
    }
    // If next event is not loaded then check why.
    else {
      // If failed because it was the last event in the whole chain to process, exit the loop.
      // You're done!
      if (allEventsInChainProcessed()) {
        //                printProgressBar();
        // Set the boolean that a new file has not opened for this event
        isnewfileopened = false;
        return false;
      }
      // If failed because it's last in the tree then load the next tree and the event
      else if (allEventsInTreeProcessed()) {
        // Load the next tree if it returns true, then proceed to next event in tree.
        while (nextTree()) {
          // If the next event in tree was successfully loaded return true, that it's good.
          if (nextEventInTree()) {
            // Set the boolean that a new file has opened for this event
            isnewfileopened = true;
            return true;
          }
          // If the first event in this tree was not good, continue to the next tree
          else
            continue;
        }

        // If looping over all trees, we fail to find first event that's good,
        // return false and call it quits.
        // Again you're done!
        //                printProgressBar();
        // Set the boolean that a new file has not opened for this event
        isnewfileopened = false;
        return false;
      } else {
        // Why are you even here?
        // spit error and return false to avoid warnings
        error("You should not be here! Please contact philip@physics.ucsd.edu", __FUNCTION__);
        // Set the boolean that a new file has not opened for this event
        isnewfileopened = false;
        return false;
      }
    }
  }
}

//_________________________________________________________________________________________________
bool RooUtil::Looper::isNewFileInChain() { return isnewfileopened; }

//_________________________________________________________________________________________________
void RooUtil::Looper::setFileList() {
  if (!fileIter) {
    listOfFiles = tchain->GetListOfFiles();
    fileIter = new TObjArrayIter(listOfFiles);
  }
}

//_________________________________________________________________________________________________
void RooUtil::Looper::setNEventsToProcess() {
  if (tchain) {
    nEventsTotalInChain = tchain->GetEntries();

    if (nEventsToProcess < 0)
      nEventsToProcess = nEventsTotalInChain;

    if (nEventsToProcess > (int)nEventsTotalInChain) {
      print(TString::Format("Asked to process %d events, but there aren't that many events", nEventsToProcess));
      nEventsToProcess = nEventsTotalInChain;
    }

    print(TString::Format("Total Events in this Chain to process = %d", nEventsToProcess));
  }
}

//_________________________________________________________________________________________________
void RooUtil::Looper::initProgressBar() {
  /// Init progress bar
  my_timer.Start();
  bar_id = 0;
}

//_________________________________________________________________________________________________
void RooUtil::Looper::printProgressBar(bool force) {
  if (silent)
    return;

  /// Print progress bar

  int entry = nEventsProcessed;
  int totalN = nEventsToProcess;

  //    if (use_tqdm_progress_bar)
  //    {
  //        if (force) return;  // N.B. If i am not using my own scheme i shouldn't force it.
  //        bar.progress(nEventsProcessed-1, nEventsToProcess); // tqdm expects 0 to N-1 index not 1 to N
  //        return;
  //    }

  if (use_treeclass_progress) {
    if (force)
      return;  // N.B. If i am not using my own scheme i shouldn't force it.
    // treeclass->progress(nEventsProcessed, nEventsToProcess);
    return;
  }

  if (totalN < 20)
    totalN = 20;

  // Progress bar
  if (entry > totalN) {
    //        printf( "Why are you here?\n" );
  } else if (entry == totalN) {
    Double_t elapsed = my_timer.RealTime();
    Double_t rate;

    if (elapsed != 0)
      rate = entry / elapsed;
    else
      rate = -999;

    const int mins_in_hour = 60;
    const int secs_to_min = 60;
    Int_t input_seconds = elapsed;
    Int_t seconds = input_seconds % secs_to_min;
    Int_t minutes = input_seconds / secs_to_min % mins_in_hour;
    Int_t hours = input_seconds / secs_to_min / mins_in_hour;

    printf("\rRooUtil::");
    printf("+");
    printf("|====================");

    //for ( int nb = 0; nb < 20; ++nb )
    //{
    //  printf("=");
    //}

    printf("| %.1f %% (%d/%d) with  [avg. %d Hz]   Total Time: %.2d:%.2d:%.2d         \n",
           100.0,
           entry,
           totalN,
           (int)rate,
           hours,
           minutes,
           seconds);
    fflush(stdout);
  }
  //else if ( entry % ( ( ( int ) print_rate ) ) < (0.3) * print_rate || force )
  else if (entry % (((int)print_rate)) == 0 || force) {
    // sanity check
    if (entry >=
        totalN +
            10)  // +2 instead of +1 since, the loop might be a while loop where to check I got a bad event the index may go over 1.
    {
      TString msg = TString::Format("%d %d", entry, totalN);
      RooUtil::print(msg, __FUNCTION__);
      RooUtil::error("Total number of events processed went over max allowed! Check your loop boundary conditions!!",
                     __FUNCTION__);
    }

    int nbars = entry / (totalN / 20);
    Double_t elapsed = my_timer.RealTime();
    Double_t rate;

    if (elapsed != 0)
      rate = entry / elapsed;
    else
      rate = -999;

    Double_t percentage = entry / (totalN * 1.) * 100;
    const int mins_in_hour = 60;
    const int secs_to_min = 60;
    Int_t input_seconds = (totalN - entry) / rate;
    Int_t seconds = input_seconds % secs_to_min;
    Int_t minutes = input_seconds / secs_to_min % mins_in_hour;
    Int_t hours = input_seconds / secs_to_min / mins_in_hour;

    print_rate = (int)(rate / 5) + 1;

    printf("RooUtil:: ");

    if (bar_id % 4 == 3)
      printf("-");

    if (bar_id % 4 == 2)
      printf("/");

    if (bar_id % 4 == 1)
      printf("|");

    if (bar_id % 4 == 0)
      printf("\\");

    printf("|");
    bar_id++;

    for (int nb = 0; nb < 20; ++nb) {
      if (nb < nbars)
        printf("=");
      else
        printf(".");
    }

    printf("| %.1f %% (%d/%d) with  [%d Hz]   ETA %.2d:%.2d:%.2d         \r",
           percentage,
           entry + 1,
           totalN,
           (int)rate,
           hours,
           minutes,
           seconds);
    fflush(stdout);
  }

  my_timer.Start(kFALSE);
}

//_________________________________________________________________________________________________
void RooUtil::Looper::setSkim(TString ofilename) {
  skimfilename = ofilename;
  doskim = true;
}

//_________________________________________________________________________________________________
void RooUtil::Looper::createSkimTree() {
  skimfile = new TFile(skimfilename, "recreate");
  TObjArray* toa = ttree->GetListOfBranches();

  if (skimbrfiltpttn.size() > 0) {
    ttree->SetBranchStatus("*", 0);

    for (auto& pttn : skimbrfiltpttn) {
      for (const auto& brobj : *toa) {
        TString brname = brobj->GetName();

        if (pttn.Contains("*")) {
          TString modpttn = pttn;
          modpttn.ReplaceAll("*", "");
          if (brname.Contains(modpttn) && brname.BeginsWith(modpttn)) {
            // std::cout << brname << std::endl;
            ttree->SetBranchStatus(brname + "*", 1);
          }
        } else {
          if (brname.EqualTo(pttn)) {
            // std::cout << brname << std::endl;
            ttree->SetBranchStatus(brname, 1);
          }
        }
      }
    }
  }

  skimtree = ttree->CloneTree(0);
}

//_________________________________________________________________________________________________
void RooUtil::Looper::copyAddressesToSkimTree() { ttree->CopyAddresses(skimtree); }

//_________________________________________________________________________________________________
void RooUtil::Looper::fillSkim() {
  treeclass->loadAllBranches();
  skimtree->Fill();
  nEventsSkimmed++;
}

//_________________________________________________________________________________________________
void RooUtil::Looper::saveSkim() {
  double frac_skimmed = (double)nEventsSkimmed / (double)nEventsProcessed * 100;
  RooUtil::print(Form("Skimmed events %d out of %d. [%f%%]", nEventsSkimmed, nEventsProcessed, frac_skimmed));
  skimtree->GetCurrentFile()->cd();
  skimtree->Write();
  //    skimfile->Close();
}

//_________________________________________________________________________________________________
bool RooUtil::Looper::handleBadEvent() {
  using namespace std;
  cout << endl;
  cout << "RooUtil::Looper [CheckCorrupt] Caught an I/O failure in the ROOT file." << endl;
  cout << "RooUtil::Looper [CheckCorrupt] Possibly corrupted hadoop file." << endl;
  cout << "RooUtil::Looper [CheckCorrupt] event index = " << getCurrentEventIndex() << " out of "
       << tchain->GetEntries() << endl;
  cout << endl;

  // If the total nskip reaches a threshold just fail the whole thing...
  if (nskipped >= nskipped_threshold) {
    nskipped += tchain->GetEntries() - getCurrentEventIndex() - 1;
    return false;
  }

  nskipped_batch++;

  // If the nskipped is quite large than skip the entire file
  if (nskipped_batch > nbatch_skip_threshold) {
    nskipped += nskipped_batch;
    nskipped_batch = 0;
    for (unsigned int i = 0; i < nbatch_to_skip; ++i) {
      if (!nextEvent())
        return false;
      nskipped++;
    }
  }

  return true;
}

//_________________________________________________________________________________________________
void RooUtil::Looper::printStatus() {
  getTree()->PrintCacheStats();
  printSkippedBadEventStatus();
}

//_________________________________________________________________________________________________
TString RooUtil::Looper::getListOfFileNames() {
  TString rtnstring = "";
  TObjArray* filepaths = tchain->GetListOfFiles();
  TObjArrayIter* iter = new TObjArrayIter(listOfFiles);
  for (Int_t ifile = 0; ifile < filepaths->GetEntries(); ++ifile) {
    TChainElement* chainelement = (TChainElement*)iter->Next();
    if (chainelement) {
      TString filepath = chainelement->GetTitle();
      if (rtnstring.IsNull())
        rtnstring = filepath;
      else
        rtnstring += "," + filepath;
    }
  }
  return rtnstring;
}

//_________________________________________________________________________________________________
void RooUtil::Looper::printSkippedBadEventStatus() {
  using namespace std;
  nskipped += nskipped_batch;

  if (nskipped) {
    cout << "RooUtil:Looper [CheckCorrupt] Skipped " << nskipped << " events out of " << tchain->GetEntries() << " ["
         << float(nskipped) / float(tchain->GetEntries()) * 100 << "% loss]"
         << " POSSIBLE BADFILES = " << getListOfFileNames() << endl;
  }
}

//_________________________________________________________________________________________________
bool RooUtil::Looper::doesBranchExist(TString bname) {
  if (ttree->GetBranch(bname))
    return true;
  if (ttree->GetBranch(ttree->GetAlias(bname)))
    return true;
  return false;
}
