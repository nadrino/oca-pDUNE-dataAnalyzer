// Plot total number of hits vs sigma threshold
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TH1.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <regex>
#include <filesystem>

#include "CmdLineParser.h"
#include "Logger.h"
#include "ocaEvent.h"
#include "constants.h"

LoggerInit([]{ Logger::getUserHeader() << "[" << FILENAME << "]"; });

namespace {
  const int kNDet = oca::N_DETECTORS;
  const int kNChan = oca::N_CHANNELS;

  // Edge channel helper: mask first and last channel on each ASIC (assume 64-chan ASICs if 384 channels)
  inline bool isEdgeChannel(int ch) {
    if (kNChan % 64 != 0) return (ch == 0 || ch == (kNChan-1));
    int mod = ch % 64;
    return (mod == 0 || mod == 63);
  }
}

int main(int argc, char* argv[]) {
  CmdLineParser clp;
  clp.getDescription() << "> hits vs sigma plotter" << std::endl;
  clp.addDummyOption("Options");
  clp.addOption("inputRootFile",  {"-r", "--root-file"}, "Input ROOT file");
  clp.addOption("inputCalFile",   {"-c", "--cal-file"},  "Calibration file (.cal). If omitted, auto-select nearest CAL .cal from ROOT folder (prefer previous CAL run)");
  clp.addOption("runNumber",      {"-n", "--run-number"}, "Run number (optional; inferred from root filename if omitted)");
  clp.addOption("sigmaMin",       {"--smin"},             "Min sigma (default 1)");
  clp.addOption("sigmaMax",       {"--smax"},             "Max sigma (default 15)");
  clp.addOption("sigmaStep",      {"--sstep"},            "Step (default 1)");
  clp.addOption("outputPdf",      {"-o", "--output"},    "Output PDF path");
  clp.addDummyOption();

  LogInfo << clp.getDescription().str() << std::endl;
  LogInfo << clp.getConfigSummary() << std::endl;
  clp.parseCmdLine(argc, argv);
  LogThrowIf(clp.isNoOptionTriggered(), "No option was provided.");

  std::string calPath  = clp.isOptionTriggered("inputCalFile") ? clp.getOptionVal<std::string>("inputCalFile") : std::string("");
  std::string rootPath = clp.getOptionVal<std::string>("inputRootFile");
  double smin = clp.isOptionTriggered("sigmaMin")  ? clp.getOptionVal<double>("sigmaMin")  : 1.0;
  double smax = clp.isOptionTriggered("sigmaMax")  ? clp.getOptionVal<double>("sigmaMax")  : 15.0;
  double sstep= clp.isOptionTriggered("sigmaStep") ? clp.getOptionVal<double>("sigmaStep") : 1.0;
  std::string outPdf   = clp.isOptionTriggered("outputPdf") ? clp.getOptionVal<std::string>("outputPdf") : std::string("hits_vs_sigma.pdf");
  // Determine run number (for title)
  int runNumForTitle = -1;
  if (clp.isOptionTriggered("runNumber")) {
    runNumForTitle = clp.getOptionVal<int>("runNumber");
  } else {
    std::smatch m; std::regex reTitle("SCD_RUN([0-9]{5})_");
    std::string baseForTitle = std::filesystem::path(rootPath).filename().string();
    if (std::regex_search(baseForTitle, m, reTitle)) {
      runNumForTitle = std::stoi(m[1]);
    }
  }

  // If calibration file not provided, auto-select from the ROOT file's folder using nearest CAL policy
  if (calPath.empty()) {
    namespace fs = std::filesystem;
    fs::path rpath(rootPath);
    LogThrowIf(!fs::exists(rpath), "Input ROOT file does not exist");
    fs::path dir = rpath.parent_path();
    std::string base = rpath.filename().string();
    // Extract 5-digit run number from filename pattern SCD_RUNNNNNN_...
    int runNum = -1;
    if (clp.isOptionTriggered("runNumber")) {
      runNum = clp.getOptionVal<int>("runNumber");
    } else {
      std::smatch m; std::regex re("SCD_RUN([0-9]{5})_");
      if (std::regex_search(base, m, re)) {
        runNum = std::stoi(m[1]);
      }
    }
    LogThrowIf(runNum < 0, "Could not infer run number; pass --run-number or --cal-file");

    // Prefer same-run CAL first
    char buf[64]; snprintf(buf, sizeof(buf), "SCD_RUN%05d_", runNum);
    std::string sameRunPrefix(buf);
    std::string sameRunCalFound;
    for (auto &de : fs::directory_iterator(dir)) {
      if (!de.is_regular_file()) continue;
      auto name = de.path().filename().string();
      if (name.rfind(sameRunPrefix, 0) == 0 && name.find("CAL") != std::string::npos && de.path().extension() == ".cal") {
        sameRunCalFound = de.path().string();
        break;
      }
    }

    // Otherwise search previous CAL; if none, next CAL
    int bestPrev = -1; std::string bestPrevCal;
    int bestNext = -1; std::string bestNextCal;
    std::regex rre("SCD_RUN([0-9]{5}).*CAL.*\\.cal$");
    for (auto &de : fs::directory_iterator(dir)) {
      if (!de.is_regular_file()) continue;
      auto name = de.path().filename().string();
      std::smatch mm;
      if (std::regex_match(name, mm, rre)) {
        int rn = std::stoi(mm[1]);
        if (rn < runNum) {
          if (bestPrev < 0 || rn > bestPrev) { bestPrev = rn; bestPrevCal = de.path().string(); }
        } else if (rn > runNum) {
          if (bestNext < 0 || rn < bestNext) { bestNext = rn; bestNextCal = de.path().string(); }
        }
      }
    }

    if (!sameRunCalFound.empty()) {
      calPath = sameRunCalFound;
      LogInfo << "Auto-selected same-run CAL: " << calPath << std::endl;
    } else if (!bestPrevCal.empty()) {
      calPath = bestPrevCal;
      LogInfo << "Auto-selected previous CAL: " << calPath << std::endl;
    } else if (!bestNextCal.empty()) {
      calPath = bestNextCal;
      LogInfo << "Auto-selected next CAL: " << calPath << std::endl;
    } else {
      LogThrow("Could not locate a suitable CAL file in directory " + dir.string());
    }
  }

  LogInfo << "Reading calibration file: " << calPath << std::endl;
  std::ifstream calFile(calPath);
  LogThrowIf(!calFile.is_open(), "Calibration file could not be opened");

  std::vector<std::vector<float>> baseline(kNDet, std::vector<float>(kNChan, 0.0f));
  std::vector<std::vector<float>> sigma(kNDet,    std::vector<float>(kNChan, 0.0f));
  std::vector<std::vector<bool>>  mask(kNDet,     std::vector<bool>(kNChan, false));

  std::string line;
  for (int det=0; det<kNDet; ++det) {
    for (int i=0;i<18;i++) std::getline(calFile, line); // skip header
    for (int ch=0; ch<kNChan; ++ch) {
      if (!std::getline(calFile, line)) { LogError << "Unexpected end of calibration file" << std::endl; return 1; }
      std::istringstream iss(line);
      std::vector<float> vals; vals.reserve(8);
      std::string tok; float v;
      while (std::getline(iss, tok, ',')) { std::istringstream iv(tok); if (iv>>v) vals.push_back(v); }
      if (vals.size()!=8) { LogError << "Bad calib line: size=" << vals.size() << std::endl; return 1; }
      int channel = (int)vals[0];
      baseline[det][channel] = vals[3];
      sigma[det][channel]    = vals[5];
      int badflag = (int)vals[6]; // 0 good, 1 bad
      bool calibBad = (badflag != 0);
      mask[det][channel] = calibBad || isEdgeChannel(channel);
    }
  }

  LogInfo << "Opening ROOT file: " << rootPath << std::endl;
  TFile* f = TFile::Open(rootPath.c_str(), "READ");
  LogThrowIf(!f || !f->IsOpen(), "ROOT file not open");

  // Use 3 trees as in event_display
  std::vector<TTree*> trees;
  trees.emplace_back( (TTree*)f->Get("raw_events") );
  trees.emplace_back( (TTree*)f->Get("raw_events_B") );
  trees.emplace_back( (TTree*)f->Get("raw_events_C") );
  for (int i=0;i<kNDet;i++) { LogThrowIf(!trees[i], Form("Missing tree for detector %d", i)); }

  // Branch holders
  std::vector<std::vector<float>*> data(kNDet);
  for (int i=0;i<kNDet;i++) data[i] = new std::vector<float>();
  trees[0]->SetBranchAddress("RAW Event",   &data[0]);
  trees[1]->SetBranchAddress("RAW Event B", &data[1]);
  trees[2]->SetBranchAddress("RAW Event C", &data[2]);

  Long64_t nEntries = trees[0]->GetEntries();
  LogThrowIf(nEntries<=0, "No entries found in raw_events tree");

  // Prepare sweep of sigma values
  std::vector<double> xs; std::vector<double> ys; xs.reserve(256); ys.reserve(256);
  for (double s = smin; s <= smax + 1e-9; s += sstep) {
    xs.push_back(s);
    ys.push_back(0.0);
  }

  // Count hits per sigma threshold
  for (Long64_t ievt=0; ievt<nEntries; ++ievt) {
    for (int det=0; det<kNDet; ++det) { data[det]->clear(); trees[det]->GetEntry(ievt); }
    for (int det=0; det<kNDet; ++det) {
      const auto &vec = *data[det];
      int n = std::min((int)vec.size(), kNChan);
      for (int ch=0; ch<n; ++ch) {
        if (mask[det][ch]) continue; // calib-bad or edge channel
        float amp = vec[ch] - baseline[det][ch];
        float sg  = sigma[det][ch];
        if (sg <= 0) continue;
        double sval = std::abs(amp) / sg;
        // For each sigma threshold bucket, if sval >= threshold, it counts as a hit
        for (size_t i=0; i<xs.size(); ++i) {
          if (sval >= xs[i]) ys[i] += 1.0; else break; // thresholds are increasing
        }
      }
    }
  }

  // Plot
  gStyle->SetOptStat(0);
  TCanvas c("c_hvS", "Hits vs Sigma", 900, 700);
  TGraph g(xs.size());
  c.SetLogy();
  for (size_t i=0; i<xs.size(); ++i) g.SetPoint(i, xs[i], ys[i]);
  if (runNumForTitle >= 0) {
    g.SetTitle(Form("Total hits vs Sigma (Run %d);Sigma threshold;Total hits", runNumForTitle));
  } else {
    g.SetTitle("Total hits vs Sigma;Sigma threshold;Total hits");
  }
  g.SetMarkerStyle(20); g.SetMarkerSize(1.0);
  g.Draw("ALP");
  c.SetGrid();
  c.SaveAs(outPdf.c_str());

  LogInfo << "Saved plot to: " << outPdf << std::endl;
  return 0;
}
