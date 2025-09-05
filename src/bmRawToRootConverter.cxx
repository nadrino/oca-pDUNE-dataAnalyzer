//
// Created by Nadrino on 04/07/2025.
//

#include "bmRawToRootConverter.h"
#include "PAPERO.h"

#include "CmdLineParser.h"
#include "Logger.h"

#include "TFile.h"
#include "TTree.h"

#include <filesystem>
#include <cstdlib>

#include "GenericToolbox.Root.h"
#include "GenericToolbox.Utils.h"


template <typename T>
std::vector<T> reorder_DUNE(std::vector<T> const &v)
{
  std::vector<T> reordered_vec(v.size());
  int j = 0;
  constexpr int order[] = {1, 0, 3, 2, 4, 8, 6, 5, 9, 7};
  for (int ch = 0; ch < 192; ch++)
  {
    for (int adc : order)
    {
      reordered_vec.at(adc * 192 + ch) = v.at(j);
      j++;
    }
  }
  return reordered_vec;
}



int main(int argc, char **argv){
  LogInfo << "Running raw to ROOT converter..." << std::endl;

  CmdLineParser parser = CmdLineParser();

  parser.addDummyOption("Required");
  parser.addOption("inputFile", {"-i"}, "Input .dat file to be converted.");

  parser.addDummyOption("Optional");
  parser.addOption("outputFolder", {"-o"}, "Output folder where the ROOT file will be writen.");
  parser.addOption("calibrationFile", {"-c"}, "Use calibration ROOT file (calib data processed by this app).");
  parser.addOption("threshold", {"-t"}, "Skip event with all channel peaks < baseline + threshold*std-dev (requires calib file)");
  parser.addOption("maxNbEvents", {"-me"}, "Stop after reading N events");
  parser.addTriggerOption("writeCalibData", {"--is-calib"}, "Will compute the baseline, the std-dev and create the calib TTree");
  parser.addTriggerOption("calcCovCalib", {"--cov"}, "Calculate covariance matrix of the std-dev");
  parser.addTriggerOption("skipEventTree", {"--skip-event"}, "Don't write the events TTree (useful if calib)");
  parser.addTriggerOption("enableIfBeamOutput", {"-if"}, "Also generate .txt files for filling up the IFBeam database");
  parser.addTriggerOption("verbose", {"-v"}, "Enable verbode");
  parser.addTriggerOption("enableZeroSuppr", {"-z"}, "Enable zero suppression");

  LogInfo << parser.getDescription().str() << std::endl;
  LogInfo << "Usage: " << std::endl;
  LogInfo << parser.getConfigSummary() << std::endl << std::endl;

  parser.parseCmdLine(argc, argv);

  LogThrowIf( parser.isNoOptionTriggered(), "No option was provided." );

  LogInfo << "Provided arguments: " << std::endl;
  LogInfo << parser.getValueSummary() << std::endl << std::endl;

  auto inputDatFilePath = parser.getOptionVal<std::string>("inputFile");
  std::string calibFilePath = {};
  if(parser.isOptionTriggered("calibrationFile")) {
    calibFilePath = parser.getOptionVal<std::string>("calibrationFile");
  }
  int nMaxEvents{-1};
  if(parser.isOptionTriggered("maxNbEvents") ){ nMaxEvents = parser.getOptionVal<int>("maxNbEvents"); }
  bool verbose = parser.isOptionTriggered("verbose");
  bool skipEventTree = parser.isOptionTriggered("skipEventTree");
  bool writeCalibData = parser.isOptionTriggered("writeCalibData");
  bool useCalibThreshold = parser.isOptionTriggered("threshold");
  bool calcCovCalib = parser.isOptionTriggered("calcCovCalib");
  bool enableZeroSuppr = parser.isOptionTriggered("enableZeroSuppr");
  double threshold{std::nan("unset")};
  if(useCalibThreshold) threshold = parser.getOptionVal<double>("threshold");

  LogExitIf(
    writeCalibData and not calibFilePath.empty(),
    "Can't write calib data if using external calib file."
    );

  std::string outFolderPath{};
  if( parser.isOptionTriggered("outputFolder") ){ outFolderPath = parser.getOptionVal<std::string>("outputFolder"); }

  std::filesystem::path outputRootFilePath(inputDatFilePath);
  outputRootFilePath.replace_extension(".root");
  outputRootFilePath = outFolderPath / outputRootFilePath.filename();
  LogWarning << inputDatFilePath << " -> " << outputRootFilePath << std::endl;

  // calib (write or read)
  double peakBaseline[N_DETECTORS][N_CHANNELS]{};
  double peakStdDev[N_DETECTORS][N_CHANNELS]{};

  std::unique_ptr<TFile> calibFile{nullptr};
  TTree *calibTree{nullptr};
  if( not calibFilePath.empty() ){
    LogInfo << "Reading calibration data from " << calibFilePath << std::endl;
    int detectorIdx; int channelIdx; double baseline; double stddev;
    calibFile = std::make_unique<TFile>(calibFilePath.c_str(), "READ");
    LogThrowIf(calibFile==nullptr, "Can't open calibration file.");
    calibTree = calibFile->Get<TTree>("calibration");
    LogThrowIf(calibTree==nullptr, "Can't open calibration tree.");

    calibTree->SetBranchAddress("detectorIdx", &detectorIdx);
    calibTree->SetBranchAddress("channelIdx", &channelIdx);
    calibTree->SetBranchAddress("baseline", &baseline);
    calibTree->SetBranchAddress("stddev", &stddev);

    int nCalibEntries = int(calibTree->GetEntries());
    LogThrowIf(nCalibEntries != N_DETECTORS*N_CHANNELS, "Wrong number of channels in calibration file.");
    for( int iEntry = 0; iEntry < nCalibEntries; ++iEntry ){
      calibTree->GetEntry(iEntry);
      peakBaseline[detectorIdx][channelIdx] = baseline;
      peakStdDev[detectorIdx][channelIdx] = stddev;
    }
  }

  auto inputDatFile = std::fstream(inputDatFilePath, std::ios::in | std::ios::binary);
  LogThrowIf(inputDatFile.fail());

  BeamMonitorEventBuffer bmEvent{};

  uint64_t offset = 0;
  uint64_t padding_offset = 0;
  uint64_t nEntries{0};

  LogInfo << "Checking the nb of entries..." << std::endl;
  offset = seek_first_evt_header(inputDatFile, 0, verbose);
  while( not inputDatFile.eof() ) {
    nEntries++;

    if( nMaxEvents != -1 and nEntries == nMaxEvents ){ break; }

    // check for event header if this is the first board
    if( not read_evt_header(inputDatFile, offset, verbose) ){ break; }

    bmEvent.readTuple( read_de10_header(inputDatFile, offset, verbose) );

    if( not bmEvent.isGood ){ continue; }

    offset = bmEvent.offset;
    read_event(inputDatFile, offset, int(bmEvent.eventSize), verbose, false);
    offset = (uint64_t) inputDatFile.tellg() + padding_offset + 8;

  }

  std::unique_ptr<TFile> outputRootFile = std::make_unique<TFile>(outputRootFilePath.c_str(), "RECREATE");

  if( not calibFilePath.empty() ) {
    outputRootFile->cd();
    GenericToolbox::writeInTFile(outputRootFile.get(), TNamed("calibFilePath", calibFilePath.c_str()));
    calibTree->CloneTree()->Write();
  }


  std::ofstream ifBeamOutput;
  size_t fileTimestamp;
  if( parser.isOptionTriggered("enableIfBeamOutput") ) {
    std::filesystem::path ifBeamOutputPath(outputRootFilePath);
    ifBeamOutputPath.replace_extension(".txt");
    LogInfo << "Opening output IFBeam output txt file: " << ifBeamOutputPath << std::endl;
    ifBeamOutput = std::ofstream(ifBeamOutputPath);
    LogThrowIf(ifBeamOutput.fail(), "Could not open IFBeam output file: " << ifBeamOutputPath);
  }

  // SCD_RUN00528_BEAM_20250904_065657.dat
  std::string dateTimeStr = GenericToolbox::splitString(inputDatFilePath, "/").back().substr(18, 15);
  DEBUG_VAR(dateTimeStr);
  std::tm tmbuff = {};
  std::stringstream ss(dateTimeStr);
  ss >> std::get_time(&tmbuff, "%Y%m%d_%H%M%S");
  std::time_t timestamp = std::mktime(&tmbuff);
  fileTimestamp = static_cast<size_t>(timestamp) * 1E9; // Convert to ns
  DEBUG_VAR(fileTimestamp);

  outputRootFile->cd();
  TTree* tree{nullptr};

  if( not skipEventTree ) {
    tree = new TTree("events", "events");

    // tree->Branch("size", &bmEvent.eventSize);
    // tree->Branch("fwVersion", &bmEvent.fwVersion);
    tree->Branch("triggerNumber", &bmEvent.triggerNumber);
    // tree->Branch("boardId", &bmEvent.boardId); // always board 0
    tree->Branch("timestamp", &bmEvent.timestamp);
    tree->Branch("timestampNs", &bmEvent.timestampNs);
    tree->Branch("timestampUtcNs", &bmEvent.timestampUtcNs);
    tree->Branch("extTimestamp", &bmEvent.extTimestamp);
    tree->Branch("triggerId", &bmEvent.triggerId);
    tree->Branch("peakAdc", &bmEvent.peakAdc, Form("peakAdc[%d][%d]/i", N_DETECTORS, N_CHANNELS));
    tree->Branch("deltaTimeNsLastEvent", &bmEvent.deltaTimeNsLastEvent);
    // tree->Branch("peakAdcSum", &bmEvent.peakAdcSum, Form("peakAdcSum[%d]/i", N_DETECTORS));

    if( not calibFilePath.empty() ){
      tree->Branch("peak", &bmEvent.peak, Form("peak[%d][%d]/D", N_DETECTORS, N_CHANNELS));
      tree->Branch("peakSum", &bmEvent.peakSum, Form("peakSum[%d]/D", N_DETECTORS));
      if( useCalibThreshold ) {
        tree->Branch("peakZeroSuppr", &bmEvent.peakZeroSuppr, Form("peakZeroSuppr[%d][%d]/D", N_DETECTORS, N_CHANNELS));
        tree->Branch("peakZeroSupprSum", &bmEvent.peakZeroSupprSum, Form("peakZeroSupprSum[%d]/D", N_DETECTORS));
        tree->Branch("deltaTimeNsLastTrigEvent", &bmEvent.deltaTimeNsLastTrigEvent);
        tree->Branch("nClusters", &bmEvent.nClusters, Form("nClusters[%d]/i", N_DETECTORS));
        tree->Branch("xBarycenter", &bmEvent.xBarycenter, Form("xBarycenter[%d]/D", N_DETECTORS));
        tree->Branch("yBarycenter", &bmEvent.yBarycenter, Form("yBarycenter/D"));
      }
    }
  }

  std::unique_ptr<TMatrixD> covMatrix{nullptr};
  if( calcCovCalib ){
    covMatrix = std::make_unique<TMatrixD>(
      N_DETECTORS*N_CHANNELS,
      N_DETECTORS*N_CHANNELS
      );
    for(int iFlat = 0; iFlat < N_DETECTORS*N_CHANNELS; ++iFlat ) {
      for(int jFlat = 0; jFlat < N_DETECTORS*N_CHANNELS; ++jFlat ) {
        (*covMatrix)[iFlat][jFlat] = 0;
      }
    }
  }

  LogInfo << "Reading " << nEntries << " entries..." << std::endl;
  inputDatFile = std::fstream(inputDatFilePath, std::ios::in | std::ios::binary);
  offset = seek_first_evt_header(inputDatFile, 0, verbose);
  auto iEvent{nEntries}; iEvent = 0;
  auto nWriten{iEvent};
  while( not inputDatFile.eof() ) {
    iEvent++;
    GenericToolbox::displayProgressBar(iEvent, nEntries, "Writing events...");

    if( iEvent == nEntries ){ break; }

    // check for event header if this is the first board
    if( not read_evt_header(inputDatFile, offset, verbose) ){ break; }

    // reading header
    bmEvent.readTuple( read_de10_header(inputDatFile, offset, verbose) );

    if( not bmEvent.isGood ){ continue; }

    bmEvent.timestampNs = bmEvent.timestamp * 20;
    // DEBUG_VAR(fileTimestamp);
    // DEBUG_VAR(bmEvent.timestampNs);
    // DEBUG_VAR(bmEvent.timestampUtc);
    bmEvent.timestampUtcNs = bmEvent.timestampNs + fileTimestamp;
    if(bmEvent.lastTimestampNs != 0) {
      bmEvent.deltaTimeNsLastEvent = long(bmEvent.timestampNs) - long(bmEvent.lastTimestampNs);
    }
    // DEBUG_VAR(bmEvent.timestamp);
    bmEvent.lastTimestampNs = bmEvent.timestampNs;

    offset = bmEvent.offset;
    auto data = read_event(inputDatFile, offset, int(bmEvent.eventSize), verbose, false);
    data = reorder_DUNE(data);

    // should be the amount of channel
    size_t nbOfValuesPerDet = 2 * data.size() / ADC_N;
    LogThrowIf(nbOfValuesPerDet - N_CHANNELS != 0, "Invalid data size: " << nbOfValuesPerDet - N_CHANNELS);

    bool skip{true}; // if zeroSuppress and at no signal is over the threshold
    bmEvent.yBarycenter=0;
    for (size_t iDet = 0; iDet < N_DETECTORS; ++iDet) {
      memcpy(&bmEvent.peakAdc[iDet][0], &data[iDet * N_CHANNELS], N_CHANNELS * sizeof(unsigned int));
      bmEvent.peakAdcSum[iDet] = std::accumulate(
              &bmEvent.peakAdc[iDet][0],
              &bmEvent.peakAdc[iDet][N_CHANNELS],
              static_cast<uint32_t>(0)
            );

      if( not calibFilePath.empty() ) {
        bmEvent.xBarycenter[iDet] = 0;
        for( size_t iCh = 0; iCh < N_CHANNELS; ++iCh ) {
          if( iCh%64 == 0 or iCh%64 == 63 or iCh%64 == 1 or iCh%64 == 62 ){ continue; }
          bmEvent.peak[iDet][iCh] = static_cast<double>(bmEvent.peakAdc[iDet][iCh]) - peakBaseline[iDet][iCh];

          if( useCalibThreshold and bmEvent.peak[iDet][iCh] >= peakStdDev[iDet][iCh]*threshold ) {
            bmEvent.peakZeroSuppr[iDet][iCh] = bmEvent.peak[iDet][iCh];
            bmEvent.xBarycenter[iDet] += double(iCh)*bmEvent.peakZeroSuppr[iDet][iCh];
            skip = false;
          }
        }

        bmEvent.peakSum[iDet] = std::accumulate(
              &bmEvent.peak[iDet][0],
              &bmEvent.peak[iDet][N_CHANNELS],
              0.0
            );
        if( useCalibThreshold ) {
          bmEvent.nClusters[iDet] = 0;

          bool isLastChOn = false;
          for( int iCh = 0; iCh < N_CHANNELS; ++iCh ) {
            if( bmEvent.peakZeroSuppr[iDet][iCh] > 0. ) {
              // count once
              if( not isLastChOn ) { bmEvent.nClusters[iDet]++; }
              isLastChOn = true;
            }
            else{ isLastChOn = false; }
          }

          bmEvent.peakZeroSupprSum[iDet] = std::accumulate(
              &bmEvent.peakZeroSuppr[iDet][0],
              &bmEvent.peakZeroSuppr[iDet][N_CHANNELS],
              0.0
            );

          bmEvent.xBarycenter[iDet] /= bmEvent.peakZeroSupprSum[iDet];
          bmEvent.xBarycenter[iDet] -= double(N_CHANNELS)/2.;

          // bmEvent.yBarycenter += std::sin(double(iDet)*15.*M_PI/180.)*bmEvent.xBarycenter[iDet]/2.; // 2 detectors will give this info
          if( iDet == 2 ){ bmEvent.yBarycenter += std::sin(15.*M_PI/180.)*bmEvent.xBarycenter[iDet]/2.; }
          if( iDet == 1 ){ bmEvent.yBarycenter += std::sin(30.*M_PI/180.)*bmEvent.xBarycenter[iDet]/2.; }
        }
      }
    }
    if( useCalibThreshold and enableZeroSuppr and skip ){ continue; } // skip TTree::Fill();

    if( useCalibThreshold ){
      if(bmEvent.lastTriggeredTimestampNs != 0) {
        bmEvent.deltaTimeNsLastTrigEvent = long(bmEvent.timestampNs) - long(bmEvent.lastTriggeredTimestampNs);
      }
      bmEvent.lastTriggeredTimestampNs = bmEvent.timestampNs;
    }

    if( ifBeamOutput.is_open() and not (useCalibThreshold and skip) ) {
      ifBeamOutput << "z,pdune " << bmEvent.timestampUtcNs << " " << 0 << std::endl;
      // <device_name>\t<timestamp_in_ms>\t<optional_unit_name>\t<scalar_value|string|null>\t<array_value|null>
      ifBeamOutput << "dip/acc/NORTH/NP02/BPM/timestampNs\t" << bmEvent.timestampUtcNs << "\tnull\t" << bmEvent.timestampNs << "\tnull" << std::endl;
      ifBeamOutput << "dip/acc/NORTH/NP02/BPM/extTimestamp\t" << bmEvent.timestampUtcNs << "\tnull\t" << bmEvent.extTimestamp << "\tnull" << std::endl;
      ifBeamOutput << "dip/acc/NORTH/NP02/BPM/triggerNumber\t" << bmEvent.timestampUtcNs << "\tnull\t" << bmEvent.triggerNumber << "\tnull" << std::endl;
      ifBeamOutput << "dip/acc/NORTH/NP02/BPM/nClusters[]\t" << bmEvent.timestampUtcNs << "\tnull\tnull\t{" << bmEvent.nClusters[0] << "," << bmEvent.nClusters[1] << "," << bmEvent.nClusters[2] << "}" << std::endl;
      ifBeamOutput << "dip/acc/NORTH/NP02/BPM/xBarycenter[]\t" << bmEvent.timestampUtcNs << "\tnull\tnull\t{" << bmEvent.xBarycenter[0] << "," << bmEvent.xBarycenter[1] << "," << bmEvent.xBarycenter[2] << "}" << std::endl;
      ifBeamOutput << "dip/acc/NORTH/NP02/BPM/yBarycenter\t" << bmEvent.timestampUtcNs << "\tnull\t" << bmEvent.yBarycenter << "\tnull" << std::endl;
    }

    if( not skipEventTree ){ tree->Fill(); nWriten++; }

    if( writeCalibData ){
      for( int iDet = 0 ; iDet < N_DETECTORS ; iDet++ ) {
        for( int iCh = 0 ; iCh < N_CHANNELS ; iCh++ ) {
          auto val_i = static_cast<double>(bmEvent.peakAdc[iDet][iCh]);
          peakBaseline[iDet][iCh] += val_i;
          peakStdDev[iDet][iCh] += val_i * val_i;

          if( covMatrix != nullptr ) {
            int iFlat = iDet*N_CHANNELS+iCh;
            for (int jFlat = iFlat; jFlat < N_DETECTORS * N_CHANNELS; ++jFlat) {
              auto val_j = static_cast<double>(bmEvent.peakAdc[jFlat / N_CHANNELS][jFlat % N_CHANNELS]);
              (*covMatrix)[iFlat][jFlat] += val_i * val_j;
            }
          }
        }
      }
    }

    // next offset
    offset = static_cast<uint64_t>(inputDatFile.tellg()) + padding_offset + 8;

  }
  if( not skipEventTree ) {
    LogInfo << nWriten << " events have been writen." << std::endl;
    tree->Write(tree->GetName(), TObject::kOverwrite);
  }

  if( writeCalibData ){
    LogInfo << "Writing calibration data..." << std::endl;
    for (int iDet = 0; iDet < N_DETECTORS; ++iDet) {
      for (int iCh = 0; iCh < N_CHANNELS; ++iCh) {
        peakBaseline[iDet][iCh] /= double(nEntries);
        peakStdDev[iDet][iCh] = std::sqrt(
          peakStdDev[iDet][iCh] / double(nEntries) - peakBaseline[iDet][iCh] * peakBaseline[iDet][iCh]
        );
      }
    }

    if( covMatrix != nullptr ) {
      LogInfo << "Writing correlation matrix..." << std::endl;
      for (int i = 0; i < N_DETECTORS * N_CHANNELS; ++i) {
        for (int j = i; j < N_DETECTORS * N_CHANNELS; ++j) {
          double mean_i = peakBaseline[i / N_CHANNELS][i % N_CHANNELS];
          double mean_j = peakBaseline[j / N_CHANNELS][j % N_CHANNELS];
          (*covMatrix)[i][j] = ((*covMatrix)[i][j] / double(nEntries)) - mean_i * mean_j;
          (*covMatrix)[j][i] = (*covMatrix)[i][j];
        }
      }

      outputRootFile->cd();
      auto* corr = GenericToolbox::convertToCorrelationMatrix(covMatrix.get());
      auto* corrHist = GenericToolbox::convertTMatrixDtoTH2D(
        corr,
        "Calibration covariance matrix",
        "correlation",
        "Channel #", "Channel #");

      corrHist->SetMinimum(-1);
      corrHist->SetMaximum(1);
      corrHist->SetDrawOption("COLZ");
      GenericToolbox::fixTH2display(corrHist);
      corrHist->Write("correlationMatrix");
      delete corrHist;
    }

    outputRootFile->cd();
    auto* outCalibTree = new TTree("calibration", "calibration");
    int detectorIdx;
    int channelIdx;
    double baseline;
    double stddev;
    outCalibTree->Branch("detectorIdx", &detectorIdx);
    outCalibTree->Branch("channelIdx", &channelIdx);
    outCalibTree->Branch("baseline", &baseline);
    outCalibTree->Branch("stddev", &stddev);

    for (detectorIdx = 0; detectorIdx < N_DETECTORS; ++detectorIdx) {
      for (channelIdx = 0; channelIdx < N_CHANNELS; ++channelIdx) {
        baseline = peakBaseline[detectorIdx][channelIdx];
        stddev = std::sqrt(peakStdDev[detectorIdx][channelIdx]);
        outCalibTree->Fill();
      }
    }
    outCalibTree->Write(outCalibTree->GetName(), TObject::kOverwrite);
  }

  LogInfo << "Written " << outputRootFile->GetPath() << std::endl;
  return EXIT_SUCCESS;
}
