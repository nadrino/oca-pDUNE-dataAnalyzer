///////////////////////////////////////
// Created by E. Villa on 2024-08-27.//
// Code to look at the data of the   //
// NP02 beam monitor.                //
///////////////////////////////////////

// #include <iostream>
// #include <string>
// #include <vector>
// #include <fstream>
// #include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TSystem.h"
#include <cmath>
#include <algorithm>

#include <nlohmann/json.hpp>

#include "CmdLineParser.h"
#include "Logger.h"
#include "ocaEvent.h"

LoggerInit([]{
  Logger::getUserHeader() << "[" << FILENAME << "]";
});

int main(int argc, char* argv[]) {

    CmdLineParser clp;

    clp.getDescription() << "> This program takes a root file and a calibration file and analyzes it." << std::endl;

    clp.addDummyOption("Main options");
    clp.addOption("runNumber", {"-n", "--run-number"}, "Specify run number.");
    clp.addOption("appSettings",    {"-j", "--json-settings"},   "Specify application settings file path.");
    clp.addOption("inputRootFile",  {"-r", "--root-file"},      "Root converted data");
    clp.addOption("inputCalFile",   {"-c", "--cal-file"},       "Calibration file.");
    clp.addOption("outputDir",      {"-o", "--output"},         "Specify output directory path");
    clp.addOption("nSigma",         {"-s", "--n-sigma"},        "Number of sigmas above pedestal to consider signal");

    clp.addDummyOption("Triggers");
    clp.addTriggerOption("verboseMode",     {"-v"},             "RunVerboseMode, bool");
    clp.addTriggerOption("debugMode",       {"-d"},             "RunDebugMode, bool");

    clp.addDummyOption();

    // usage always displayed
    LogInfo << clp.getDescription().str() << std::endl;

    LogInfo << "Usage: " << std::endl;
    LogInfo << clp.getConfigSummary() << std::endl << std::endl;

    clp.parseCmdLine(argc, argv);

    LogThrowIf( clp.isNoOptionTriggered(), "No option was provided." );

    LogInfo << "Provided arguments: " << std::endl;
    LogInfo << clp.getValueSummary() << std::endl << std::endl;

    bool verbose = clp.isOptionTriggered("verboseMode");
    bool debug = clp.isOptionTriggered("debugMode");
    
    // read json settings
    std::string jsonSettingsFile = "";
    std::ifstream i(jsonSettingsFile);
    nlohmann::json jsonSettings;
    if (clp.isOptionTriggered("appSettings")) {
        jsonSettingsFile = clp.getOptionVal<std::string>("appSettings");
        LogInfo << "Reading JSON settings from: " << jsonSettingsFile << std::endl;
        i.open(jsonSettingsFile);
        if (i.is_open()) {
            i >> jsonSettings;
            LogInfo << "JSON settings loaded successfully." << std::endl;
        } else {
            LogError << "Failed to read JSON settings from: " << jsonSettingsFile << std::endl;
            return 1;
        }
    } else {
        LogInfo << "No JSON settings file provided, using default values." << std::endl;
    }

    ///////////////////////////
    // Some variables

    int nDetectors = 4;
    int nChannels = 384;

    ///////////////////////////

    // get calibration file
    std::string inputCalFile = clp.getOptionVal<std::string>("inputCalFile");
    LogInfo << "Calibration file: " << inputCalFile << std::endl;
    std::ifstream calFile(inputCalFile);
    // check if the file is open correctly
    if (!calFile.is_open()) {
        std::cout << "Error: calibration file not open" << std::endl;
        return 1;
    }

    // the format is // TODO move to some docs
    // 0: channel
    // 1: channel / 64      // I think this is the chip
    // 2: va_chan           // I think this is the channel on the chuip
    // 3: pedestals->at(ch)
    // 4: rsigma->at(ch)
    // 5: sigma_value
    // 6: badchan
    // 7: 0.000

    /// Calib file
    // read the calibration file and store the values in a vector
    std::vector <std::vector <float>> baseline (nDetectors, std::vector<float>(nChannels));
    
    std::vector<std::vector<float>> baseline_sigma(nDetectors, std::vector<float>(nChannels));

    // skip header lines, starting with #

    LogInfo << "Reading calibration file..." << std::endl;
    std::string line;
    for (int detit = 0; detit < nDetectors ; detit++){
        // skip the first 18 lines, header
        for (int i = 0; i < 18; i++) std::getline(calFile, line); // skipping, header

        for (int i = 0; i < nChannels; i++) {
            std::getline(calFile, line);
            // values are separated by commas, read all and store line by line
            std::istringstream iss(line);
            std::vector <float> values;
            float value;
            std::string token;
            while (std::getline(iss, token, ',')) {
                std::istringstream iss_value(token);
                float value;
                if (iss_value >> value) {
                values.push_back(value);
                }
            }

            // check if the values are correct
            if (values.size() != 8) {
                LogError << "Error: wrong number of values in the calibration file" << std::endl;
                LogError << "Values size: " << values.size() << std::endl;
                return 1;
            }

            // store the values
            int this_channel = values.at(0);
            float this_baseline = values.at(3);
            float this_baseline_sigma = values.at(5);

            if (verbose) LogInfo << "Channel: " << this_channel << " Detector: " << detit << " Baseline: " << this_baseline << " Baseline sigma: " << this_baseline_sigma << std::endl;

            baseline.at(detit).emplace(baseline.at(detit).begin() + this_channel, this_baseline);
            baseline_sigma.at(detit).emplace(baseline_sigma.at(detit).begin() + this_channel, this_baseline_sigma);

        }
    }

    LogInfo << "Stored baseline and sigma for all channels locally" << std::endl;

    /// ROOT file

    // Get root file name
    std::string input_root_filename = clp.getOptionVal<std::string>("inputRootFile");
    LogInfo << "Root file: " << input_root_filename << std::endl;   
    TFile *input_root_file = new TFile(input_root_filename.c_str(), "READ");
    // check if the file is open correctly
    if( !input_root_file->IsOpen() ){
        LogError << "Error: file not open" << std::endl;
        return 1;
    }

    // get trees
    // TODO change at converter level, I don't like the names
    std::vector<TTree*> raw_events_trees = std::vector <TTree*>();
    raw_events_trees.reserve(nDetectors);
    TTree *these_raw_events = (TTree*)input_root_file->Get("raw_events");
    raw_events_trees.emplace_back(these_raw_events);
    these_raw_events = (TTree*)input_root_file->Get("raw_events_B");
    raw_events_trees.emplace_back(these_raw_events);
    these_raw_events = (TTree*)input_root_file->Get("raw_events_C");
    raw_events_trees.emplace_back(these_raw_events);
    these_raw_events = (TTree*)input_root_file->Get("raw_events_D");
    raw_events_trees.emplace_back(these_raw_events);

    LogInfo << "Got the trees" << std::endl;

    // get the number of entries
    std::vector <int> nEntries = std::vector <int>();
    nEntries.reserve(nDetectors);
    for (int detit = 0; detit < nDetectors; detit++) {
        nEntries.emplace_back(raw_events_trees.at(detit)->GetEntries());
        LogInfo << "Detector " << detit << " has " << nEntries.at(detit) << " entries" << std::endl;
    }

    // Entries should always be the same for all detectors
    // if (nEntries.at(0) != nEntries.at(1) || nEntries.at(0) != nEntries.at(2) || nEntries.at(0) != nEntries.at(3)) {
    //     LogError << "Error: number of entries is different for the detectors! Something went wrong" << std::endl;
    //     return 1;
    // }

    ///////////////////////////
    
    /// Create some objects to plot results

    // Root app
    TApplication *app = new TApplication("app", &argc, argv);

    // set root to displat overflow and underflow in stat box
    gStyle->SetOptStat("emruo"); // e: entries, m:

    // Create a vector of TF1 objects to show the channels that fire, one for each detector
    std::vector <TH1F*> *h_firingChannels = new std::vector <TH1F*>;
    h_firingChannels->reserve(nDetectors);
    for (int i = 0; i < nDetectors; i++) {
        TH1F *this_h_firingChannels = new TH1F(Form("Firing channels (Detector %d)", i), Form("Firing channels (Detector %d)", i), nChannels, 0, nChannels);
        this_h_firingChannels->GetXaxis()->SetTitle("Channel");
        this_h_firingChannels->GetYaxis()->SetTitle("Counts");
        h_firingChannels->emplace_back(this_h_firingChannels);
    }

    // plot values of sigma per channel in a tgraph
    std::vector <TGraph*> *g_sigma = new std::vector <TGraph*>;
    g_sigma->reserve(nDetectors);
    for (int i = 0; i < nDetectors; i++) {
        TGraph *this_g_sigma = new TGraph(nChannels);
        this_g_sigma->SetTitle(Form("Sigma (Detector %d)", i));
        this_g_sigma->GetXaxis()->SetTitle("Channel");
        this_g_sigma->GetYaxis()->SetTitle("Sigma");
        this_g_sigma->SetMarkerStyle(20);
        this_g_sigma->SetMarkerSize(0.8);
        g_sigma->emplace_back(this_g_sigma);
    }

    // plot baseline
    std::vector <TGraph*> *g_baseline = new std::vector <TGraph*>;
    g_baseline->reserve(nDetectors);
    for (int i = 0; i < nDetectors; i++) {
        TGraph *this_g_baseline = new TGraph(nChannels);
        this_g_baseline->SetTitle(Form("Baseline (Detector %d)", i));
        this_g_baseline->GetXaxis()->SetTitle("Channel");
        this_g_baseline->GetYaxis()->SetTitle("Baseline");
        this_g_baseline->SetMarkerStyle(20);
        this_g_baseline->SetMarkerSize(0.8);
        g_baseline->emplace_back(this_g_baseline);
    }

    // raw peak for each channel
    std::vector <std::vector <TH1F*>*> *h_rawPeak = new std::vector <std::vector <TH1F*>*>;
    h_rawPeak->reserve(nDetectors);
    for (int i = 0; i < nDetectors; i++) {
        LogInfo << "Detector " << i << std::endl;
        std::vector  <TH1F*> *this_h_rawPeak_vector = new std::vector <TH1F*>;
        for (int j = 0; j < nChannels; j++) {
            // LogInfo << "Channel " << j << std::endl;
            TH1F *this_h_rawPeak = new TH1F(Form("Raw peak (Detector %d, Channel %d)", i, j), Form("Raw peak (Detector %d, Channel %d)", i, j), 1000, 0, 2000);
            this_h_rawPeak->GetXaxis()->SetTitle("Peak");
            this_h_rawPeak->GetYaxis()->SetTitle("Counts");
            this_h_rawPeak->GetYaxis()->SetRangeUser(0,5);
            // fill color blue
            this_h_rawPeak->SetFillColor(kBlue);
            this_h_rawPeak_vector->emplace_back(this_h_rawPeak);
        }
        h_rawPeak->emplace_back(this_h_rawPeak_vector);
    }

    std::vector <TH1F*> *h_amplitude = new std::vector <TH1F*>;
    h_amplitude->reserve(nDetectors);
    for (int i = 0; i < nDetectors; i++) {
    // Reduce amplitude range to 0–200 as requested
    TH1F *this_h_amplitude = new TH1F(Form("Amplitude (Detector %d)", i), Form("Amplitude (Detector %d)", i), 100, 0, 200);
        this_h_amplitude->GetXaxis()->SetTitle("Amplitude");
        this_h_amplitude->GetYaxis()->SetTitle("Counts");
        h_amplitude->emplace_back(this_h_amplitude);
    }

    // Add 2D histogram for amplitude vs channel
    std::vector <TH2F*> *h_amplitudeVsChannel = new std::vector <TH2F*>;
    h_amplitudeVsChannel->reserve(nDetectors);
    for (int i = 0; i < nDetectors; i++) {
        TH2F *this_h_amplitudeVsChannel = new TH2F(Form("Amplitude vs Channel (Detector %d)", i), 
                                                   Form("Amplitude vs Channel (Detector %d)", i), 
                                                   nChannels, 0, nChannels, 100, 0, 1000);
        this_h_amplitudeVsChannel->GetXaxis()->SetTitle("Channel");
        this_h_amplitudeVsChannel->GetYaxis()->SetTitle("Amplitude");
        h_amplitudeVsChannel->emplace_back(this_h_amplitudeVsChannel);
    }

    // Extend hits-per-event range to cover up to 20 hits
    TH1F *h_hitsInEvent = new TH1F("Hits in event", "Hits in event", 21, -0.5, 20.5);
    h_hitsInEvent->GetXaxis()->SetTitle("Hits");
    h_hitsInEvent->GetYaxis()->SetTitle("Counts");

    // loop over the entries to get the peak. For each entry, store the values in the event
    // TODO temporary, do the real thing

    LogInfo << "Reading the entries and storing in a vector of Events" << std::endl;
    std::vector <Event> * events = new std::vector <Event>();
    events->reserve(nEntries.at(0)); // should all be the same

    LogInfo << "Size of events: " << events->size() << std::endl;
    std::vector<std::vector<float>*>* data = new std::vector<std::vector<float>*>();
    data->reserve(nDetectors);
    for (int i = 0; i < nDetectors; i++) {
        std::vector<float> *this_data = new std::vector<float>;
        this_data->reserve(nChannels);
        data->emplace_back(this_data);
    }

    // TODO consider changing branch names?
    raw_events_trees.at(0)->SetBranchAddress("RAW Event", &data->at(0));
    raw_events_trees.at(1)->SetBranchAddress("RAW Event B", &data->at(1));
    raw_events_trees.at(2)->SetBranchAddress("RAW Event C", &data->at(2));
    raw_events_trees.at(3)->SetBranchAddress("RAW Event D", &data->at(3));
    
    int limit = nEntries.at(0);
    int maxEvents = 1e9; // default value
    // Overwrite maxEvents if present in json settings file
    if (jsonSettings.contains("maxEvents")) {
        // Check if the value is a number and handle it properly
        if (jsonSettings["maxEvents"].is_number_integer()) {
            maxEvents = jsonSettings["maxEvents"];
            LogInfo << "Overwriting maxEvents from JSON: " << maxEvents << std::endl;
        } else {
            LogError << "maxEvents in JSON is not a valid integer, using default: " << maxEvents << std::endl;
        }
    } else {
        LogInfo << "Using default maxEvents: " << maxEvents << std::endl;
    }
    if (limit > maxEvents) limit = maxEvents;

    // Check if edge channels should be masked
    auto parseBool = [](const nlohmann::json& j, const std::string& key, bool defVal)->bool{
        if (!j.contains(key)) return defVal;
        const auto& v = j.at(key);
        if (v.is_boolean()) return v.get<bool>();
        if (v.is_string()) {
            std::string s = v.get<std::string>();
            std::transform(s.begin(), s.end(), s.begin(), ::tolower);
            return (s == "true" || s == "1" || s == "yes" || s == "y");
        }
        if (v.is_number_integer()) return v.get<int>() != 0;
        return defVal;
    };

    bool maskEdgeChannels = parseBool(jsonSettings, "maskEdgeChannels", false);
    LogInfo << "Mask edge channels: " << (maskEdgeChannels ? "true" : "false") << std::endl;

    // Function to check if a channel is an edge channel (chip boundaries)
    auto isEdgeChannel = [](int channel) -> bool {
        // Edge channels are at chip boundaries: 0, 63, 64, 127, 128, 191, 192, 255, 256, 319, 320, 383
        // Each chip has 64 channels, so edge channels are: 0, 63 of each chip
        int chipNumber = channel / 64;
        int channelInChip = channel % 64;
        return (channelInChip == 0 || channelInChip == 1 || channelInChip == 62 || channelInChip == 63); // it seems also the ones next are noisy
    };

    int hitsInEvent = 0;
    int triggeredEvents = 0;
    
    // Geometry / interpolation setup for detectors 0 (-15°), 1 (+15°), 2 (0°)
    const int geomDetN = 3;
    const double detAnglesDeg[geomDetN] = {-15.0, +15.0, 0.0};
    double detAnglesRad[geomDetN];
    for (int i=0;i<geomDetN;i++) detAnglesRad[i] = detAnglesDeg[i]*M_PI/180.0;
    // Physical pitch: 10 cm across strip-normal over nChannels
    double channelPitch = 100.0 / nChannels; // mm per strip
    
    // Histograms to accumulate reconstructed centers (in mm)
    // Further increase bin size (coarser bins): 120 mm range with 40 bins => 3 mm/bin
    TH2F *h_recoCenter_3 = new TH2F("h_recoCenter_3", "Interpolated centers (exactly 3 clusters, 1 per detector);X [mm];Y [mm]", 40, -60.0, 60.0, 40, -60.0, 60.0);
    TH2F *h_recoCenter_2to3 = new TH2F("h_recoCenter_2to3", "Interpolated centers (2 or 3 clusters);X [mm];Y [mm]", 40, -60.0, 60.0, 40, -60.0, 60.0);
    TGraph *g_recoCenters_3 = new TGraph();
    g_recoCenters_3->SetName("g_recoCenters_3");
    g_recoCenters_3->SetTitle("Interpolated centers (exactly 3 clusters, 1 per detector);X [mm];Y [mm]");
    g_recoCenters_3->SetMarkerStyle(20);
    g_recoCenters_3->SetMarkerSize(0.5);
    TGraph *g_recoCenters_2to3 = new TGraph();
    g_recoCenters_2to3->SetName("g_recoCenters_2to3");
    g_recoCenters_2to3->SetTitle("Interpolated centers (2 or 3 clusters);X [mm];Y [mm]");
    g_recoCenters_2to3->SetMarkerStyle(20);
    g_recoCenters_2to3->SetMarkerSize(0.5);

    for (int entryit = 0; entryit < limit; entryit++) {

        hitsInEvent = 0;
        
        Event * this_event; // across detectors
        this_event = new Event();

        if (entryit % 10000 == 0 ) LogInfo << "Entry " << entryit << std::endl;

        this_event->SetBaseline(baseline);
        this_event->SetSigma(baseline_sigma); 
        this_event->SetNSigma(clp.getOptionVal<int>("nSigma"));

        // clear data
        for (int detit = 0; detit < nDetectors; detit++)   data->at(detit)->clear();

        for (int detit = 0; detit < nDetectors; detit++) {
            raw_events_trees.at(detit)->GetEntry(entryit);
            this_event->AddPeak(detit, *data->at(detit));
            for (int chit = 0; chit < nChannels; chit++) {
                // Skip edge channels if masking is enabled
                if (maskEdgeChannels && isEdgeChannel(chit)) {
                    continue;
                }
                // LogInfo << "DetId " << detit << ", channel " << chit << ", peak: " << this_event->GetPeak(detit, chit) << ", baseline: " << this_event->GetBaseline(detit, chit) << ", sigma: " << this_event->GetSigma(detit, chit) << "\t";
                h_rawPeak->at(detit)->at(chit)->Fill(this_event->GetPeak(detit, chit));
            }
        }

        this_event->ExtractTriggeredHits();

    if (verbose) this_event->PrintOverview();        
        if (debug) this_event->PrintInfo(); // this should rather be debug

        // Apply event-level filter for 2D centers: accept only events with exactly
        // one cluster per detector (0,1,2). A cluster is one or more hits in contiguous
        // channels within the same detector (edge-masked if enabled).
    std::vector <std::pair<int, int>> triggeredHits = this_event->GetTriggeredHits();
    int eventHitsMasked = 0;
    if (!triggeredHits.empty()){
            // Always fill per-hit histograms for all triggered (masked) hits
            for (auto &hit : triggeredHits) {
                int det = hit.first;
                int ch  = hit.second;
                if (det < 0 || det >= nDetectors) continue;
                if (maskEdgeChannels && isEdgeChannel(ch)) continue;
        eventHitsMasked++;
                float amplitude = this_event->GetPeak(det, ch) - this_event->GetBaseline(det, ch);
                h_firingChannels->at(det)->Fill(ch);
                h_amplitude->at(det)->Fill(amplitude);
                h_amplitudeVsChannel->at(det)->Fill(ch, amplitude);
            }

            // Gather masked hits per detector
            std::vector<int> detHits[4];
            for (auto &hit : triggeredHits) {
                int det = hit.first;
                int ch  = hit.second;
                if (det < 0 || det >= 4) continue;
                if (maskEdgeChannels && isEdgeChannel(ch)) continue;
                detHits[det].push_back(ch);
            }

            // Compute clusters per detector (contiguous runs)
            auto countClusters = [](std::vector<int>& chans)->std::pair<int, std::pair<int,int>>{
                if (chans.empty()) return {0, { -1, -1 }};
                std::sort(chans.begin(), chans.end());
                chans.erase(std::unique(chans.begin(), chans.end()), chans.end());
                int clusters = 0;
                int runStart = chans.front();
                int prev = chans.front();
                int cMin = runStart, cMax = prev; // track the only cluster's min/max when clusters==1
                for (size_t i=1; i<chans.size(); ++i) {
                    if (chans[i] != prev + 1) {
                        clusters++;
                        if (clusters == 1) { cMin = runStart; cMax = prev; }
                        runStart = chans[i];
                    }
                    prev = chans[i];
                }
                // end last run
                clusters++;
                if (clusters == 1) { cMin = runStart; cMax = prev; }
                return {clusters, {cMin, cMax}};
            };

            auto c0 = countClusters(detHits[0]);
            auto c1 = countClusters(detHits[1]);
            auto c2 = countClusters(detHits[2]);
            // Accept events with exactly 2 or 3 clusters total, each on different detectors
            bool perDetOk = (c0.first <= 1 && c1.first <= 1 && c2.first <= 1);
            int totalClusters = (c0.first > 0) + (c1.first > 0) + (c2.first > 0);
            bool acceptTwoOrThree = perDetOk && (totalClusters == 2 || totalClusters == 3);
            bool acceptExactlyThree = (c0.first == 1 && c1.first == 1 && c2.first == 1);
            if (acceptTwoOrThree) {
                // Use cluster center channel per detector for (x,y) interpolation
                int firstTrigChan[geomDetN];
                firstTrigChan[0] = (c0.second.first >= 0 ? (c0.second.first + c0.second.second)/2 : -1);
                firstTrigChan[1] = (c1.second.first >= 0 ? (c1.second.first + c1.second.second)/2 : -1);
                firstTrigChan[2] = (c2.second.first >= 0 ? (c2.second.first + c2.second.second)/2 : -1);
                triggeredEvents++;
                hitsInEvent = (int)(detHits[0].size() + detHits[1].size() + detHits[2].size());
            // Reconstruct (x,y) from projections u_i = x cosθ_i + y sinθ_i, with u_i centered around detector origin.
            // Center channel indices so that strip centers map to u = 0 at detector center.
            double centerOffset = (nChannels * channelPitch) / 2.0; // 50 mm
            // Accumulate for least squares if >=2 detectors
            double Scc=0, Sss=0, Scs=0, Su_c=0, Su_s=0; int used=0;
            // Temporary storage for exactly-2-detector direct solve
            double c_a=0,s_a=0,u_a=0,c_b=0,s_b=0,u_b=0; int pairCount=0;
            for (int i=0;i<geomDetN;i++) {
                if (firstTrigChan[i] >= 0) {
                    double u = ( (firstTrigChan[i] + 0.5) * channelPitch ) - centerOffset; // mm
                    double c = cos(detAnglesRad[i]);
                    double s = sin(detAnglesRad[i]);
                    Scc += c*c; Sss += s*s; Scs += c*s; Su_c += u*c; Su_s += u*s; used++;
                    if (pairCount==0) {c_a=c; s_a=s; u_a=u; pairCount=1;} else if (pairCount==1){c_b=c; s_b=s; u_b=u; pairCount=2;}
                }
            }
            double cx=0, cy=0; bool haveXY=false;
            if (used >= 2) {
                if (used == 2) {
                    // Direct 2x2 solve
                    double det = c_a*s_b - s_a*c_b; // = sin(theta_b - theta_a)
                    if (fabs(det) > 1e-9) {
                        cx = ( u_a*s_b - s_a*u_b)/det;
                        cy = ( c_a*u_b - u_a*c_b)/det;
                        haveXY = true;
                    }
                } else { // used ==3 (or more if extended)
                    double det = Scc*Sss - Scs*Scs;
                    if (fabs(det) > 1e-12) {
                        cx = ( Sss*Su_c - Scs*Su_s)/det;
                        cy = ( Scc*Su_s - Scs*Su_c)/det;
                        haveXY = true;
                    }
                }
            }
            if (haveXY) {
                // Fill 2-or-3 clusters map always
                h_recoCenter_2to3->Fill(cx, cy);
                g_recoCenters_2to3->SetPoint(g_recoCenters_2to3->GetN(), cx, cy);
                // Additionally fill the exactly-3-clusters map when applicable
                if (acceptExactlyThree) {
                    h_recoCenter_3->Fill(cx, cy);
                    g_recoCenters_3->SetPoint(g_recoCenters_3->GetN(), cx, cy);
                }
            }
            }
    }
    // Fill hits-per-event with the number of masked triggered hits (not cluster-gated)
    h_hitsInEvent->Fill(eventHitsMasked);
    // Do not constrain other histograms; only centers are gated by the cluster criterion

        if (debug) this_event->PrintValidHits();      
        
        if (verbose) LogInfo << "Stored entry " << entryit << " in the vector of Events" << std::endl;

        delete this_event;

    }

    LogInfo << "Read all entries" << std::endl;

    LogInfo << "Number of triggered events: " << (double) triggeredEvents/limit *100 << "%" << std::endl;
    
    ///////////////////////////
    
    // plots

    // create a canvas
    LogInfo << "Creating canvas" << std::endl;
    std::string runNumber = "99999"; // fallback value
    if (clp.isOptionTriggered("runNumber")) {
        runNumber = clp.getOptionVal<std::string>("runNumber");
        LogInfo << "Using provided run number: " << runNumber << std::endl;
    } else {
        LogInfo << "No run number provided, using fallback: " << runNumber << std::endl;
    }

    // Populate sigma and baseline TGraphs with the correct values per channel
    LogInfo << "Populating sigma and baseline graphs" << std::endl;
    for (int detit = 0; detit < nDetectors; detit++) {
        for (int chit = 0; chit < nChannels; chit++) {
            g_sigma->at(detit)->SetPoint(chit, chit, baseline_sigma.at(detit).at(chit));
            g_baseline->at(detit)->SetPoint(chit, chit, baseline.at(detit).at(chit));
        }
        // Update graph titles with run number
        g_sigma->at(detit)->SetTitle(Form("Sigma (Detector %d) - Run %s", detit, runNumber.c_str()));
        g_baseline->at(detit)->SetTitle(Form("Baseline (Detector %d) - Run %s", detit, runNumber.c_str()));
    }

    // Update histogram titles with run number
    LogInfo << "Updating histogram titles with run number" << std::endl;
    for (int i = 0; i < nDetectors; i++) {
        h_firingChannels->at(i)->SetTitle(Form("Firing channels (Detector %d) - Run %s", i, runNumber.c_str()));
        h_amplitude->at(i)->SetTitle(Form("Amplitude (Detector %d) - Run %s", i, runNumber.c_str()));
        h_amplitudeVsChannel->at(i)->SetTitle(Form("Amplitude vs Channel (Detector %d) - Run %s", i, runNumber.c_str()));
        
        // Update raw peak histogram titles if in verbose mode
        if (verbose) {
            for (int j = 0; j < nChannels; j++) {
                h_rawPeak->at(i)->at(j)->SetTitle(Form("Raw peak (Detector %d, Channel %d) - Run %s", i, j, runNumber.c_str()));
            }
        }
    }
    h_hitsInEvent->SetTitle(Form("Hits in event - Run %s", runNumber.c_str()));

    // Update canvas titles
    TCanvas *c_channelsFiring = new TCanvas(Form("c_channelsFiring_Run%s", runNumber.c_str()), Form("Channels Firing (Run %s)", runNumber.c_str()), 800, 600);
    c_channelsFiring->Divide(2, 2);

    TCanvas *c_sigma = new TCanvas(Form("c_sigma_Run%s", runNumber.c_str()), Form("Sigma (Run %s)", runNumber.c_str()), 800, 600);
    c_sigma->Divide(2, 2);

    TCanvas *c_baseline = new TCanvas(Form("c_baseline_Run%s", runNumber.c_str()), Form("Baseline (Run %s)", runNumber.c_str()), 800, 600);
    c_baseline->Divide(2, 2);

    // Update raw peak canvas titles
    std::vector<TCanvas*> *c_rawPeak = new std::vector<TCanvas*>;
    c_rawPeak->reserve(6);
    if (verbose) {
        for (int i = 0; i < 6; i++) {
            TCanvas *this_c_rawPeak = new TCanvas(Form("c_rawPeak%d_Run%s", i, runNumber.c_str()), Form("Raw Peak %d (Run %s)", i, runNumber.c_str()), 800, 600);
            this_c_rawPeak->Divide(8, 8);
            c_rawPeak->emplace_back(this_c_rawPeak);
        }
    }

    TCanvas *c_amplitude = new TCanvas(Form("c_amplitude_Run%s", runNumber.c_str()), Form("Amplitude (Run %s)", runNumber.c_str()), 800, 600);
    c_amplitude->Divide(2, 2);

    TCanvas *c_amplitudeVsChannel = new TCanvas(Form("c_amplitudeVsChannel_Run%s", runNumber.c_str()), Form("Amplitude vs Channel (Run %s)", runNumber.c_str()), 800, 600);
    c_amplitudeVsChannel->Divide(2, 2);

    TCanvas *c_hitsInEvent = new TCanvas(Form("c_hitsInEvent_Run%s", runNumber.c_str()), Form("Hits in Event (Run %s)", runNumber.c_str()), 800, 600);

    TCanvas *c_recoCenter_3 = new TCanvas("c_recoCenter_3", "Interpolated Centers (exactly 3)", 700, 500);
    c_recoCenter_3->cd();
    h_recoCenter_3->Draw("COLZ");
    g_recoCenters_3->SetMarkerColor(kBlack);
    g_recoCenters_3->Draw("P SAME");
    c_recoCenter_3->Update();

    TCanvas *c_recoCenter_2to3 = new TCanvas("c_recoCenter_2to3", "Interpolated Centers (2 or 3)", 700, 500);
    c_recoCenter_2to3->cd();
    h_recoCenter_2to3->Draw("COLZ");
    g_recoCenters_2to3->SetMarkerColor(kBlack);
    g_recoCenters_2to3->Draw("P SAME");
    c_recoCenter_2to3->Update();

    LogInfo << "Drawing histograms" << std::endl;
    for (int i = 0; i < nDetectors; i++) {
        c_channelsFiring->cd(i+1);
    // Use log scale for counts and ensure a positive minimum
    gPad->SetLogy();
    h_firingChannels->at(i)->SetMinimum(0.5);
    h_firingChannels->at(i)->Draw();
    }
    c_channelsFiring->Update();

    // Print top-N spike channels per detector for quick inspection
    const int topN = 10;
    LogInfo << "Top firing channels per detector (channel:count)" << std::endl;
    for (int i = 0; i < nDetectors; i++) {
        std::vector<std::pair<int,double>> counts; counts.reserve(nChannels);
        for (int b = 1; b <= nChannels; ++b) {
            double c = h_firingChannels->at(i)->GetBinContent(b);
            if (c > 0) counts.emplace_back(b-1, c); // channel index = bin-1
        }
        std::sort(counts.begin(), counts.end(), [](const auto& a, const auto& b){return a.second > b.second;});
        LogInfo << "  Detector " << i << ": ";
        int limitN = std::min(topN, (int)counts.size());
        for (int k = 0; k < limitN; ++k) {
            LogInfo << counts[k].first << ":" << (long long)counts[k].second << (k==limitN-1?"":"  ");
        }
        LogInfo << std::endl;
    }

    for (int i = 0; i < nDetectors; i++) {
        c_sigma->cd(i+1);
        g_sigma->at(i)->Draw("AP");
    }

    for (int i = 0; i < nDetectors; i++) {
        c_baseline->cd(i+1);
        g_baseline->at(i)->Draw("AP");
    }

    // note that only raw peaks of detector 0 are being plotted
    if (verbose){
        for (int ch = 0; ch < nChannels; ch++) {
            c_rawPeak->at(ch/64)->cd(ch%64+1);
            h_rawPeak->at(0)->at(ch)->Draw();
        }
    }


    for (int i = 0; i < nDetectors; i++) {
        c_amplitude->cd(i+1);
        h_amplitude->at(i)->Draw();
    }

    for (int i = 0; i < nDetectors; i++) {
        c_amplitudeVsChannel->cd(i+1);
        h_amplitudeVsChannel->at(i)->Draw("COLZ");
    }

    c_hitsInEvent->cd();
    gPad->SetLogy();
    h_hitsInEvent->Draw();

    // create a pdf report containing firing channels, sigma and amplitude
    std::string outputDir = clp.getOptionVal<std::string>("outputDir");
    // LogInfo << "Output directory: " << outputDir << std::endl;
    // output filename same as input root file, but remove .root and create single PDF report
    std::string input_file_base = input_root_filename.substr(input_root_filename.find_last_of("/\\") + 1);
    size_t lastdot = input_file_base.find_last_of(".");
    if (lastdot != std::string::npos) {
        input_file_base = input_file_base.substr(0, lastdot);
    }
    // Single PDF report filename
    std::string output_filename_report = outputDir + "/" + input_file_base + "_" + std::to_string(clp.getOptionVal<int>("nSigma")) + "sigma_report.pdf";

    LogInfo << "Output PDF report: " << output_filename_report << std::endl;
    
    // Check output directory permissions
    if (system(("test -w " + outputDir).c_str()) != 0) {
        LogError << "Output directory " << outputDir << " is not writable!" << std::endl;
        return 1;
    }

    LogInfo << "Saving canvases to multi-page pdf report" << std::endl;
    
    // Ensure output directory exists
    system(("mkdir -p " + outputDir).c_str());
    
    // save the canvases to a multi-page pdf report with error handling
    try {
        // Check that all canvases are valid before attempting to save
        if (!c_channelsFiring || !c_sigma || !c_baseline || !c_amplitude || !c_amplitudeVsChannel || !c_hitsInEvent) {
            LogError << "One or more canvases are null, cannot save plots" << std::endl;
            return 1;
        }
        
        LogInfo << "Creating multi-page PDF report..." << std::endl;
        
        // Page 1: Channels firing plot
        c_channelsFiring->Update();
        c_channelsFiring->SaveAs((output_filename_report + "(").c_str());
        
        // Page 2: Sigma plot
        c_sigma->Update();
        c_sigma->SaveAs(output_filename_report.c_str());
        
        // Page 3: Baseline plot
        c_baseline->Update();
        c_baseline->SaveAs(output_filename_report.c_str());
        
        // Page 4: Amplitude plot
        c_amplitude->Update();
        c_amplitude->SaveAs(output_filename_report.c_str());
        
        // Page 5: Reconstructed centers heatmap (exactly 3)
        if (c_recoCenter_3) {
            c_recoCenter_3->Update();
            c_recoCenter_3->SaveAs(output_filename_report.c_str());
        }
        // Page 6: Reconstructed centers heatmap (2 or 3)
        if (c_recoCenter_2to3) {
            c_recoCenter_2to3->Update();
            c_recoCenter_2to3->SaveAs(output_filename_report.c_str());
        }
        
        // Page 7: Hits in event plot (final page)
        c_hitsInEvent->Update();
        c_hitsInEvent->SaveAs((output_filename_report + ")").c_str());
        
        LogInfo << "Multi-page PDF report saved successfully: " << output_filename_report << std::endl;
    } catch (const std::exception& e) {
        LogError << "Error saving PDF report: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        LogError << "Unknown error occurred while saving PDF report" << std::endl;
        return 1;
    }

    // Check if showPlots is enabled in JSON settings
    bool showPlots = parseBool(jsonSettings, "showPlots", false);
    LogInfo << "Show plots: " << (showPlots ? "true" : "false") << std::endl;

    if (showPlots) {
        LogInfo << "Running root viewer..." << std::endl;
        app->Run();
    }
    else { LogInfo << "If you wish to see the plots, you need to set showPlots option to true.\n";}

    return 0;
}
