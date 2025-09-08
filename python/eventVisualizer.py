import os
import argparse
import ROOT
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='Event Visualizer for .root files')
    parser.add_argument('-f', '--file', type=str, required=True,
                        help='Input .root file path')

    args = parser.parse_args()
    input_file = Path(args.file)

    if not input_file.exists():
        print(f"Error: Input file '{input_file}' does not exist.")
        return

    if not str(input_file).endswith('.root'):
        print(f"Error: Input file '{input_file}' is not a ROOT file.")
        return

    root_file = ROOT.TFile(str(input_file), "READ")

    tree = root_file.Get("events")
    if not tree:
        print("Error: Could not find 'events' tree in the ROOT file.")
        root_file.Close()
        return

    n_entries = tree.GetEntries()
    print(f"Found {n_entries} events in the ROOT file")

    # Create canvas and histograms
    canvasPeak = ROOT.TCanvas("canvasPeak", "Peak Values by Detector", 1200, 800)
    canvasPeak.Divide(3, 2)

    hists = []
    histsSuppr = []
    for det in range(3):
        hist = ROOT.TH1D(f"hist{det}", f"Detector {det}", 384, 0, 384)
        hist.SetLineColor(det + 1)
        hists.append(hist)
        histsSuppr.append(hist.Clone())

    for entry in range(n_entries):
        tree.GetEntry(entry)
        print(f"Visualizing event {entry}/{n_entries}")

        # Fill histograms
        for det in range(3):
            hists[det].Reset()
            histsSuppr[det].Reset()
            for ch in range(384):
                hists[det].Fill(ch, tree.peak[det][ch])
                histsSuppr[det].Fill(ch, tree.peakZeroSuppr[det][ch])

            canvasPeak.cd(det + 1)
            hists[det].Draw()
            canvasPeak.cd(3 + det + 1)

            histsSuppr[det].SetTitle(f"nClusters = {tree.nClusters[det]}")
            histsSuppr[det].Draw()

        canvasPeak.Update()
        ROOT.gPad.Update()

        # Process user input between events
        input("Press Enter to continue to next event...")
        # ROOT.gPad.WaitPrimitive()

    root_file.Close()


if __name__ == '__main__':
    main()
