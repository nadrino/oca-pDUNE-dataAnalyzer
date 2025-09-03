#!/usr/bin/env python3

import argparse
import ROOT
from ROOT import TFile, TTree


def main():
    parser = argparse.ArgumentParser(description='Draw beam monitor data from ROOT file')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input ROOT file')
    args = parser.parse_args()

    # Open ROOT file
    input_file = TFile(args.input, "READ")
    if not input_file or input_file.IsZombie():
        print(f"Error: Could not open file {args.input}")
        return 1

    event_tree = input_file.Get("events")
    histograms = []
    for i in range(3):
        h = ROOT.TH1F(f"h_xBarycenter_{i}", f"X Barycenter Detector {i}", 100, -50, 50)
        event_tree.Draw(f"xBarycenter[{i}]>>h_xBarycenter_{i}")
        if i == 0:
            h2d = ROOT.TH2F("h_xyBarycenter", "Y vs X Barycenter", 100, -50, 50, 100, -50, 50)
            event_tree.Draw("yBarycenter:xBarycenter[0]>>h_xyBarycenter", "", "colz")
        histograms.append(h)



    # Create canvas with 3x2 grid
    canvas = ROOT.TCanvas("canvas", "Beam Monitor Data", 1200, 800)
    canvas.Divide(3, 2)



    input_file.Close()
    return 0


if __name__ == "__main__":
    exit(main())
