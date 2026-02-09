// -------------------------------------------------------------------------
// Plot PDF Variations (BJ Method v3 - Dynamic Range & 0.01 Binning)
// File: plot_pdf_variations_BJ_v3.cpp
//
// [Logic]
// 1. Event Loop (Collection):
//    - Collect all sys_pdf values into memory.
//    - Simultaneously find the Global Min and Max of sys_pdf values.
// 2. Histogram Setup:
//    - Set Y-axis Range: [Global Min, Global Max]
//    - Calculate N_bins: (Max - Min) / 0.01
// 3. Bin Loop (Filling):
//    - Fill the dynamically created Heatmap.
//
// compile: g++ -o plot_pdf_variations_BJ_v3.exe plot_pdf_variations_BJ_v3.cpp $(root-config --cflags --glibs)
// run: ./plot_pdf_variations_BJ_v3.exe final_output.root
// -------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm> // for min/max

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TLine.h"

using namespace std;

// --- Binning Definition ---
const int nBins = 14;
const int binNumbers[nBins] = {
    22, 23, 24, // Nb=0
    25, 26, 27, // Nb=1
    28, 29, 30, // Nb=2
    31, 32, 33, // Nb=3
    35, 36      // Nb>=4 (Bin 34 skipped)
};

// Helper: Get Array Index (0~13) from Bin Number
int getIdx(int binNum) {
    for(int i=0; i<nBins; ++i) {
        if(binNumbers[i] == binNum) return i;
    }
    return -1;
}

// Helper: Determine Bin Number
int getBinNumber(int njets, int nbm) {
    int j_cat = -1;
    if (njets >= 4 && njets <= 5) j_cat = 0;
    else if (njets >= 6 && njets <= 7) j_cat = 1;
    else if (njets >= 8) j_cat = 2;

    if (j_cat == -1) return -1;

    if (nbm == 0) return (22 + j_cat);
    if (nbm == 1) return (25 + j_cat);
    if (nbm == 2) return (28 + j_cat);
    if (nbm == 3) return (31 + j_cat);
    if (nbm >= 4) {
        if (j_cat == 0) return 31; // Merged
        if (j_cat == 1) return 35;
        if (j_cat == 2) return 36;
    }
    return -1;
}

int main(int argc, char* argv[]) {
    // Style Settings
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(255);

    if (argc < 2) {
        cout << "Usage: ./plot_pdf_variations_BJ_v3.exe [root_file]" << endl;
        return 1;
    }

    TString filename = argv[1];
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        cout << "Error opening file: " << filename << endl;
        return 1;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        cout << "Tree 'tree' not found!" << endl;
        return 1;
    }

    // --- Branch Setup ---
    vector<float> *sys_pdf = nullptr;
    int nleps, njets, nbm;

    tree->SetBranchAddress("nleps", &nleps);
    tree->SetBranchAddress("njets", &njets);
    tree->SetBranchAddress("nbm", &nbm);

    if (tree->GetBranch("sys_pdf")) {
        tree->SetBranchAddress("sys_pdf", &sys_pdf);
    } else {
        cout << "[Error] 'sys_pdf' branch is required!" << endl;
        return 1;
    }

    // --- Data Storage ---
    // [BinIndex] -> List of sys_pdf values
    vector<vector<double>> bin_data(nBins);

    // Variables for Dynamic Range
    double y_min = 1.0e9;  // Initialize with large number
    double y_max = -1.0e9; // Initialize with small number

    // --- Step 1: Event Loop (Collect & Find Min/Max) ---
    Long64_t nentries = tree->GetEntries();
    cout << "Step 1: Collecting events from " << nentries << " entries..." << endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        if (!sys_pdf || sys_pdf->size() < 2) continue;

        if (nleps != 1) continue;

        int binNum = getBinNumber(njets, nbm);
        if (binNum == -1) continue;
        int bIdx = getIdx(binNum);
        if (bIdx == -1) continue;

        double val_up   = sys_pdf->at(0);
        double val_down = sys_pdf->at(1);

        // Store data
        bin_data[bIdx].push_back(val_up);
        bin_data[bIdx].push_back(val_down);

        // Check Min/Max
        if (val_up < y_min) y_min = val_up;
        if (val_up > y_max) y_max = val_up;

        if (val_down < y_min) y_min = val_down;
        if (val_down > y_max) y_max = val_down;
    }

    // Safety check if no data found
    if (y_min > y_max) {
        cout << "No valid events found within cuts. Setting default range." << endl;
        y_min = 0.0; y_max = 2.0;
    }

    // --- Histogram Setup (Dynamic) ---
    // Calculate Number of Bins for 0.01 width
    double range = y_max - y_min;
    int nYbins = (int)ceil(range / 0.01);

    // Safety for single bin
    if (nYbins < 1) nYbins = 1;

    cout << "Dynamic Range: [" << y_min << ", " << y_max << "]" << endl;
    cout << "Bin Width: 0.01 -> Total Bins: " << nYbins << endl;

    TH2D* h_map = new TH2D("h_map", "", nBins, 0, nBins, nYbins, y_min, y_max);

    // --- Step 2: Bin Loop (Fill) ---
    cout << "Step 2: Filling Heatmap..." << endl;

    for (int b = 0; b < nBins; ++b) {
        vector<double>& values = bin_data[b];

        for (double val : values) {
            h_map->Fill(b + 0.5, val);
        }
    }

    // --- Drawing ---
    TCanvas* c1 = new TCanvas("c1", "PDF Variations BJ v3 Dynamic", 1200, 700);
    c1->SetRightMargin(0.15);
    c1->SetGridx();
    c1->SetGridy();

    // Axis Settings
    h_map->GetXaxis()->SetLabelSize(0.04);
    for(int i=0; i<nBins; ++i) {
        h_map->GetXaxis()->SetBinLabel(i+1, Form("Bin %d", binNumbers[i]));
    }

    h_map->GetYaxis()->SetTitle("sys_pdf Value (0.01 binning)");
    h_map->GetYaxis()->SetTitleOffset(1.2);

    h_map->Draw("COLZ");

    // Nominal Reference Line (at 1.0)
    TLine* line = new TLine(0, 1.0, nBins, 1.0);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw();

    TLegend* leg = new TLegend(0.65, 0.82, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->AddEntry(h_map, "Event Density", "f");
    leg->AddEntry(line, "Nominal (1.0)", "l");
    leg->Draw();

    c1->SaveAs("pdf_variations_BJ_v3.png");
    c1->SaveAs("pdf_variations_BJ_v3.pdf");

    cout << "Plot saved as pdf_variations_BJ_v3.png" << endl;

    delete h_map;
    delete line;
    delete c1;
    file->Close();
    return 0;
}
