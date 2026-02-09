// -------------------------------------------------------------------------
// Plot PDF Variations (BJ Method v4 - TLine Fixed)
// File: plot_pdf_variations_BJ_v4.cpp
//
// [Fix]
// - Replaced TH1D with TLine for drawing events.
// - Removes unwanted vertical lines connecting to y=0.
//
// compile: g++ -o plot_pdf_variations_BJ_v4.exe plot_pdf_variations_BJ_v4.cpp $(root-config --cflags --glibs)
// run: ./plot_pdf_variations_BJ_v4.exe output_nominal_newnt_UL2018.root
// -------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
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
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    if (argc < 2) {
        cout << "Usage: ./plot_pdf_variations_BJ_v4.exe [root_file]" << endl;
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
    vector<float> *weight_vec = nullptr;
    int nleps, njets, nbm;

    tree->SetBranchAddress("weight", &weight_vec);
    tree->SetBranchAddress("nleps", &nleps);
    tree->SetBranchAddress("njets", &njets);
    tree->SetBranchAddress("nbm", &nbm);

    // --- Data Storage ---
    vector<vector<double>> bin_data(nBins);
    double y_min = 1.0e9;
    double y_max = -1.0e9;

    // --- Step 1: Event Loop (Collect) ---
    Long64_t nentries = tree->GetEntries();
    cout << "Processing " << nentries << " events..." << endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        if (!weight_vec || weight_vec->empty()) continue;
        if (nleps != 1) continue;

        int binNum = getBinNumber(njets, nbm);
        if (binNum == -1) continue;
        int bIdx = getIdx(binNum);
        if (bIdx == -1) continue;

        double sum = 0.0;
        int limit = (weight_vec->size() < 100) ? weight_vec->size() : 100;
        for(int k=0; k<limit; ++k) sum += weight_vec->at(k);

        double avg_val = sum / 100.0;
        bin_data[bIdx].push_back(avg_val);

        if (avg_val < y_min) y_min = avg_val;
        if (avg_val > y_max) y_max = avg_val;
    }

    if (y_min > y_max) { y_min = 0; y_max = 1; }

    // --- Step 2: Sorting & Drawing with TLine ---
    TCanvas* c1 = new TCanvas("c1", "PDF Variations BJ v4", 1200, 600);
    c1->SetGridy();

    // Draw Axis Frame first
    TH1D* h_frame = new TH1D("h_frame", "", nBins, 0, nBins);
    double y_range = y_max - y_min;
    h_frame->GetYaxis()->SetRangeUser(y_min - y_range*0.1, y_max + y_range*0.1);
    h_frame->GetYaxis()->SetTitle("Average Weight (Sum / 100)");
    for(int i=0; i<nBins; ++i) h_frame->GetXaxis()->SetBinLabel(i+1, Form("Bin %d", binNumbers[i]));
    h_frame->GetXaxis()->SetLabelSize(0.04);
    h_frame->Draw("AXIS"); // Draw axes only

    // Store lines to delete later
    vector<TLine*> lines_cyan;
    vector<TLine*> lines_blue;

    cout << "Creating lines..." << endl;

    for (int b = 0; b < nBins; ++b) {
        vector<double>& values = bin_data[b];
        int N = values.size();
        if (N == 0) continue;

        // 1. Sort
        std::sort(values.begin(), values.end());

        // 2. Identify Indices
        int idx_16 = (int)(N * 0.16);
        int idx_84 = (int)(N * 0.84);
        if (idx_84 >= N) idx_84 = N - 1;

        // 3. Create Lines
        for (int i = 0; i < N; ++i) {
            double val = values[i];

            // Line from x=b to x=b+1
            TLine* line = new TLine(b, val, b+1, val);

            if (i == idx_16 || i == idx_84) {
                line->SetLineColor(kBlue);
                line->SetLineWidth(2);
                lines_blue.push_back(line);
            } else {
                line->SetLineColor(kCyan);
                line->SetLineWidth(1);
                lines_cyan.push_back(line);
            }
        }
    }

    // Draw Cyan lines first (Background)
    for (auto l : lines_cyan) l->Draw();

    // Draw Blue lines on top (Foreground)
    for (auto l : lines_blue) l->Draw();

    // Redraw Axis on top
    h_frame->Draw("AXIS SAME");

    // Legend (Use dummy lines for correct style)
    TLine* dummy_cyan = new TLine(0,0,1,1); dummy_cyan->SetLineColor(kCyan); dummy_cyan->SetLineWidth(1);
    TLine* dummy_blue = new TLine(0,0,1,1); dummy_blue->SetLineColor(kBlue); dummy_blue->SetLineWidth(2);

    TLegend* leg = new TLegend(0.65, 0.78, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(dummy_cyan, "Individual Events", "l");
    leg->AddEntry(dummy_blue, "16th/84th Percentile", "l");
    leg->Draw();

    c1->SaveAs("pdf_variations_BJ_v4.png");
    c1->SaveAs("pdf_variations_BJ_v4.pdf");

    cout << "Plot saved as pdf_variations_BJ_v4.png" << endl;

    // Cleanup
    for(auto l : lines_cyan) delete l;
    for(auto l : lines_blue) delete l;
    delete dummy_cyan; delete dummy_blue;
    delete h_frame;
    delete c1;
    file->Close();
    return 0;
}
