// -------------------------------------------------------------------------
// Plot PDF Variations (BJ Method v3 - Sum of Weights Logic)
// Author : Gemini
//
// [Logic]
// 1. Loop Events -> Get 'weight' vector.
// 2. Calculate Sum = Sum(weight[0]...weight[100]) for the event.
// 3. Normalize: Ratio = Sum / weight[0].
// 4. Collect these Ratios per Bin.
// 5. Sort & Plot:
//    - Cyan: Lines for ALL events.
//    - Blue: Lines for 16th/84th percentile events.
//
// compile: g++ -o plot_pdf_variations_BJ_v3.exe plot_pdf_variations_BJ_v3.cpp $(root-config --cflags --glibs)
// run: ./plot_pdf_variations_BJ_v3.exe final_output.root
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
    vector<float> *weight_vec = nullptr;
    int nleps, njets, nbm;

    tree->SetBranchAddress("weight", &weight_vec);
    tree->SetBranchAddress("nleps", &nleps);
    tree->SetBranchAddress("njets", &njets);
    tree->SetBranchAddress("nbm", &nbm);

    // --- Data Storage ---
    // Stores the calculated Ratio for each event in each bin
    // [Bin][EventIndex]
    vector<vector<double>> bin_event_ratios(nBins);

    // --- Step 1: Single Event Loop ---
    Long64_t nentries = tree->GetEntries();
    cout << "Processing " << nentries << " events (Single Loop)..." << endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        if (!weight_vec || weight_vec->empty()) continue;

        // Cuts
        if (nleps != 1) continue;

        // Determine Bin
        int binNum = getBinNumber(njets, nbm);
        if (binNum == -1) continue;

        int bIdx = getIdx(binNum);
        if (bIdx == -1) continue;

        // [Logic] Calculate Ratio per Event
        double w_nom = weight_vec->at(0);
        if (w_nom == 0) w_nom = 1.0;

        double w_sum = 0.0;
        // Sum weights from index 0 to 99 (Total 100 items)
        int count = 0;
        for (int k = 0; k < 100; ++k) {
            if (k < weight_vec->size()) {
                w_sum += weight_vec->at(k);
                count++;
            }
        }

        // Ratio = Sum / (Nominal * 100)
        // Note: Using 100.0 exactly as requested, even if vector size < 100 (though unlikely)
        double ratio = w_sum / (w_nom * 100.0);

        bin_event_ratios[bIdx].push_back(ratio);
    }

    // --- Prepare Histograms ---
    vector<TH1D*> cyan_lines;
    vector<TH1D*> blue_lines;

    // Use a frame histogram for axes
    TH1D* h_frame = new TH1D("h_frame", "", nBins, 0, nBins);

    // --- Step 2: Sort and Create Lines per Bin ---
    cout << "Analyzing event distributions per bin..." << endl;

    for (int b = 0; b < nBins; ++b) {
        h_frame->SetBinContent(b+1, 1.0); // Reference Line at 1.0

        vector<double>& events = bin_event_ratios[b];
        int N = events.size();

        if (N > 0) {
            // 1. Sort the events based on Ratio
            std::sort(events.begin(), events.end());

            // 2. Identify Indices for Blue Lines (16th, 84th percentile)
            int idx_16 = (int)(N * 0.16);
            int idx_84 = (int)(N * 0.84);
            if(idx_84 >= N) idx_84 = N - 1;

            // 3. Create Lines
            for (int i = 0; i < N; ++i) {
                double val = events[i];

                TH1D* h_line = new TH1D("", "", nBins, 0, nBins);
                h_line->SetBinContent(b+1, val);
                h_line->SetLineWidth(1);

                if (i == idx_16 || i == idx_84) {
                    h_line->SetLineColor(kBlue);
                    h_line->SetLineWidth(2);
                    blue_lines.push_back(h_line);
                } else {
                    h_line->SetLineColor(kCyan);
                    cyan_lines.push_back(h_line);
                }
            }
        }
    }

    // --- Drawing ---
    TCanvas* c1 = new TCanvas("c1", "PDF Variations BJ v3 Ratio", 1000, 600);
    c1->SetGridy();

    // Axis Settings
    // Values should be around 1.0
    h_frame->GetYaxis()->SetRangeUser(0.5, 1.15); // Adjust range as needed!
    h_frame->GetYaxis()->SetTitle("#Sigma(w_{0}..w_{99}) / (w_{0} * 100)");
    for(int i=0; i<nBins; ++i) h_frame->GetXaxis()->SetBinLabel(i+1, Form("Bin %d", binNumbers[i]));
    h_frame->GetXaxis()->SetLabelSize(0.04);
    h_frame->Draw("HIST");

    // 1. Draw Cyan Lines (All Events)
    cout << "Drawing " << cyan_lines.size() << " cyan event lines..." << endl;
    for (auto h : cyan_lines) {
        h->Draw("HIST SAME");
    }

    // 2. Draw Blue Lines (16th/84th Percentile Events)
    cout << "Drawing " << blue_lines.size() << " blue event lines..." << endl;
    for (auto h : blue_lines) {
        h->Draw("HIST SAME");
    }

    // 3. Draw Nominal Reference (Black)
    h_frame->SetLineColor(kBlack);
    h_frame->SetLineWidth(2);
    h_frame->SetLineStyle(2);
    h_frame->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.55, 0.78, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h_frame, "Nominal (1.0)", "l");
    if(blue_lines.size() > 0) leg->AddEntry(blue_lines[0], "16th/84th Percentile", "l");
    if(cyan_lines.size() > 0) leg->AddEntry(cyan_lines[0], "Individual Events", "l");
    leg->Draw();

    c1->SaveAs("pdf_variations_BJ_v3.png");
//    c1->SaveAs("pdf_variations_BJ_v3.pdf");

    cout << "Plot saved as pdf_variations_BJ_v3.png" << endl;

    // Cleanup
    for(auto h : cyan_lines) delete h;
    for(auto h : blue_lines) delete h;
    file->Close();
    return 0;
}
