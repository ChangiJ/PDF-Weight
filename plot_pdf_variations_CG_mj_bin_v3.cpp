// -------------------------------------------------------------------------
// Plot PDF Variations (CG Method with Mj Binning - Grid Layout)
// File: plot_pdf_variations_CG_mj_bin_v3.cpp
//
// [Logic Flow]
// 1. Accumulate:
//    - Loop over Events.
//    - Identify (Physical Bin, Mj Bin).
//    - Sum weights for all 100 replicas:
//      sum[bin][mj][k] += weight[k]
//
// 2. Process (Per Physical Bin & Per Mj Bin):
//    - Calculate Nominal Sum (k=0).
//    - Calculate Ratios for all k: Ratio[k] = Sum[k] / Sum[0].
//    - Fill Cyan Lines (All 100 Ratios).
//    - Sort Ratios.
//    - Identify 16th/84th percentile Ratios.
//    - Fill Blue Lines (Envelope).
//
// compile: g++ -o plot_pdf_variations_CG_mj_bin_v3.exe plot_pdf_variations_CG_mj_bin_v3.cpp $(root-config --cflags --glibs)
// run: ./plot_pdf_variations_CG_mj_bin_v3.exe final_output.root
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
#include "TLatex.h"
#include "TPad.h"

using namespace std;

// --- Physical Binning Definition ---
const int nBins = 14;
const int binNumbers[nBins] = {
    22, 23, 24, // Nb=0
    25, 26, 27, // Nb=1
    28, 29, 30, // Nb=2
    31, 32, 33, // Nb=3
    35, 36      // Nb>=4 (Bin 34 skipped)
};

// --- Mj12 Binning Definition ---
const int nMjBins = 3;
const string mjLabels[nMjBins] = {
    "500-800",
    "800-1100",
    "1100+"
};

// Helper: Get Array Index (0~13) from Bin Number
int getIdx(int binNum) {
    for(int i=0; i<nBins; ++i) {
        if(binNumbers[i] == binNum) return i;
    }
    return -1;
}

// Helper: Determine Physical Bin Number
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

// Helper: Map Bin Number to Pad Number (1~15)
// Grid: 3 Cols (Njets) x 5 Rows (Nb)
int getPadNumber(int binNum) {
    if (binNum == 22) return 1; if (binNum == 23) return 2; if (binNum == 24) return 3;
    if (binNum == 25) return 4; if (binNum == 26) return 5; if (binNum == 27) return 6;
    if (binNum == 28) return 7; if (binNum == 29) return 8; if (binNum == 30) return 9;
    if (binNum == 31) return 10; if (binNum == 32) return 11; if (binNum == 33) return 12;
    if (binNum == 35) return 14; if (binNum == 36) return 15;
    return -1;
}

// Helper: Determine Mj12 Bin Index
int getMjBinIndex(float mj12) {
    if (mj12 >= 500 && mj12 < 800) return 0;
    else if (mj12 >= 800 && mj12 < 1100) return 1;
    else if (mj12 >= 1100) return 2;
    return -1;
}

int main(int argc, char* argv[]) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    if (argc < 2) {
        cout << "Usage: ./plot_pdf_variations_CG_mj_bin_v3.exe [root_file]" << endl;
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
    float mj12;

    tree->SetBranchAddress("weight", &weight_vec);
    tree->SetBranchAddress("nleps", &nleps);
    tree->SetBranchAddress("njets", &njets);
    tree->SetBranchAddress("nbm", &nbm);
    tree->SetBranchAddress("mj12", &mj12);

    // --- Data Storage (Accumulator) ---
    // [PhysicalBin][MjBin][Replica]
    // 14 Bins, 3 Mj Bins, 100 Replicas
    vector<vector<vector<double>>> bin_mj_replica_sums(nBins,
        vector<vector<double>>(nMjBins,
            vector<double>(100, 0.0)
        )
    );

    // --- Step 1: Event Loop (Accumulate Sums) ---
    Long64_t nentries = tree->GetEntries();
    cout << "Step 1: Accumulating weights from " << nentries << " events..." << endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        if (!weight_vec || weight_vec->empty()) continue;
        if (nleps != 1) continue;

        // 1. Identify Bins
        int binNum = getBinNumber(njets, nbm);
        if (binNum == -1) continue;
        int bIdx = getIdx(binNum);
        if (bIdx == -1) continue;

        int mIdx = getMjBinIndex(mj12);
        if (mIdx == -1) continue;

        // 2. Accumulate Weights (Summing w_pdf[evt][k])
        if (weight_vec->size() >= 100) {
            for(int k=0; k<100; ++k) {
                bin_mj_replica_sums[bIdx][mIdx][k] += weight_vec->at(k);
            }
        }
    }

    // --- Step 2: Drawing on Grid Canvas ---
    cout << "Step 2: Processing and Drawing..." << endl;

    TCanvas* c1 = new TCanvas("c1", "PDF Variations Grid v3", 1200, 1600);
    c1->Divide(3, 5, 0.01, 0.01);

    vector<TH1D*> trash_bin; // To keep histograms alive

    // Loop over Physical Bins (The 15 Pads)
    for (int b = 0; b < nBins; ++b) {
        int binNum = binNumbers[b];
        int padNum = getPadNumber(binNum);
        if (padNum < 1 || padNum > 15) continue;

        c1->cd(padNum);
        gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);

        // Prepare Histograms for this Pad
        TH1D* h_nom = new TH1D(Form("h_nom_%d", b), "", nMjBins, 0, nMjBins);
        TH1D* h_up  = new TH1D(Form("h_up_%d", b), "", nMjBins, 0, nMjBins);
        TH1D* h_down = new TH1D(Form("h_down_%d", b), "", nMjBins, 0, nMjBins);
        trash_bin.push_back(h_nom); trash_bin.push_back(h_up); trash_bin.push_back(h_down);

        // 100 Histograms for Cyan Lines
        vector<TH1D*> h_reps(100);
        for(int k=0; k<100; ++k) {
            h_reps[k] = new TH1D(Form("h_rep_%d_%d", b, k), "", nMjBins, 0, nMjBins);
            trash_bin.push_back(h_reps[k]);
        }

        // --- Inner Loop: Mj Bins (X-axis) ---
        // For each Mj bin, calculate ratios, sort, and fill histograms
        for (int m = 0; m < nMjBins; ++m) {

            // 1. Get the accumulated sums for this (Bin, Mj)
            // This vector contains [Sum_k0, Sum_k1, ..., Sum_k99]
            vector<double> current_sums = bin_mj_replica_sums[b][m];

            // 2. Get Nominal Sum (Index 0)
            double nom_sum = current_sums[0];
            if (nom_sum == 0) nom_sum = 1.0; // Safety

            // 3. Calculate Ratios & Fill Cyan Lines (Before Sorting)
            // We fill h_reps[k] with the ratio of the k-th universe
            for(int k=0; k<100; ++k) {
                double ratio = current_sums[k] / nom_sum;
                h_reps[k]->SetBinContent(m+1, ratio);
            }

            // 4. Sort to find Envelope (CG Method Logic)
            // We sort the raw sums (or ratios) to find 16th/84th percentiles
            std::sort(current_sums.begin(), current_sums.end());

            double val_16 = current_sums[15]; // 16th value
            double val_84 = current_sums[83]; // 84th value

            // 5. Fill Envelope Lines (Blue)
            h_nom->SetBinContent(m+1, 1.0);
            h_down->SetBinContent(m+1, val_16 / nom_sum); // Ratio of 16th sum
            h_up->SetBinContent(m+1, val_84 / nom_sum);   // Ratio of 84th sum
        }

        // --- Drawing ---
        h_nom->GetYaxis()->SetRangeUser(0.85, 1.15);
        h_nom->GetXaxis()->SetLabelSize(0.08);
        h_nom->GetYaxis()->SetLabelSize(0.08);
        h_nom->GetYaxis()->SetNdivisions(505);
        for(int m=0; m<nMjBins; ++m) h_nom->GetXaxis()->SetBinLabel(m+1, mjLabels[m].c_str());

        h_nom->Draw("HIST"); // Axis Frame

        // Draw all 100 Replicas (Cyan)
        for(auto h : h_reps) {
            h->SetLineColor(kCyan);
            h->SetLineWidth(1);
            h->Draw("HIST SAME");
        }

        // Draw Envelope (Blue)
        h_up->SetLineColor(kBlue); h_up->SetLineWidth(2); h_up->Draw("HIST SAME");
        h_down->SetLineColor(kBlue); h_down->SetLineWidth(2); h_down->Draw("HIST SAME");

        // Draw Nominal (Black Dashed)
        h_nom->SetLineColor(kBlack); h_nom->SetLineWidth(2); h_nom->SetLineStyle(2);
        h_nom->Draw("HIST SAME");

        // Label
        TLatex latex;
        latex.SetNDC(); latex.SetTextSize(0.12);
        latex.DrawLatex(0.20, 0.82, Form("Bin %d", binNum));
    }

    // Info for empty pad
    c1->cd(13);
    TLatex info; info.SetNDC(); info.SetTextSize(0.08);
//    info.DrawLatex(0.1, 0.5, "Bin 31 Merged");

    c1->SaveAs("plot_pdf_variations_CG_mj_bin_v3.png");
    c1->SaveAs("plot_pdf_variations_CG_mj_bin_v3.pdf");

    cout << "Saved grid plots to plot_pdf_variations_CG_mj_bin_v3.png" << endl;

    for(auto h : trash_bin) delete h;
    delete c1;
    file->Close();
    return 0;
}
