// -------------------------------------------------------------------------
// Plot PDF Variations (CG Method v3 - Optimized Structure)
// Author : Gemini
//
// [Optimization]
// - Adopts the efficient "Single Event Loop" structure from v4.
// - Uses the correct Physical Binning logic from v4.
// - Maintains the "CG Method" logic (Summing Yields per Replica).
//
//  compile: g++ -o plot_pdf_variations_CG_v3.exe plot_pdf_variations_CG_v3.cpp $(root-config --cflags --glibs)
//  run: ./plot_pdf_variations_CG_v3.exe final_output.root
// -------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm> // for sort

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

// --- Binning Definition (Same as v4) ---
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
        cout << "Usage: ./plot_pdf_variations_CG_v3.exe [root_file]" << endl;
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

    // --- Data Storage for CG Method ---
    // Instead of looping bins, we store sums for ALL bins at once.
    // bin_replica_sums[binIdx][replicaIdx]
    // replicaIdx 0: Nominal Sum
    // replicaIdx 1~99: Replica Sums
    vector<vector<double>> bin_replica_sums(nBins, vector<double>(100, 0.0));

    // --- Step 1: Single Event Loop (Efficient) ---
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

        int bIdx = getIdx(binNum); // Map to 0~13
        if (bIdx == -1) continue;

        // [CG Method Logic] Accumulate Weights Directly
        // weight_vec index k corresponds to Replica k (0 is Nominal)
        if (weight_vec->size() >= 100) {
            for(int k=0; k<=100; ++k) {
                bin_replica_sums[bIdx][k] += weight_vec->at(k);
            }
        }
    }

    // --- Prepare Histograms ---
    TH1D* h_nom = new TH1D("h_nom", "", nBins, 0, nBins);
    TH1D* h_up  = new TH1D("h_up", "", nBins, 0, nBins);
    TH1D* h_down = new TH1D("h_down", "", nBins, 0, nBins);

    vector<TH1D*> h_reps_plot;
    for(int k=0; k<100; ++k) {
        h_reps_plot.push_back(new TH1D(Form("h_rep_%d", k), "", nBins, 0, nBins));
    }

    // --- Step 2: Process Accumulated Data per Bin ---
    cout << "Calculating systematic uncertainties per bin..." << endl;

    for (int b = 0; b < nBins; ++b) {
        double nom_sum = bin_replica_sums[b][0];

        // 1. Collect sums for the 100 replicas
        vector<double> replica_yields;
        for(int k=1; k<=100; ++k) {
            double rep_sum = bin_replica_sums[b][k];
            replica_yields.push_back(rep_sum);

            // Fill Cyan Lines (Raw Yield)
            h_reps_plot[k-1]->SetBinContent(b+1, rep_sum);
        }

        // 2. Sort to find Envelope (CG Method)
        std::sort(replica_yields.begin(), replica_yields.end());

        double val_16 = replica_yields[15]; // 16th percentile
        double val_84 = replica_yields[83]; // 84th percentile

        // 3. Fill Result Histograms (Yields)
        h_nom->SetBinContent(b+1, nom_sum);
        h_down->SetBinContent(b+1, val_16);
        h_up->SetBinContent(b+1, val_84);
    }

    // --- Step 3: Convert to Ratios & Draw ---
    for (int b = 1; b <= nBins; ++b) {
        double nom_val = h_nom->GetBinContent(b);
        if(nom_val == 0) nom_val = 1.0;

        h_nom->SetBinContent(b, 1.0);
        h_up->SetBinContent(b, h_up->GetBinContent(b) / nom_val);
        h_down->SetBinContent(b, h_down->GetBinContent(b) / nom_val);

        for(int k=0; k<100; ++k) {
            h_reps_plot[k]->SetBinContent(b, h_reps_plot[k]->GetBinContent(b) / nom_val);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "PDF Variations CG v3 Optimized", 1000, 600);
    c1->SetGridy();

    // Axis Settings
    h_nom->GetYaxis()->SetRangeUser(0.85, 1.15);
    h_nom->GetYaxis()->SetTitle("Ratio to Nominal");
    for(int i=0; i<nBins; ++i) h_nom->GetXaxis()->SetBinLabel(i+1, Form("Bin %d", binNumbers[i]));
    h_nom->GetXaxis()->SetLabelSize(0.04);
    h_nom->Draw("HIST");

    // Draw Replicas (Cyan)
    for(int k=0; k<100; ++k) {
        h_reps_plot[k]->SetLineColor(kCyan);
        h_reps_plot[k]->SetLineWidth(1);
        h_reps_plot[k]->Draw("HIST SAME");
    }

    // Draw Envelope (Blue)
    h_up->SetLineColor(kBlue);
    h_up->SetLineWidth(2);
    h_up->Draw("HIST SAME");

    h_down->SetLineColor(kBlue);
    h_down->SetLineWidth(2);
    h_down->Draw("HIST SAME");

    // Draw Nominal (Black)
    h_nom->SetLineColor(kBlack);
    h_nom->SetLineWidth(2);
    h_nom->SetLineStyle(2);
    h_nom->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.78, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h_nom, "Nominal", "l");
    leg->AddEntry(h_up, "PDF 68% CL (Total Yield)", "l");
    leg->AddEntry(h_reps_plot[0], "Replica Yields", "l");
    leg->Draw();

    c1->SaveAs("pdf_variations_CG_v3.png");
    c1->SaveAs("pdf_variations_CG_v3.pdf");

    cout << "Plot saved as pdf_variations_CG_v3.png" << endl;
    cout << "Used optimized single-loop structure with correct binning." << endl;

    file->Close();
    return 0;
}
