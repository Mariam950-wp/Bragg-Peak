/// \file analyze_pg_spectrum.C
/// \brief ROOT macro: prompt-gamma intensity vs proton range (depth).
///
/// Reads every PG_Spectrum_VS_Angle_<depth>.root file in `dataDir`,
/// integrates each angular histogram, and writes results to one output file:
///
///   graphs/4p4MeV/         – TGraph per detector angle for the 4.4 MeV line
///   graphs/9p6MeV/         – TGraph per detector angle for the 9.6 MeV line
///   canvases/              – two styled TCanvas objects
///
/// Usage (interactive):
///   root -l 'analyze_pg_spectrum.C("./", "PG_analysis.root")'
/// Usage (batch):
///   root -l -b -q 'analyze_pg_spectrum.C("./", "PG_analysis.root")'

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cstdio>

// ================================================================
//  Configuration: histogram names, labels, colours
// ================================================================

static const Int_t kNAngles = 9;
static const Int_t kSpecificAngles[kNAngles] = {30, 50, 60, 65, 90, 115, 120, 130, 150};
static const Color_t kAngleColors[kNAngles] = {
    kBlue+1, kRed+1, kGreen+2, kMagenta+1, kCyan+2,
    kOrange+7, kViolet+2, kTeal+3, kGray+1
};
static const Style_t kAngleMarkers[kNAngles] = {20, 21, 22, 23, 29, 34, 33, 47, 43};

// One entry per input file
struct FileEntry {
    Double_t    depth;
    std::string path;
    bool operator<(const FileEntry& o) const { return depth < o.depth; }
};

// ================================================================
//  1. ApplyGlobalStyle
//     Set ROOT style options used by all canvases.
// ================================================================
void ApplyGlobalStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(kGray);
}

// ================================================================
//  2. CollectInputFiles
//     Scan `dataDir` for PG_Spectrum_VS_Angle_*.root files,
//     parse the depth from each filename, and return them sorted.
// ================================================================
Double_t DepthFromFilename(const std::string& path)
{
    std::string base = path.substr(path.find_last_of("/\\") + 1);
    std::size_t pos  = base.rfind('_');
    if (pos == std::string::npos) return -1.;
    std::string num  = base.substr(pos + 1);
    num = num.substr(0, num.find(".root"));
    return std::stod(num);
}

std::vector<FileEntry> CollectInputFiles(const char* dataDir)
{
    std::vector<FileEntry> files;
    void* dh = gSystem->OpenDirectory(dataDir);
    if (!dh) {
        ::Error("CollectInputFiles", "Cannot open directory: %s", dataDir);
        return files;
    }
    const char* entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dh)) != nullptr) {
        std::string name(entry);
        if (name.find("PG_Spectrum_VS_Angle_") == std::string::npos) continue;
        if (name.find(".root")                  == std::string::npos) continue;
        FileEntry fe;
        fe.path  = std::string(dataDir) + "/" + name;
        fe.depth = DepthFromFilename(fe.path);
        if (fe.depth < 0) continue;
        files.push_back(fe);
    }
    gSystem->FreeDirectory(dh);
    std::sort(files.begin(), files.end());
    return files;
}

// ================================================================
//  3. FillGammaLineGraphs
//     For each depth file, integrate the specific gamma-line spectra
//     (one histogram per detector angle) and fill one TGraph per angle.
//     `energyTag` selects the histogram family, e.g. "4.400000" or "9.600000".
//     `graphTag`  is used to build object names, e.g. "4p4MeV" or "9p6MeV".
//     GetSumOfWeights() returns the raw sum of bin contents with no
//     per-entry normalisation; the histograms are already normalised
//     per primary particle by the simulation.
// ================================================================
void FillGammaLineGraphs(const std::vector<FileEntry>& files,
                          TGraph* graphs[],
                          const char* energyTag,
                          const char* graphTag)
{
    const Int_t nFiles = (Int_t)files.size();
    for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
        TFile* f = TFile::Open(files[iFile].path.c_str(), "READ");
        if (!f || f->IsZombie()) {
            ::Warning("FillGammaLineGraphs", "Cannot open %s",
                      files[iFile].path.c_str());
            continue;
        }
        for (Int_t j = 0; j < kNAngles; ++j) {
            char hname[128];
            std::snprintf(hname, sizeof(hname), "%s_MeV_Gamma_%d_deg",
                          energyTag, kSpecificAngles[j]);
            TH1D* h = (TH1D*)f->Get(hname);
            graphs[j]->SetPoint(iFile, files[iFile].depth,
                                h ? h->GetSumOfWeights() : 0.);
        }
        f->Close();
        delete f;
    }
    for (Int_t j = 0; j < kNAngles; ++j) {
        graphs[j]->SetName(Form("%s_Gamma_%ddeg", graphTag, kSpecificAngles[j]));
        graphs[j]->SetTitle(Form("%s gamma %d deg;Depth (mm);Intensity (counts / primary)",
                                 graphTag, kSpecificAngles[j]));
    }
}

// ================================================================
//  4. DrawGammaLineCanvas
//     One canvas with intensity-vs-depth curves for one gamma line
//     (e.g. 4.4 MeV or 9.6 MeV), one curve per detector angle.
//     TMultiGraph is used so the Y-axis range is derived from the
//     actual data maximum and all content stays within the frame.
// ================================================================
TCanvas* DrawGammaLineCanvas(TGraph* graphs[],
                              const std::vector<FileEntry>& files,
                              const char* energyStr,
                              const char* canvasName)
{
    TCanvas* c = new TCanvas(canvasName,
                             Form("%s MeV gamma line", energyStr), 900, 600);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);

    TMultiGraph* mg = new TMultiGraph(canvasName, "");

    TLegend* leg = new TLegend(0.55, 0.50, 0.88, 0.88);
    leg->SetHeader("Detector angle", "C");
    leg->SetBorderSize(1);
    leg->SetTextSize(0.035);
    leg->SetNColumns(2);

    Double_t ymax = 0.;
    for (Int_t j = 0; j < kNAngles; ++j) {
        Bool_t nonzero = kFALSE;
        for (Int_t p = 0; p < (Int_t)files.size(); ++p) {
            ymax = std::max(ymax, graphs[j]->GetY()[p]);
            if (graphs[j]->GetY()[p] > 0.) nonzero = kTRUE;
        }
        if (!nonzero) continue;

        graphs[j]->SetLineColor(kAngleColors[j]);
        graphs[j]->SetMarkerColor(kAngleColors[j]);
        graphs[j]->SetMarkerStyle(kAngleMarkers[j]);
        graphs[j]->SetMarkerSize(1.2);
        graphs[j]->SetLineWidth(2);
        mg->Add(graphs[j], "LP");
        leg->AddEntry(graphs[j], Form("%d#circ", kSpecificAngles[j]), "lp");
    }

    if (ymax == 0.) ymax = 1e-8;
    mg->SetMinimum(0.);
    mg->SetMaximum(ymax * 1.15);
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("Depth (mm)");
    mg->GetYaxis()->SetTitle("Intensity (counts / primary)");
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetYaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetLabelSize(0.04);
    mg->GetYaxis()->SetLabelSize(0.04);
    mg->GetYaxis()->SetTitleOffset(1.3);

    leg->Draw();

    TLatex title;
    title.SetNDC();
    title.SetTextSize(0.045);
    title.SetTextAlign(22);
    title.DrawLatex(0.50, 0.96,
        Form("Prompt #gamma Intensity vs Depth (%s MeV line)", energyStr));

    // Repaint the frame border on top of all primitives so that any graph
    // line segment that reaches the axis boundary is visually clipped.
    gPad->RedrawAxis();

    return c;
}

// ================================================================
//  5. WriteResultsToFile
//     Save all TGraphs and TCanvas objects into the output ROOT file
//     under named sub-directories.
// ================================================================
void WriteResultsToFile(TFile* fOut,
                         TGraph* g44[], TGraph* g96[],
                         TCanvas* c44, TCanvas* c96)
{
    TDirectory* dGraphs = fOut->mkdir("graphs");

    TDirectory* d44 = dGraphs->mkdir("4p4MeV");
    d44->cd();
    for (Int_t j = 0; j < kNAngles; ++j) g44[j]->Write();

    TDirectory* d96 = dGraphs->mkdir("9p6MeV");
    d96->cd();
    for (Int_t j = 0; j < kNAngles; ++j) g96[j]->Write();

    TDirectory* dCan = fOut->mkdir("canvases");
    dCan->cd();
    c44->Write();
    c96->Write();

    Printf("Output structure:");
    Printf("  graphs/4p4MeV/        – %d TGraphs", kNAngles);
    Printf("  graphs/9p6MeV/        – %d TGraphs", kNAngles);
    Printf("  canvases/             – 2 TCanvas objects");
}

// ================================================================
//  Main entry point
//  Orchestrates the steps above; each step is a single function call.
// ================================================================
void analyze_pg_spectrum(const char* dataDir = "./",
                          const char* outFile = "PG_analysis.root")
{
    ApplyGlobalStyle();

    // --- collect input files ---
    std::vector<FileEntry> files = CollectInputFiles(dataDir);
    if (files.empty()) return;

    Printf("Found %d file(s):", (Int_t)files.size());
    for (auto& fe : files)
        Printf("  depth = %.1f mm  ->  %s", fe.depth, fe.path.c_str());

    // --- fill intensity-vs-depth graphs ---
    TGraph* g44[kNAngles], *g96[kNAngles];
    for (Int_t j = 0; j < kNAngles; ++j) {
        g44[j] = new TGraph((Int_t)files.size());
        g96[j] = new TGraph((Int_t)files.size());
    }

    FillGammaLineGraphs(files, g44, "4.400000", "4p4MeV");
    FillGammaLineGraphs(files, g96, "9.600000", "9p6MeV");

    // --- build canvases ---
    TCanvas* c44 = DrawGammaLineCanvas(g44, files, "4.4", "c_4p4MeV");
    TCanvas* c96 = DrawGammaLineCanvas(g96, files, "9.6", "c_9p6MeV");

    // --- write everything to the output file ---
    TFile* fOut = TFile::Open(outFile, "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        ::Error("analyze_pg_spectrum", "Cannot create: %s", outFile);
        return;
    }
    WriteResultsToFile(fOut, g44, g96, c44, c96);
    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    Printf("\nResults written to: %s", outFile);
}
