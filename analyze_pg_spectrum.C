/// \file analyze_pg_spectrum.C
/// \brief ROOT macro: prompt-gamma intensity vs proton range (depth) at 90 degrees.
///
/// Reads every PG_Spectrum_VS_Angle_<depth>.root file in `dataDir`,
/// integrates the 90-degree detector histogram for each gamma line,
/// and writes results to one output file:
///
///   graphs/4p4MeV/         – TGraph for the 4.4 MeV line at 90 deg
///   graphs/6p13MeV/        – TGraph for the 6.13 MeV line at 90 deg
///   graphs/9p6MeV/         – TGraph for the 9.6 MeV line at 90 deg
///   canvases/              – one canvas: 4.4 + 9.6 MeV overlaid
///                            + one dedicated canvas for the 6.13 MeV line
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
#include "TMath.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cstdio>

// ================================================================
//  Configuration: histogram names, labels, colours
// ================================================================

// Bragg peak position for 130 MeV protons in PMMA (mm)
static const Double_t kBraggPeakMm  = 107.0;
static const Int_t    kNAngles       = 1;
static const Int_t    kSpecificAngles[kNAngles] = {90};

// Detector geometry at 90 degrees
// FOV: target slice thickness (mm)
// kDeltaTheta: polar angular range covered by detector (rad)
// kDeltaOmega = 2*pi*sin(90 deg)*DeltaTheta  (sr)
static const Double_t kFOV        = 2.0;
static const Double_t kDeltaTheta = 0.018865;
static const Double_t kDeltaOmega = 2.0 * TMath::Pi() * kDeltaTheta; // sin(90)=1

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
//     `eLo` / `eHi` define the energy integration range in MeV.
//     The integral is energy-weighted (bin content * bin centre) so the
//     result is total photon energy in MeV per primary, not a photon count.
// ================================================================
void FillGammaLineGraphs(const std::vector<FileEntry>& files,
                          TGraph* graphs[],
                          const char* energyTag,
                          const char* graphTag,
                          Double_t eLo,
                          Double_t eHi)
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
            Double_t nGamma = 0.;
            if (h) {
                TAxis* ax  = h->GetXaxis();
                Int_t  bLo = ax->FindBin(eLo);
                Int_t  bHi = ax->FindBin(eHi);
                for (Int_t ib = bLo; ib <= bHi; ++ib)
                    nGamma += h->GetBinContent(ib);
            }
            // Yield per proton per (FOV * DeltaOmega)
            Double_t yield = nGamma / (kFOV * kDeltaOmega);
            graphs[j]->SetPoint(iFile, files[iFile].depth / kBraggPeakMm, yield);
        }
        f->Close();
        delete f;
    }
    for (Int_t j = 0; j < kNAngles; ++j) {
        graphs[j]->SetName(Form("%s_Gamma_%ddeg", graphTag, kSpecificAngles[j]));
        graphs[j]->SetTitle(Form("%s gamma %d deg;Depth / d_{BP};N_{#gamma} / (FOV #cdot #Delta#Omega)  [proton^{-1} mm^{-1} sr^{-1}]",
                                 graphTag, kSpecificAngles[j]));
    }
}

// ================================================================
//  3b. GraphYMax – return the maximum Y value of a TGraph.
// ================================================================
Double_t GraphYMax(TGraph* g)
{
    Double_t yMax = -1e300;
    for (Int_t i = 0; i < g->GetN(); ++i)
        if (g->GetY()[i] > yMax) yMax = g->GetY()[i];
    return yMax;
}

// ================================================================
//  4. Draw44_96Canvas
//     Draw 4.4 MeV and 9.6 MeV normalised intensity vs depth
//     on the same canvas at 90 degrees.
// ================================================================
TCanvas* Draw44_96Canvas(TGraph* g44, TGraph* g96, Int_t angle)
{
    TCanvas* c = new TCanvas(Form("c_44_96_%ddeg", angle),
                             Form("4.4 + 9.6 MeV – %d deg", angle), 800, 600);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);

    g44->SetLineColor(kBlue + 1);
    g44->SetMarkerColor(kBlue + 1);
    g44->SetMarkerStyle(20);
    g44->SetMarkerSize(0.9);
    g44->SetLineWidth(2);

    g96->SetLineColor(kGreen + 2);
    g96->SetMarkerColor(kGreen + 2);
    g96->SetMarkerStyle(22);
    g96->SetMarkerSize(0.9);
    g96->SetLineWidth(2);

    Double_t yMax44_96 = std::max(GraphYMax(g44), GraphYMax(g96));
    g44->Draw("APL");
    g44->GetXaxis()->SetTitle("Depth / d_{BP}");
    g44->GetYaxis()->SetTitle("N_{#gamma} / (FOV #cdot #Delta#Omega)  [proton^{-1} mm^{-1} sr^{-1}]");
    g44->GetXaxis()->SetTitleSize(0.05);
    g44->GetYaxis()->SetTitleSize(0.05);
    g44->GetYaxis()->SetRangeUser(0., yMax44_96 * 1.1);
    g96->Draw("PL SAME");

    gPad->RedrawAxis();

    TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(g44, "4.4 MeV",  "lp");
    leg->AddEntry(g96, "9.6 MeV",  "lp");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.14, 0.92, Form("Prompt gammas – 4.4 + 9.6 MeV – %d#circ detector", angle));

    return c;
}

// ================================================================
//  5. Draw6p13Canvas
//     Plot the 6.13 MeV normalised prompt-gamma line at 90 degrees vs depth.
// ================================================================
TCanvas* Draw6p13Canvas(TGraph* g613[], Int_t angle)
{
    TCanvas* c = new TCanvas(Form("c_6p13MeV_%ddeg", angle),
                             Form("6.13 MeV – %d deg", angle), 900, 650);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);

    g613[0]->SetLineColor(kRed + 1);
    g613[0]->SetMarkerColor(kRed + 1);
    g613[0]->SetMarkerStyle(21);
    g613[0]->SetMarkerSize(0.9);
    g613[0]->SetLineWidth(2);
    Double_t yMax613 = GraphYMax(g613[0]);
    g613[0]->Draw("APL");
    g613[0]->GetXaxis()->SetTitle("Depth / d_{BP}");
    g613[0]->GetYaxis()->SetTitle("N_{#gamma} / (FOV #cdot #Delta#Omega)  [proton^{-1} mm^{-1} sr^{-1}]");
    g613[0]->GetXaxis()->SetTitleSize(0.05);
    g613[0]->GetYaxis()->SetTitleSize(0.05);
    g613[0]->GetYaxis()->SetRangeUser(0., yMax613 * 1.1);

    gPad->RedrawAxis();

    TLegend* leg = new TLegend(0.65, 0.78, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->AddEntry(g613[0], "6.13 MeV", "lp");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.14, 0.92, Form("Prompt gammas – 6.13 MeV – %d#circ detector", angle));

    return c;
}

// ================================================================
//  6. WriteResultsToFile
//     Save all TGraphs and canvases into the output ROOT file
//     under named sub-directories.
// ================================================================
void WriteResultsToFile(TFile* fOut,
                         TGraph* g44[], TGraph* g613[], TGraph* g96[],
                         TCanvas* c44_96, TCanvas* c6p13)
{
    TDirectory* dGraphs = fOut->mkdir("graphs");

    TDirectory* d44 = dGraphs->mkdir("4p4MeV");
    d44->cd();
    for (Int_t j = 0; j < kNAngles; ++j) g44[j]->Write();

    TDirectory* d613 = dGraphs->mkdir("6p13MeV");
    d613->cd();
    for (Int_t j = 0; j < kNAngles; ++j) g613[j]->Write();

    TDirectory* d96 = dGraphs->mkdir("9p6MeV");
    d96->cd();
    for (Int_t j = 0; j < kNAngles; ++j) g96[j]->Write();

    TDirectory* dCanv = fOut->mkdir("canvases");
    dCanv->cd();
    c44_96->Write();
    c6p13->Write();

    Printf("Output structure (90 deg only):");
    Printf("  graphs/4p4MeV/        – 1 TGraph (90 deg)");
    Printf("  graphs/6p13MeV/       – 1 TGraph (90 deg)");
    Printf("  graphs/9p6MeV/        – 1 TGraph (90 deg)");
    Printf("  canvases/             – 4.4+9.6 MeV overlay canvas + 6.13 MeV canvas");
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
    TGraph* g44[kNAngles], *g613[kNAngles], *g96[kNAngles];
    for (Int_t j = 0; j < kNAngles; ++j) {
        g44[j]  = new TGraph((Int_t)files.size());
        g613[j] = new TGraph((Int_t)files.size());
        g96[j]  = new TGraph((Int_t)files.size());
    }

    FillGammaLineGraphs(files, g44,  "4.400000",  "4p4MeV",   4.0, 50.0);
    FillGammaLineGraphs(files, g613, "6.130000",  "6p13MeV",  5.6,  6.3);
    FillGammaLineGraphs(files, g96,  "9.600000",  "9p6MeV",   9.0, 10.0);

    // --- canvas 1: 4.4 MeV + 9.6 MeV overlaid ---
    TCanvas* c44_96 = Draw44_96Canvas(g44[0], g96[0], kSpecificAngles[0]);

    // --- canvas 2: 6.13 MeV dedicated ---
    TCanvas* c6p13 = Draw6p13Canvas(g613, kSpecificAngles[0]);

    // --- write everything to the output file ---
    TFile* fOut = TFile::Open(outFile, "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        ::Error("analyze_pg_spectrum", "Cannot create: %s", outFile);
        return;
    }
    WriteResultsToFile(fOut, g44, g613, g96, c44_96, c6p13);
    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    Printf("\nResults written to: %s", outFile);
}
