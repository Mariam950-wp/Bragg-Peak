/// \file compare_pg_profiles.C
/// \brief ROOT macro: overlay prompt-gamma depth profiles from several
/// PG_analysis.root files (one per subdirectory of a "mother" directory,
/// e.g. one per physics list) on shared canvases.
///
/// Expected layout (subdirectories for several beam energies can coexist
/// in the same motherDir; `energyTag` selects which ones to compare):
///   motherDir/
///     prompt_gamma_spectra_130MeV_QGSP_BIC_HP_EMZ/PG_analysis.root
///     prompt_gamma_spectra_130MeV_QBBC/PG_analysis.root
///     prompt_gamma_spectra_100MeV_QGSP_BIC_HP_EMZ/PG_analysis.root
///     ...
///
/// Each PG_analysis.root must contain (as produced by analyze_pg_spectrum.C):
///   graphs/4p4MeV/4p4MeV_Gamma_90deg
///   graphs/6p13MeV/6p13MeV_Gamma_90deg
///   graphs/9p6MeV/9p6MeV_Gamma_90deg
///
/// Only subdirectories whose name contains `energyTag` (e.g. "130MeV") are
/// used. For each gamma line, the matching subdirectories' graphs are drawn
/// on one canvas with a legend (labelled by subdirectory name, common prefix
/// stripped). Output (canvases as PNG + one ROOT file with all
/// canvases/graphs), named after `energyTag`, is written to `motherDir`.
///
/// Usage:
///   root -l -b -q 'compare_pg_profiles.C("/path/to/motherDir", "130MeV")'

#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"

#include <vector>
#include <string>
#include <algorithm>

// ================================================================
//  Configuration: one entry per gamma line to overlay
// ================================================================
struct LineCfg {
    const char* dirTag;    // e.g. "4p4MeV"
    const char* label;     // e.g. "4.4 MeV"
    const char* fileTag;   // output file suffix
};

static const LineCfg kLines[] = {
    { "4p4MeV",  "4.4 MeV",  "4p4MeV"  },
    { "6p13MeV", "6.13 MeV", "6p13MeV" },
    { "9p6MeV",  "9.6 MeV",  "9p6MeV"  },
};
static const Int_t kNLines = sizeof(kLines) / sizeof(kLines[0]);

static const Int_t kColors[] = {
    kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
    kOrange + 7, kCyan + 2, kAzure + 1, kViolet + 1
};
static const Int_t kMarkers[] = { 20, 21, 22, 23, 33, 34, 47, 43 };
static const Int_t kNStyles = sizeof(kColors) / sizeof(kColors[0]);

struct RunEntry {
    std::string label;    // legend label (subdirectory name, prefix stripped)
    std::string path;     // full path to PG_analysis.root
};

// ================================================================
//  1. ApplyGlobalStyle
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
//  2. CollectRunDirs
//     Every immediate subdirectory of motherDir whose name contains
//     `energyTag` and that holds a PG_analysis.root becomes one overlay
//     entry. Pass an empty energyTag to match every subdirectory.
// ================================================================
std::vector<RunEntry> CollectRunDirs(const char* motherDir, const char* energyTag)
{
    std::vector<RunEntry> runs;
    void* dh = gSystem->OpenDirectory(motherDir);
    if (!dh) {
        ::Error("CollectRunDirs", "Cannot open directory: %s", motherDir);
        return runs;
    }
    const char* entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dh)) != nullptr) {
        std::string name(entry);
        if (name == "." || name == "..") continue;
        if (energyTag[0] != '\0' && name.find(energyTag) == std::string::npos) continue;

        std::string subDir  = std::string(motherDir) + "/" + name;
        std::string rootFile = subDir + "/PG_analysis.root";
        if (gSystem->AccessPathName(rootFile.c_str()) != 0) continue; // not found

        RunEntry re;
        re.label = name;
        re.path  = rootFile;
        runs.push_back(re);
    }
    gSystem->FreeDirectory(dh);
    std::sort(runs.begin(), runs.end(),
              [](const RunEntry& a, const RunEntry& b) { return a.label < b.label; });
    return runs;
}

// ================================================================
//  3. StripCommonPrefix
//     Shortens legend labels by removing the prefix shared by all runs
//     (e.g. "prompt_gamma_spectra_130MeV_"), so only the distinguishing
//     part (e.g. "QGSP_BIC_HP_EMZ") is shown.
// ================================================================
void StripCommonPrefix(std::vector<RunEntry>& runs)
{
    if (runs.size() < 2) return;

    std::string prefix = runs.front().label;
    for (const auto& r : runs) {
        std::size_t n = 0;
        std::size_t maxN = std::min(prefix.size(), r.label.size());
        while (n < maxN && prefix[n] == r.label[n]) ++n;
        prefix.resize(n);
        if (prefix.empty()) return;
    }
    for (auto& r : runs) {
        std::string stripped = r.label.substr(prefix.size());
        if (!stripped.empty()) r.label = stripped;
    }
}

// ================================================================
//  4. DrawOverlayCanvas
//     Build one canvas overlaying `graphs` (already loaded/cloned) with
//     a legend, for a given gamma line.
// ================================================================
TCanvas* DrawOverlayCanvas(const LineCfg& cfg,
                            const std::vector<TGraph*>& graphs,
                            const std::vector<std::string>& labels)
{
    TCanvas* c = new TCanvas(Form("c_compare_%s", cfg.fileTag),
                              Form("%s comparison", cfg.label), 900, 650);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);

    Double_t yMax = -1e300;
    for (TGraph* g : graphs) {
        for (Int_t i = 0; i < g->GetN(); ++i)
            if (g->GetY()[i] > yMax) yMax = g->GetY()[i];
    }

    TLegend* leg = new TLegend(0.55, 0.72, 0.90, 0.90);
    leg->SetBorderSize(0);

    for (std::size_t i = 0; i < graphs.size(); ++i) {
        TGraph* g = graphs[i];
        Int_t style = (Int_t)(i % kNStyles);
        g->SetLineColor(kColors[style]);
        g->SetMarkerColor(kColors[style]);
        g->SetMarkerStyle(kMarkers[style]);
        g->SetMarkerSize(0.9);
        g->SetLineWidth(2);

        g->Draw(i == 0 ? "APL" : "PL SAME");
        if (i == 0) {
            g->GetXaxis()->SetTitle("Depth / d_{BP}");
            g->GetYaxis()->SetTitle("#varepsilon #cdot N_{#gamma} / (FOV #cdot #Delta#Omega)  [proton^{-1} mm^{-1} sr^{-1}]");
            g->GetXaxis()->SetTitleSize(0.05);
            g->GetYaxis()->SetTitleSize(0.05);
            g->GetYaxis()->SetRangeUser(0., yMax * 1.1);
        }
        leg->AddEntry(g, labels[i].c_str(), "lp");
    }
    leg->Draw();
    gPad->RedrawAxis();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.14, 0.92, Form("Prompt gammas – %s – 90#circ detector", cfg.label));

    return c;
}

// ================================================================
//  Main entry point
// ================================================================
void compare_pg_profiles(const char* motherDir = "./", const char* energyTag = "130MeV")
{
    ApplyGlobalStyle();

    std::vector<RunEntry> runs = CollectRunDirs(motherDir, energyTag);
    if (runs.empty()) {
        ::Error("compare_pg_profiles", "No PG_analysis.root found under %s matching \"%s\"",
                 motherDir, energyTag);
        return;
    }
    StripCommonPrefix(runs);

    Printf("Found %d run(s) for %s:", (Int_t)runs.size(), energyTag);
    for (auto& r : runs) Printf("  %-30s  ->  %s", r.label.c_str(), r.path.c_str());

    std::string tagSuffix = (energyTag[0] != '\0') ? (std::string("_") + energyTag) : "";
    std::string outPath = std::string(motherDir) + "/PG_comparison" + tagSuffix + ".root";
    TFile* fOut = TFile::Open(outPath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        ::Error("compare_pg_profiles", "Cannot create: %s", outPath.c_str());
        return;
    }
    TDirectory* dGraphs = fOut->mkdir("graphs");

    for (Int_t iLine = 0; iLine < kNLines; ++iLine) {
        const LineCfg& cfg = kLines[iLine];
        std::vector<TGraph*> graphs;
        std::vector<std::string> labels;

        for (const auto& run : runs) {
            TFile* f = TFile::Open(run.path.c_str(), "READ");
            if (!f || f->IsZombie()) {
                ::Warning("compare_pg_profiles", "Cannot open %s", run.path.c_str());
                continue;
            }
            std::string gname = std::string("graphs/") + cfg.dirTag + "/" +
                                 cfg.dirTag + "_Gamma_90deg";
            TGraph* g = (TGraph*)f->Get(gname.c_str());
            if (!g) {
                ::Warning("compare_pg_profiles", "Missing %s in %s",
                          gname.c_str(), run.path.c_str());
                f->Close();
                delete f;
                continue;
            }
            graphs.push_back((TGraph*)g->Clone());
            labels.push_back(run.label);
            f->Close();
            delete f;
        }
        if (graphs.empty()) continue;

        TCanvas* c = DrawOverlayCanvas(cfg, graphs, labels);

        std::string pngPath = std::string(motherDir) + "/PG_comparison" + tagSuffix +
                              "_" + cfg.fileTag + ".png";
        c->SaveAs(pngPath.c_str());

        c->Write();

        TDirectory* dLine = dGraphs->mkdir(cfg.dirTag);
        dLine->cd();
        for (std::size_t i = 0; i < graphs.size(); ++i) {
            graphs[i]->SetName(Form("%s_%s", cfg.dirTag, labels[i].c_str()));
            graphs[i]->Write();
        }
        fOut->cd();
    }

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    Printf("\nComparison canvases and graphs written to: %s", outPath.c_str());
}
