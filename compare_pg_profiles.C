/// \file compare_pg_profiles.C
/// \brief ROOT macro: overlay prompt-gamma depth profiles from several
/// PG_analysis.root files (one per subdirectory of a "mother" directory,
/// e.g. one per physics list and/or beam energy) on shared canvases.
///
/// Expected layout (subdirectories for several beam energies can coexist
/// in the same motherDir; `energyTags` selects which ones to compare):
///   motherDir/
///     prompt_gamma_spectra_130MeV_QGSP_BIC_HP_EMZ/PG_analysis.root
///     prompt_gamma_spectra_130MeV_QBBC/PG_analysis.root
///     prompt_gamma_spectra_70MeV_QGSP_BIC_HP_EMZ/PG_analysis.root
///     ...
///
/// Each PG_analysis.root must contain (as produced by analyze_pg_spectrum.C):
///   graphs/4p4MeV/4p4MeV_Gamma_90deg   graphs/4p4MeV/4p4MeV_Gamma_120deg
///   graphs/9p6MeV/9p6MeV_Gamma_90deg   graphs/9p6MeV/9p6MeV_Gamma_120deg
///
/// The runs to overlay come from one of two sources:
///   * auto-discovery (default): every subdirectory of motherDir whose name
///     contains ANY tag in `energyTags` (comma-separated, e.g.
///     "130MeV,70MeV") and that holds a PG_analysis.root; or
///   * an explicit `dirList` (comma-separated subdirectory names, in order)
///     — pass this when a directory does NOT follow the "<energy>MeV"
///     naming convention (e.g. a 70 MeV run named just
///     "prompt_gamma_spectra_FTFP_BERT_HP"), so auto-discovery would miss
///     it, or to control exactly which runs appear.
///
/// Styling follows the experimental reference plots: COLOUR encodes the beam
/// energy (130 MeV blue, 70 MeV black, in energyTags order) while MARKER
/// shape + line style encode the physics list, so every (energy, physics
/// list) pair is a unique combination and the two energies read at a glance.
/// The x-axis (z - d_{BP}, "effective target thickness - proton range") is
/// fixed to -35..10 mm (the 4.4 MeV / 120 deg panel zooms to -12..6 mm to
/// match the experimental reference) with a red line at 0 marking the Bragg
/// peak, x/y gridlines, and a centred nuclear-transition title (e.g.
/// "^{12}C_{4.44 -> g.s.}").
///
/// Each (gamma line, detector angle) combination is drawn on its own
/// separate square picture (one JPG each), overlaying every run — e.g. the
/// 4.4 MeV line at 90 deg with 130 MeV and 70 MeV for three physics lists
/// each = 6 histograms on one picture. With two gamma lines (4.4/9.6 MeV)
/// and two angles (90/120 deg) this yields four pictures. Each has a legend
/// (labelled by subdirectory name, common prefix stripped). Output (the
/// JPGs + one ROOT file with all canvases/graphs), named after
/// `energyTags`, is written to `motherDir`.
///
/// Usage:
///   root -l -b -q 'compare_pg_profiles.C("/path/to/motherDir", "130MeV,70MeV")'
///   // explicit list (include a dir that lacks the "70MeV" tag):
///   root -l -b -q 'compare_pg_profiles.C("/motherDir", "130MeV,70MeV", "prompt_gamma_spectra_130MeV_FTFP_BERT_HP,prompt_gamma_spectra_130MeV_QBBC,prompt_gamma_spectra_130MeV_QGSP_BIC_HP_EMZ,prompt_gamma_spectra_70MeV_QBBC,prompt_gamma_spectra_70MeV_QGSP_BIC_HP_EMZ,prompt_gamma_spectra_FTFP_BERT_HP")'

#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"

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
    const char* title;     // physics title shown centred on top of the picture
};

static const LineCfg kLines[] = {
    { "4p4MeV", "4.4 MeV", "4p4MeV", "^{12}C_{4.44 #rightarrow g.s.}" },
    { "9p6MeV", "9.6 MeV", "9p6MeV", "9.6 MeV prompt #gamma" },
};
static const Int_t kNLines = sizeof(kLines) / sizeof(kLines[0]);

static const Int_t kAngles[] = { 90, 120 };
static const Int_t kNAngles = sizeof(kAngles) / sizeof(kAngles[0]);

// x-axis window (z - d_{BP}, mm), matching the experimental reference plot.
static const Double_t kXmin = -35.;
static const Double_t kXmax =  10.;

// Single-number rescale applied to every plotted y-value. analyze_pg_spectrum.C
// stores the yield per single proton (proton^-1); the experimental reference
// (Kelleter et al.) plots it per 10^9 protons, so 1e9 overlays the two on the
// same axis. Set to 1.0 to keep the raw per-proton scale. This is a pure
// display factor: the profile shape and the Bragg-peak fall-off are unchanged.
static const Double_t kYScale = 1.0e9;

// Colour encodes the beam ENERGY (experimental-reference style: 130 MeV
// blue, 70 MeV black), assigned in the order the tags appear in energyTags.
static const Int_t kEnergyColors[] = { kBlue, kBlack, kRed + 1, kGreen + 2 };
static const Int_t kNEnergyColors = sizeof(kEnergyColors) / sizeof(kEnergyColors[0]);

// Marker shape + line style encode the PHYSICS LIST (open markers for the
// clean reference look); the two arrays are indexed together.
static const Int_t kListMarkers[]    = { 24, 25, 26, 32, 27, 28, 30, 5 };
static const Int_t kListLineStyles[] = {  1,  2,  9,  7,  3,  5,  6, 8 };
static const Int_t kNListStyles = sizeof(kListMarkers) / sizeof(kListMarkers[0]);

struct RunEntry {
    std::string label;     // legend label (subdirectory name, prefix stripped)
    std::string path;      // full path to PG_analysis.root
    std::string energyTag; // which energyTags entry matched (selects colour)
    std::string physList;  // label with energyTag removed (selects marker/line)
    Int_t       color     = kBlack; // by beam energy
    Int_t       marker    = 24;     // by physics list
    Int_t       lineStyle = 1;      // by physics list
};

// ================================================================
//  1. ApplyGlobalStyle
// ================================================================
void ApplyGlobalStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(kTRUE);   // x and y gridlines
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetGridStyle(3);      // dotted, light grey
    gStyle->SetGridColor(kGray);
    gStyle->SetPadTickX(1);       // ticks on all four frame sides
    gStyle->SetPadTickY(1);
    gStyle->SetEndErrorSize(3);
}

// ================================================================
//  2. CollectRunDirs
//     Every immediate subdirectory of motherDir whose name contains ANY
//     of `energyTags` and that holds a PG_analysis.root becomes one
//     overlay entry. Pass an empty energyTags to match every subdirectory.
// ================================================================
std::vector<std::string> SplitTags(const std::string& tagList)
{
    std::vector<std::string> tags;
    std::size_t start = 0;
    while (start <= tagList.size()) {
        std::size_t comma = tagList.find(',', start);
        std::string tag = tagList.substr(start, comma - start);
        std::size_t b = tag.find_first_not_of(" \t"); // trim surrounding spaces
        std::size_t e = tag.find_last_not_of(" \t");
        tag = (b == std::string::npos) ? std::string() : tag.substr(b, e - b + 1);
        if (!tag.empty()) tags.push_back(tag);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return tags;
}

std::vector<RunEntry> CollectRunDirs(const char* motherDir, const std::vector<std::string>& energyTags)
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
        std::string matchedTag;
        if (!energyTags.empty()) {
            bool matches = false;
            for (const auto& tag : energyTags)
                if (name.find(tag) != std::string::npos) { matches = true; matchedTag = tag; break; }
            if (!matches) continue;
        }

        std::string subDir  = std::string(motherDir) + "/" + name;
        std::string rootFile = subDir + "/PG_analysis.root";
        if (gSystem->AccessPathName(rootFile.c_str()) != 0) continue; // not found

        RunEntry re;
        re.label     = name;
        re.path      = rootFile;
        re.energyTag = matchedTag;
        runs.push_back(re);
    }
    gSystem->FreeDirectory(dh);
    std::sort(runs.begin(), runs.end(),
              [](const RunEntry& a, const RunEntry& b) { return a.label < b.label; });
    return runs;
}

// Build runs from an explicit, ordered list of subdirectory names (relative
// to motherDir), keeping the given order. Use this when a directory does not
// carry an "<energy>MeV" tag and so would be skipped by CollectRunDirs. The
// energy tag is still detected from each name (when present) for line style.
std::vector<RunEntry> CollectRunsFromList(const char* motherDir,
                                          const std::vector<std::string>& dirNames,
                                          const std::vector<std::string>& energyTags)
{
    std::vector<RunEntry> runs;
    for (const auto& name : dirNames) {
        std::string rootFile = std::string(motherDir) + "/" + name + "/PG_analysis.root";
        if (gSystem->AccessPathName(rootFile.c_str()) != 0) {
            ::Warning("CollectRunsFromList", "Skipping \"%s\": no PG_analysis.root", name.c_str());
            continue;
        }
        RunEntry re;
        re.label = name;
        re.path  = rootFile;
        for (const auto& tag : energyTags)
            if (name.find(tag) != std::string::npos) { re.energyTag = tag; break; }
        runs.push_back(re);
    }
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
//  4. AssignPlotStyles
//     Colour encodes the beam energy (energyTags order: 130 MeV blue,
//     70 MeV black) so the two energies read at a glance like the
//     experimental reference. Marker shape + line style encode the physics
//     list (label with the energy tag removed), so e.g. QBBC at 130 and
//     70 MeV share a marker/line but differ in colour. Every (energy,
//     physics list) pair is therefore a unique colour+marker combination.
// ================================================================
void AssignPlotStyles(std::vector<RunEntry>& runs, const std::vector<std::string>& tags)
{
    // Derive each run's physics list and collect the distinct ones.
    std::vector<std::string> lists;
    for (auto& r : runs) {
        std::string key = r.label;
        if (!r.energyTag.empty()) {
            std::size_t pos = key.find(r.energyTag);
            if (pos != std::string::npos) key.erase(pos, r.energyTag.size());
        }
        while (!key.empty() && key.front() == '_') key.erase(key.begin());
        while (!key.empty() && key.back()  == '_') key.pop_back();
        r.physList = key.empty() ? r.label : key;
        if (std::find(lists.begin(), lists.end(), r.physList) == lists.end())
            lists.push_back(r.physList);
    }
    std::sort(lists.begin(), lists.end());

    for (auto& r : runs) {
        // Colour by beam energy (position in energyTags; grey if untagged).
        Int_t eIdx = -1;
        for (std::size_t i = 0; i < tags.size(); ++i)
            if (tags[i] == r.energyTag) { eIdx = (Int_t)i; break; }
        r.color = (eIdx >= 0) ? kEnergyColors[eIdx % kNEnergyColors] : (kGray + 2);

        // Marker + line style by physics list.
        auto it = std::find(lists.begin(), lists.end(), r.physList);
        Int_t lIdx = (Int_t)std::distance(lists.begin(), it);
        r.marker    = kListMarkers[lIdx % kNListStyles];
        r.lineStyle = kListLineStyles[lIdx % kNListStyles];
    }
}

// ================================================================
//  5. DrawOverlayPad
//     Draw `graphs` (already loaded/cloned) overlaid with a legend into
//     the current pad, for a given gamma line and detector angle.
// ================================================================
void DrawOverlayPad(const LineCfg& cfg,
                     Int_t angle,
                     const std::vector<TGraph*>& graphs,
                     const std::vector<std::string>& labels,
                     const std::vector<Int_t>& colors,
                     const std::vector<Int_t>& markers,
                     const std::vector<Int_t>& lineStyles,
                     Double_t xmin,
                     Double_t xmax)
{
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.10);
    gPad->SetRightMargin(0.05);

    // y-scale from the points inside the visible x-window only, so a wide
    // energy (e.g. 130 MeV reaching far upstream) does not squash the rest.
    Double_t yMax = -1e300;
    for (TGraph* g : graphs) {
        for (Int_t i = 0; i < g->GetN(); ++i) {
            Double_t x = g->GetX()[i];
            if (x < xmin || x > xmax) continue;
            if (g->GetY()[i] > yMax) yMax = g->GetY()[i];
        }
    }
    if (yMax <= 0.) yMax = 1.;

    // Legend top-left (as in the reference); the data rises toward the peak
    // on the right, leaving this corner free. One column, thin border. The
    // header counts as one extra row (+1) when sizing the box.
    Int_t nEntries = (Int_t)graphs.size();
    TLegend* leg = new TLegend(0.18, 0.87 - 0.042 * (nEntries + 1), 0.55, 0.87);
    leg->SetTextSize(0.024);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->SetNColumns(1);
    leg->SetHeader(Form("%d#circ detector", angle));

    // Top margin above the in-window peak, generous enough that the upstream
    // (left-hand) rise stays clear of the top-left legend box.
    Double_t headroom = 1.35;

    for (std::size_t i = 0; i < graphs.size(); ++i) {
        TGraph* g = graphs[i];
        g->SetLineColor(colors[i]);
        g->SetMarkerColor(colors[i]);
        g->SetMarkerStyle(markers[i]);
        g->SetMarkerSize(1.1);
        g->SetLineWidth(2);
        g->SetLineStyle(lineStyles[i]);

        g->Draw(i == 0 ? "APL" : "PL SAME");
        if (i == 0) {
            g->GetXaxis()->SetTitle("effective target thickness - proton range (mm)");
            g->GetYaxis()->SetTitle(kYScale == 1.0e9
                ? "#varepsilon #cdot N_{#gamma} / (FOV #cdot #Delta#Omega) / 10^{9} protons  (mm^{-1} sr^{-1})"
                : "#varepsilon #cdot N_{#gamma} / (FOV #cdot #Delta#Omega)  (proton^{-1} mm^{-1} sr^{-1})");
            g->GetXaxis()->SetTitleSize(0.042);
            g->GetYaxis()->SetTitleSize(0.042);
            g->GetXaxis()->SetLimits(xmin, xmax);
            g->GetYaxis()->SetRangeUser(0., yMax * headroom);
        }
        leg->AddEntry(g, labels[i].c_str(), "lp");
    }

    // Red vertical line at x = 0 marking the Bragg peak (z - d_{BP} = 0).
    TLine* peakLine = new TLine(0., 0., 0., yMax * headroom);
    peakLine->SetLineColor(kRed);
    peakLine->SetLineWidth(2);
    peakLine->Draw();

    leg->Draw();
    gPad->RedrawAxis();

    // Physics title centred on top, e.g. "^{12}C_{4.44 -> g.s.}".
    TLatex lat;
    lat.SetNDC();
    lat.SetTextAlign(23); // centred, top
    lat.SetTextSize(0.052);
    lat.DrawLatex(0.55, 0.98, cfg.title);
}

// ================================================================
//  Main entry point
// ================================================================
void compare_pg_profiles(const char* motherDir = "./",
                         const char* energyTags = "130MeV,70MeV",
                         const char* dirList = "")
{
    ApplyGlobalStyle();

    std::vector<std::string> tags = SplitTags(energyTags);

    // Two ways to pick the runs to overlay: an explicit comma-separated list
    // of subdirectory names (dirList), or auto-discovery by energyTags. Use
    // dirList when a directory lacks an "<energy>MeV" tag (so auto-discovery
    // would miss it), or to control exactly which runs / what order appear.
    std::vector<std::string> dirNames = SplitTags(dirList);
    std::vector<RunEntry> runs = dirNames.empty()
        ? CollectRunDirs(motherDir, tags)
        : CollectRunsFromList(motherDir, dirNames, tags);
    if (runs.empty()) {
        ::Error("compare_pg_profiles", "No PG_analysis.root found under %s (%s \"%s\")",
                 motherDir, dirNames.empty() ? "energyTags" : "dirList",
                 dirNames.empty() ? energyTags : dirList);
        return;
    }
    StripCommonPrefix(runs);
    AssignPlotStyles(runs, tags);

    Printf("Found %d run(s) for %s:", (Int_t)runs.size(), energyTags);
    for (auto& r : runs) Printf("  %-30s  ->  %s", r.label.c_str(), r.path.c_str());

    std::string tagForFile = energyTags;
    std::replace(tagForFile.begin(), tagForFile.end(), ',', '_');
    std::string tagSuffix = !tagForFile.empty() ? (std::string("_") + tagForFile) : "";
    std::string outPath = std::string(motherDir) + "/PG_comparison" + tagSuffix + ".root";
    TFile* fOut = TFile::Open(outPath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        ::Error("compare_pg_profiles", "Cannot create: %s", outPath.c_str());
        return;
    }
    TDirectory* dGraphs = fOut->mkdir("graphs");
    TDirectory* dLines[kNLines];
    for (Int_t iLine = 0; iLine < kNLines; ++iLine)
        dLines[iLine] = dGraphs->mkdir(kLines[iLine].dirTag);

    // One separate square picture per (gamma line, detector angle) combo,
    // each overlaying every matching run (e.g. 130 and 70 MeV x several
    // physics lists = up to 6 curves).
    const Int_t kPadSize = 650; // square picture
    for (Int_t iLine = 0; iLine < kNLines; ++iLine) {
        const LineCfg& cfg = kLines[iLine];

        for (Int_t iAngle = 0; iAngle < kNAngles; ++iAngle) {
            Int_t angle = kAngles[iAngle];

            std::vector<TGraph*> graphs;
            std::vector<std::string> labels;
            std::vector<Int_t> colors;
            std::vector<Int_t> markers;
            std::vector<Int_t> lineStyles;

            for (const auto& run : runs) {
                TFile* f = TFile::Open(run.path.c_str(), "READ");
                if (!f || f->IsZombie()) {
                    ::Warning("compare_pg_profiles", "Cannot open %s", run.path.c_str());
                    continue;
                }
                std::string gname = std::string("graphs/") + cfg.dirTag + "/" +
                                     cfg.dirTag + Form("_Gamma_%ddeg", angle);
                TGraph* g = (TGraph*)f->Get(gname.c_str());
                if (!g) {
                    ::Warning("compare_pg_profiles", "Missing %s in %s",
                              gname.c_str(), run.path.c_str());
                    f->Close();
                    delete f;
                    continue;
                }
                TGraph* gClone = (TGraph*)g->Clone();
                // Rescale the per-proton yield onto the experimental "per 10^9
                // protons" axis (kYScale). Shape is untouched.
                for (Int_t ip = 0; ip < gClone->GetN(); ++ip)
                    gClone->GetY()[ip] *= kYScale;
                graphs.push_back(gClone);
                labels.push_back(run.label);
                colors.push_back(run.color);
                markers.push_back(run.marker);
                lineStyles.push_back(run.lineStyle);
                f->Close();
                delete f;
            }

            if (graphs.empty()) continue;

            TCanvas* c = new TCanvas(Form("c_compare_%s_%ddeg", cfg.fileTag, angle),
                                      Form("PG comparison – %s, %d deg", cfg.label, angle),
                                      kPadSize, kPadSize);
            // The experimental reference zooms the 4.4 MeV / 120 deg panel to
            // [-12, 6] mm; every other (gamma line, angle) keeps the full
            // -35..10 mm window (kXmin/kXmax).
            Double_t padXmin = kXmin, padXmax = kXmax;
            if (std::string(cfg.dirTag) == "4p4MeV" && angle == 120) {
                padXmin = -12.;
                padXmax =   6.;
            }

            c->cd();
            DrawOverlayPad(cfg, angle, graphs, labels, colors, markers, lineStyles,
                           padXmin, padXmax);

            TDirectory* dAngle = dLines[iLine]->mkdir(Form("%ddeg", angle));
            dAngle->cd();
            for (std::size_t i = 0; i < graphs.size(); ++i) {
                graphs[i]->SetName(Form("%s_%s", cfg.dirTag, labels[i].c_str()));
                graphs[i]->Write();
            }
            fOut->cd();

            std::string jpgPath = std::string(motherDir) + "/PG_comparison" + tagSuffix +
                                  "_" + cfg.fileTag + Form("_%ddeg", angle) + ".jpg";
            c->SaveAs(jpgPath.c_str());
            c->Write();
        }
    }

    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    Printf("\nComparison canvases and graphs written to: %s", outPath.c_str());
}
