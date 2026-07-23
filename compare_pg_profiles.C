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
/// `energyTags` is a comma-separated list (e.g. "130MeV,70MeV"); every
/// subdirectory whose name contains ANY of the listed tags is included, so
/// runs from different beam energies land in the same legend/picture (the
/// x-axis, z - d_{BP}, already lines their Bragg peaks up at 0). Pass a
/// single tag (e.g. "130MeV") to compare only that energy across physics
/// lists, as before. Colour/marker are assigned per physics list (the
/// label with its energyTag removed) and shared across energies; the
/// energyTag instead selects the line style (1st tag solid, 2nd dashed,
/// ...), so e.g. QBBC at 130 and 70 MeV are drawn in the same colour, one
/// solid and one dashed. Each (gamma line, detector angle) combination is
/// drawn on its own separate square picture (one JPG each), overlaying
/// every matching run — e.g. the 4.4 MeV line at 90 deg with 130 MeV and
/// 70 MeV for three physics lists each = 6 histograms on one picture. With
/// two gamma lines (4.4/9.6 MeV) and two angles (90/120 deg) this yields
/// four pictures. Each has a legend (labelled by subdirectory name, common
/// prefix stripped). Output (the JPGs + one ROOT file with all
/// canvases/graphs), named after `energyTags`, is written to `motherDir`.
///
/// Usage:
///   root -l -b -q 'compare_pg_profiles.C("/path/to/motherDir", "130MeV,70MeV")'

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
    { "4p4MeV", "4.4 MeV", "4p4MeV" },
    { "9p6MeV", "9.6 MeV", "9p6MeV" },
};
static const Int_t kNLines = sizeof(kLines) / sizeof(kLines[0]);

static const Int_t kAngles[] = { 90, 120 };
static const Int_t kNAngles = sizeof(kAngles) / sizeof(kAngles[0]);

static const Int_t kColors[] = {
    kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
    kOrange + 7, kCyan + 2, kAzure + 1, kViolet + 1
};
static const Int_t kMarkers[] = { 20, 21, 22, 23, 33, 34, 47, 43 };
static const Int_t kNStyles = sizeof(kColors) / sizeof(kColors[0]);

// Line style per energyTags entry (1=solid, 2=dashed, 3=dotted, 4=dash-dot),
// so e.g. the 1st tag ("130MeV") is always drawn solid and the 2nd
// ("70MeV") always dashed, regardless of colour.
static const Int_t kLineStyles[] = { 1, 2, 3, 4 };
static const Int_t kNLineStyles = sizeof(kLineStyles) / sizeof(kLineStyles[0]);

struct RunEntry {
    std::string label;     // legend label (subdirectory name, prefix stripped)
    std::string path;      // full path to PG_analysis.root
    std::string energyTag; // which energyTags entry matched (selects line style)
    std::string styleKey;  // label with energyTag removed (selects colour/marker)
    Int_t       colorIdx  = 0; // index into kColors/kMarkers, shared across energies
    Int_t       lineStyle = 1; // solid by default
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
//     Group runs by "styleKey" (their label with the energyTag substring
//     removed, e.g. "130MeV_QBBC" and "70MeV_QBBC" both key to "QBBC") and
//     give each distinct key its own colour/marker. The energyTag instead
//     selects the line style, so e.g. QBBC at 130 MeV and 70 MeV are drawn
//     in the same colour, one solid and one dashed.
// ================================================================
void AssignPlotStyles(std::vector<RunEntry>& runs, const std::vector<std::string>& tags)
{
    for (auto& r : runs) {
        std::string key = r.label;
        if (!r.energyTag.empty()) {
            std::size_t pos = key.find(r.energyTag);
            if (pos != std::string::npos) key.erase(pos, r.energyTag.size());
        }
        while (!key.empty() && key.front() == '_') key.erase(key.begin());
        while (!key.empty() && key.back()  == '_') key.pop_back();
        r.styleKey = key.empty() ? r.label : key;
    }

    std::vector<std::string> styleKeys;
    for (const auto& r : runs) {
        if (std::find(styleKeys.begin(), styleKeys.end(), r.styleKey) == styleKeys.end())
            styleKeys.push_back(r.styleKey);
    }
    std::sort(styleKeys.begin(), styleKeys.end());

    for (auto& r : runs) {
        auto it = std::find(styleKeys.begin(), styleKeys.end(), r.styleKey);
        r.colorIdx = (Int_t)std::distance(styleKeys.begin(), it) % kNStyles;

        r.lineStyle = 1; // solid fallback (e.g. no energyTags filter was used)
        for (std::size_t i = 0; i < tags.size(); ++i) {
            if (tags[i] == r.energyTag) {
                r.lineStyle = kLineStyles[i % kNLineStyles];
                break;
            }
        }
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
                     const std::vector<Int_t>& colorIdx,
                     const std::vector<Int_t>& lineStyles)
{
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.10);

    Double_t yMax = -1e300;
    for (TGraph* g : graphs) {
        for (Int_t i = 0; i < g->GetN(); ++i)
            if (g->GetY()[i] > yMax) yMax = g->GetY()[i];
    }

    Int_t nCols = (Int_t)graphs.size() > 2 ? 2 : 1;
    Int_t nRows = ((Int_t)graphs.size() + nCols - 1) / nCols;
    TLegend* leg = new TLegend(0.55, 0.89 - 0.05 * nRows, 0.89, 0.89);
    leg->SetTextSize(0.028);
    leg->SetBorderSize(0);
    leg->SetNColumns(nCols);

    for (std::size_t i = 0; i < graphs.size(); ++i) {
        TGraph* g = graphs[i];
        Int_t style = colorIdx[i];
        g->SetLineColor(kColors[style]);
        g->SetMarkerColor(kColors[style]);
        g->SetMarkerStyle(kMarkers[style]);
        g->SetMarkerSize(0.9);
        g->SetLineWidth(2);
        g->SetLineStyle(lineStyles[i]);

        g->Draw(i == 0 ? "APL" : "PL SAME");
        if (i == 0) {
            g->GetXaxis()->SetTitle("z - d_{BP}  [mm]");
            g->GetYaxis()->SetTitle("#varepsilon #cdot N_{#gamma} / (FOV #cdot #Delta#Omega)  [proton^{-1} mm^{-1} sr^{-1}]");
            g->GetXaxis()->SetTitleSize(0.045);
            g->GetYaxis()->SetTitleSize(0.045);
            g->GetYaxis()->SetRangeUser(0., yMax * 1.15);
        }
        leg->AddEntry(g, labels[i].c_str(), "lp");
    }
    leg->Draw();
    gPad->RedrawAxis();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.DrawLatex(0.15, 0.93, Form("%s – %d#circ detector", cfg.label, angle));
}

// ================================================================
//  Main entry point
// ================================================================
void compare_pg_profiles(const char* motherDir = "./", const char* energyTags = "130MeV,70MeV")
{
    ApplyGlobalStyle();

    std::vector<std::string> tags = SplitTags(energyTags);
    std::vector<RunEntry> runs = CollectRunDirs(motherDir, tags);
    if (runs.empty()) {
        ::Error("compare_pg_profiles", "No PG_analysis.root found under %s matching \"%s\"",
                 motherDir, energyTags);
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
            std::vector<Int_t> colorIdx;
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
                graphs.push_back((TGraph*)g->Clone());
                labels.push_back(run.label);
                colorIdx.push_back(run.colorIdx);
                lineStyles.push_back(run.lineStyle);
                f->Close();
                delete f;
            }

            if (graphs.empty()) continue;

            TCanvas* c = new TCanvas(Form("c_compare_%s_%ddeg", cfg.fileTag, angle),
                                      Form("PG comparison – %s, %d deg", cfg.label, angle),
                                      kPadSize, kPadSize);
            c->cd();
            DrawOverlayPad(cfg, angle, graphs, labels, colorIdx, lineStyles);

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
