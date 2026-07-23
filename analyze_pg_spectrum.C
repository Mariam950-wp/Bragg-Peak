/// \file analyze_pg_spectrum.C
/// \brief ROOT macro: prompt-gamma intensity vs proton range (depth) at 90 and 120 degrees.
///
/// Reads every PG_Spectrum_VS_Angle_<depth>.root file in `dataDir`,
/// integrates the 90 and 120-degree detector histograms for each gamma line,
/// and writes results to one output file:
///
///   graphs/4p4MeV/         – TGraph per angle for the 4.4 MeV line
///   graphs/6p13MeV/        – TGraph per angle for the 6.13 MeV line
///   graphs/9p6MeV/         – TGraph per angle for the 9.6 MeV line
///
/// Each graph's x-axis is z - d_{BP} (mm), so the Bragg peak sits at 0.
/// d_{BP} is picked from the "130MeV"/"70MeV" tag found in `dataDir`
/// (BraggPeakDepthMm); pass braggPeakMm explicitly to override.
///
/// Usage (interactive):
///   root -l 'analyze_pg_spectrum.C("./prompt_gamma_spectra_70MeV_QBBC", "PG_analysis.root")'
/// Usage (batch):
///   root -l -b -q 'analyze_pg_spectrum.C("./prompt_gamma_spectra_70MeV_QBBC", "PG_analysis.root")'

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cstdio>

// ================================================================
//  Configuration: histogram names, labels, colours
// ================================================================

// Bragg peak depth in PMMA (mm), per beam energy used in this study.
// Looked up from `dataDir`'s name (see BraggPeakDepthMm), or overridden
// explicitly via the braggPeakMm argument.
static const Double_t kDbpMm130 = 107.0; // 130.87 MeV protons
static const Double_t kDbpMm70  = 33.0;  // 70.54 MeV protons
static const Int_t    kNAngles       = 2;
static const Int_t    kSpecificAngles[kNAngles] = {90, 120};

// Detector geometry (same FOV/solid angle used for all detector angles)
static const Double_t kFOV        = 2.0;
static const Double_t kDeltaTheta = 0.018865;
static const Double_t kDeltaOmega = 2.0 * TMath::Pi() * kDeltaTheta; // assumes sin(theta)=1, exact at 90 deg only

// HPGe detector efficiency from Kelleter et al. 2017 (Monte-Carlo determined).
static const Double_t kEff44  = 0.039;
static const Double_t kEff613 = 0.025;
static const Double_t kEff96  = 1.0;  // unknown — no experimental counterpart

// One entry per input file
struct FileEntry {
    Double_t    depth;
    std::string path;
    bool operator<(const FileEntry& o) const { return depth < o.depth; }
};

// ================================================================
//  1. CollectInputFiles
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

// Infer the Bragg peak depth from the "130MeV"/"70MeV" tag in dataDir's
// path, matching the naming convention CollectRunDirs relies on in
// compare_pg_profiles.C. Falls back to the 130 MeV value with a warning
// if neither tag is found.
Double_t BraggPeakDepthMm(const char* dataDir)
{
    std::string dir(dataDir);
    if (dir.find("130MeV") != std::string::npos) return kDbpMm130;
    if (dir.find("70MeV")  != std::string::npos) return kDbpMm70;
    ::Warning("BraggPeakDepthMm",
              "Cannot infer beam energy from \"%s\" (expected \"130MeV\" or \"70MeV\" "
              "in the path); defaulting to d_BP = %.1f mm. Pass braggPeakMm explicitly to override.",
              dataDir, kDbpMm130);
    return kDbpMm130;
}

// ================================================================
//  2. FillGammaLineGraphs
// ================================================================
void FillGammaLineGraphs(const std::vector<FileEntry>& files,
                          TGraph* graphs[],
                          const char* energyTag,
                          const char* graphTag,
                          Double_t eLo,
                          Double_t eHi,
                          Double_t efficiency,
                          Double_t braggPeakMm)
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
            Double_t yield = nGamma * efficiency / (kFOV * kDeltaOmega);
            graphs[j]->SetPoint(iFile, files[iFile].depth - braggPeakMm, yield);
        }
        f->Close();
        delete f;
    }
    for (Int_t j = 0; j < kNAngles; ++j) {
        graphs[j]->SetName(Form("%s_Gamma_%ddeg", graphTag, kSpecificAngles[j]));
        graphs[j]->SetTitle(Form("%s gamma %d deg;z - d_{BP}  [mm];#varepsilon #cdot N_{#gamma} / (FOV #cdot #Delta#Omega)  [proton^{-1} mm^{-1} sr^{-1}]",
                                 graphTag, kSpecificAngles[j]));
    }
}

// ================================================================
//  3. WriteResultsToFile
// ================================================================
void WriteResultsToFile(TFile* fOut, TGraph* g44[], TGraph* g613[], TGraph* g96[])
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

    Printf("Output structure:");
    Printf("  graphs/4p4MeV/        – %d TGraph(s) (90, 120 deg)", kNAngles);
    Printf("  graphs/6p13MeV/       – %d TGraph(s) (90, 120 deg)", kNAngles);
    Printf("  graphs/9p6MeV/        – %d TGraph(s) (90, 120 deg)", kNAngles);
}

// ================================================================
//  Main entry point
// ================================================================
void analyze_pg_spectrum(const char* dataDir = "./",
                          const char* outFile = "PG_analysis.root",
                          Double_t braggPeakMm = -1.)
{
    std::vector<FileEntry> files = CollectInputFiles(dataDir);
    if (files.empty()) return;

    Double_t dBP = (braggPeakMm > 0.) ? braggPeakMm : BraggPeakDepthMm(dataDir);
    Printf("Using d_BP = %.1f mm", dBP);

    Printf("Found %d file(s):", (Int_t)files.size());
    for (auto& fe : files)
        Printf("  depth = %.1f mm  ->  %s", fe.depth, fe.path.c_str());

    TGraph* g44[kNAngles], *g613[kNAngles], *g96[kNAngles];
    for (Int_t j = 0; j < kNAngles; ++j) {
        g44[j]  = new TGraph((Int_t)files.size());
        g613[j] = new TGraph((Int_t)files.size());
        g96[j]  = new TGraph((Int_t)files.size());
    }

    FillGammaLineGraphs(files, g44,  "4.400000",  "4p4MeV",   4.0, 50.0, kEff44,  dBP);
    FillGammaLineGraphs(files, g613, "6.130000",  "6p13MeV",  5.6,  6.3, kEff613, dBP);
    FillGammaLineGraphs(files, g96,  "9.600000",  "9p6MeV",   9.0, 10.0, kEff96,  dBP);

    std::string outPath = std::string(dataDir) + "/" + outFile;
    TFile* fOut = TFile::Open(outPath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        ::Error("analyze_pg_spectrum", "Cannot create: %s", outPath.c_str());
        return;
    }
    WriteResultsToFile(fOut, g44, g613, g96);
    fOut->Write("", TObject::kOverwrite);
    fOut->Close();

    Printf("\nResults written to: %s", outPath.c_str());
}
