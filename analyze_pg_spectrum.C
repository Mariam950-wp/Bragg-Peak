/// \file analyze_pg_spectrum.C
/// \brief ROOT macro: plot prompt-gamma intensities vs proton range (depth).
///
/// Reads every PG_Spectrum_VS_Angle_<depth>.root file found in `dataDir`,
/// integrates each angular histogram, and writes results to a single output
/// ROOT file with the following structure:
///
///   graphs/broad_angles/   – one TGraph (intensity vs depth) per angular bin
///   graphs/4p4MeV/         – one TGraph per detector angle for the 4.4 MeV line
///   graphs/9p6MeV/         – one TGraph per detector angle for the 9.6 MeV line
///   canvases/              – four TCanvas objects (fully styled, ready to view)
///
/// Usage (interactive):
///   root -l 'analyze_pg_spectrum.C("./", "PG_analysis.root")'
///
/// Usage (batch):
///   root -l -b -q 'analyze_pg_spectrum.C("./", "PG_analysis.root")'

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cstdio>

// ============================================================
//  Helper: extract depth (mm) from filename
//  e.g. "PG_Spectrum_VS_Angle_330.0.root" -> 330.0
// ============================================================
Double_t DepthFromFile(const std::string& path)
{
    std::string base = path.substr(path.find_last_of("/\\") + 1);
    // find last '_' before the number
    std::size_t pos = base.rfind('_');
    if (pos == std::string::npos) return -1.;
    std::string numStr = base.substr(pos + 1);         // e.g. "330.0.root"
    numStr = numStr.substr(0, numStr.find(".root"));   // e.g. "330.0"
    return std::stod(numStr);
}

// ============================================================
//  Colour palette – 9 distinguishable colours
// ============================================================
static const Int_t kNColours = 9;
static const Color_t kPalette[kNColours] = {
    kBlue+1, kRed+1, kGreen+2, kMagenta+1, kCyan+2,
    kOrange+7, kViolet+2, kTeal+3, kGray+1
};

// ============================================================
//  Main macro
// ============================================================
void analyze_pg_spectrum(const char* dataDir  = "./",
                         const char* outFile  = "PG_analysis.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(kGray);

    // --------------------------------------------------------
    //  1. Collect ROOT files and sort by depth
    // --------------------------------------------------------
    void* dirHandle = gSystem->OpenDirectory(dataDir);
    if (!dirHandle) {
        ::Error("analyze_pg_spectrum", "Cannot open directory: %s", dataDir);
        return;
    }

    std::vector<std::pair<Double_t, std::string>> files; // (depth, fullpath)
    const char* entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dirHandle)) != nullptr) {
        std::string name(entry);
        if (name.find("PG_Spectrum_VS_Angle_") == std::string::npos) continue;
        if (name.find(".root") == std::string::npos)                 continue;
        std::string full = std::string(dataDir) + "/" + name;
        Double_t depth = DepthFromFile(full);
        if (depth < 0) continue;
        files.emplace_back(depth, full);
    }
    gSystem->FreeDirectory(dirHandle);
    std::sort(files.begin(), files.end());

    if (files.empty()) {
        ::Error("analyze_pg_spectrum",
                "No PG_Spectrum_VS_Angle_*.root files found in %s", dataDir);
        return;
    }

    const Int_t nFiles = (Int_t)files.size();
    Printf("Found %d ROOT file(s):", nFiles);
    for (auto& kv : files)
        Printf("  depth = %.1f mm  ->  %s", kv.first, kv.second.c_str());

    // --------------------------------------------------------
    //  Open output ROOT file and create sub-directories
    // --------------------------------------------------------
    TFile* fOut = TFile::Open(outFile, "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        ::Error("analyze_pg_spectrum", "Cannot create output file: %s", outFile);
        return;
    }
    TDirectory* dGraphs       = fOut->mkdir("graphs");
    TDirectory* dBroadGraphs  = dGraphs->mkdir("broad_angles");
    TDirectory* d44Graphs     = dGraphs->mkdir("4p4MeV");
    TDirectory* d96Graphs     = dGraphs->mkdir("9p6MeV");
    TDirectory* dCanvases     = fOut->mkdir("canvases");

    // --------------------------------------------------------
    //  2. Define histograms to read
    // --------------------------------------------------------

    // --- broad angular bins ---
    const Int_t nBroad = 6;
    const char* broadNames[nBroad] = {
        "PG_spectra_0_to_30_deg",
        "PG_spectra_30_to_60_deg",
        "PG_spectra_60_to_90_deg",
        "PG_spectra_90_to_120_deg",
        "PG_spectra_120_to_150_deg",
        "PG_spectra_150_to_180_deg"
    };
    const char* broadLabels[nBroad] = {
        "0#circ#font[122]{-}30#circ",
        "30#circ#font[122]{-}60#circ",
        "60#circ#font[122]{-}90#circ",
        "90#circ#font[122]{-}120#circ",
        "120#circ#font[122]{-}150#circ",
        "150#circ#font[122]{-}180#circ"
    };
    static const Color_t broadColors[nBroad] = {
        kBlue+1, kRed+1, kGreen+2, kMagenta+1, kCyan+2, kOrange+7
    };

    // --- specific angles for the gamma-line histograms ---
    const Int_t nAngles = 9;
    const Int_t specificAngles[nAngles] = {30, 50, 60, 65, 90, 115, 120, 130, 150};

    // --------------------------------------------------------
    //  3. Fill TGraphs: integral vs depth
    // --------------------------------------------------------

    // broad bins: one TGraph per histogram name
    TGraph* gBroad[nBroad];
    for (Int_t i = 0; i < nBroad; ++i) gBroad[i] = new TGraph(nFiles);

    // specific gamma lines: one TGraph per (energy, angle)
    TGraph* g44[nAngles];   // 4.4 MeV
    TGraph* g96[nAngles];   // 9.6 MeV
    for (Int_t j = 0; j < nAngles; ++j) {
        g44[j] = new TGraph(nFiles);
        g96[j] = new TGraph(nFiles);
    }

    for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
        Double_t depth = files[iFile].first;
        TFile* f = TFile::Open(files[iFile].second.c_str(), "READ");
        if (!f || f->IsZombie()) {
            ::Warning("analyze_pg_spectrum", "Cannot open %s",
                      files[iFile].second.c_str());
            continue;
        }

        // broad angular bins
        for (Int_t i = 0; i < nBroad; ++i) {
            TH1D* h = (TH1D*)f->Get(broadNames[i]);
            Double_t integral = h ? h->Integral() : 0.;
            gBroad[i]->SetPoint(iFile, depth, integral);
        }

        // specific gamma-line histograms
        for (Int_t j = 0; j < nAngles; ++j) {
            char hname[128];

            std::snprintf(hname, sizeof(hname),
                          "4.400000_MeV_Gamma_%d_deg", specificAngles[j]);
            TH1D* h44 = (TH1D*)f->Get(hname);
            g44[j]->SetPoint(iFile, depth, h44 ? h44->Integral() : 0.);

            std::snprintf(hname, sizeof(hname),
                          "9.600000_MeV_Gamma_%d_deg", specificAngles[j]);
            TH1D* h96 = (TH1D*)f->Get(hname);
            g96[j]->SetPoint(iFile, depth, h96 ? h96->Integral() : 0.);
        }

        f->Close();
        delete f;
    }

    // --------------------------------------------------------
    //  Name TGraphs and save them to the output file
    // --------------------------------------------------------
    for (Int_t i = 0; i < nBroad; ++i) {
        gBroad[i]->SetName(broadNames[i]);
        gBroad[i]->SetTitle(Form("%s;Depth (mm);Intensity (counts / primary)",
                                 broadNames[i]));
        dBroadGraphs->cd();
        gBroad[i]->Write();
    }
    for (Int_t j = 0; j < nAngles; ++j) {
        g44[j]->SetName(Form("4p4MeV_Gamma_%ddeg", specificAngles[j]));
        g44[j]->SetTitle(Form("4.4 MeV gamma %d deg;Depth (mm);Intensity (counts / primary)",
                               specificAngles[j]));
        d44Graphs->cd();
        g44[j]->Write();

        g96[j]->SetName(Form("9p6MeV_Gamma_%ddeg", specificAngles[j]));
        g96[j]->SetTitle(Form("9.6 MeV gamma %d deg;Depth (mm);Intensity (counts / primary)",
                               specificAngles[j]));
        d96Graphs->cd();
        g96[j]->Write();
    }

    // --------------------------------------------------------
    //  Helper lambda: style one TGraph
    // --------------------------------------------------------
    auto StyleGraph = [](TGraph* g, Color_t col, Style_t marker, Int_t idx) {
        g->SetLineColor(col);
        g->SetMarkerColor(col);
        g->SetMarkerStyle(marker);
        g->SetMarkerSize(1.2);
        g->SetLineWidth(2);
        (void)idx;
    };

    // --------------------------------------------------------
    //  4. Plot: broad angular bins
    // --------------------------------------------------------
    {
        TCanvas* c = new TCanvas("cBroad", "Broad angles", 900, 600);
        c->SetLeftMargin(0.13);
        c->SetBottomMargin(0.13);

        // find axis range
        Double_t ymax = 0.;
        for (Int_t i = 0; i < nBroad; ++i)
            for (Int_t p = 0; p < nFiles; ++p)
                ymax = std::max(ymax, gBroad[i]->GetY()[p]);

        TGraph* frame = new TGraph(2);
        frame->SetPoint(0, files.front().first - 5, 0.);
        frame->SetPoint(1, files.back().first  + 5, ymax * 1.25);
        frame->SetMarkerStyle(1);
        frame->Draw("AP");
        frame->GetXaxis()->SetTitle("Depth (mm)");
        frame->GetYaxis()->SetTitle("Intensity (counts / primary)");
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetTitleOffset(1.3);

        TLegend* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
        leg->SetHeader("Detection angle", "C");
        leg->SetBorderSize(1);
        leg->SetTextSize(0.038);

        Style_t markers[nBroad] = {20, 21, 22, 23, 29, 34};
        for (Int_t i = 0; i < nBroad; ++i) {
            StyleGraph(gBroad[i], broadColors[i], markers[i], i);
            gBroad[i]->Draw("LP SAME");
            leg->AddEntry(gBroad[i], broadLabels[i], "lp");
        }
        leg->Draw();

        TLatex title;
        title.SetNDC();
        title.SetTextSize(0.045);
        title.SetTextAlign(22);
        title.DrawLatex(0.50, 0.96,
            "Prompt #gamma Intensity vs Depth (full spectrum 1.5#font[122]{-}12 MeV)");

        dCanvases->cd();
        c->Write("cBroadAngles");
        Printf("  Written canvas: canvases/cBroadAngles");
    }

    // --------------------------------------------------------
    //  Helper lambda: draw one "specific gamma-line" canvas
    // --------------------------------------------------------
    auto DrawGammaLine = [&](TGraph* grArr[], const char* energyStr,
                              const char* tag) {
        TCanvas* c = new TCanvas(Form("c_%s", tag),
                                 Form("%s MeV gamma", energyStr), 900, 600);
        c->SetLeftMargin(0.13);
        c->SetBottomMargin(0.13);

        Double_t ymax = 0.;
        for (Int_t j = 0; j < nAngles; ++j)
            for (Int_t p = 0; p < nFiles; ++p)
                ymax = std::max(ymax, grArr[j]->GetY()[p]);

        if (ymax == 0.) ymax = 1e-8;   // avoid empty frame

        TGraph* frame = new TGraph(2);
        frame->SetPoint(0, files.front().first - 5, 0.);
        frame->SetPoint(1, files.back().first  + 5, ymax * 1.30);
        frame->SetMarkerStyle(1);
        frame->Draw("AP");
        frame->GetXaxis()->SetTitle("Depth (mm)");
        frame->GetYaxis()->SetTitle("Intensity (counts / primary)");
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetTitleOffset(1.3);

        TLegend* leg = new TLegend(0.55, 0.50, 0.88, 0.88);
        leg->SetHeader("Detector angle", "C");
        leg->SetBorderSize(1);
        leg->SetTextSize(0.035);
        leg->SetNColumns(2);

        Style_t mkrs[nAngles] = {20, 21, 22, 23, 29, 34, 33, 47, 43};
        Bool_t hasAny = kFALSE;
        for (Int_t j = 0; j < nAngles; ++j) {
            Bool_t nonzero = kFALSE;
            for (Int_t p = 0; p < nFiles; ++p)
                if (grArr[j]->GetY()[p] > 0.) { nonzero = kTRUE; break; }
            if (!nonzero) continue;
            hasAny = kTRUE;
            StyleGraph(grArr[j], kPalette[j % kNColours], mkrs[j], j);
            grArr[j]->Draw("LP SAME");
            leg->AddEntry(grArr[j], Form("%d#circ", specificAngles[j]), "lp");
        }
        if (hasAny) leg->Draw();

        TLatex title;
        title.SetNDC();
        title.SetTextSize(0.045);
        title.SetTextAlign(22);
        title.DrawLatex(0.50, 0.96,
                        Form("Prompt #gamma Intensity vs Depth (%s MeV line)",
                             energyStr));

        dCanvases->cd();
        c->Write(Form("c_%s", tag));
        Printf("  Written canvas: canvases/c_%s", tag);
    };

    // --------------------------------------------------------
    //  5. Plot: 4.4 MeV and 9.6 MeV gamma lines
    // --------------------------------------------------------
    DrawGammaLine(g44, "4.4", "4p4MeV");
    DrawGammaLine(g96, "9.6", "9p6MeV");

    // --------------------------------------------------------
    //  6. Plot: spectrum overlay – one pad per angular bin
    // --------------------------------------------------------
    {
        const Int_t ncols = 3;
        const Int_t nrows = (nBroad + ncols - 1) / ncols;
        TCanvas* c = new TCanvas("cOverlay", "Spectra overlay",
                                 500 * ncols, 380 * nrows);
        c->Divide(ncols, nrows, 0.001, 0.001);

        // depth → colour
        Int_t dColours[20];
        Int_t palette[] = {kBlue+1, kCyan+1, kGreen+2, kOrange+7, kRed+1};
        for (Int_t k = 0; k < nFiles; ++k)
            dColours[k] = palette[k % 5];

        for (Int_t i = 0; i < nBroad; ++i) {
            TPad* pad = (TPad*)c->cd(i + 1);
            pad->SetLeftMargin(0.18);
            pad->SetBottomMargin(0.18);

            TH1D* hFrame = nullptr;
            Bool_t first = kTRUE;

            for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
                TFile* f = TFile::Open(files[iFile].second.c_str(), "READ");
                if (!f || f->IsZombie()) continue;
                TH1D* h = (TH1D*)f->Get(broadNames[i]);
                if (!h) { f->Close(); delete f; continue; }

                // clone so we can close the file
                TH1D* hc = (TH1D*)h->Clone(Form("h_%d_%d", i, iFile));
                hc->SetDirectory(nullptr);
                f->Close();
                delete f;

                hc->SetLineColor(dColours[iFile]);
                hc->SetLineWidth(2);
                hc->SetFillStyle(0);

                if (first) {
                    hc->Draw("HIST");
                    hc->GetXaxis()->SetTitle("Energy (MeV)");
                    hc->GetYaxis()->SetTitle("dN/dE (/ primary / MeV)");
                    hc->GetXaxis()->SetTitleSize(0.06);
                    hc->GetYaxis()->SetTitleSize(0.055);
                    hc->GetXaxis()->SetLabelSize(0.05);
                    hc->GetYaxis()->SetLabelSize(0.05);
                    hc->GetYaxis()->SetTitleOffset(1.5);
                    hFrame = hc;
                    first = kFALSE;
                } else {
                    // rescale frame if needed
                    if (hc->GetMaximum() > hFrame->GetMaximum())
                        hFrame->SetMaximum(hc->GetMaximum() * 1.2);
                    hc->Draw("HIST SAME");
                }
            }

            // pad title
            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(0.07);
            lab.SetTextAlign(22);
            lab.DrawLatex(0.55, 0.93, broadLabels[i]);
        }

        // hide unused pads
        for (Int_t i = nBroad + 1; i <= ncols * nrows; ++i)
            c->cd(i)->SetFillStyle(4000);

        // shared legend on first pad
        c->cd(1);
        TLegend* leg = new TLegend(0.22, 0.60, 0.70, 0.90);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.055);
        for (Int_t k = 0; k < nFiles; ++k) {
            TGraph* dummy = new TGraph(1);
            dummy->SetLineColor(dColours[k]);
            dummy->SetLineWidth(2);
            leg->AddEntry(dummy, Form("%.0f mm", files[k].first), "l");
        }
        leg->Draw();

        // super-title
        c->cd(0);
        TLatex suptitle;
        suptitle.SetNDC();
        suptitle.SetTextSize(0.025);
        suptitle.SetTextAlign(22);
        suptitle.DrawLatex(0.50, 0.995,
            "Prompt #gamma Spectra at Different Depths");

        dCanvases->cd();
        c->Write("cSpectraOverlay");
        Printf("  Written canvas: canvases/cSpectraOverlay");
    }

    // --------------------------------------------------------
    //  Close output file
    // --------------------------------------------------------
    fOut->Write("", TObject::kOverwrite);
    fOut->Close();
    Printf("\nResults written to: %s", outFile);
    Printf("  graphs/broad_angles/  – %d TGraphs (intensity vs depth)", nBroad);
    Printf("  graphs/4p4MeV/        – %d TGraphs", nAngles);
    Printf("  graphs/9p6MeV/        – %d TGraphs", nAngles);
    Printf("  canvases/             – 4 TCanvas objects");
    Printf("\nDone.");
}
