////////////////////////////////////////////////////////////////////////////////
//   analyze_angular_homogeneity.C  for Hadron Theraphy (Bragg-Peak project)   //
//                                                                            //
//   Depth-resolved angular-homogeneity (isotropy) analysis of prompt-gamma   //
//   (PG) emission, following the method of the manuscript                    //
//   (main14.tex, Sec. "Angular homogeneity of the prompt-gamma emission").   //
//                                                                            //
//   INPUT                                                                     //
//   -----                                                                     //
//   A "mother" directory is given as the single argument.  Each simulation   //
//   configuration (beam energy x physics list) lives in its own              //
//   sub-directory, e.g.                                                       //
//        <mother>/prompt_gamma_spectra_130MeV_FTFP_BERT_HP/                   //
//                     PG_Spectrum_VS_Angle_67.root                            //
//                     PG_Spectrum_VS_Angle_69.root ...                        //
//   The <depth> token in each file name is the degrader thickness [mm].       //
//   Every sub-directory that contains such files is analysed independently;   //
//   if the mother directory itself holds the files it is treated as one       //
//   configuration (so a single folder of example files also works).          //
//                                                                            //
//   Inside every PG_Spectrum_VS_Angle_<depth>.root the polar emission is      //
//   segmented into 30-deg rings stored as full PG energy spectra             //
//        PG_spectra_0_to_30_deg, PG_spectra_30_to_60_deg, ... 150_to_180.     //
//   The rings are discovered automatically from their names.                  //
//                                                                            //
//   METHOD (per depth, per line)                                             //
//   ---------------------------                                              //
//     1. Differential angular yield.  The line intensity N_gamma is obtained //
//        by integrating each ring spectrum over the line energy window and    //
//        divided by the TRUE ring solid angle                                 //
//              dOmega = 2*pi*(cos theta_lo - cos theta_hi),                   //
//        giving  I(theta) = dN/dOmega , flat in theta for an isotropic        //
//        source.  Poisson errors sqrt(N) are propagated to each ring.         //
//     2. Even-order Legendre fit                                              //
//              W(theta) = A0 [ 1 + a2 P2(cos theta) + a4 P4(cos theta) ]      //
//        -> anisotropy coefficients a2(z), a4(z) with fit uncertainties.      //
//     3. Model-free uniformity indices                                        //
//              IU = (I_max - I_min)/(I_max + I_min)                           //
//              CV = sigma_I / <I>                                             //
//        plus the reduced chi2/ndf of the flat (isotropic) hypothesis.        //
//        IU, CV and chi2/ndf remain defined in the low-yield distal region    //
//        where the Legendre fit becomes unstable.                            //
//                                                                            //
//   For isotropic emission a2 = a4 = 0, IU = CV = 0 and chi2/ndf = 1.         //
//                                                                            //
//   OUTPUT (per configuration, under AngularHomogeneity_out/<config>/)        //
//   -----                                                                     //
//     * a2/a4/IU/CV/chi2ndf_vs_depth_<line>.jpg   (trends vs depth)           //
//     * angular_dist_<line>_z<depth>.jpg          (I(theta) + Legendre fit)   //
//     * angular_homogeneity.root                  (TGraphErrors)              //
//     * angular_homogeneity.csv                   (table)                     //
//                                                                            //
//   Manuscript scope: the 4.44 MeV and 9.6 MeV lines are analysed; the        //
//   6.13 MeV line is intentionally excluded (Geant4 cross-section issue).     //
//                                                                            //
//   Run:  root -l -b -q 'analyze_angular_homogeneity.C("/path/to/mother")'    //
//         root -l -b -q  analyze_angular_homogeneity.C     (uses ".")         //
//                                                                            //
//              - 29. Jul. 2026.  Bragg-Peak / HadronTheraphy1                 //
////////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TList.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TMath.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TLegend.h"

#include <vector>
#include <string>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <functional>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//                              CONFIGURATION                                  //
//   Everything a user may need to adjust lives in this block.                 //
////////////////////////////////////////////////////////////////////////////////
namespace AH {

// ---- input depth files ------------------------------------------------------
// Files are auto-discovered as  <cfgDir>/<kFilePrefix><depth><kFileSuffix>
// and the <depth> part of the name is parsed into the numeric depth axis.
const TString kFilePrefix = "PG_Spectrum_VS_Angle_";
const TString kFileSuffix = ".root";

// ---- polar-ring spectra inside each file ------------------------------------
// Ring histograms are named  PG_spectra_<lo>_to_<hi>_deg  and are discovered by
// this scanf pattern; <lo>,<hi> are the ring edges in degrees.
const char*  kRingScanf = "PG_spectra_%d_to_%d_deg";

// ---- prompt-gamma lines analysed --------------------------------------------
// name : file/label-safe tag        E : line energy [MeV]
// halfWin : half-width of the yield-integration window [MeV]
struct Line { const char* name; const char* label; double E; double halfWin; };
const std::vector<Line> kLines = {
    { "4p44MeV", "4.44 MeV", 4.44, 0.20 },
    { "9p6MeV",  "9.6 MeV",  9.60, 0.30 },
    // 6.13 MeV (16-O) intentionally excluded (see manuscript / benchmark).
};

// ---- misc -------------------------------------------------------------------
const TString kOutDir    = "AngularHomogeneity_out";
const TString kDepthAxis = "Degrader thickness z  [mm]";
const int     kMinAngles = 3;   // minimum rings required for the a2/a4 fit

} // namespace AH

////////////////////////////////////////////////////////////////////////////////
//                          SMALL MATH HELPERS                                 //
////////////////////////////////////////////////////////////////////////////////

// Legendre polynomials in x = cos(theta).
static inline double P2(double c) { return 0.5 * (3.0 * c * c - 1.0); }
static inline double P4(double c) { return 0.125 * (35.0 * c * c * c * c - 30.0 * c * c + 3.0); }

// Even Legendre expansion used for the fit:  p0 (1 + p1 P2 + p2 P4).
static double LegendreW(double* x, double* p)
{
    const double c = TMath::Cos(x[0] * TMath::DegToRad());
    return p[0] * (1.0 + p[1] * P2(c) + p[2] * P4(c));
}

// Integral (max-min) uniformity.
static double IntegralUniformity(const std::vector<double>& y)
{
    if (y.size() < 2) return 0.0;
    const double ymax = *std::max_element(y.begin(), y.end());
    const double ymin = *std::min_element(y.begin(), y.end());
    const double den  = ymax + ymin;
    return (den > 0.0) ? (ymax - ymin) / den : 0.0;
}

// Coefficient of variation (sample standard deviation / mean).
static double CoeffOfVariation(const std::vector<double>& y)
{
    const int n = static_cast<int>(y.size());
    if (n < 2) return 0.0;
    double mean = 0.0;
    for (double v : y) mean += v;
    mean /= n;
    if (mean == 0.0) return 0.0;
    double s2 = 0.0;
    for (double v : y) s2 += (v - mean) * (v - mean);
    s2 /= (n - 1);
    return std::sqrt(s2) / mean;
}

// Linear (finite-difference) error propagation of a scalar functional f(y)
// through independent Poisson uncertainties sig[i] on the yields y[i]:
//   var(f) = sum_i ( df/dy_i )^2 sig_i^2 ,  df/dy_i via central differences.
static double PropagateError(const std::function<double(const std::vector<double>&)>& f,
                             const std::vector<double>& y,
                             const std::vector<double>& sig)
{
    double var = 0.0;
    for (size_t i = 0; i < y.size(); ++i) {
        const double h = sig[i];
        if (h <= 0.0) continue;
        std::vector<double> yp = y, ym = y;
        yp[i] += h;
        ym[i] -= h;
        const double deriv = (f(yp) - f(ym)) / (2.0 * h);
        var += deriv * deriv * h * h;
    }
    return std::sqrt(var);
}

// Parse the first signed decimal number found in a string; returns false if none.
static bool ParseLeadingNumber(const TString& s, double& out)
{
    const char* p = s.Data();
    while (*p && !(std::isdigit((unsigned char)*p) ||
                   ((*p == '+' || *p == '-' || *p == '.') && std::isdigit((unsigned char)*(p + 1)))))
        ++p;
    if (!*p) return false;
    char* end = nullptr;
    out = std::strtod(p, &end);
    return end != p;
}

////////////////////////////////////////////////////////////////////////////////
//                    PER-DEPTH DATA EXTRACTION (I/O ADAPTER)                   //
//   Reads the polar-ring spectra of one depth file and returns, for one       //
//   analysed line, the ring angles [deg] and the solid-angle-normalised       //
//   differential yield  I(theta) = dN/dOmega  with Poisson uncertainties.     //
////////////////////////////////////////////////////////////////////////////////

struct AngleYields {
    std::vector<double> theta;   // ring-centre emission angle [deg]
    std::vector<double> yield;   // differential line yield  dN/dOmega [1/sr]
    std::vector<double> err;     // Poisson uncertainty on the differential yield
};

// Integrate a 1-D energy spectrum in [E-halfWin, E+halfWin] with its error.
static void IntegrateLine(TH1* spec, double E, double halfWin, double& yield, double& err)
{
    const int b1 = spec->GetXaxis()->FindBin(E - halfWin);
    const int b2 = spec->GetXaxis()->FindBin(E + halfWin);
    double e = 0.0;
    const double y = spec->IntegralAndError(b1, b2, e);
    yield = y;
    // Guard against a null error on a non-empty window (unweighted fills).
    err = (e > 0.0) ? e : (y > 0.0 ? std::sqrt(y) : 0.0);
}

// One discovered polar ring.
struct Ring { int lo; int hi; TString name; };

static bool ExtractRingYields(const TString& path, const AH::Line& L, AngleYields& out)
{
    TFile f(path, "READ");
    if (f.IsZombie()) {
        std::cerr << "  [warn] cannot open " << path << std::endl;
        return false;
    }

    // ---- discover the polar-ring spectra by name -----------------------------
    std::vector<Ring> rings;
    if (TList* keys = f.GetListOfKeys()) {
        TIter next(keys);
        while (TKey* k = static_cast<TKey*>(next())) {
            if (!TString(k->GetClassName()).BeginsWith("TH1")) continue;
            const TString hname = k->GetName();
            int lo = 0, hi = 0;
            if (std::sscanf(hname.Data(), AH::kRingScanf, &lo, &hi) == 2 && hi > lo)
                rings.push_back({ lo, hi, hname });
        }
    }
    if (rings.empty()) {
        std::cerr << "  [warn] no 'PG_spectra_<lo>_to_<hi>_deg' rings in " << path << "\n";
        f.Close();
        return false;
    }
    std::sort(rings.begin(), rings.end(),
              [](const Ring& a, const Ring& b) { return a.lo < b.lo; });

    // ---- differential yield per ring -----------------------------------------
    const double twoPi = TMath::TwoPi();
    for (const Ring& r : rings) {
        TH1* h = dynamic_cast<TH1*>(f.Get(r.name));
        if (!h) continue;
        double N = 0.0, eN = 0.0;
        IntegrateLine(h, L.E, L.halfWin, N, eN);

        const double dOmega = twoPi * (TMath::Cos(r.lo * TMath::DegToRad()) -
                                       TMath::Cos(r.hi * TMath::DegToRad()));
        if (dOmega <= 0.0) continue;

        out.theta.push_back(0.5 * (r.lo + r.hi));   // ring centre [deg]
        out.yield.push_back(N  / dOmega);           // dN/dOmega [1/sr]
        out.err.push_back(eN / dOmega);
    }
    f.Close();
    return !out.theta.empty();
}

////////////////////////////////////////////////////////////////////////////////
//                          ANGULAR-HOMOGENEITY METRICS                        //
////////////////////////////////////////////////////////////////////////////////

struct DepthResult {
    double z = 0.0;
    double a0 = 0.0, a0e = 0.0;   // Legendre normalisation (needed to draw W(theta))
    double a2 = 0.0, a2e = 0.0;
    double a4 = 0.0, a4e = 0.0;
    double iu = 0.0, iue = 0.0;
    double cv = 0.0, cve = 0.0;
    double chi2ndf = 0.0;         // reduced chi2 of the flat (isotropic) hypothesis
    int    nAngles = 0;
    bool   fitOK  = false;
    bool   chiOK  = false;
};

static DepthResult ComputeMetrics(double z, const AngleYields& d)
{
    DepthResult r;
    r.z = z;
    r.nAngles = static_cast<int>(d.theta.size());

    // ---- even Legendre fit -> a2, a4 (rings with a valid error only) ---------
    std::vector<double> ft, fy, fe;
    for (int i = 0; i < r.nAngles; ++i)
        if (d.yield[i] > 0.0 && d.err[i] > 0.0) {
            ft.push_back(d.theta[i]);
            fy.push_back(d.yield[i]);
            fe.push_back(d.err[i]);
        }
    const int nFit = static_cast<int>(ft.size());

    if (nFit >= AH::kMinAngles) {
        double amin = *std::min_element(ft.begin(), ft.end());
        double amax = *std::max_element(ft.begin(), ft.end());
        TGraphErrors g(nFit);
        double y0 = 0.0;
        for (int i = 0; i < nFit; ++i) {
            g.SetPoint(i, ft[i], fy[i]);
            g.SetPointError(i, 0.0, fe[i]);
            y0 += fy[i];
        }
        y0 /= nFit;

        TF1 fW("fW", LegendreW, amin, amax, 3);
        fW.SetParameters(y0, 0.0, 0.0);
        fW.SetParNames("A0", "a2", "a4");
        // Poisson weighting comes from the per-point 1/err^2 (default chi2 fit).
        TFitResultPtr res = g.Fit(&fW, "Q S N");
        if (res.Get() && res->IsValid()) {
            r.a0  = res->Parameter(0);  r.a0e = res->ParError(0);
            r.a2  = res->Parameter(1);  r.a2e = res->ParError(1);
            r.a4  = res->Parameter(2);  r.a4e = res->ParError(2);
            r.fitOK = true;
        }
    }

    // ---- model-free uniformity indices (all rings) ---------------------------
    if (r.nAngles >= 2) {
        r.iu  = IntegralUniformity(d.yield);
        r.iue = PropagateError(IntegralUniformity, d.yield, d.err);
        r.cv  = CoeffOfVariation(d.yield);
        r.cve = PropagateError(CoeffOfVariation, d.yield, d.err);
    }

    // ---- chi2/ndf of the flat (isotropic) hypothesis -------------------------
    // Weighted mean over rings with a valid error, then reduced chi2 (ndf=n-1).
    if (nFit >= 2) {
        double sw = 0.0, swy = 0.0;
        for (int i = 0; i < nFit; ++i) {
            const double w = 1.0 / (fe[i] * fe[i]);
            sw += w;  swy += w * fy[i];
        }
        const double mean = swy / sw;
        double chi2 = 0.0;
        for (int i = 0; i < nFit; ++i) {
            const double dlt = (fy[i] - mean) / fe[i];
            chi2 += dlt * dlt;
        }
        r.chi2ndf = chi2 / (nFit - 1);
        r.chiOK   = true;
    }
    return r;
}

////////////////////////////////////////////////////////////////////////////////
//                              PLOTTING                                       //
////////////////////////////////////////////////////////////////////////////////

static void DrawAndSave(TGraphErrors* g, const TString& title, const TString& ytitle,
                        const TString& outfile, int color)
{
    TCanvas c(Form("c_%s", outfile.Data()), title, 900, 650);
    c.SetGrid();
    c.SetLeftMargin(0.14);
    c.SetBottomMargin(0.12);

    g->SetTitle(title);
    g->GetXaxis()->SetTitle(AH::kDepthAxis);
    g->GetYaxis()->SetTitle(ytitle);
    g->GetYaxis()->SetTitleOffset(1.5);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.1);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
    g->Draw("AP");

    c.SaveAs(outfile);
}

// Per-depth angular distribution:  differential yield I(theta)=dN/dOmega vs the
// polar emission angle, with the even-Legendre fit W(theta)=A0[1+a2 P2+a4 P4]
// overlaid (drawn only when the fit succeeded).  One JPG per line and depth.
static void DrawAngularDistribution(const AH::Line& L, const DepthResult& r,
                                    const AngleYields& d, const TString& outfile)
{
    const int n = static_cast<int>(d.theta.size());
    if (n < 1) return;

    TCanvas c(Form("cad_%s", outfile.Data()), "AngularDistribution", 900, 650);
    c.SetGrid();
    c.SetLeftMargin(0.14);
    c.SetBottomMargin(0.12);

    TGraphErrors g(n);
    double ymax = 0.0, amin = d.theta[0], amax = d.theta[0];
    for (int i = 0; i < n; ++i) {
        g.SetPoint(i, d.theta[i], d.yield[i]);
        g.SetPointError(i, 0.0, d.err[i]);
        ymax = std::max(ymax, d.yield[i] + d.err[i]);
        amin = std::min(amin, d.theta[i]);
        amax = std::max(amax, d.theta[i]);
    }
    g.SetTitle(Form("Differential PG yield vs emission angle  -  %s   (z = %g mm)", L.label, r.z));
    g.GetXaxis()->SetTitle("Emission angle  #theta  [deg]");
    g.GetYaxis()->SetTitle("dN/d#Omega  [sr^{-1}]");
    g.GetYaxis()->SetTitleOffset(1.5);
    if (ymax > 0.0) g.GetYaxis()->SetRangeUser(0.0, ymax * 1.25);
    g.SetMarkerStyle(20);
    g.SetMarkerSize(1.1);
    g.SetMarkerColor(kAzure + 2);
    g.SetLineColor(kAzure + 2);
    g.SetLineWidth(2);
    g.Draw("AP");

    TLegend leg(0.58, 0.74, 0.88, 0.88);
    leg.SetBorderSize(1);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.030);
    leg.AddEntry(&g, "dN/d#Omega #pm stat.", "lp");

    // Even-Legendre fit overlay (same model used to extract a2/a4).
    TF1* fit = nullptr;
    if (r.fitOK && amax > amin) {
        fit = new TF1(Form("fWdraw_%s", outfile.Data()), LegendreW, amin, amax, 3);
        fit->SetParameters(r.a0, r.a2, r.a4);
        fit->SetLineColor(kRed + 1);
        fit->SetLineWidth(2);
        fit->SetNpx(400);
        fit->Draw("L SAME");
        leg.AddEntry(fit, "A_{0}[1+a_{2}P_{2}+a_{4}P_{4}]", "l");
    }
    leg.Draw();

    // Numeric summary of the homogeneity metrics for this depth.
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.033);
    double yline = 0.86;
    if (r.fitOK) {
        tx.DrawLatex(0.17, yline, Form("a_{2} = %.3f #pm %.3f", r.a2, r.a2e)); yline -= 0.05;
        tx.DrawLatex(0.17, yline, Form("a_{4} = %.3f #pm %.3f", r.a4, r.a4e)); yline -= 0.05;
    }
    tx.DrawLatex(0.17, yline, Form("IU = %.3f   CV = %.3f", r.iu, r.cv));        yline -= 0.05;
    if (r.chiOK)
        tx.DrawLatex(0.17, yline, Form("#chi^{2}/ndf = %.2f", r.chi2ndf));

    gPad->RedrawAxis();
    c.SaveAs(outfile);
    delete fit;
}

////////////////////////////////////////////////////////////////////////////////
//                        PER-CONFIGURATION PROCESSING                         //
////////////////////////////////////////////////////////////////////////////////

// Discover the PG_Spectrum_VS_Angle_<depth>.root files in one directory.
static std::vector<std::pair<double, TString>> ListDepthFiles(const TString& dirPath)
{
    std::vector<std::pair<double, TString>> files;
    TSystemDirectory dir(dirPath, dirPath);
    TList* list = dir.GetListOfFiles();
    if (!list) return files;
    TIter next(list);
    while (TSystemFile* sf = static_cast<TSystemFile*>(next())) {
        if (sf->IsDirectory()) continue;
        TString name = sf->GetName();
        if (!name.BeginsWith(AH::kFilePrefix) || !name.EndsWith(AH::kFileSuffix)) continue;
        TString tag = name;
        tag.Remove(0, AH::kFilePrefix.Length());
        tag.Remove(tag.Length() - AH::kFileSuffix.Length(), AH::kFileSuffix.Length());
        double z = 0.0;
        if (!ParseLeadingNumber(tag, z)) {
            std::cerr << "  [warn] cannot parse depth from '" << name << "', skipped\n";
            continue;
        }
        files.emplace_back(z, dirPath + "/" + name);
    }
    std::sort(files.begin(), files.end(),
              [](const std::pair<double, TString>& a, const std::pair<double, TString>& b) {
                  return a.first < b.first;
              });
    return files;
}

// Analyse one configuration (one directory holding a depth scan).
static void ProcessConfig(const TString& dirPath, const TString& label)
{
    std::vector<std::pair<double, TString>> files = ListDepthFiles(dirPath);
    if (files.empty()) return;

    const TString outDir = AH::kOutDir + "/" + label;
    gSystem->mkdir(outDir, kTRUE);

    std::cout << "\n############################################################\n"
              << "# Configuration: " << label << "   (" << files.size() << " depths)\n"
              << "############################################################\n";

    TFile fout(outDir + "/angular_homogeneity.root", "RECREATE");
    std::ofstream csv((outDir + "/angular_homogeneity.csv").Data());
    csv << "line,depth,n_rings,a2,a2_err,a4,a4_err,IU,IU_err,CV,CV_err,chi2ndf\n";

    for (const auto& L : AH::kLines) {
        std::cout << "\n=== Line " << L.label << " ("
                  << L.E << " +/- " << L.halfWin << " MeV) ===\n";

        std::vector<DepthResult> R;
        for (const auto& fp : files) {
            AngleYields d;
            if (!ExtractRingYields(fp.second, L, d)) {
                std::cerr << "  [warn] no ring yields for " << L.label
                          << " in " << fp.second << "\n";
                continue;
            }
            DepthResult r = ComputeMetrics(fp.first, d);
            R.push_back(r);
            std::cout << "  z=" << r.z << "  nRings=" << r.nAngles
                      << "  a2=" << r.a2 << "+/-" << r.a2e
                      << "  a4=" << r.a4 << "+/-" << r.a4e
                      << "  IU=" << r.iu
                      << "  CV=" << r.cv
                      << "  chi2/ndf=" << r.chi2ndf
                      << (r.fitOK ? "" : "  [fit skipped]") << "\n";

            csv << L.name << "," << r.z << "," << r.nAngles << ","
                << r.a2 << "," << r.a2e << "," << r.a4 << "," << r.a4e << ","
                << r.iu << "," << r.iue << "," << r.cv << "," << r.cve << ","
                << r.chi2ndf << "\n";

            // JPG of the differential angular distribution + Legendre fit.
            DrawAngularDistribution(L, r, d,
                outDir + "/angular_dist_" + L.name + Form("_z%g", r.z) + ".jpg");
        }
        if (R.empty()) continue;

        // Build the depth-trend graphs.
        TGraphErrors gA2, gA4, gIU, gCV, gChi;
        int nFit = 0, nAll = 0, nChi = 0;
        for (const auto& r : R) {
            if (r.fitOK) {
                gA2.SetPoint(nFit, r.z, r.a2);  gA2.SetPointError(nFit, 0.0, r.a2e);
                gA4.SetPoint(nFit, r.z, r.a4);  gA4.SetPointError(nFit, 0.0, r.a4e);
                ++nFit;
            }
            gIU.SetPoint(nAll, r.z, r.iu);  gIU.SetPointError(nAll, 0.0, r.iue);
            gCV.SetPoint(nAll, r.z, r.cv);  gCV.SetPointError(nAll, 0.0, r.cve);
            ++nAll;
            if (r.chiOK) { gChi.SetPoint(nChi, r.z, r.chi2ndf); gChi.SetPointError(nChi, 0.0, 0.0); ++nChi; }
        }

        const TString tag = L.name;
        gA2.SetName(Form("a2_vs_depth_%s", tag.Data()));
        gA4.SetName(Form("a4_vs_depth_%s", tag.Data()));
        gIU.SetName(Form("IU_vs_depth_%s", tag.Data()));
        gCV.SetName(Form("CV_vs_depth_%s", tag.Data()));
        gChi.SetName(Form("chi2ndf_vs_depth_%s", tag.Data()));

        DrawAndSave(&gA2, Form("Legendre a_{2}(z)  -  %s  [%s]", L.label, label.Data()), "a_{2}",
                    outDir + "/a2_vs_depth_" + tag + ".jpg", kAzure + 2);
        DrawAndSave(&gA4, Form("Legendre a_{4}(z)  -  %s  [%s]", L.label, label.Data()), "a_{4}",
                    outDir + "/a4_vs_depth_" + tag + ".jpg", kViolet + 1);
        DrawAndSave(&gIU, Form("Integral uniformity IU(z)  -  %s  [%s]", L.label, label.Data()),
                    "IU = (I_{max}-I_{min})/(I_{max}+I_{min})",
                    outDir + "/IU_vs_depth_" + tag + ".jpg", kOrange + 7);
        DrawAndSave(&gCV, Form("Coefficient of variation CV(z)  -  %s  [%s]", L.label, label.Data()),
                    "CV = #sigma_{I} / #LTI#GT",
                    outDir + "/CV_vs_depth_" + tag + ".jpg", kTeal + 2);
        DrawAndSave(&gChi, Form("Flat-hypothesis #chi^{2}/ndf(z)  -  %s  [%s]", L.label, label.Data()),
                    "#chi^{2} / ndf",
                    outDir + "/chi2ndf_vs_depth_" + tag + ".jpg", kGray + 2);

        fout.cd();
        gA2.Write(); gA4.Write(); gIU.Write(); gCV.Write(); gChi.Write();
    }

    csv.close();
    fout.Close();

    std::cout << "\nWrote outputs to '" << outDir << "/'\n";
}

////////////////////////////////////////////////////////////////////////////////
//                                 MAIN                                        //
////////////////////////////////////////////////////////////////////////////////

void analyze_angular_homogeneity(const char* motherDir = ".")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFontSize(0.045);

    gSystem->mkdir(AH::kOutDir, kTRUE);

    const TString mother = motherDir;

    // ---- collect configurations ----------------------------------------------
    // A "configuration" is any directory that directly contains depth files:
    // each sub-directory of <mother>, and <mother> itself as a fallback.
    std::vector<std::pair<TString, TString>> configs;   // (label, dirPath)

    if (!ListDepthFiles(mother).empty()) {
        TString base = gSystem->BaseName(mother);
        if (base == "" || base == ".") base = "root";
        configs.emplace_back(base, mother);
    }

    {
        TSystemDirectory dir(mother, mother);
        if (TList* list = dir.GetListOfFiles()) {
            TIter next(list);
            while (TSystemFile* sf = static_cast<TSystemFile*>(next())) {
                if (!sf->IsDirectory()) continue;
                TString name = sf->GetName();
                if (name == "." || name == "..") continue;
                TString sub = mother + "/" + name;
                if (!ListDepthFiles(sub).empty())
                    configs.emplace_back(name, sub);
            }
        }
    }

    if (configs.empty()) {
        std::cerr << "\nNo '" << AH::kFilePrefix << "*<depth>" << AH::kFileSuffix
                  << "' files found in '" << mother << "' or its sub-directories.\n"
                  << "Usage:  root -l -b -q 'analyze_angular_homogeneity.C(\"/path/to/mother\")'\n"
                  << "where <mother> holds the  prompt_gamma_spectra_*  sub-directories.\n"
                  << std::endl;
        return;
    }

    std::sort(configs.begin(), configs.end(),
              [](const std::pair<TString, TString>& a, const std::pair<TString, TString>& b) {
                  return a.first < b.first;
              });

    std::cout << "\nFound " << configs.size() << " configuration(s) under '" << mother << "':\n";
    for (const auto& c : configs) std::cout << "   " << c.first << "   (" << c.second << ")\n";

    for (const auto& c : configs) ProcessConfig(c.second, c.first);

    std::cout << "\nDone. All outputs under '" << AH::kOutDir << "/<config>/':\n"
              << "   a2/a4/IU/CV/chi2ndf_vs_depth_<line>.jpg  (trends vs depth)\n"
              << "   angular_dist_<line>_z<depth>.jpg         (dN/dOmega + Legendre fit)\n"
              << "   angular_homogeneity.root / .csv\n" << std::endl;
}
