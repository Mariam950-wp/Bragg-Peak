////////////////////////////////////////////////////////////////////////////////
//   estimate_pg_growth_vs_peak.C  for Hadron Theraphy (Bragg-Peak project)    //
//                                                                            //
//   Angular-homogeneity (isotropy) analysis of prompt-gamma (PG) emission     //
//   in the INTENSITY-GROWTH region on the proximal (rising) side of the       //
//   distal fall-off, compared against the region right at the Bragg peak.     //
//   Companion of  estimate_pg_homogeneity.C  (which treats the Bragg-peak     //
//   region alone); it re-uses the very same method - solid-angle-normalised   //
//   differential yield  I(theta)=dN/dOmega  and even-order Legendre fit        //
//   W(theta)=A0[1+a2 P2(cos theta)+a4 P4(cos theta)]  (manuscript main14.tex,  //
//   Sec. "Angular homogeneity of the prompt-gamma emission").                 //
//                                                                            //
//   MOTIVATION                                                                //
//   ----------                                                                //
//   The manuscript notes that the yield of both lines "rises towards the end   //
//   of range and drops sharply at the distal fall-off".  estimate_pg_homo-    //
//   geneity.C characterises the emission AT the Bragg peak (the kNearPeak      //
//   depths closest to d_BP).  This macro does the same for the depth window   //
//   where the intensity STARTS TO GROW, on the near side of the peak, and     //
//   overlays the two regions so the change of the angular distribution        //
//   between "far from" and "close to" the Bragg peak can be read directly.    //
//                                                                            //
//   The growth window is line-dependent, measured as an interval of depth     //
//   BEFORE (upstream of) the Bragg peak:                                       //
//        9.6 MeV  line :  10 mm ... 5 mm  before the Bragg peak               //
//        4.44 MeV line :   4 mm ... 2 mm  before the Bragg peak               //
//   i.e. in the shifted coordinate  z - d_BP  the windows are [-10,-5] mm and //
//   [-4,-2] mm respectively.  "Near the Bragg peak" keeps the definition of   //
//   estimate_pg_homogeneity.C: the kNearPeak (=3) depths closest to d_BP.     //
//                                                                            //
//   INPUT  (identical layout to estimate_pg_homogeneity.C)                     //
//   -----                                                                      //
//   A "mother" directory is given as the single argument.  Each simulation    //
//   configuration (beam energy x physics list) lives in its own               //
//   sub-directory whose name encodes both, e.g.                               //
//        <mother>/prompt_gamma_spectra_130MeV_FTFP_BERT_HP/                    //
//                     PG_Spectrum_VS_Angle_105.root                           //
//                     PG_Spectrum_VS_Angle_107.root ...                       //
//   The <depth> token in each file name is the degrader thickness [mm].       //
//   Inside every file the polar emission is segmented into 30-deg rings        //
//        PG_spectra_0_to_30_deg ... PG_spectra_150_to_180_deg.                //
//                                                                            //
//   WHAT THIS MACRO PRODUCES                                                  //
//   ------------------------                                                  //
//   For every (beam energy x line) it draws ONE picture overlaying the        //
//   prompt-gamma intensity  dN/dOmega  versus the polar emission angle theta   //
//   for the SIX cases  {3 physics lists} x {near Bragg peak, growth region}.  //
//   Colour encodes the physics list; FILLED markers are the near-peak region  //
//   and OPEN markers the growth region.  The even-Legendre fit is NOT drawn    //
//   (to keep the crowded 6-curve canvas readable): instead the fitted         //
//   polynomial W(theta) and its a2,a4 are reported in the legend, each line    //
//   printed in the colour of its physics list.  Exactly four JPGs are written //
//        peak_vs_growth_130MeV_4p44MeV.jpg   peak_vs_growth_130MeV_9p6MeV.jpg //
//        peak_vs_growth_70MeV_4p44MeV.jpg    peak_vs_growth_70MeV_9p6MeV.jpg  //
//   plus a summary table peak_vs_growth_summary.csv (region, a2, a4, IU, CV,  //
//   chi2/ndf per physics list).  All are written into the mother directory.   //
//                                                                            //
//   Run:  root -l -b -q 'estimate_pg_growth_vs_peak.C("/path/to/mother")'     //
//         root -l -b -q  estimate_pg_growth_vs_peak.C     (uses ".")          //
//                                                                            //
//              - 02. Aug. 2026.  Bragg-Peak / HadronTheraphy1                 //
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
#include "TLegend.h"
#include "TLegendEntry.h"

#include <vector>
#include <string>
#include <map>
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
namespace GP {

// ---- input depth files ------------------------------------------------------
// Files are auto-discovered as  <cfgDir>/<kFilePrefix><depth><kFileSuffix>
// and the <depth> part of the name is parsed into the numeric depth axis.
const TString kFilePrefix = "PG_Spectrum_VS_Angle_";
const TString kFileSuffix = ".root";

// ---- configuration sub-directory names --------------------------------------
// Parsed as  prompt_gamma_spectra_<energy>MeV_<physicsList>  to recover the
// beam energy [MeV] and the physics-list label.
const char*  kConfigScanf = "prompt_gamma_spectra_%dMeV_%127s";

// ---- polar-ring spectra inside each file ------------------------------------
// Ring histograms are named  PG_spectra_<lo>_to_<hi>_deg  and are discovered by
// this scanf pattern; <lo>,<hi> are the ring edges in degrees.
const char*  kRingScanf = "PG_spectra_%d_to_%d_deg";

// ---- prompt-gamma lines analysed --------------------------------------------
// name : file/label-safe tag        E : line energy [MeV]
// halfWin : half-width of the yield-integration window [MeV]
// growInner/growOuter : the intensity-growth depth window, expressed as the
//   near/far edges measured in [mm] BEFORE (upstream of) the Bragg peak, so the
//   window in the shifted coordinate z-d_BP is  [-growOuter, -growInner].
struct Line {
    const char* name;
    const char* label;
    double E;
    double halfWin;
    double growInner;   // closest-to-peak edge of the growth window [mm before peak]
    double growOuter;   // furthest-from-peak edge of the growth window [mm before peak]
};
const std::vector<Line> kLines = {
    { "4p44MeV", "4.44 MeV", 4.44, 0.20, 2.0,  4.0  },   // growth window: 4..2 mm before peak
    { "9p6MeV",  "9.6 MeV",  9.60, 0.30, 5.0, 10.0  },   // growth window: 10..5 mm before peak
    // 6.13 MeV (16-O) intentionally excluded (see manuscript / benchmark).
};

// ---- Bragg-peak depth in PMMA per beam energy [mm] --------------------------
// From the manuscript: 130.87 MeV -> ~107 mm, 70.54 MeV -> ~33 mm.
static double BraggPeakDepth(int energyMeV)
{
    if (energyMeV >= 115 && energyMeV <= 145) return 107.0;   // ~130 MeV
    if (energyMeV >=  55 && energyMeV <=  85) return  33.0;   // ~70  MeV
    return -1.0;                                              // unknown
}

// ---- misc -------------------------------------------------------------------
const int kNearPeak  = 3;   // number of near-Bragg-peak depths ("close" region)
const int kMinAngles = 3;   // minimum rings required for the a2/a4 fit

} // namespace GP

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
// through independent Poisson uncertainties sig[i] on the yields y[i].
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
    err = (e > 0.0) ? e : (y > 0.0 ? std::sqrt(y) : 0.0);
}

// One discovered polar ring.
struct Ring { int lo; int hi; TString name; };

static bool ExtractRingYields(const TString& path, const GP::Line& L, AngleYields& out)
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
            if (std::sscanf(hname.Data(), GP::kRingScanf, &lo, &hi) == 2 && hi > lo)
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
    double a0 = 0.0, a0e = 0.0;
    double a2 = 0.0, a2e = 0.0;
    double a4 = 0.0, a4e = 0.0;
    double iu = 0.0, iue = 0.0;
    double cv = 0.0, cve = 0.0;
    double chi2ndf = 0.0;
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

    if (nFit >= GP::kMinAngles) {
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

// Sum the ring line-yields over a set of depth files and return the depth-
// averaged differential angular distribution I(theta) = dN/dOmega, so the
// homogeneity is characterised over the whole depth window at once.
static AngleYields AggregateDepths(const std::vector<std::pair<double, TString>>& files,
                                   const GP::Line& L)
{
    std::map<double, double> sumI, sumE2;
    std::map<double, int>    cnt;
    for (const auto& fp : files) {
        AngleYields d;
        if (!ExtractRingYields(fp.second, L, d)) continue;
        for (size_t i = 0; i < d.theta.size(); ++i) {
            sumI [d.theta[i]] += d.yield[i];
            sumE2[d.theta[i]] += d.err[i] * d.err[i];
            cnt  [d.theta[i]] += 1;
        }
    }
    AngleYields out;                       // std::map keeps the rings sorted in theta
    for (const auto& kv : sumI) {
        const double th = kv.first;
        const int    n  = cnt[th];
        if (n <= 0) continue;
        out.theta.push_back(th);
        out.yield.push_back(kv.second / n);              // depth-averaged dN/dOmega
        out.err.push_back(std::sqrt(sumE2[th]) / n);
    }
    return out;
}

////////////////////////////////////////////////////////////////////////////////
//                     DIRECTORY / DEPTH-FILE DISCOVERY                         //
////////////////////////////////////////////////////////////////////////////////

// Discover the PG_Spectrum_VS_Angle_<depth>.root files in one directory,
// sorted by depth.
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
        if (!name.BeginsWith(GP::kFilePrefix) || !name.EndsWith(GP::kFileSuffix)) continue;
        TString tag = name;
        tag.Remove(0, GP::kFilePrefix.Length());
        tag.Remove(tag.Length() - GP::kFileSuffix.Length(), GP::kFileSuffix.Length());
        double z = 0.0;
        if (!ParseLeadingNumber(tag, z)) continue;
        files.emplace_back(z, dirPath + "/" + name);
    }
    std::sort(files.begin(), files.end(),
              [](const std::pair<double, TString>& a, const std::pair<double, TString>& b) {
                  return a.first < b.first;
              });
    return files;
}

// Keep the kNearPeak depths closest to dRef ("near Bragg peak"), sorted by depth.
static std::vector<std::pair<double, TString>>
SelectNearPeak(std::vector<std::pair<double, TString>> files, double dRef, int nWant)
{
    std::sort(files.begin(), files.end(),
              [dRef](const std::pair<double, TString>& a, const std::pair<double, TString>& b) {
                  return std::fabs(a.first - dRef) < std::fabs(b.first - dRef);
              });
    if (static_cast<int>(files.size()) > nWant) files.resize(nWant);
    std::sort(files.begin(), files.end(),
              [](const std::pair<double, TString>& a, const std::pair<double, TString>& b) {
                  return a.first < b.first;
              });
    return files;
}

// Keep the depth files whose depth lies in the intensity-growth window on the
// proximal side of the peak:  dRef-outerMM <= z <= dRef-innerMM  (i.e. between
// innerMM and outerMM before the Bragg peak).  Returned sorted by depth.
static std::vector<std::pair<double, TString>>
SelectGrowthWindow(const std::vector<std::pair<double, TString>>& files,
                   double dRef, double innerMM, double outerMM)
{
    const double zLo = dRef - outerMM;   // furthest-upstream edge
    const double zHi = dRef - innerMM;   // closest-to-peak edge
    const double eps = 1e-6;             // tolerance so window edges are inclusive
    std::vector<std::pair<double, TString>> sel;
    for (const auto& f : files)
        if (f.first >= zLo - eps && f.first <= zHi + eps)
            sel.push_back(f);
    std::sort(sel.begin(), sel.end(),
              [](const std::pair<double, TString>& a, const std::pair<double, TString>& b) {
                  return a.first < b.first;
              });
    return sel;
}

////////////////////////////////////////////////////////////////////////////////
//                                 PLOTTING                                    //
////////////////////////////////////////////////////////////////////////////////

// One drawn case = one physics list in one depth region.
struct Series {
    TString     phys;      // physics-list label
    TString     region;    // "near peak" | "growth"
    AngleYields d;         // averaged differential angular distribution I(theta)
    DepthResult r;         // Legendre A0/a2/a4 (+ IU, CV, chi2) of that distribution
    int         nDepths;   // number of depth files averaged
    int         color;     // colour (== physics list)
    int         marker;    // marker style (filled == near peak, open == growth)
};

// Overlay the six cases {3 physics lists} x {near peak, growth} of one
// (beam energy, line): prompt-gamma intensity dN/dOmega versus the polar
// emission angle.  Colour encodes the physics list, FILLED markers the near-
// peak region and OPEN markers the growth region.  The even-Legendre fit is
// deliberately NOT drawn to keep the six-case canvas readable; instead its
// polynomial and (a2,a4) are reported in the legend, coloured per physics list.
static void DrawPeakVsGrowth(int energy, const GP::Line& L,
                             const std::vector<Series>& series, const TString& outfile)
{
    if (series.empty()) return;

    TCanvas c(Form("c_%dMeV_%s", energy, L.name), "peak_vs_growth", 950, 700);
    c.SetGrid();
    c.SetLeftMargin(0.15);
    c.SetBottomMargin(0.12);

    double ymax = 0.0;
    for (const Series& s : series)
        for (size_t i = 0; i < s.d.theta.size(); ++i)
            ymax = std::max(ymax, s.d.yield[i] + s.d.err[i]);
    if (ymax <= 0.0) ymax = 1.0;

    // Legend: header states the Legendre polynomial form; one line per case,
    // each printed in its physics-list colour, giving the fitted a2, a4.
    TLegend leg(0.14, 0.60, 0.93, 0.90);
    leg.SetBorderSize(1);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.0235);
    leg.SetHeader(Form("%d MeV,  %s        "
                       "W(#theta)=A_{0}[1+a_{2}P_{2}(cos#theta)+a_{4}P_{4}(cos#theta)]",
                       energy, L.label));

    for (size_t is = 0; is < series.size(); ++is) {
        const Series& s = series[is];
        const int n = static_cast<int>(s.d.theta.size());
        if (n <= 0) continue;

        TGraphErrors* g = new TGraphErrors(n);
        for (int i = 0; i < n; ++i) {
            g->SetPoint(i, s.d.theta[i], s.d.yield[i]);
            g->SetPointError(i, 0.0, s.d.err[i]);
        }
        g->SetMarkerStyle(s.marker);
        g->SetMarkerSize(1.5);
        g->SetMarkerColor(s.color);
        g->SetLineColor(s.color);

        if (is == 0) {
            g->SetTitle(Form("Prompt-#gamma angular distribution: near Bragg peak vs "
                             "intensity-growth region  -  %d MeV,  %s", energy, L.label));
            g->Draw("AP");
            g->GetXaxis()->SetTitle("Emission angle  #theta  [deg]");
            g->GetYaxis()->SetTitle("dN/d#Omega   [sr^{-1}]");
            g->GetYaxis()->SetTitleOffset(1.7);
            g->GetXaxis()->SetLimits(0.0, 180.0);
            g->GetYaxis()->SetRangeUser(0.0, ymax * 1.9);
        } else {
            g->Draw("P SAME");
        }

        // NOTE: the even-Legendre fit is intentionally NOT drawn here.
        // Only its coefficients are reported, in the legend, below.
        TString disp = s.phys;
        disp.ReplaceAll("_", " ");          // underscores are TLatex subscripts
        TString txt;
        if (s.r.fitOK)
            txt = Form("%-16s %-9s:  a_{2}=%+.3f#pm%.3f,  a_{4}=%+.3f#pm%.3f",
                       disp.Data(), s.region.Data(), s.r.a2, s.r.a2e, s.r.a4, s.r.a4e);
        else
            txt = Form("%-16s %-9s:  (fit unavailable)", disp.Data(), s.region.Data());

        TLegendEntry* e = leg.AddEntry(g, txt.Data(), "p");
        if (e) e->SetTextColor(s.color);    // "corresponding colour for fit"
    }

    leg.Draw();
    gPad->RedrawAxis();
    c.SaveAs(outfile);
}

////////////////////////////////////////////////////////////////////////////////
//                                 MAIN                                        //
////////////////////////////////////////////////////////////////////////////////

void estimate_pg_growth_vs_peak(const char* motherDir = ".")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFontSize(0.04);

    const TString mother = motherDir;

    // ---- group configuration sub-directories by beam energy ------------------
    //   energy [MeV]  ->  list of (physicsList, dirPath)
    std::map<int, std::vector<std::pair<TString, TString>>> byEnergy;
    {
        TSystemDirectory dir(mother, mother);
        if (TList* list = dir.GetListOfFiles()) {
            TIter next(list);
            while (TSystemFile* sf = static_cast<TSystemFile*>(next())) {
                if (!sf->IsDirectory()) continue;
                TString name = sf->GetName();
                if (name == "." || name == "..") continue;
                int energy = 0;
                char phys[128] = { 0 };
                if (std::sscanf(name.Data(), GP::kConfigScanf, &energy, phys) != 2) continue;
                const TString sub = mother + "/" + name;
                if (ListDepthFiles(sub).empty()) continue;
                byEnergy[energy].emplace_back(TString(phys), sub);
            }
        }
    }

    if (byEnergy.empty()) {
        std::cerr << "\nNo 'prompt_gamma_spectra_<energy>MeV_<physics>' sub-directories with "
                  << GP::kFilePrefix << "*" << GP::kFileSuffix << " files found under '"
                  << mother << "'.\n"
                  << "Usage:  root -l -b -q 'estimate_pg_growth_vs_peak.C(\"/path/to/mother\")'\n"
                  << std::endl;
        return;
    }

    std::ofstream csv((mother + "/peak_vs_growth_summary.csv").Data());
    csv << "energy_MeV,line,region,physics_list,n_depths,"
        << "a2,a2_err,a4,a4_err,IU,IU_err,CV,CV_err,chi2ndf\n";

    // Colour == physics list; filled marker == near peak, open marker == growth.
    const int physColor[]      = { kAzure + 2, kRed + 1, kGreen + 2, kViolet + 1, kOrange + 7 };
    const int markerNear[]     = { 20, 21, 22, 33, 29 };   // filled circle/square/triangle/...
    const int markerGrowth[]   = { 24, 25, 26, 27, 30 };   // matching OPEN markers
    const int nCol             = 5;

    int nJpg = 0;
    for (auto& en : byEnergy) {
        const int energy = en.first;
        std::vector<std::pair<TString, TString>>& models = en.second;
        std::sort(models.begin(), models.end(),
                  [](const std::pair<TString, TString>& a, const std::pair<TString, TString>& b) {
                      return a.first < b.first;
                  });

        // Bragg-peak depth for this energy (fallback: median of the first scan).
        double dRef = GP::BraggPeakDepth(energy);
        if (dRef <= 0.0 && !models.empty()) {
            std::vector<std::pair<double, TString>> f0 = ListDepthFiles(models[0].second);
            if (!f0.empty()) dRef = f0[f0.size() / 2].first;
            std::cerr << "  [warn] unknown Bragg-peak depth for " << energy
                      << " MeV; using z_ref = " << dRef << " mm\n";
        }

        std::cout << "\n=== " << energy << " MeV   (Bragg peak z = " << dRef << " mm,  "
                  << models.size() << " physics lists) ===\n";

        for (const auto& L : GP::kLines) {
            // Build, for every physics list, the near-peak case then the growth
            // case, pushed in that order so the legend groups the two regions of
            // each physics list (each colour) together.
            std::vector<Series> series;

            for (size_t im = 0; im < models.size(); ++im) {
                const TString& phys = models[im].first;
                const TString& dirp = models[im].second;
                std::vector<std::pair<double, TString>> all = ListDepthFiles(dirp);

                // ---- near-peak ("close") region --------------------------------
                std::vector<std::pair<double, TString>> near =
                    SelectNearPeak(all, dRef, GP::kNearPeak);
                AngleYields aggN = AggregateDepths(near, L);
                if (!aggN.theta.empty()) {
                    Series s;
                    s.phys = phys;  s.region = "near peak";
                    s.d = aggN;     s.r = ComputeMetrics(dRef, aggN);
                    s.nDepths = static_cast<int>(near.size());
                    s.color = physColor[im % nCol];
                    s.marker = markerNear[im % nCol];
                    series.push_back(s);

                    csv << energy << "," << L.name << ",near_peak," << phys << ","
                        << near.size() << ","
                        << s.r.a2 << "," << s.r.a2e << "," << s.r.a4 << "," << s.r.a4e << ","
                        << s.r.iu << "," << s.r.iue << "," << s.r.cv << "," << s.r.cve << ","
                        << s.r.chi2ndf << "\n";
                } else {
                    std::cerr << "   [warn] no near-peak ring yields for " << phys
                              << " / " << L.label << "\n";
                }

                // ---- intensity-growth ("far") region ---------------------------
                std::vector<std::pair<double, TString>> grow =
                    SelectGrowthWindow(all, dRef, L.growInner, L.growOuter);
                AngleYields aggG = AggregateDepths(grow, L);
                if (!aggG.theta.empty()) {
                    Series s;
                    s.phys = phys;  s.region = "growth";
                    s.d = aggG;     s.r = ComputeMetrics(dRef - 0.5 * (L.growInner + L.growOuter), aggG);
                    s.nDepths = static_cast<int>(grow.size());
                    s.color = physColor[im % nCol];
                    s.marker = markerGrowth[im % nCol];
                    series.push_back(s);

                    csv << energy << "," << L.name << ",growth," << phys << ","
                        << grow.size() << ","
                        << s.r.a2 << "," << s.r.a2e << "," << s.r.a4 << "," << s.r.a4e << ","
                        << s.r.iu << "," << s.r.iue << "," << s.r.cv << "," << s.r.cve << ","
                        << s.r.chi2ndf << "\n";
                } else {
                    std::cerr << "   [warn] no growth-region ring yields for " << phys
                              << " / " << L.label << " (window "
                              << L.growOuter << "-" << L.growInner << " mm before peak, z in ["
                              << dRef - L.growOuter << ", " << dRef - L.growInner << "] mm)\n";
                }

                std::cout << "   " << phys << " / " << L.label << ": near="
                          << near.size() << " depth(s), growth=" << grow.size() << " depth(s)\n";
            }

            if (series.empty()) continue;

            const TString out =
                mother + "/" + Form("peak_vs_growth_%dMeV_%s.jpg", energy, L.name);
            DrawPeakVsGrowth(energy, L, series, out);
            std::cout << "   -> " << out << "   (" << series.size() << " cases)\n";
            ++nJpg;
        }
    }

    csv.close();

    std::cout << "\nDone. Wrote " << nJpg << " JPG(s) and peak_vs_growth_summary.csv into '"
              << mother << "/'.\n"
              << "Each JPG overlays dN/dOmega vs emission angle for the three physics lists in "
              << "the near-Bragg-peak region (filled markers) and the intensity-growth region "
              << "(open markers); the Legendre a2,a4 are given in the legend, no fit curve drawn.\n"
              << std::endl;
}
