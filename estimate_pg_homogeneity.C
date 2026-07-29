////////////////////////////////////////////////////////////////////////////////
//   estimate_pg_homogeneity.C  for Hadron Theraphy (Bragg-Peak project)      //
//                                                                            //
//   Angular-homogeneity (isotropy) analysis of prompt-gamma (PG) emission     //
//   near the Bragg peak, following the method of the manuscript               //
//   (main14.tex, Sec. "Angular homogeneity of the prompt-gamma emission").   //
//                                                                            //
//   INPUT                                                                     //
//   -----                                                                     //
//   A "mother" directory is given as the single argument.  Each simulation    //
//   configuration (beam energy x physics list) lives in its own               //
//   sub-directory whose name encodes both, e.g.                               //
//        <mother>/prompt_gamma_spectra_130MeV_FTFP_BERT_HP/                    //
//                     PG_Spectrum_VS_Angle_105.root                           //
//                     PG_Spectrum_VS_Angle_107.root ...                       //
//   The <depth> token in each file name is the degrader thickness [mm].       //
//   Inside every file the polar emission is segmented into 30-deg rings        //
//   stored as full PG energy spectra                                          //
//        PG_spectra_0_to_30_deg ... PG_spectra_150_to_180_deg,                //
//   discovered automatically from their names.                                //
//                                                                            //
//   METHOD (per depth, per line) - manuscript technique                       //
//   --------------------------------------------------                        //
//     * differential angular yield  I(theta) = dN/dOmega , with the TRUE      //
//       ring solid angle  dOmega = 2*pi*(cos theta_lo - cos theta_hi)  and    //
//       propagated Poisson errors (flat in theta for an isotropic source);    //
//     * even-order Legendre fit                                               //
//           W(theta) = A0 [ 1 + a2 P2(cos theta) + a4 P4(cos theta) ]         //
//       -> anisotropy coefficients a2, a4;                                     //
//     * model-free indices IU, CV and the reduced chi2/ndf of the flat        //
//       (isotropic) hypothesis.  For isotropy a2=a4=0, IU=CV=0, chi2/ndf=1.   //
//                                                                            //
//   WHAT THIS MACRO PRODUCES                                                  //
//   ------------------------                                                  //
//   For every (beam energy x line) it draws ONE picture: the prompt-gamma     //
//   emission intensity  dN/dOmega  versus the polar emission angle theta,     //
//   averaged over the kNearPeak (=3) depths closest to the Bragg peak, the    //
//   three physics lists overlaid, each distribution fitted with the even      //
//   Legendre expansion W(theta).  Exactly four JPGs are written:              //
//        inhomogeneity_130MeV_4p44MeV.jpg   inhomogeneity_130MeV_9p6MeV.jpg   //
//        inhomogeneity_70MeV_4p44MeV.jpg    inhomogeneity_70MeV_9p6MeV.jpg    //
//   plus a summary table angular_homogeneity_summary.csv (a2, a4, IU, CV,     //
//   chi2/ndf per physics list).  All are written into the mother (input)      //
//   directory.                                                                //
//                                                                            //
//   Run:  root -l -b -q 'estimate_pg_homogeneity.C("/path/to/mother")'       //
//         root -l -b -q  estimate_pg_homogeneity.C        (uses ".")         //
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
#include "TLegend.h"

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
namespace AH {

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
struct Line { const char* name; const char* label; double E; double halfWin; };
const std::vector<Line> kLines = {
    { "4p44MeV", "4.44 MeV", 4.44, 0.20 },
    { "9p6MeV",  "9.6 MeV",  9.60, 0.30 },
    // 6.13 MeV (16-O) intentionally excluded (see manuscript / benchmark).
};

// ---- Bragg-peak depth in PMMA per beam energy [mm] --------------------------
// From the manuscript: 130.87 MeV -> ~107 mm, 70.54 MeV -> ~33 mm.  Only the
// 3 depths closest to this value are analysed.
static double BraggPeakDepth(int energyMeV)
{
    if (energyMeV >= 115 && energyMeV <= 145) return 107.0;   // ~130 MeV
    if (energyMeV >=  55 && energyMeV <=  85) return  33.0;   // ~70  MeV
    return -1.0;                                              // unknown
}

// ---- misc -------------------------------------------------------------------
const int kNearPeak = 3;   // number of near-Bragg-peak depths analysed
const int kMinAngles = 3;  // minimum rings required for the a2/a4 fit

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

// Sum the ring line-yields over the near-peak depth files and return the
// depth-averaged differential angular distribution I(theta) = dN/dOmega, so the
// homogeneity is characterised over the whole near-Bragg-peak region at once.
static AngleYields AggregateNearPeak(const std::vector<std::pair<double, TString>>& files,
                                     const AH::Line& L)
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
        if (!name.BeginsWith(AH::kFilePrefix) || !name.EndsWith(AH::kFileSuffix)) continue;
        TString tag = name;
        tag.Remove(0, AH::kFilePrefix.Length());
        tag.Remove(tag.Length() - AH::kFileSuffix.Length(), AH::kFileSuffix.Length());
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

// Keep the kNearPeak depths closest to dRef, returned sorted by depth.
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

////////////////////////////////////////////////////////////////////////////////
//                                 PLOTTING                                    //
////////////////////////////////////////////////////////////////////////////////

// One physics list: its near-peak angular distribution and the fitted metrics.
struct ModelPlot {
    TString     phys;   // physics-list label
    AngleYields d;      // averaged differential angular distribution I(theta)
    DepthResult r;      // Legendre A0/a2/a4 (+ IU, CV, chi2) of that distribution
};

// Overlay the physics lists of one (beam energy, line) case: the prompt-gamma
// intensity dN/dOmega versus the polar emission angle, each angular distribution
// fitted with the even Legendre expansion W(theta)=A0[1+a2 P2+a4 P4].  One JPG.
static void DrawAngularComparison(int energy, const AH::Line& L,
                                  const std::vector<ModelPlot>& models, const TString& outfile)
{
    if (models.empty()) return;

    TCanvas c(Form("c_%dMeV_%s", energy, L.name), "angular", 900, 650);
    c.SetGrid();
    c.SetLeftMargin(0.15);
    c.SetBottomMargin(0.12);

    double ymax = 0.0;
    for (const ModelPlot& m : models)
        for (size_t i = 0; i < m.d.theta.size(); ++i)
            ymax = std::max(ymax, m.d.yield[i] + m.d.err[i]);
    if (ymax <= 0.0) ymax = 1.0;

    const int cols[] = { kAzure + 2, kRed + 1, kTeal + 2, kViolet + 1, kOrange + 7 };
    const int mks[]  = { 20, 21, 22, 33, 29 };
    const int nCol   = 5;

    TLegend leg(0.15, 0.72, 0.72, 0.88);
    leg.SetBorderSize(1);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.028);
    leg.SetHeader(Form("%d MeV,  %s   (near Bragg peak)", energy, L.label));

    for (size_t im = 0; im < models.size(); ++im) {
        const ModelPlot& m = models[im];
        const int n   = static_cast<int>(m.d.theta.size());
        const int col = cols[im % nCol];

        TGraphErrors* g = new TGraphErrors(n);
        double tmin = 180.0, tmax = 0.0;
        for (int i = 0; i < n; ++i) {
            g->SetPoint(i, m.d.theta[i], m.d.yield[i]);
            g->SetPointError(i, 0.0, m.d.err[i]);
            tmin = std::min(tmin, m.d.theta[i]);
            tmax = std::max(tmax, m.d.theta[i]);
        }
        g->SetMarkerStyle(mks[im % nCol]);
        g->SetMarkerSize(1.4);
        g->SetMarkerColor(col);
        g->SetLineColor(col);

        if (im == 0) {
            g->SetTitle(Form("Prompt-#gamma angular distribution near Bragg peak  -  %d MeV,  %s",
                             energy, L.label));
            g->Draw("AP");
            g->GetXaxis()->SetTitle("Emission angle  #theta  [deg]");
            g->GetYaxis()->SetTitle("dN/d#Omega   [sr^{-1}]");
            g->GetYaxis()->SetTitleOffset(1.6);
            g->GetXaxis()->SetLimits(0.0, 180.0);
            g->GetYaxis()->SetRangeUser(0.0, ymax * 1.25);
        } else {
            g->Draw("P SAME");
        }

        // Even-Legendre fit overlay in the same colour.
        if (m.r.fitOK && tmax > tmin) {
            TF1* fit = new TF1(Form("fit_%dMeV_%s_%d", energy, L.name, static_cast<int>(im)),
                               LegendreW, tmin, tmax, 3);
            fit->SetParameters(m.r.a0, m.r.a2, m.r.a4);
            fit->SetLineColor(col);
            fit->SetLineWidth(2);
            fit->SetNpx(300);
            fit->Draw("L SAME");
        }

        TString disp = m.phys;
        disp.ReplaceAll("_", " ");          // underscores are TLatex subscripts
        leg.AddEntry(g, Form("%s   (a_{2}=%.2f, a_{4}=%.2f)", disp.Data(), m.r.a2, m.r.a4), "lp");
    }

    leg.Draw();
    gPad->RedrawAxis();
    c.SaveAs(outfile);
}

////////////////////////////////////////////////////////////////////////////////
//                                 MAIN                                        //
////////////////////////////////////////////////////////////////////////////////

void estimate_pg_homogeneity(const char* motherDir = ".")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFontSize(0.045);

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
                if (std::sscanf(name.Data(), AH::kConfigScanf, &energy, phys) != 2) continue;
                const TString sub = mother + "/" + name;
                if (ListDepthFiles(sub).empty()) continue;
                byEnergy[energy].emplace_back(TString(phys), sub);
            }
        }
    }

    if (byEnergy.empty()) {
        std::cerr << "\nNo 'prompt_gamma_spectra_<energy>MeV_<physics>' sub-directories with "
                  << AH::kFilePrefix << "*" << AH::kFileSuffix << " files found under '"
                  << mother << "'.\n"
                  << "Usage:  root -l -b -q 'estimate_pg_homogeneity.C(\"/path/to/mother\")'\n"
                  << std::endl;
        return;
    }

    std::ofstream csv((mother + "/angular_homogeneity_summary.csv").Data());
    csv << "energy_MeV,line,physics_list,n_depths,"
        << "a2,a2_err,a4,a4_err,IU,IU_err,CV,CV_err,chi2ndf\n";

    int nJpg = 0;
    for (auto& en : byEnergy) {
        const int energy = en.first;
        std::vector<std::pair<TString, TString>>& models = en.second;
        std::sort(models.begin(), models.end(),
                  [](const std::pair<TString, TString>& a, const std::pair<TString, TString>& b) {
                      return a.first < b.first;
                  });

        // Bragg-peak depth for this energy (fallback: median of the first scan).
        double dRef = AH::BraggPeakDepth(energy);
        if (dRef <= 0.0 && !models.empty()) {
            std::vector<std::pair<double, TString>> f0 = ListDepthFiles(models[0].second);
            if (!f0.empty()) dRef = f0[f0.size() / 2].first;
            std::cerr << "  [warn] unknown Bragg-peak depth for " << energy
                      << " MeV; using z_ref = " << dRef << " mm\n";
        }

        std::cout << "\n=== " << energy << " MeV   (Bragg peak z = " << dRef << " mm,  "
                  << models.size() << " physics lists) ===\n";

        for (const auto& L : AH::kLines) {
            std::vector<ModelPlot> plots;
            for (const auto& mdl : models) {
                const TString& phys = mdl.first;
                std::vector<std::pair<double, TString>> near =
                    SelectNearPeak(ListDepthFiles(mdl.second), dRef, AH::kNearPeak);

                AngleYields agg = AggregateNearPeak(near, L);
                if (agg.theta.empty()) {
                    std::cerr << "   [warn] no ring yields for " << phys
                              << " / " << L.label << "\n";
                    continue;
                }
                DepthResult r = ComputeMetrics(dRef, agg);

                std::cout << "   " << phys << ": " << near.size() << " depth(s), "
                          << agg.theta.size() << " rings,  a2=" << r.a2
                          << "  a4=" << r.a4 << (r.fitOK ? "" : "  [fit skipped]") << "\n";

                csv << energy << "," << L.name << "," << phys << "," << near.size() << ","
                    << r.a2 << "," << r.a2e << "," << r.a4 << "," << r.a4e << ","
                    << r.iu << "," << r.iue << "," << r.cv << "," << r.cve << ","
                    << r.chi2ndf << "\n";

                ModelPlot mp; mp.phys = phys; mp.d = agg; mp.r = r;
                plots.push_back(mp);
            }
            if (plots.empty()) continue;

            const TString out = mother + "/" + Form("inhomogeneity_%dMeV_%s.jpg", energy, L.name);
            DrawAngularComparison(energy, L, plots, out);
            std::cout << "   -> " << out << "\n";
            ++nJpg;
        }
    }

    csv.close();

    std::cout << "\nDone. Wrote " << nJpg << " JPG(s) and angular_homogeneity_summary.csv into '"
              << mother << "/'.\n"
              << "Each JPG: dN/dOmega vs emission angle near the Bragg peak, three "
              << "physics lists, with even-Legendre fits.\n" << std::endl;
}
