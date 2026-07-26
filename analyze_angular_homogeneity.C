////////////////////////////////////////////////////////////////////////////////
//   analyze_angular_homogeneity.C  for Hadron Theraphy (Bragg-Peak project)   //
//                                                                            //
//   Depth-resolved angular-homogeneity analysis of prompt-gamma (PG)         //
//   emission.  For every phantom depth z the macro reads the per-angle PG    //
//   energy spectra stored in  PG_Spectrum_VS_Angle_<depth>.root  and, for    //
//   each analysed gamma line, quantifies how uniform the PG yield is over    //
//   the polar emission angle theta by:                                       //
//                                                                            //
//     * a Poisson-weighted even Legendre expansion                           //
//           W(theta) = a0 [ 1 + a2 P2(cos theta) + a4 P4(cos theta) ]        //
//       -> anisotropy coefficients a2(z), a4(z) with fit uncertainties;      //
//     * two model-free uniformity indices                                    //
//           IU(z) = (Ymax - Ymin) / (Ymax + Ymin)   (integral uniformity)    //
//           CV(z) = sigma_Y / <Y>                    (coefficient of var.)    //
//       with finite-difference (linear) error propagation from the           //
//       Poisson yield uncertainties.                                         //
//                                                                            //
//   Every angle present in each input file is discovered automatically       //
//   (not just the 90 deg / 120 deg detectors used by the depth-profile       //
//   macros).  Results are written as TGraphErrors vs depth, JPGs, a ROOT     //
//   summary file and a CSV table.                                            //
//                                                                            //
//   Manuscript scope: the 4.44 MeV and 9.6 MeV lines are analysed; the       //
//   6.13 MeV line is intentionally excluded (same rationale as the           //
//   depth-profile benchmark).                                                //
//                                                                            //
//   Run:   root -l -b -q analyze_angular_homogeneity.C                       //
//                                                                            //
//              - 26. Jul. 2026.  Bragg-Peak / HadronTheraphy1                //
////////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
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

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <functional>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//                              CONFIGURATION                                  //
//   Everything a user may need to adjust lives in this block.                 //
////////////////////////////////////////////////////////////////////////////////
namespace AH {

// ---- input files ------------------------------------------------------------
// Files are auto-discovered as  <kInputDir>/<kFilePrefix><depth><kFileSuffix>
// and the <depth> part of the name is parsed into the numeric depth axis.
const TString kInputDir   = ".";
const TString kFilePrefix = "PG_Spectrum_VS_Angle_";
const TString kFileSuffix = ".root";

// ---- how the per-angle PG spectra are stored inside each file ---------------
// Primary layout: a single TH2 with the emission-angle bin on X and the PG
// energy [MeV] on Y (the natural generalisation of the simulation histogram
// "Angular_dist_of_Prompt_Gamma_EDep_on_SD").  A per-angle TH1 layout is used
// as an automatic fallback (see kTH1AngleToken).
const TString kAngleEnergyTH2 = "Angular_dist_of_Prompt_Gamma_EDep_on_SD";

// X-axis meaning of the primary TH2:
//   kXAxisInDegrees == false : X is the detector-division index (HT1 spherical
//                              detector).  theta = binCentre * kThetaSpanDeg / Nx
//   kXAxisInDegrees == true  : X-axis bin centres already are angles in degrees.
const bool   kXAxisInDegrees = false;
const double kThetaSpanDeg   = 180.0;   // polar coverage of the detector [deg]

// Fallback per-angle TH1 layout: any TH1 whose name contains this token is
// treated as one angle's spectrum, with the number just before the token read
// as the angle in degrees, e.g. "PG_Spectrum_90deg", "spectrum_120_deg".
const TString kTH1AngleToken = "deg";

// ---- prompt-gamma lines analysed --------------------------------------------
// name : file/label-safe tag        E : line energy [MeV]
// halfWin : half-width of the yield-integration window [MeV]
struct Line { const char* name; const char* label; double E; double halfWin; };
const std::vector<Line> kLines = {
    { "4p44MeV", "4.44 MeV", 4.44, 0.20 },
    { "9p6MeV",  "9.6 MeV",  9.60, 0.30 },
    // 6.13 MeV (16-O) intentionally excluded (see depth-profile benchmark).
};

// ---- misc -------------------------------------------------------------------
const TString kOutDir      = "AngularHomogeneity_out";
const TString kDepthAxis   = "Depth z  [file tag]";
const int     kMinAngles   = 3;   // minimum angles required for the a2/a4 fit

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
//   This is the only section tied to how the PG-vs-angle files are laid out.  //
//   It returns, for one analysed line, the angles [deg] and their yields with //
//   Poisson uncertainties.                                                    //
////////////////////////////////////////////////////////////////////////////////

struct AngleYields {
    std::vector<double> theta;   // emission angle [deg]
    std::vector<double> yield;   // integrated line yield
    std::vector<double> err;     // Poisson uncertainty on the yield
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

// Try the primary TH2(angle, energy) layout.
static bool ExtractFromTH2(TFile* f, const AH::Line& L, AngleYields& out)
{
    TH2* h2 = dynamic_cast<TH2*>(f->Get(AH::kAngleEnergyTH2));
    if (!h2) return false;

    const int nx = h2->GetNbinsX();
    const double scale = AH::kXAxisInDegrees ? 1.0 : (AH::kThetaSpanDeg / nx);

    for (int ix = 1; ix <= nx; ++ix) {
        TH1* proj = h2->ProjectionY(Form("_py_%s_%d", L.name, ix), ix, ix);
        double y = 0.0, e = 0.0;
        IntegrateLine(proj, L.E, L.halfWin, y, e);
        delete proj;
        if (y <= 0.0) continue;                       // no signal at this angle
        out.theta.push_back(h2->GetXaxis()->GetBinCenter(ix) * scale);
        out.yield.push_back(y);
        out.err.push_back(e > 0.0 ? e : std::sqrt(y));
    }
    return !out.theta.empty();
}

// Fallback: one TH1 spectrum per angle, angle encoded in the histogram name.
static bool ExtractFromTH1s(TFile* f, const AH::Line& L, AngleYields& out)
{
    TList* keys = f->GetListOfKeys();
    if (!keys) return false;

    std::map<double, std::pair<double, double>> byAngle;   // theta -> (yield, err)
    TIter next(keys);
    while (TKey* k = static_cast<TKey*>(next())) {
        const TString cname = k->GetClassName();
        if (!cname.BeginsWith("TH1")) continue;
        TString hname = k->GetName();
        const Ssiz_t tok = hname.Index(AH::kTH1AngleToken, 0, TString::kIgnoreCase);
        if (tok == kNPOS) continue;
        double theta = 0.0;
        if (!ParseLeadingNumber(hname, theta)) continue;   // needs an angle number
        TH1* spec = dynamic_cast<TH1*>(k->ReadObj());
        if (!spec) continue;
        double y = 0.0, e = 0.0;
        IntegrateLine(spec, L.E, L.halfWin, y, e);
        delete spec;
        if (y > 0.0) byAngle[theta] = { y, e > 0.0 ? e : std::sqrt(y) };
    }
    for (const auto& kv : byAngle) {
        out.theta.push_back(kv.first);
        out.yield.push_back(kv.second.first);
        out.err.push_back(kv.second.second);
    }
    return !out.theta.empty();
}

static bool ExtractAngleYields(const TString& path, const AH::Line& L, AngleYields& out)
{
    TFile f(path, "READ");
    if (f.IsZombie()) {
        std::cerr << "  [warn] cannot open " << path << std::endl;
        return false;
    }
    bool ok = ExtractFromTH2(&f, L, out);
    if (!ok) ok = ExtractFromTH1s(&f, L, out);
    f.Close();
    return ok;
}

////////////////////////////////////////////////////////////////////////////////
//                          ANGULAR-HOMOGENEITY METRICS                        //
////////////////////////////////////////////////////////////////////////////////

struct DepthResult {
    double z = 0.0;
    double a2 = 0.0, a2e = 0.0;
    double a4 = 0.0, a4e = 0.0;
    double iu = 0.0, iue = 0.0;
    double cv = 0.0, cve = 0.0;
    int    nAngles = 0;
    bool   fitOK = false;
};

static DepthResult ComputeMetrics(double z, const AngleYields& d)
{
    DepthResult r;
    r.z = z;
    r.nAngles = static_cast<int>(d.theta.size());

    // ---- even Legendre fit -> a2, a4 -----------------------------------------
    if (r.nAngles >= AH::kMinAngles) {
        double amin = *std::min_element(d.theta.begin(), d.theta.end());
        double amax = *std::max_element(d.theta.begin(), d.theta.end());
        TGraphErrors g(r.nAngles);
        for (int i = 0; i < r.nAngles; ++i) {
            g.SetPoint(i, d.theta[i], d.yield[i]);
            g.SetPointError(i, 0.0, d.err[i]);
        }
        double y0 = 0.0;
        for (double v : d.yield) y0 += v;
        y0 /= r.nAngles;

        TF1 fW("fW", LegendreW, amin, amax, 3);
        fW.SetParameters(y0, 0.0, 0.0);
        fW.SetParNames("a0", "a2", "a4");
        // Poisson weighting comes from the per-point 1/err^2 (default chi2 fit).
        TFitResultPtr res = g.Fit(&fW, "Q S N");
        if (res.Get() && res->IsValid()) {
            r.a2  = res->Parameter(1);
            r.a2e = res->ParError(1);
            r.a4  = res->Parameter(2);
            r.a4e = res->ParError(2);
            r.fitOK = true;
        }
    }

    // ---- model-free uniformity indices ---------------------------------------
    if (r.nAngles >= 2) {
        r.iu  = IntegralUniformity(d.yield);
        r.iue = PropagateError(IntegralUniformity, d.yield, d.err);
        r.cv  = CoeffOfVariation(d.yield);
        r.cve = PropagateError(CoeffOfVariation, d.yield, d.err);
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

////////////////////////////////////////////////////////////////////////////////
//                                 MAIN                                        //
////////////////////////////////////////////////////////////////////////////////

void analyze_angular_homogeneity()
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFontSize(0.045);

    gSystem->mkdir(AH::kOutDir, kTRUE);

    // ---- discover the PG-vs-angle files and their depths ---------------------
    std::vector<std::pair<double, TString>> files;   // (depth, path)
    {
        TSystemDirectory dir(AH::kInputDir, AH::kInputDir);
        TList* list = dir.GetListOfFiles();
        if (list) {
            TIter next(list);
            while (TSystemFile* sf = static_cast<TSystemFile*>(next())) {
                if (sf->IsDirectory()) continue;
                TString name = sf->GetName();
                if (!name.BeginsWith(AH::kFilePrefix) || !name.EndsWith(AH::kFileSuffix))
                    continue;
                TString tag = name;
                tag.Remove(0, AH::kFilePrefix.Length());
                tag.Remove(tag.Length() - AH::kFileSuffix.Length(), AH::kFileSuffix.Length());
                double z = 0.0;
                if (!ParseLeadingNumber(tag, z)) {
                    std::cerr << "  [warn] cannot parse depth from '" << name << "', skipped\n";
                    continue;
                }
                TString path = AH::kInputDir + "/" + name;
                files.emplace_back(z, path);
            }
        }
    }
    std::sort(files.begin(), files.end(),
              [](const std::pair<double, TString>& a, const std::pair<double, TString>& b) {
                  return a.first < b.first;
              });

    if (files.empty()) {
        std::cerr << "\nNo input files matched  " << AH::kInputDir << "/"
                  << AH::kFilePrefix << "*<depth>*" << AH::kFileSuffix << "\n"
                  << "Set AH::kInputDir / AH::kFilePrefix to point at your "
                  << "PG_Spectrum_VS_Angle_<depth>.root files.\n" << std::endl;
        return;
    }

    std::cout << "\nFound " << files.size() << " depth file(s):\n";
    for (const auto& fp : files)
        std::cout << "   z = " << fp.first << "   " << fp.second << "\n";

    // ---- per-line analysis over depth ----------------------------------------
    TFile fout(AH::kOutDir + "/angular_homogeneity.root", "RECREATE");
    std::ofstream csv((AH::kOutDir + "/angular_homogeneity.csv").Data());
    csv << "line,depth,n_angles,a2,a2_err,a4,a4_err,IU,IU_err,CV,CV_err\n";

    for (const auto& L : AH::kLines) {
        std::cout << "\n=== Line " << L.label << " ("
                  << L.E << " +/- " << L.halfWin << " MeV) ===\n";

        std::vector<DepthResult> R;
        for (const auto& fp : files) {
            AngleYields d;
            if (!ExtractAngleYields(fp.second, L, d)) {
                std::cerr << "  [warn] no angular yields for " << L.label
                          << " in " << fp.second << "\n";
                continue;
            }
            DepthResult r = ComputeMetrics(fp.first, d);
            R.push_back(r);
            std::cout << "  z=" << r.z << "  nAng=" << r.nAngles
                      << "  a2=" << r.a2 << "+/-" << r.a2e
                      << "  a4=" << r.a4 << "+/-" << r.a4e
                      << "  IU=" << r.iu << "+/-" << r.iue
                      << "  CV=" << r.cv << "+/-" << r.cve
                      << (r.fitOK ? "" : "  [fit skipped]") << "\n";

            csv << L.name << "," << r.z << "," << r.nAngles << ","
                << r.a2 << "," << r.a2e << "," << r.a4 << "," << r.a4e << ","
                << r.iu << "," << r.iue << "," << r.cv << "," << r.cve << "\n";
        }
        if (R.empty()) continue;

        // Build TGraphErrors (a2/a4 only from valid fits).
        TGraphErrors gA2, gA4, gIU, gCV;
        int nFit = 0, nAll = 0;
        for (const auto& r : R) {
            if (r.fitOK) {
                gA2.SetPoint(nFit, r.z, r.a2);  gA2.SetPointError(nFit, 0.0, r.a2e);
                gA4.SetPoint(nFit, r.z, r.a4);  gA4.SetPointError(nFit, 0.0, r.a4e);
                ++nFit;
            }
            gIU.SetPoint(nAll, r.z, r.iu);  gIU.SetPointError(nAll, 0.0, r.iue);
            gCV.SetPoint(nAll, r.z, r.cv);  gCV.SetPointError(nAll, 0.0, r.cve);
            ++nAll;
        }

        gA2.SetName(Form("a2_vs_depth_%s", L.name));
        gA4.SetName(Form("a4_vs_depth_%s", L.name));
        gIU.SetName(Form("IU_vs_depth_%s", L.name));
        gCV.SetName(Form("CV_vs_depth_%s", L.name));

        const TString tag = L.name;
        DrawAndSave(&gA2, Form("Legendre a_{2}(z)  -  %s", L.label), "a_{2}",
                    AH::kOutDir + "/a2_vs_depth_" + tag + ".jpg", kAzure + 2);
        DrawAndSave(&gA4, Form("Legendre a_{4}(z)  -  %s", L.label), "a_{4}",
                    AH::kOutDir + "/a4_vs_depth_" + tag + ".jpg", kViolet + 1);
        DrawAndSave(&gIU, Form("Integral uniformity IU(z)  -  %s", L.label),
                    "IU = (Y_{max}-Y_{min})/(Y_{max}+Y_{min})",
                    AH::kOutDir + "/IU_vs_depth_" + tag + ".jpg", kOrange + 7);
        DrawAndSave(&gCV, Form("Coefficient of variation CV(z)  -  %s", L.label),
                    "CV = #sigma_{Y} / #LTY#GT",
                    AH::kOutDir + "/CV_vs_depth_" + tag + ".jpg", kTeal + 2);

        fout.cd();
        gA2.Write(); gA4.Write(); gIU.Write(); gCV.Write();
    }

    csv.close();
    fout.Close();

    std::cout << "\nDone. Outputs written to '" << AH::kOutDir << "/':\n"
              << "   angular_homogeneity.root  (TGraphErrors)\n"
              << "   angular_homogeneity.csv   (table)\n"
              << "   a2/a4/IU/CV_vs_depth_<line>.jpg\n" << std::endl;
}

// Allow  root -l -b -q analyze_angular_homogeneity.C  to run main directly.
void analyze_angular_homogeneity_C() { analyze_angular_homogeneity(); }
