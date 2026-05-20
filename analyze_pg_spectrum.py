"""
Analyze prompt gamma spectra from Geant4 simulations.

Reads PG_Spectrum_VS_Angle_<depth>.root files, extracts the integrated
intensity of each angular/direction histogram, and plots intensity vs depth.

Usage:
    python3 analyze_pg_spectrum.py [--data-dir DIR] [--output-dir DIR]

By default, looks for ROOT files in the current directory.
"""

import argparse
import os
import re
import glob

import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Broad angular-bin histograms present in each file
BROAD_ANGLE_HISTS = [
    "PG_spectra_0_to_30_deg",
    "PG_spectra_30_to_60_deg",
    "PG_spectra_60_to_90_deg",
    "PG_spectra_90_to_120_deg",
    "PG_spectra_120_to_150_deg",
    "PG_spectra_150_to_180_deg",
]

# Pretty labels for the broad bins
BROAD_ANGLE_LABELS = {
    "PG_spectra_0_to_30_deg":    "0°–30°",
    "PG_spectra_30_to_60_deg":   "30°–60°",
    "PG_spectra_60_to_90_deg":   "60°–90°",
    "PG_spectra_90_to_120_deg":  "90°–120°",
    "PG_spectra_120_to_150_deg": "120°–150°",
    "PG_spectra_150_to_180_deg": "150°–180°",
}

# Specific gamma-line energies
GAMMA_ENERGIES = ["4.400000", "9.600000"]

# Specific angles present for the gamma-line histograms
SPECIFIC_ANGLES = [30, 50, 60, 65, 90, 115, 120, 130, 150]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def extract_depth(filename):
    """Return the depth value (float) encoded in the ROOT filename."""
    match = re.search(r"PG_Spectrum_VS_Angle_([\d.]+)\.root", os.path.basename(filename))
    if match:
        return float(match.group(1))
    raise ValueError(f"Cannot extract depth from filename: {filename}")


def integral(hist):
    """Return the sum of bin contents (total counts) of a TH1 histogram."""
    values, _ = hist.to_numpy()
    return float(np.sum(values))


def load_data(root_files):
    """
    Load all ROOT files and return a dict with structure:
        data[hist_name][depth] = intensity
    """
    data = {}

    for path in root_files:
        depth = extract_depth(path)
        f = uproot.open(path)

        # --- broad angular bins ---
        for hname in BROAD_ANGLE_HISTS:
            if hname not in data:
                data[hname] = {}
            if hname in f:
                data[hname][depth] = integral(f[hname])
            else:
                data[hname][depth] = 0.0

        # --- specific gamma lines ---
        for energy in GAMMA_ENERGIES:
            for angle in SPECIFIC_ANGLES:
                hname = f"{energy}_MeV_Gamma_{angle}_deg"
                if hname not in data:
                    data[hname] = {}
                if hname in f:
                    data[hname][depth] = integral(f[hname])
                else:
                    data[hname][depth] = 0.0

    return data


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_broad_angles(data, depths, output_dir):
    """One figure: intensity vs depth for each broad angular bin."""
    fig, ax = plt.subplots(figsize=(9, 6))

    colors = cm.tab10(np.linspace(0, 0.9, len(BROAD_ANGLE_HISTS)))

    for hname, color in zip(BROAD_ANGLE_HISTS, colors):
        intensities = [data[hname].get(d, 0.0) for d in depths]
        ax.plot(depths, intensities, marker="o", label=BROAD_ANGLE_LABELS[hname],
                color=color, linewidth=1.8, markersize=7)

    ax.set_xlabel("Depth (mm)", fontsize=13)
    ax.set_ylabel("Intensity (counts / primary)", fontsize=13)
    ax.set_title("Prompt Gamma Intensity vs Depth\n(broad angular bins, full spectrum 1.5–12 MeV)", fontsize=13)
    ax.legend(title="Detection angle", fontsize=10, title_fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.tight_layout()

    out = os.path.join(output_dir, "PG_intensity_vs_depth_broad_angles.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


def plot_specific_gamma(data, depths, energy_label, output_dir):
    """One figure per gamma energy: intensity vs depth at each specific angle."""
    hnames = [f"{energy_label}_MeV_Gamma_{a}_deg" for a in SPECIFIC_ANGLES]
    colors = cm.tab10(np.linspace(0, 0.9, len(SPECIFIC_ANGLES)))

    fig, ax = plt.subplots(figsize=(9, 6))

    for hname, angle, color in zip(hnames, SPECIFIC_ANGLES, colors):
        intensities = [data[hname].get(d, 0.0) for d in depths]
        if any(v > 0 for v in intensities):   # skip empty histograms
            ax.plot(depths, intensities, marker="o", label=f"{angle}°",
                    color=color, linewidth=1.8, markersize=7)

    energy_mev = float(energy_label)
    ax.set_xlabel("Depth (mm)", fontsize=13)
    ax.set_ylabel("Intensity (counts / primary)", fontsize=13)
    ax.set_title(
        f"Prompt Gamma Intensity vs Depth\n"
        f"({energy_mev:.1f} MeV gamma line, specific detector angles)",
        fontsize=13,
    )
    ax.legend(title="Detector angle", fontsize=10, title_fontsize=10, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.tight_layout()

    tag = energy_label.replace(".", "p")
    out = os.path.join(output_dir, f"PG_intensity_vs_depth_{tag}MeV.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


def plot_combined_spectra(root_files, output_dir):
    """
    One combined figure: overlay spectra from all depths for each broad angular bin.
    Creates one subplot per angular bin.
    """
    depths_sorted = sorted([(extract_depth(p), p) for p in root_files])
    cmap = cm.viridis(np.linspace(0.1, 0.9, len(depths_sorted)))

    n_bins = len(BROAD_ANGLE_HISTS)
    ncols = 3
    nrows = (n_bins + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4 * nrows), sharex=True)
    axes = axes.flatten()

    for ax_idx, hname in enumerate(BROAD_ANGLE_HISTS):
        ax = axes[ax_idx]
        for (depth, path), color in zip(depths_sorted, cmap):
            f = uproot.open(path)
            if hname in f:
                vals, edges = f[hname].to_numpy()
                centers = 0.5 * (edges[:-1] + edges[1:])
                ax.step(centers, vals, where="mid", color=color,
                        linewidth=1.2, label=f"{depth:.0f} mm")
        ax.set_title(BROAD_ANGLE_LABELS[hname], fontsize=11)
        ax.set_ylabel("dN/dE (/ primary / MeV)", fontsize=9)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.grid(True, alpha=0.25)

    for ax_idx in range(n_bins, len(axes)):
        axes[ax_idx].set_visible(False)

    # shared x label
    for ax in axes[-ncols:]:
        ax.set_xlabel("Energy (MeV)", fontsize=10)

    # legend on the last visible subplot
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Depth", loc="lower right",
               bbox_to_anchor=(1.0, 0.02), fontsize=9, title_fontsize=9)

    fig.suptitle("Prompt Gamma Spectra at Different Depths", fontsize=14, y=1.01)
    fig.tight_layout()

    out = os.path.join(output_dir, "PG_spectra_overlay_all_depths.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


def plot_summary_panel(data, depths, output_dir):
    """
    Summary panel: broad-angle intensity + 4.4 MeV + 9.6 MeV side by side.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # --- left: broad angular bins ---
    ax = axes[0]
    colors = cm.tab10(np.linspace(0, 0.9, len(BROAD_ANGLE_HISTS)))
    for hname, color in zip(BROAD_ANGLE_HISTS, colors):
        intensities = [data[hname].get(d, 0.0) for d in depths]
        ax.plot(depths, intensities, marker="o", label=BROAD_ANGLE_LABELS[hname],
                color=color, linewidth=1.8, markersize=6)
    ax.set_title("Full spectrum (1.5–12 MeV)", fontsize=11)
    ax.legend(title="Angle bin", fontsize=8, title_fontsize=8)
    ax.set_xlabel("Depth (mm)", fontsize=11)
    ax.set_ylabel("Intensity (counts / primary)", fontsize=11)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.grid(True, alpha=0.3)

    # --- middle: 4.4 MeV ---
    ax = axes[1]
    colors2 = cm.plasma(np.linspace(0.1, 0.9, len(SPECIFIC_ANGLES)))
    for angle, color in zip(SPECIFIC_ANGLES, colors2):
        hname = f"4.400000_MeV_Gamma_{angle}_deg"
        intensities = [data[hname].get(d, 0.0) for d in depths]
        if any(v > 0 for v in intensities):
            ax.plot(depths, intensities, marker="s", label=f"{angle}°",
                    color=color, linewidth=1.8, markersize=6)
    ax.set_title("4.4 MeV gamma line", fontsize=11)
    ax.legend(title="Angle", fontsize=8, title_fontsize=8, ncol=2)
    ax.set_xlabel("Depth (mm)", fontsize=11)
    ax.set_ylabel("Intensity (counts / primary)", fontsize=11)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.grid(True, alpha=0.3)

    # --- right: 9.6 MeV ---
    ax = axes[2]
    for angle, color in zip(SPECIFIC_ANGLES, colors2):
        hname = f"9.600000_MeV_Gamma_{angle}_deg"
        intensities = [data[hname].get(d, 0.0) for d in depths]
        if any(v > 0 for v in intensities):
            ax.plot(depths, intensities, marker="^", label=f"{angle}°",
                    color=color, linewidth=1.8, markersize=6)
    ax.set_title("9.6 MeV gamma line", fontsize=11)
    ax.legend(title="Angle", fontsize=8, title_fontsize=8, ncol=2)
    ax.set_xlabel("Depth (mm)", fontsize=11)
    ax.set_ylabel("Intensity (counts / primary)", fontsize=11)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.grid(True, alpha=0.3)

    fig.suptitle("Prompt Gamma Intensity vs Proton Range (Depth)", fontsize=14)
    fig.tight_layout()

    out = os.path.join(output_dir, "PG_intensity_vs_depth_summary.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default=".",
                        help="Directory containing PG_Spectrum_VS_Angle_*.root files")
    parser.add_argument("--output-dir", default=".",
                        help="Directory where output PNG files are written")
    args = parser.parse_args()

    # Allow passing explicit ROOT-file paths as positional arguments
    pattern = os.path.join(args.data_dir, "PG_Spectrum_VS_Angle_*.root")
    root_files = sorted(glob.glob(pattern), key=lambda p: extract_depth(p))

    if not root_files:
        print(f"No ROOT files found matching: {pattern}")
        return

    print(f"Found {len(root_files)} ROOT file(s):")
    for p in root_files:
        print(f"  depth={extract_depth(p):.1f} mm  ->  {os.path.basename(p)}")

    os.makedirs(args.output_dir, exist_ok=True)

    depths = [extract_depth(p) for p in root_files]
    data = load_data(root_files)

    print("\nGenerating plots ...")
    plot_broad_angles(data, depths, args.output_dir)
    for energy in GAMMA_ENERGIES:
        plot_specific_gamma(data, depths, energy, args.output_dir)
    plot_combined_spectra(root_files, args.output_dir)
    plot_summary_panel(data, depths, args.output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
