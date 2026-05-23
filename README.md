# RTSPG: 3D Crustal Velocity Model Evaluation via Reverse-Time Source-Point Gathers

Analysis scripts supporting the manuscript:

> *Evaluation of 3D crustal velocity models in the Sichuan–Yunnan region using
> reverse-time source-point gathers of ambient-noise surface waves.*
> (To be submitted to JGR: Solid Earth, 2026.)

This repository contains the MATLAB code for picking surface-wave traveltime
residuals (dt) from reverse-time source-point gathers (RTSPG) and for the
associated per-source, per-azimuth statistics used to assess four 3D crustal
velocity models of southwestern China:

- **SWChinaCVM-2.0** (`yao`)
- **USTClitho2.0** (`zhang`)
- **Bao20** (`bao`)
- **Feng20** (`feng`)

---

## Method overview

The RTSPG framework constructs source-point gathers by convolving 3D synthetic
Green's functions (from CGFD3D simulations with each candidate velocity model)
with time-reversed empirical Green's functions (EGFs) from ChinArray ambient-noise
cross-correlations. The peak of each RTSPG trace near zero time gives the dt
between observed and synthetic surface-wave arrivals.

Pipeline:

```
ChinArray ambient noise  -->  EGFs (CCFs)
                                  |
                                  v   (time-reverse)
CGFD3D synthetic GFs  --conv-->  RTSPG traces  --pick-->  dt values
                                                              |
                                                              v
                                          per-sector statistics, radar plot
```

This repository covers steps **from RTSPG dt picking onward**. The earlier
steps (CGFD3D forward simulation and ambient-noise EGF preparation) are
documented in an accompanying configuration script (`cgfd3d_config.sh`) and
rely on external software (CGFD3D, see https://github.com/zw-vis/gpu-CGFD3D).

---

## Repository layout

```
rtspg-3d-cvm-evaluation/
├── README.md                          ← this file
├── LICENSE                            ← GPL v3
├── scripts/
│   ├── rtspg_pick_dt.m                ← main dt picking routine
│   ├── analyze_dt_statistics.m        ← scatter + linear fit + radar chart
│   ├── run_example.m                  ← driver: how to call rtspg_pick_dt
│   └── run_analysis_example.m         ← driver: how to call analyze_dt_statistics
├── third_party/                       ← bundled external utilities
│   ├── README.md
│   ├── polarPcolor/                   ← Cheynet (BSD)
│   ├── radarChart/                    ← slandarer (BSD)
│   ├── wiggle/                        ← Hansen (GPL v2+)
│   └── jetwr/                         ← Katzman
└── data_example/                      ← sample data for source 179
    ├── StaX1.list                     ← 263 ChinArray station coordinates
    └── 179prefil_con/                 ← filtered RTSPG SAC traces
        ├── yaocon5-10/    *.sac       ← SWChinaCVM-2.0,  5–10 s
        ├── yaocon8-18/    *.sac       ← SWChinaCVM-2.0,  8–18 s
        ├── yaocon15-35/   *.sac       ← SWChinaCVM-2.0, 15–35 s
        ├── yaocon20-45/   *.sac       ← SWChinaCVM-2.0, 20–45 s
        ├── zhangcon5-10/  *.sac       ← USTClitho2.0,    5–10 s
        ├── zhangcon8-18/  *.sac       ← USTClitho2.0,    8–18 s
        ├── zhangcon15-35/ *.sac       ← USTClitho2.0,   15–35 s
        ├── zhangcon20-45/ *.sac       ← USTClitho2.0,   20–45 s
        ├── baocon5-10/    *.sac       ← Bao20,           5–10 s
        ├── baocon8-18/    *.sac       ← Bao20,           8–18 s
        ├── baocon15-35/   *.sac       ← Bao20,          15–35 s
        ├── baocon20-45/   *.sac       ← Bao20,          20–45 s
        ├── fengcon5-10/   *.sac       ← Feng20,          5–10 s
        ├── fengcon8-18/   *.sac       ← Feng20,          8–18 s
        ├── fengcon15-35/  *.sac       ← Feng20,         15–35 s
        └── fengcon20-45/  *.sac       ← Feng20,         20–45 s
```

The repository ships with a complete sample data subset for source
point 179 (all 4 velocity models × 4 period bands, ~4200 SAC files,
~10 MB total). This is sufficient to fully reproduce all figures
related to source 179. Data for the other five source points are
available on Zenodo (see "Data availability" below).

(The repository ships with one example source point; full data are available
on Zenodo, see "Data availability" below.)

---

## Dependencies

- **MATLAB R2019b or later** (uses `omitnan`, `exportgraphics`)
- **readsac** for reading SAC files (e.g. from the SAC Toolbox, not
  redistributed here; users should install it separately)

The following plotting utilities are bundled in `third_party/` (see
`third_party/README.md` for authors and licenses):

- `polarPcolor.m` — Etienne Cheynet (BSD)
- `radarChart.m` — slandarer / Zhaoxu Liu (BSD)
- `wiggle.m` — Thomas Mejer Hansen (GPL v2-or-later)
- `jetwr.m` — Rafi Katzman

---

## Quick start

1. Clone the repository and `cd` into it:
   ```bash
   git clone <repo-url>
   cd rtspg-3d-cvm-evaluation
   ```

2. In MATLAB, edit the paths at the top of `scripts/run_example.m` to point
   to your data, then run:
   ```matlab
   cd rtspg-3d-cvm-evaluation
   run scripts/run_example.m
   ```
   This processes one source point (#179) + one velocity model (`zhang`,
   i.e., USTClitho2.0) across four period bands (5–10, 8–18, 15–35,
   20–45 s), producing dt text files, gather diagnostic plots, and
   circle plots.

3. After running `rtspg_pick_dt` for all 4 models, generate per-source
   statistics with:
   ```matlab
   run scripts/run_analysis_example.m
   ```

---

## Input/output specification

### dt output text file (8 columns)

Each row corresponds to one station; missing or rejected picks are stored as
`NaN` in column 6.

| Column | Quantity            | Unit |
|-------:|---------------------|------|
| 1–5    | Station metadata (from `StaX1.list`) | — |
| 6      | dt = peak − t₀      | s    |
| 7      | Epicentral distance | m    |
| 8      | Sector azimuth      | °    |

### `StaX1.list` format

ASCII table with one row per station and 5 columns:

| Column | Quantity                       | Unit       |
|-------:|--------------------------------|------------|
| 1      | Longitude                      | degrees    |
| 2      | Latitude                       | degrees    |
| 3      | Station code (numeric)         | —          |
| 4      | Projected x (Easting)          | meters     |
| 5      | Projected y (Northing)         | meters     |

The dt picking uses columns 4–5 for sector geometry; the full row is
copied into columns 1–5 of the dt output file. See
`data_example/StaX1.list` for the 263-station ChinArray Phase II
deployment used in this study.

---

## Parameters

### Globally fixed parameters (in `rtspg_pick_dt.m`)

| Parameter       | Value   | Description |
|-----------------|---------|-------------|
| `t0`            | 299     | Zero-time sample index of RTSPG output |
| `nt`            | 597     | Total samples per RTSPG trace |
| `sc_ncan`       | 4       | Peak candidates retained per trace |
| `sc_ref_hw`     | 6       | Spatial-consistency reference window (traces) |
| `min_peak_amp`  | 0.2     | Normalized amplitude threshold for peak acceptance |

### Per-band parameters

| Band   | cHW (s) | dtmax (s) | scTol (s) |
|--------|--------:|----------:|----------:|
| 5–10   | 5       | 4         | 4.0       |
| 8–18   | 9       | 6         | 6.0       |
| 15–35  | 18      | 11        | 11.0      |
| 20–45  | 23      | 15        | 15.0      |

`cHW` = half-width of coarse peak search window;
`dtmax` = maximum allowed |dt|;
`scTol` = spatial consistency tolerance.

### Per-source parameters

A few parameters are set per source point inside `get_source_config()` in
`rtspg_pick_dt.m`:

- **Azimuthal sector widths.** The 4 cardinal sectors (N/S/W/E) use
  perpendicular half-widths that vary between 30 and 50 km among source
  points; the 8 diagonal sectors share a constant 30 km half-width. The
  cardinal half-widths are manually set per source to compensate for the
  uneven ChinArray station density in the Sichuan–Yunnan region: sectors
  with sparser local coverage are given wider bands to ensure enough
  stations per sector for the consistency-based pick refinement. This is
  an empirical data-driven choice and does not affect the picking
  algorithm; it only controls which stations are grouped into each
  azimuthal bin.

- **Circle-plot radial grid.** The polar pseudocolor plot uses a radial
  grid spanning the distance range of available stations around each
  source, with 23–30 radial nodes per source. Settings are listed in
  `get_source_config()`.

- **Applying the code to a new source point.** For any `source_id` not
  in the manuscript set (60, 104, 142, 179, 203, 246), `get_source_config()`
  falls through to a default branch that uses 40 km half-widths for all
  4 cardinal sectors and triggers automatic detection of the circle-plot
  radial grid from the observed data range. Users wanting tighter control
  can simply add a new `case` clause in `get_source_config()`.

---

## Picking algorithm (summary)

For each azimuthal sector, traces are sorted by epicentral distance and
processed sequentially:

1. **Candidate peaks.** In the window `[t0 − cHW, t0 + cHW]`, find up to
   `sc_ncan = 4` positive peaks above 0.05 normalized amplitude, ordered
   by amplitude descending.

2. **Initial coarse pick.** Take the largest-amplitude candidate.

3. **Spatial consistency check.** For each trace `i`, compute a reference
   dt as the median of dt values from the previous (up to `sc_ref_hw = 6`)
   valid traces in the same sector. If `|dt_i − ref_dt| > scTol`, scan
   candidates in amplitude-descending order and keep the first that
   satisfies the tolerance; if none does, set NaN.

4. **Amplitude check.** Reject (set NaN) picks whose normalized amplitude
   at the peak sample is below `min_peak_amp = 0.2`.

5. **`dtmax` rejection.** Reject picks with `|dt| > dtmax` for the band.

---

## Source-time function and time-delay correction

The CGFD3D forward simulations use a Gaussian source-time function with
`t0 = 1 s`, which imposes a 1 s onset delay on synthetic Green's functions.
This delay is corrected during RTSPG construction by interpolating the
convolved trace to a 1 s sampling interval, which aligns the zero-time of
the gather with the EGF arrival rather than the synthetic source onset.
See `cgfd3d_config.sh` for the CGFD3D source configuration and the
manuscript Methods section for a derivation.

---

## Data availability

This repository ships with a complete sample for one source point:

- `data_example/StaX1.list` — 263 ChinArray Phase II station coordinates
- `data_example/179prefil_con/` — filtered RTSPG SAC traces for source
  point 179, covering all 4 candidate velocity models (SWChinaCVM-2.0,
  USTClitho2.0, Bao20, Feng20) and all 4 period bands (5–10, 8–18,
  15–35, 20–45 s) — approximately 4,200 SAC files, ~10 MB.

Running `scripts/run_example.m` with this sample reproduces all
RTSPG dt picking outputs for source 179; running
`scripts/run_analysis_example.m` reproduces the per-source statistics
and radar chart for source 179.

Data for the other five source points reported in the manuscript
(60, 104, 142, 203, 246) are archived on Zenodo with DOI [TBD upon
submission]. After downloading from Zenodo, extract the archive so the
directory layout matches `<source_id>prefil_con/<model>con<band>/`,
update the paths in `run_example.m`, and re-run.

Raw ChinArray Phase II waveform data are available from the China
Seismic Array Data Management Centre at the Institute of Geophysics,
China Earthquake Administration.

---

## Citation

If you use this code, please cite:

> [Manuscript citation pending]

---

## Contact

Tong Li · Southern University of Science and Technology · lit7@sustech.edu.cn

---

## License

This repository is released under the **GNU General Public License v3.0**
(see `LICENSE`).

The choice of GPL v3 is required because the bundled `third_party/wiggle/wiggle.m`
is GPL v2-or-later: any project that incorporates GPL-licensed code must
itself be GPL-compatible. The other bundled utilities (`polarPcolor`,
`radarChart`, `jetwr`) are under BSD-compatible terms. See
`third_party/README.md` for per-component details.
