# Third-Party Code

This directory contains external MATLAB utilities used by the main
analysis scripts. Each subdirectory preserves the original code along
with the original author's name, source, and license notice.

| Subdirectory   | Author                          | Original License | Used by                          |
|----------------|---------------------------------|------------------|----------------------------------|
| `polarPcolor/` | Etienne Cheynet (Univ. Stavanger) | BSD (FEX default) | `draw_circle` in `rtspg_pick_dt.m` |
| `radarChart/`  | slandarer (Zhaoxu Liu)          | BSD (FEX default) | `analyze_dt_statistics.m`        |
| `wiggle/`      | Thomas Mejer Hansen             | **GPL v2 or later** | `plot_daoji` in `rtspg_pick_dt.m` |
| `jetwr/`       | Rafi Katzman (modified from MATLAB `jet`) | Algorithm reuse | `draw_circle` in `rtspg_pick_dt.m` |

## Why the whole repository is GPL v3

`wiggle.m` is licensed under GPL v2 or later. GPL is a "copyleft"
license: any project that incorporates GPL-licensed code must itself
be released under a compatible GPL license. Since this repository
includes `wiggle.m`, the repository as a whole is therefore released
under **GPL v3** (compatible with GPL v2-or-later).

The other three utilities (`polarPcolor`, `radarChart`, `jetwr`) are
under BSD-compatible terms and are fully compatible with GPL v3.

## No modifications

The files in this directory are redistributed verbatim from their
original sources, with the original author attribution preserved in
each file's header comments. No functional modifications have been
made by the repository maintainers.

## How to add new third-party utilities

1. Create a new subdirectory under `third_party/`.
2. Place the file(s) inside, keeping the original header comments.
3. Add a short `README.md` in the subdirectory noting the author,
   source URL, and original license.
4. Verify the license is compatible with GPL v3 before adding.
