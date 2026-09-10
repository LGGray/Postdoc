# Rebuild existing_figures/snRNA_Xist_by_LOX_VCM_adult_aged.png from the cluster PDF.
#
# The slide panel is panel b of core_escape_block_new_AR_Xist_umap_panel_VCM.pdf
# (04_core_escape.R), cropped to the 9w and 78w facets: Xist expression in
# LOX-like vs other ventricular cardiomyocytes, Wilcoxon p in the strip.
# There is no standalone PDF for it, hence the crop.
#
#   python3 crop_xist_lox_panel.py [RESULTS_ROOT]
#
# RESULTS_ROOT defaults to the doublet-free tree; pass Allelic_ratio_results
# for the version that still contains doublets.
import os, subprocess, sys, tempfile
import numpy as np
from PIL import Image

CL = os.environ.get("CLUSTER_MOUNT", "/Users/graylachlan/cluster")
RES = sys.argv[1] if len(sys.argv) > 1 else "Allelic_ratio_results_nodoublet"
SRC = os.path.join(CL, "OCM", RES, "core_escape_cutoff_5",
                   "core_escape_block_new_AR_Xist_umap_panel_VCM.pdf")
OUT = os.path.join(os.environ.get("FIG_EXISTING",
    "/Users/graylachlan/LRZ Sync+Share/LGray/Presentations/IPT_Seminar_2026-09-15/existing_figures"),
    "snRNA_Xist_by_LOX_VCM_adult_aged.png")
KEEP_FACETS = 2          # 9w, 78w (the deck is adult vs aged)
DPI = 300

with tempfile.TemporaryDirectory() as td:
    subprocess.run(["pdftoppm", "-png", "-r", str(DPI), SRC, os.path.join(td, "p")], check=True)
    page = os.path.join(td, "p-1.png")
    im = Image.open(page)
    a = np.array(im.convert("L")); H, W = a.shape

    # panel b begins after the tallest all-white band in the middle of the page
    blank = (a < 200).sum(axis=1) == 0
    runs, s = [], None
    for i, b in enumerate(blank):
        if b and s is None: s = i
        elif not b and s is not None: runs.append((s, i - s)); s = None
    if s is not None: runs.append((s, len(blank) - s))
    mid = [r for r in runs if H * 0.25 < r[0] < H * 0.75]
    if not mid:
        sys.exit("could not find the panel a/b separator; layout changed?")
    gap = max(mid, key=lambda r: r[1])
    top = gap[0] + gap[1] // 2

    panel = a[top:, :] < 200
    cols = panel.any(axis=0)
    runs, s = [], None
    for i, v in enumerate(~cols):
        if v and s is None: s = i
        elif not v and s is not None: runs.append((s, i - s)); s = None
    if s is not None: runs.append((s, len(cols) - s))
    gutters = [r for r in runs if r[1] >= 8 and 200 < r[0] < W - 100]
    # gutters[0] is the gap between the y axis and the first facet; the facet
    # separators are the ones after it.
    facet_gaps = gutters[1:]
    if len(facet_gaps) < KEEP_FACETS:
        sys.exit(f"expected >= {KEEP_FACETS} facet separators, found {len(facet_gaps)}")
    cut = facet_gaps[KEEP_FACETS - 1]
    right = cut[0] + cut[1] // 2

    # drop the "b" panel tag: first ink run at the far left, before the axis title
    tag = np.where(panel[:120, :300].any(axis=0))[0]
    left = 0
    if len(tag):
        brk = tag[0]
        for c in tag[1:]:
            if c - brk > 5: break
            brk = c
        left = brk + 6

    rows = np.where(panel.any(axis=1))[0]
    im.crop((left, top + rows.min() - 8, right, top + rows.max() + 16)).save(OUT)
    print(f"{RES} -> {OUT}")
