# Append the results slides to the intro deck, using the deck's own layouts.
#   python3 presentations/ipt_seminar_2026-09-15/build_deck.py
import copy
import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.oxml.ns import qn
from lxml import etree
from PIL import Image

SRC = "/Users/graylachlan/LRZ Sync+Share/LGray/Presentations/IPT Seminar 15.9.26.pptx"
BASE = "/Users/graylachlan/LRZ Sync+Share/LGray/Presentations/IPT_Seminar_2026-09-15"
OUT = os.path.join(BASE, "IPT Seminar 15.9.26 - v2.pptx")
FIG = os.path.join(BASE, "figures")
EX = os.path.join(BASE, "existing_figures")

prs = Presentation(SRC)
L_TITLE_ONLY = prs.slide_layouts[5]
L_SECTION = prs.slide_layouts[2]
GREY = RGBColor(0x40, 0x40, 0x40)
ACCENT = RGBColor(0xC0, 0x39, 0x2B)

# content area below the title placeholder (title: top 0.4in, height 1.45in)
TOP, BOTTOM, LEFT, RIGHT = 1.9, 7.2, 0.5, 12.85
H = BOTTOM - TOP  # 5.3 in of content height


def fig(name):
    for d in (FIG, EX):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    return None


MAX_PX = 2600  # longest side of the embedded copy; originals stay at 300 dpi in figures/
EMBED = os.path.join(BASE, "data", "embedded"); os.makedirs(EMBED, exist_ok=True)


def embed_copy(path):
    """Downscaled copy for the deck so the .pptx stays a sensible size."""
    im = Image.open(path)
    if max(im.size) <= MAX_PX:
        return path
    scale = MAX_PX / max(im.size)
    out = os.path.join(EMBED, os.path.basename(path))
    im.convert("RGB").resize((round(im.size[0] * scale), round(im.size[1] * scale)), Image.LANCZOS).save(out, optimize=True)
    return out


def fit_picture(slide, path, left, top, width, height, align="center"):
    """Place an image inside a box, preserving aspect ratio."""
    path = embed_copy(path)
    im = Image.open(path)
    ar = im.size[0] / im.size[1]
    if width / height > ar:
        ph, pw = height, height * ar
    else:
        pw, ph = width, width / ar
    l = left + ((width - pw) / 2 if align == "center" else 0)
    t = top + ((height - ph) / 2 if align == "center" else 0)
    return slide.shapes.add_picture(path, Inches(l), Inches(t), Inches(pw), Inches(ph))


def placeholder_box(slide, text, left, top, width, height):
    box = slide.shapes.add_shape(1, Inches(left), Inches(top), Inches(width), Inches(height))
    box.fill.solid(); box.fill.fore_color.rgb = RGBColor(0xF2, 0xF2, 0xF2)
    box.line.color.rgb = RGBColor(0xBF, 0xBF, 0xBF)
    tf = box.text_frame; tf.word_wrap = True
    tf.paragraphs[0].text = text
    for r in tf.paragraphs[0].runs:
        r.font.size = Pt(14); r.font.color.rgb = GREY
    return box


def picture_or_placeholder(slide, name, left, top, width, height, note=None):
    p = fig(name)
    if p:
        return fit_picture(slide, p, left, top, width, height)
    return placeholder_box(slide, note or f"[{name} - not generated yet]", left, top, width, height)


def add_text(slide, items, left, top, width, height, size=15, bullets=True, space=5):
    """items: list of (text, level). level -1 = bold header without bullet, 0/1 = bullet levels."""
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = tb.text_frame; tf.word_wrap = True
    for i, (text, level) in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        run = p.add_run(); run.text = text
        p.space_after = Pt(space)
        if level < 0:
            run.font.bold = True; run.font.size = Pt(size)
            run.font.color.rgb = GREY
            if i > 0:
                p.space_before = Pt(8)
            continue
        run.font.size = Pt(size - 2 * level)
        if bullets:
            pPr = p._p.get_or_add_pPr()
            indent = 0.22
            pPr.set("marL", str(int(Inches(indent + 0.3 * level))))
            pPr.set("indent", str(-int(Inches(indent))))
            bu = etree.SubElement(pPr, qn("a:buChar")); bu.set("char", "•" if level == 0 else "–")
    return tb


def add_caption(slide, text, left, top, width, size=11):
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(0.4))
    tf = tb.text_frame; tf.word_wrap = True
    r = tf.paragraphs[0].add_run(); r.text = text
    r.font.size = Pt(size); r.font.italic = True; r.font.color.rgb = GREY
    return tb


def add_table(slide, rows, left, top, width, col_widths=None, size=12, header=True):
    nrows, ncols = len(rows), len(rows[0])
    shp = slide.shapes.add_table(nrows, ncols, Inches(left), Inches(top), Inches(width), Inches(0.3 * nrows))
    tbl = shp.table
    if col_widths:
        for j, w in enumerate(col_widths):
            tbl.columns[j].width = Inches(w)
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = tbl.cell(i, j)
            cell.text = str(val)
            for p in cell.text_frame.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(size)
                    r.font.bold = (i == 0 and header)
            cell.margin_top = cell.margin_bottom = Inches(0.03)
    return shp


def new_slide(title, notes=None, layout=None):
    s = prs.slides.add_slide(layout or L_TITLE_ONLY)
    s.shapes.title.text = title
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    return s


def move_slide(index_from, index_to):
    """python-pptx can only append; reorder the sldIdLst to place a slide."""
    lst = prs.slides._sldIdLst
    ids = list(lst)
    lst.remove(ids[index_from])
    lst.insert(index_to, ids[index_from])


def section(title, subtitle, notes=None):
    s = prs.slides.add_slide(L_SECTION)
    s.shapes.title.text = title
    s.placeholders[1].text = subtitle
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    return s


# ---------------------------------------------------------------------------
# Intro slides 2 and 5 (existing): wording carried over from the hand-edited
# 2026-09-09 deck. Edit here, not in PowerPoint, or the next build reverts it.
# ---------------------------------------------------------------------------
def replace_line(shape, old, new):
    """Swap one line of an existing placeholder, keeping its run formatting."""
    for para in shape.text_frame.paragraphs:
        if "".join(r.text for r in para.runs).strip() == old:
            for r in para.runs[1:]:
                r._r.getparent().remove(r._r)
            para.runs[0].text = new
            return True
    return False


for sh in prs.slides[1].shapes:
    if sh.has_text_frame:
        replace_line(sh, "Prior to menopause, males have a higher risk of cardiovascular disease including:",
                         "Males have a higher lifetime risk of cardiovascular disease including:")

def set_lines(shape, lines):
    """Rewrite a placeholder as one paragraph of soft-broken lines.

    PowerPoint's in-paragraph break is an <a:br/> element, not a \\x0b
    character: assigning the character makes lxml escape it to a literal
    "_x000B_". An empty string in `lines` is a blank line.
    """
    tf = shape.text_frame
    para = tf.paragraphs[0]
    src = para.runs[0]._r if para.runs else None
    rPr = copy.deepcopy(src.find(qn("a:rPr"))) if src is not None and src.find(qn("a:rPr")) is not None else None
    for extra in tf.paragraphs[1:]:
        extra._p.getparent().remove(extra._p)
    pel = para._p
    for child in list(pel):
        if child.tag != qn("a:pPr"):
            pel.remove(child)
    for i, line in enumerate(lines):
        if i:
            pel.append(pel.makeelement(qn("a:br"), {}))
        if not line:
            continue
        r = pel.makeelement(qn("a:r"), {})
        if rPr is not None:
            r.append(copy.deepcopy(rPr))
        t = pel.makeelement(qn("a:t"), {})
        t.text = line
        r.append(t)
        pel.append(r)


# Slide 5: three aims, with the multiome as Aim 3.
set_lines(prs.slides[4].shapes.title, [
    "Research question",
    "How does aging influence XCI in the heart at snRNA-seq, spatial, and multiome resolution?",
    "",
    "Research Aims",
    "Aim 1: Quantify XCI escape in snRNA-seq of adult and aged hearts",
    "",
    "Aim 2: Determine whether XCI escape is spatially structured across the myocardium",
    "",
    "Aim 3: Test whether the inactive X gains distal chromosome accessibility with age, and in which cell types",
    "",
])

# ---------------------------------------------------------------------------
# Slide 6 (existing): swap the four-sample UMAP + QC panels for adult/aged versions
# ---------------------------------------------------------------------------
s6 = prs.slides[5]
if fig("F21_snRNA_UMAP_celltypes.png") and fig("F22_snRNA_QC_panels.png"):
    for sh in list(s6.shapes):
        if sh.shape_type == 13:  # picture
            sh._element.getparent().remove(sh._element)
    fit_picture(s6, fig("F21_snRNA_UMAP_celltypes.png"), 0.3, 1.6, 5.6, 5.6)
    fit_picture(s6, fig("F22_snRNA_QC_panels.png"), 6.1, 1.6, 7.0, 5.7)
    s6.notes_slide.notes_text_frame.text = ("Adult and aged nuclei only (Sham/TAC libraries processed in the same run are kept for the aside). "
        "Marker-based cell types on the shared UMAP; QC metrics after ddqcR filtering; cell-type proportions per sample.")
s6.shapes.title.text = "Overview of snRNA-seq data"

# Aim 1 divider, ahead of the overview slide (hand-edit from the 2026-09-09 deck).
_d = new_slide("Aim 1: Quantify XCI escape in snRNA-seq of adult and aged hearts")
move_slide(len(prs.slides._sldIdLst) - 1, 5)

# ===========================================================================
# AIM 1 - snRNA-seq
# ===========================================================================
s = new_slide("Allele-specific analysis per nucleus",
    "Allelome.PRO2 was run per nucleus against the B6/CAST SNP set (Xist masked). The whole-chromosome allelic ratio is B6 reads over all SNP-overlapping reads. "
    "Adult and aged only. Left: 85-88% of nuclei clear the 30-read cutoff. Right: autosomes sit at 0.53, the B6-ward mapping bias; chrX at ~0.95 because the CAST X is the inactive X in every nucleus. "
    "Any CAST read on chrX is expression from the inactive X.")
picture_or_placeholder(s, "F01_snRNA_chrX_reads_per_nucleus.png", LEFT, TOP, 6.4, 3.6)
picture_or_placeholder(s, "F02_snRNA_autosome_vs_chrX.png", 7.0, TOP, 5.85, 3.6)
add_text(s, [
    ("Allelome.PRO2 per nucleus against ~20 M B6/CAST SNPs (Xist masked); allelic ratio = B6 / (B6 + CAST)", 0),
    ("scDblFinder doublets excluded before the allelic analysis (9.1% of nuclei flagged)", 0),
    ("≥ 30 SNP-overlapping chrX reads keeps ~85-88% of nuclei; autosomes sit at 0.53 (B6 mapping bias), chrX at ~0.95", 0),
    ("CAST reads on chrX are expression from the inactive X: the Xist-deleted B6 X is active in every nucleus", 0),
], LEFT, 5.65, 12.3, 1.5, size=14)

s = new_slide("Whole-chrX allelic ratio by cell type",
    "Adult and aged only; Sham/TAC are kept for the aside. Every cell type is dominated by monoallelic nuclei (median AR 0.91-0.97) with a tail of biallelic nuclei below the 0.9 boundary. "
    "Beta-binomial dispersion tests flagged 78w vs 9w differences in endothelial cells, macrophages and pericytes, but with one animal per condition the simulated false-positive rate of that test is 8-81%, so read this as descriptive.")
picture_or_placeholder(s, "F03_snRNA_chrX_AR_violin_by_celltype.png", LEFT, TOP, 12.35, H)

s = new_slide("Biallelic nuclei per cell type",
    "Fraction of nuclei with whole-chrX AR below 0.9, scDblFinder doublets excluded. Adult to aged: endothelial 21% to 44%, pericytes/SMC 16% to 33%, "
    "fibroblasts 15% to 23%, ventricular CM 7% to 11%. Endocardium does not move (16% to 17%). Doublet removal drops 5.9% of nuclei at this cutoff and lowers "
    "every escape estimate by 1-6 points; the age direction is unchanged. One animal per age.")
picture_or_placeholder(s, "F04_snRNA_fraction_biallelic_nuclei.png", LEFT, TOP, 8.3, H)
add_text(s, [
    ("Adult → aged", -1),
    ("Endothelial 21% → 44%", 0), ("Pericytes / SMC 16% → 33%", 0), ("Fibroblasts 15% → 23%", 0), ("Ventricular CM 7% → 11%", 0),
    ("Endocardium unchanged 16% → 17%", 0),
    ("Lymphocytes: too few nuclei to read", 0),
    ("Caveat", -1),
    ("One animal per age: descriptive, not a test of age", 0),
    ("Doublets removed (scDblFinder): 5.9% of nuclei", 0),
    ("Depth-dependent", 0),
], 9.0, TOP, 3.9, H, size=14)

s = new_slide("Whole-chrX allelic ratio on the UMAP",
    "Same nuclei projected on the shared UMAP, coloured by whole-chrX allelic ratio. Biallelic nuclei are scattered through every cluster rather than forming a sub-cluster; "
    "the endothelial and pericyte clusters of the aged heart carry visibly more green.")
picture_or_placeholder(s, "F06_snRNA_AR_umap_by_sample.png", LEFT, TOP, 12.35, H)

s = new_slide("Caveat: apparent escape tracks allelic depth",
    "Adult and aged only; the same gradient is present in Sham and TAC (figures/sham_tac). Pooled B6 fraction rises from ~0.85 at 10-25 chrX reads to ~0.97 above 200 reads, within every cell type. Shallow nuclei carry proportionally more ambient RNA, "
    "which is a mosaic of both alleles. This is why results are only comparable within one read cutoff, and why WP1 of the grant proposes deep sequencing on NovaSeq X 1.5B flow cells.")
picture_or_placeholder(s, "F05_snRNA_depth_bias.png", LEFT, TOP, 12.35, 3.8)
add_text(s, [
    ("Pooled B6 fraction rises from ~0.85 at 10-25 chrX reads to ~0.97 above 200 reads, inside every cell type", 0),
    ("Shallow nuclei carry proportionally more ambient RNA, a mosaic of both alleles, so they look biallelic", 0),
    ("Comparisons are only valid within one cutoff; deeper allelic coverage (grant WP1) is the fix, not a different threshold", 0),
], LEFT, 5.85, 12.3, 1.3, size=14)

s = new_slide("Per-gene allelic ratios by cell type",
    "Reads pooled within each cell type (>= 20 SNP-overlapping reads per gene), genes ordered along chrX, restricted to genes testable in at least six cell types in both ages so the slide stays readable; "
    "the full 300-gene version is in the figures folder. Red = monoallelic (silenced), green/blue = expression from both X chromosomes; known escapees in bold. "
    "This is the same view as Figure 2 of the grant: escape can be shared across cell types (Kdm6a, Kdm5c, Eif2s3x, Pbdc1, 5530601H04Rik) or restricted. "
    "Gene-level depth is the limit: most genes are testable only in cardiomyocytes and fibroblasts, which is the WP1 argument for deep sequencing.")
picture_or_placeholder(s, "F09_snRNA_chrX_gene_heatmap_adult_aged.png", LEFT, TOP, 12.35, H - 0.6,
                       note="Per-gene pseudobulk heatmap (chrX gene x cell type, adult vs aged) - generated once the per-nucleus tables are consolidated")
add_caption(s, "Same view as grant Figure 2. Known escapees in bold; escape shared across cell types (Kdm6a, Kdm5c, Eif2s3x, Pbdc1) or restricted to a few. Genes testable in >= 6 cell types shown; full 300-gene version in the figures folder.",
            LEFT, BOTTOM - 0.55, 12.35, size=12)

s = new_slide("Core escape genes are biallelic everywhere",
    "Kdm5c, Kdm6a, Ddx3x and Eif2s3x, the canonical escapees, sit at allelic ratios of 0.5-0.75 in every cell type and every condition. "
    "Posterior means use a beta prior centred on the Hoelzl et al. 2025 bulk medians so that thin cell types are shrunk toward the bulk value. The assay validates itself here.")
picture_or_placeholder(s, "F07_snRNA_core_escape_posterior.png", LEFT, TOP, 12.35, 4.2)
add_text(s, [
    ("Canonical escapees sit at AR 0.5-0.75 in every cell type and condition: the assay recovers known biology", 0),
    ("Posterior mean with a prior centred on the Hoelzl et al. 2025 bulk medians, so thin cell types are shrunk toward bulk", 0),
], LEFT, 6.2, 12.3, 1.0, size=14)

s = new_slide("Core escape genes: adult vs aged",
    "CAST (inactive X) fraction of reads pooled over the four core escape genes. Macrophages show less inactive-X expression in the aged heart (FDR 0.003) and endothelial cells slightly more (FDR 0.09); "
    "every other cell type is unchanged. Escape at these genes is already near its ceiling, so age effects, if real, must be cell-type specific and small.")
picture_or_placeholder(s, "F08_snRNA_core_escape_adult_vs_aged.png", LEFT, TOP, 7.6, H)
add_text(s, [
    ("Inactive-X fraction pooled over Kdm5c, Kdm6a, Ddx3x, Eif2s3x", 0),
    ("Macrophages: less inactive-X expression in aged (FDR 0.003)", 0),
    ("Endothelial: slightly more (FDR 0.09); all other cell types unchanged", 0),
    ("Escape at these genes is near its ceiling (~40% CAST), so age effects can only be cell-type specific and small", 0),
    ("Fisher test on pooled reads ignores nucleus-level overdispersion; n = 1 per age", 0),
], 8.4, TOP, 4.45, H, size=14)

s = new_slide("Xist expression is lower in biallelic nuclei",
    "Core-escape-block analysis. Left: ventricular cardiomyocytes, Xist (SCT-normalised) in nuclei with a LOX-like allelic ratio versus the rest; lower in both ages (p = 0.02 adult, 1e-5 aged). "
    "Right: beta-binomial model of the core-escape-block allelic ratio on Xist per cell type; an odds ratio below 1 means more Xist goes with a lower B6 fraction, i.e. less inactive-X expression. "
    "Significant in ventricular CM in both ages; the other cell types trend the same way but are too thin. Loss of Xist coincides with inactive-X expression, which is the grant's premise.")
picture_or_placeholder(s, "snRNA_Xist_by_LOX_VCM_adult_aged.png", LEFT, TOP, 5.2, 3.6)
picture_or_placeholder(s, "F20_snRNA_Xist_vs_AR_forest.png", 5.9, TOP, 6.95, 3.6)
add_text(s, [
    ("Ventricular CM: nuclei with a LOX-like allelic ratio have lower Xist in both ages (p = 0.02 adult, 1e-5 aged)", 0),
    ("Beta-binomial per cell type: odds ratio < 1 per unit Xist, i.e. more Xist → more monoallelic; significant in ventricular CM in both ages, other cell types too thin", 0),
    ("Loss of Xist coincides with inactive-X expression, as the model predicts", 0),
], LEFT, 5.65, 12.3, 1.5, size=14)
add_caption(s, "Left: ventricular cardiomyocytes only", LEFT, 5.3, 5.2, size=11)

s = new_slide("Pressure overload influences XCI escape",
    "Not part of the grant, but the same animals and pipeline, kept separate from the aging story. TAC induces the hypertrophic programme (Nppa, Nppb, Ankrd1, Myh7) in ventricular cardiomyocytes "
    "and raises the fraction of biallelic ventricular CM nuclei from 11% to 18%; other cell types move little. Cardiac stress without aging can therefore also perturb XCI stability. One animal per condition. "
    "Sham/TAC versions of every snRNA-seq figure are in figures/sham_tac/ if anyone asks.")
picture_or_placeholder(s, "snRNA_TAC_hypertrophy_markers.png", LEFT, TOP, 5.6, H)
picture_or_placeholder(s, "sham_tac/F04_snRNA_fraction_biallelic_nuclei.png", 6.3, TOP, 6.55, 3.4)
add_text(s, [
    ("Same pipeline, one Sham and one TAC animal; not part of the aging story", 0),
    ("Nppa, Nppb, Ankrd1, Myh7 up in TAC ventricular CM: the surgery worked", 0),
    ("Biallelic ventricular CM nuclei 11% (Sham) → 18% (TAC); other cell types move little", 0),
    ("Cardiac stress without aging can also perturb XCI stability", 0),
], 6.3, 5.4, 6.55, 1.8, size=13)

# ===========================================================================
# AIM 2 - spatial
# ===========================================================================
section("Aim 2: Determine whether XCI escape is spatially structured across the myocardium",
        "Visium HD 3′ · 8 µm bins · one adult (9w) and one aged (78w) section · same F1 model",
        "Aim 2 asks whether escape is dispersed, clonal, or concentrated in niches. Reminder of the genotype: XCI is fully skewed, so there is no mosaic and no clonal 'which X' patches by construction; "
        "any spatial structure would be structure of escape itself.")

s = new_slide("Visium HD sections and cell-type annotation",
    "Left: tissue images. Right: 8 µm bins labelled by the best-matching marker set. Both sections are mid-ventricular transverse cuts; the 78w section has a larger atrial / vessel fragment attached.")
for i, smp in enumerate(["9w", "78w"]):
    top = TOP + i * 2.7
    picture_or_placeholder(s, f"spatial_tissue_hires_{smp}.png", LEFT, top, 3.0, 2.6)
    picture_or_placeholder(s, f"spatial_celltype_map_{smp}.png", 3.7, top, 5.6, 2.6)
    add_caption(s, "Adult (9w)" if smp == "9w" else "Aged (78w)", LEFT, top + 2.35, 3.0, size=12)
add_text(s, [
    ("Visium HD 3′, 2 µm capture, analysed at 8 µm bins", 0),
    ("Bins labelled by marker-set enrichment (dominance, not purity)", 0),
    ("Both sections are transverse ventricular cuts; differences are confounded with section plane and depth", 0),
], 9.6, TOP, 3.3, 5.4, size=13)

s = new_slide("Section composition",
    "Over 60% of bins are ventricular cardiomyocyte in both sections; non-myocyte labels are 2-3% each. The 78w section has slightly fewer fibroblast and endothelial bins. One section per age, so this is descriptive.")
picture_or_placeholder(s, "F10_spatial_composition.png", LEFT, TOP, 12.35, 4.9)

s = new_slide("Choosing the tile size",
    "Rather than cell segmentation (not designed for cardiac tissue), allelic counts are aggregated over square tiles. At 64 µm, 92% (9w) and 72% (78w) of tissue tiles reach 10 informative chrX UMIs, "
    "but the per-tile standard error on the escape fraction is still 0.06-0.09, so a single tile can only detect escape differences of 17-24%. 64 µm is a coverage choice, not a patch scale.")
picture_or_placeholder(s, "F11_spatial_tile_precision.png", LEFT, TOP, 12.35, 4.0)
add_text(s, [
    ("Allelic counts aggregated over square tiles rather than segmented cells (segmentation is not built for myocardium)", 0),
    ("At 64 µm, 92% (9w) / 72% (78w) of tissue tiles reach 10 informative chrX UMIs; per-tile SE on escape is still 0.06-0.09", 0),
    ("64 µm is a coverage choice, not a patch scale; the tile size sweep found no scale to read off", 0),
], LEFT, 6.0, 12.3, 1.2, size=13)

s = new_slide("Allelic ratio maps at 64 µm",
    "Top: allelic ratio (B6 fraction) on chrX per tile, on the same 11-level colour scale as the snRNA-seq UMAPs; red = monoallelic, green/blue = inactive-X expression. Bottom: the same for autosomes, which sit at 0.5 everywhere. "
    "Why the tile median is 0.75 (9w) / 0.82 (78w), pooled 0.87, rather than the 0.95 of the snRNA-seq nuclei: (1) the tile count is chromosome-wide, and 60% of the CAST molecules on chrX come from a handful of non-genic / multicopy loci that read ~100% CAST (the artefact story two slides on); gene-body-only the section is at 0.94, matching the nuclei. In the snRNA-seq the same split barely matters (0.95 chromosome-wide vs 0.98 gene-body), so these loci are far more prominent in the spatial libraries, which capture cytoplasmic RNA. "
    "(2) A tile holds ~30 informative UMIs, so the binomial spread alone is +/- 0.06, and the autosomal control shows tile ratios are 15-35x overdispersed on top of that; hence the wide violin, and the tiles at exactly 1.0 in the shallower 78w section. "
    "No visible patches. The two sections are not depth matched, so do not read this as an age comparison.")
picture_or_placeholder(s, "F12_spatial_tile_maps_64um.png", LEFT, TOP, 6.4, H)
picture_or_placeholder(s, "F13_spatial_tile_distribution.png", 7.1, TOP, 5.75, 2.9)
add_text(s, [
    ("Same colour scale as the snRNA-seq UMAPs; autosomal tiles sit at 0.5, chrX tiles are B6-dominated everywhere, no patches", 0),
    ("Tile median 0.75 / 0.82 (pooled 0.87) vs 0.95 in nuclei: the tile count is chromosome-wide and includes the artefact loci that carry 60% of chrX CAST molecules (gene-body only: 0.94)", 0),
    ("~30 UMIs per tile, 15-35x overdispersed: the wide spread is noise, not biology; sections are not depth matched", 0),
], 7.1, 4.9, 5.75, 2.3, size=12)

s = new_slide("No spatial structure of escape",
    "C(d) is the probability that two chrX UMIs at distance d carry the same allele. It is flat from 4 µm to 2 mm in both sections and equals the no-structure value p^2 + (1-p)^2 for a global escape "
    "fraction p = 0.873 (expected 0.778, observed 0.787 at 9w; 0.779 vs 0.778 at 78w). The residual against a permutation null is zero within error; Moran's I on the tile ratio is 0.006 / 0.003. "
    "Imprinted paternal loci give C near 1, so the method resolves monoallelic expression within cells. This is what non-random XCI predicts: escape is dispersed, not clonal.")
picture_or_placeholder(s, "F14_spatial_pair_correlation.png", LEFT, TOP, 8.2, H)
add_text(s, [
    ("C(d) = P(two chrX UMIs at distance d share an allele)", 0),
    ("Flat from 4 µm to 2 mm in both sections", 0),
    ("Equals the no-structure value p² + (1−p)²: expected 0.778, observed 0.787 (9w); 0.779 vs 0.778 (78w)", 0),
    ("Residual vs permutation null ≈ 0; Moran's I on tile ratios 0.006 / 0.003", 0),
    ("Imprinted paternal loci give C ≈ 1: the assay does resolve monoallelic expression", 0),
    ("As the skewed-XCI model predicts: escape is dispersed, not clonal", 0),
], 8.9, TOP, 3.95, H, size=13)

s = new_slide("Imprinted loci validate the allele calls",
    "Paternally expressed loci read CAST (Snrpn: 2 of 2,070 UMIs on the wrong allele, 0.06%); maternally expressed loci read B6 (H19, Igf2r, Rian: 0 of 194 wrong). "
    "That bounds the rate at which a B6 molecule could be miscalled as CAST, the error that would fake escape, at <= 1.5%. Cdkn1c, Mest and Impact are biallelic in adult heart and were excluded from the control set.")
picture_or_placeholder(s, "F15_spatial_imprinted_controls.png", LEFT, TOP, 12.35, 4.1)
add_text(s, [
    ("Paternal loci read CAST (Snrpn: 2 of 2,070 UMIs wrong, 0.06%); maternal loci read B6 (H19, Igf2r, Rian: 0 of 194 wrong)", 0),
    ("A B6 molecule miscalled as CAST is the error that would fake escape: its rate is ≤ 1.5% (95% upper bound)", 0),
    ("Cdkn1c, Mest, Impact are biallelic in adult heart and were dropped from the control set", 0),
], LEFT, 6.1, 12.3, 1.1, size=13)

s = new_slide("Where does the chrX CAST signal come from?",
    "CAST fraction per 100 kb window along chrX. A handful of windows carry most of the CAST molecules: chrX:11.5 Mb, with no annotated gene, is 5 B6 vs 6,201 CAST at 9w, 42% of all chrX CAST molecules. "
    "Slc16a2, Aff2 and the Gm14719 / olfactory cluster read ~100% CAST. Escape cannot exceed 50% (the inactive X cannot out-express the active one), so these are allele-call artefacts. "
    "The leading explanation: the SNP set is C57BL/6NJ x CAST while the reference is C57BL/6J, so wherever 6NJ differs from 6J a genuine B6 read scores as CAST.")
picture_or_placeholder(s, "F16_spatial_chrX_window_scan.png", LEFT, TOP, 8.4, H)
add_text(s, [
    ("A few 100 kb windows carry most CAST molecules", 0),
    ("chrX:11.5 Mb, no annotated gene: 5 B6 vs 6,201 CAST at 9w = 42% of all chrX CAST molecules", 0),
    ("Slc16a2, Aff2, Gm14719: ~100% CAST. Escape cannot exceed 50%, so these are allele-call artefacts", 0),
    ("Likely cause: SNP set is C57BL/6NJ × CAST, reference is C57BL/6J. Where 6NJ ≠ 6J a real B6 read scores as CAST", 0),
    ("Present identically in both sections, so it cancels in any comparison but inflates the absolute estimate", 0),
], 9.1, TOP, 3.75, H, size=13)

s = new_slide("Chromosome-wide escape is ~3%",
    "Splitting the same molecules by whether they fall inside an annotated gene body: non-genic chrX molecules are 16% of chrX but carry 60% of the CAST signal, sitting at the autosomal value. "
    "Inside gene bodies escape is 6.1%, and 2.8-3.0% once the seven genes above 50% CAST are removed. Identical in both sections. An artefact SNP mask is being built from the autosomal false-positive rate.")
picture_or_placeholder(s, "F17_spatial_escape_partition.png", LEFT, TOP, 7.8, 3.5)
add_text(s, [
    ("Non-genic chrX molecules (16% of chrX) sit at the autosomal value and carry 60% of the CAST signal", 0),
    ("Inside gene bodies: 6.1%; minus the seven >50% genes: 2.8% (9w), 3.0% (78w)", 0),
    ("Identical between sections: no chromosome-wide age effect at this depth", 0),
    ("Fix in progress: artefact SNP mask calibrated on the autosomes (biallelic by construction)", 0),
    ("The per-gene view is the honest readout →", 0),
], 8.5, TOP, 4.35, H, size=13)
add_text(s, [("Values from the gene-body split of the same UMIs (spatial/NEXT_ANALYSIS.md, 2026-09-03)", 0)],
         LEFT, H, 7.8, 0.5, size=11, bullets=False)

s = new_slide("Per-gene escape in the spatial data",
    "spASE single-cell model on 16 µm pixels, chrX genes with >= 100 informative UMIs in both sections. Kdm5c (28% / 25%), Utp14a (27% / 30%) and Akap17a (26% / 20%) escape at canonical magnitudes in both sections. "
    "Ndufb11 and Smpx, the two deepest chrX genes, are the only ones with a spatial signal in the spASE test (Ndufb11 7% to 3%, Smpx 7% to 17% between sections) and need replication. "
    "Kdm6a, Ddx3x, Ftx and Jpx are present but below the depth filter. At 16 µm 90% of occupied pixels hold one molecule, so the CIs are binomial and overdispersion is not estimable.")
picture_or_placeholder(s, "F18_spatial_per_gene_escape.png", LEFT, TOP, 7.2, H)
add_text(s, [
    ("Kdm5c 28% / 25%, Utp14a 27% / 30%, Akap17a 26% / 20% (9w / 78w): canonical escape, reproduced across sections", 0),
    ("Ndufb11 7% → 3% and Smpx 7% → 17%: the two deepest genes, the only ones with a spatial signal (spASE); need replication", 0),
    ("Kdm6a, Ddx3x, Ftx, Jpx detected but below the 100-UMI filter", 0),
    ("Gene-level escape is the honest estimate; the chromosome-wide number is not", 0),
    ("16 µm pixels hold ~1 molecule, so CIs are binomial; overdispersion needs larger pixels", 0),
], 7.9, TOP, 4.95, H, size=13)

# ===========================================================================
# WP2 preview - multiome
# ===========================================================================
section("Aim 3: Test whether the inactive X gains distal chromosome accessibility with age, and in which cell types",
        "10x Epi Multiome · one adult (9w) and one aged (78w) heart · cellranger-arc 2.2.0 on GRCm39",
        "This is WP2 of the grant, run on the same biobanked animals. RNA gives escape per nucleus, ATAC gives accessibility of the inactive X, which RNA alone cannot.")

s = new_slide("Multiome: sequencing and QC",
    "GEX was sequenced five to nine times deeper than the other libraries on the flowcell, deliberately, for allelic power: 687 M and 733 M read pairs at roughly 8,000 nuclei. That gives 98% duplicates and 10-16 M unique molecules per sample. "
    "ATAC is the weaker modality (TSS enrichment 6-7, FRiP 14-20%). The 78w sample is better on every axis (2.2x nuclei, better ATAC), so any age difference is confounded with data quality on top of n = 1.")
rows = [["", "Adult (9w)", "Aged (78w)"],
        ["Nuclei called / passing joint QC", "1,758 / 1,422", "3,822 / 3,125"],
        ["GEX read pairs", "687 M", "733 M"],
        ["GEX duplicates", "98.6%", "97.8%"],
        ["Median UMIs / genes per nucleus", "1,354 / 847", "863 / 574"],
        ["ATAC read pairs", "412 M", "444 M"],
        ["Median high-quality fragments per nucleus", "9,995", "7,821"],
        ["ATAC TSS enrichment", "5.96", "7.42"],
        ["Fragments in peaks (FRiP)", "0.14", "0.20"]]
add_table(s, rows, LEFT, TOP, 6.6, col_widths=[3.6, 1.5, 1.5], size=12)
picture_or_placeholder(s, "F19_multiome_nuclei_per_celltype.png", 7.4, TOP, H, 3.1)
add_text(s, [
    ("Deliberately over-sequenced (5-9x the other libraries on the flowcell) for allelic power", 0),
    ("ATAC is the weaker modality: modest TSS enrichment and FRiP", 0),
    ("78w is the better sample on every axis, so age is confounded with quality on top of n = 1", 0),
], 7.4, 4.9, H, 2.2, size=13)

s = new_slide("Joint RNA + ATAC cell typing",
    "Weighted nearest-neighbour clustering on RNA and ATAC together, cell types labelled with the same marker panels as the snRNA-seq object so the two datasets are comparable. "
    "Seven cell types, all of which appear in both samples. The per-cell-type barcode lists are what the allelic run is now using.")
picture_or_placeholder(s, "multiome_UMAP_WNN_labelled.png", LEFT, TOP, 6.4, H)
picture_or_placeholder(s, "multiome_UMAP_RNA_ATAC_WNN.png", 7.0, TOP, 5.85, 2.2)
add_text(s, [
    ("WNN on RNA + ATAC; labels from the same marker panels as the snRNA-seq object", 0),
    ("Seven cell types, present in both samples; T and B cells too sparse to resolve", 0),
    ("RNA and ATAC embeddings agree on the major populations", 0),
    ("Per-cell-type barcode lists feed the allelic run", 0),
], 7.0, 4.3, 5.85, 2.9, size=13)

s = new_slide("ATAC signal at Xist and escape genes",
    "Pseudobulk ATAC coverage per cell type. Xist has accessible chromatin across the locus in every cell type (the active B6 allele's Xist is deleted, so this is the CAST Xi allele). "
    "Kdm5c and Kdm6a show a sharp promoter peak in every cell type, as expected for escapees. The allelic split of these peaks, i.e. whether the peak sits on the inactive X, is the next step.")
for i, (name, lab) in enumerate([("multiome_coverage_p1.png", "Xist"), ("multiome_coverage_p2.png", "Kdm5c"), ("multiome_coverage_p3.png", "Kdm6a")]):
    picture_or_placeholder(s, name, LEFT + i * 4.15, TOP, 4.0, 2.6)
picture_or_placeholder(s, "multiome_escape_gene_activity_p2.png", LEFT, 4.6, 6.0, 2.55)
add_text(s, [
    ("Pseudobulk ATAC per cell type; Xist locus accessible in every cell type", 0),
    ("Kdm5c, Kdm6a: sharp promoter peaks in every cell type, as expected for escapees", 0),
    ("Gene-activity scores for the core escapees are highest in fibroblasts, pericytes and macrophages", 0),
    ("Next: allelic split of these peaks, i.e. is the accessible copy on the inactive X?", 0),
], 6.8, 4.6, 6.05, 2.55, size=13)

s = new_slide("Allele-specific multiome: running now",
    "Status as of today. The pseudobulk Allelome.PRO2 run per chromosome has completed for both samples; the gate is whether chrX reproduces the snRNA-seq / spatial escape level and whether autosomes sit at 0.5 under identical filters. "
    "The per-cell-type run was launched this evening after the sinto barcode split. Allelic ATAC needs a different MAPQ filter (BWA caps at 60) and peaks as the annotation.")
rows = [["Step", "Status"],
        ["Pseudobulk Allelome.PRO2 per chromosome (GEX, -q 255, duplicates removed)", "Done for both samples; gate: chrX escape level and autosomal baseline"],
        ["Per-cell-type Allelome.PRO2 (GEX)", "Running (barcode split finished today)"],
        ["Allelic ATAC over peaks (-q 30, BWA MAPQ)", "Next"],
        ["Peak-to-gene linkage at escape vs non-escape chrX genes", "Next"],
        ["Distal-chromosome enrichment of inactive-X accessibility per cell type (grant question 1)", "Needs the two steps above"]]
add_table(s, rows, LEFT, TOP, 12.3, col_widths=[7.3, 5.0], size=14)
add_text(s, [
    ("What multiome adds", -1),
    ("RNA: CAST fraction on chrX = escape per nucleus, ~20-25 informative chrX molecules per nucleus, so pseudobulk per cell type is the working level", 0),
    ("ATAC: CAST fraction over chrX peaks = accessibility of the inactive X, the axis RNA cannot give", 0),
    ("Two corrections carried over: Xic mask (deletion reads 100% CAST) and the autosomal mapping-bias baseline, re-estimated for ATAC read length", 0),
], LEFT, 4.5, 12.3, 2.6, size=15)

# ===========================================================================
# Summary
# ===========================================================================
s = new_slide("Summary",
    "Wrap up against the two grant questions. Question 1 (chromatin accessibility on the inactive X, per cell type) is the multiome, now running. Question 2 (dispersed vs clonal vs niche) has a first answer from one section per age: dispersed, no spatial structure, as the skewed-XCI model predicts. "
    "Everything is n = 1 per age; the grant's 3 vs 3 design with deep sequencing is what turns these observations into tests.")
add_text(s, [
    ("Aim 1: snRNA-seq", -1),
    ("Escape from the inactive X is measurable per nucleus and per gene; core escapees biallelic in every cell type", 0),
    ("Biallelic nuclei more frequent in aged endothelial, pericyte and fibroblast populations; endocardium unchanged", 0),
    ("Apparent escape depends on allelic depth: deep sequencing is the fix", 0),
    ("Xist loss coincides with inactive-X expression in cardiomyocytes", 0),
    ("Aim 2: spatial", -1),
    ("No spatial structure of escape from 4 µm to 2 mm: dispersed, not clonal, as non-random XCI predicts", 0),
    ("Per-gene escape reproducible across sections (Kdm5c, Utp14a, Akap17a); chromosome-wide 12.7% is mostly SNP artefact, gene-level ~3%", 0),
], LEFT, TOP, 6.6, H, size=14)
add_text(s, [
    ("Caveats", -1),
    ("One animal per condition throughout: every age statement is descriptive", 0),
    ("Sections and samples are not depth matched; 78w is the better multiome sample", 0),
    ("Next steps (grant work packages)", -1),
    ("WP1: NovaSeq X 1.5B deep sequencing to break the allelic-depth ceiling; corrected SNP mask (B6J vs B6NJ)", 0),
    ("WP2: 3 vs 3 multiome, allelic ATAC → is inactive-X accessibility enriched at distal chrX, and in which cell types?", 0),
    ("WP3: 3 vs 3 Visium HD of the superior sections, beta-binomial spatial GLMM, spASE benchmark", 0),
    ("Same hearts for WP2 and WP3: inferior half to multiome, superior half to spatial", 0),
], 7.3, TOP, 5.55, H, size=14)

prs.save(OUT)
print("saved", OUT, "slides:", len(prs.slides))
