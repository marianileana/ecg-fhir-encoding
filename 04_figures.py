"""Step 4: Figures 1, 3 and 4 of the article (PNG, 600 dpi) and the Fig. 4 regression statistics.

Fig. 1  system overview (gateway and FHIR server)
Fig. 3  payload versus captured ectopic beats in the hybrid mode for several thresholds
Fig. 4  S5 payload per record versus the share of 10-s windows with non-N beats
(Fig. 2 of the article is a text listing of an S5 message.)

Uses the Times-like font TeX Gyre Termes if it is installed, otherwise the default serif font.
Input: results/encode_send.json and results/extra.json. Output: fig/*.png, results/fig4_stats.json.
"""
import json
import os

import matplotlib
import numpy as np
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "fig")
os.makedirs(FIG, exist_ok=True)
TERMES = ("texgyretermes-regular.otf", "texgyretermes-bold.otf", "texgyretermes-italic.otf")
FONT_DIRS = (os.environ.get("TEXGYRE_DIR", ""), "/usr/share/texmf/fonts/opentype/public/tex-gyre",
             "/usr/share/texlive/texmf-dist/fonts/opentype/public/tex-gyre")
FONT = None
for d in FONT_DIRS:
    if d and all(os.path.exists(os.path.join(d, f)) for f in TERMES):
        for f in TERMES:
            font_manager.fontManager.addfont(os.path.join(d, f))
        FONT = "TeX Gyre Termes"
        break
if FONT:
    plt.rcParams.update({"font.family": FONT, "mathtext.fontset": "custom", "mathtext.rm": FONT,
                         "mathtext.it": FONT + ":italic", "mathtext.bf": FONT + ":bold", "mathtext.cal": FONT})
else:
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix"})
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
BLUE, ORANGE = "#2a78d6", "#eb6834"

# ------------------------------------------------------------------ Fig. 1 (system overview)
fig, ax = plt.subplots(figsize=(6.85, 2.05), dpi=300)
ax.set_xlim(0, 150)
ax.set_ylim(0, 41)
ax.axis("off")


def box(x, y, w, h, title, lines, fc="#ffffff", ls="-"):
    """Rounded box with a bold title and centered text lines."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                lw=0.8, ec=INK, fc=fc, ls=ls, zorder=2))
    ax.text(x + w / 2, y + h - 2.6, title, ha="center", va="top", fontsize=8, weight="bold", color=INK, zorder=3)
    ax.text(x + w / 2, y + h - 7.6, "\n".join(lines), ha="center", va="top", fontsize=7,
            color=INK, linespacing=1.35, zorder=3)


def arrow(x0, y0, x1, y1, label=None, dy=1.6):
    """Arrow with an optional label above it."""
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=8, lw=0.8,
                                 color=INK, zorder=1))
    if label:
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + dy, label, ha="center", va="bottom", fontsize=6.8, color=INK2)


box(1, 9, 22, 26, "ECG source", ["single lead, 360 Hz,", "11-bit samples", "(48 MIT-BIH", "records)"])
box(34, 3, 44, 35, "Edge gateway", ["streaming QRS detector", "(Pan-Tompkins, 60-s blocks)",
                                     "HR per 10-s window",
                                     "premature-beat trigger", r"(RR$_{\mathit{i}}$ < θ · median of previous 8 RR)",
                                     "FHIR encoder S1-S5", "+ optional gzip"], fc="#f4f4f2")
box(95, 3, 30, 35, "FHIR server", ["HTTP/1.1 keep-alive", "gzip decoding", "JSON parsing",
                                      "structural validation", "(Observation, Bundle)"])
box(131, 9, 18, 26, "Clinician", ["HR trend", "flagged ECG", "segments"], ls="--")
arrow(24.0, 22, 33.0, 22, "samples")
arrow(79.0, 26, 94.0, 26, "HTTP POST")
ax.text(86.5, 18.0, "Observation\nor transaction\nBundle", ha="center", va="center", fontsize=6.8, color=INK2)
arrow(126.0, 22, 130.0, 22)  # context only (dashed box)
fig.savefig(f"{FIG}/Fig1_architecture.png", bbox_inches="tight", pad_inches=0.02, dpi=600)
plt.close(fig)

# ------------------------------------------------------------------ Fig. 3 (hybrid mode trade-off)
d = json.load(open(f"{OUT}/encode_send.json"))
sw = d["sweep"]
th = sorted(sw)
x = [sw[t]["body_kB_per_hour"] for t in th]
veb = [100 * sw[t]["veb_se"] for t in th]
sveb = [100 * sw[t]["sveb_se"] for t in th]
s3 = d["strategies"]["S3"]["body_kB_per_hour"]
s4 = d["strategies"]["S4"]["body_kB_per_hour"]

fig, ax = plt.subplots(figsize=(3.3, 2.35), dpi=300)
ax.grid(True, color=GRID, lw=0.4, zorder=0)
ax.plot(x, veb, "-o", color=BLUE, lw=1.4, ms=4, mec="white", mew=0.6, zorder=3, label="VEB")
ax.plot(x, sveb, "-s", color=ORANGE, lw=1.4, ms=4, mec="white", mew=0.6, zorder=3, label="SVEB")
POS = {"0.70": ((0, -5), "center", "top"), "0.95": ((-7, 3), "right", "bottom")}
for t, xi, yi in zip(th, x, sveb):
    off, ha, va = POS.get(f"{float(t):.2f}", ((5, -2), "left", "top"))
    ax.annotate(f"θ = {float(t):.2f}", (xi, yi), xytext=off,
                textcoords="offset points", ha=ha, va=va, fontsize=6.5, color=INK2)
ax.axvline(s3, color=INK2, lw=0.8, ls="--", zorder=2)
ax.text(s3 - 20, 3, "S3: full raw\nstream (gzip)", ha="right", va="bottom", fontsize=6.5, color=INK2)
ax.axvline(s4, color=INK2, lw=0.8, ls=":", zorder=2)
ax.text(s4 + 20, 3, "S4: HR only", ha="left", va="bottom", fontsize=6.5, color=INK2)
ax.set_xlim(0, 1600)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_ylim(0, 112)
ax.set_xlabel("Payload (kB per hour of ECG)", color=INK)
ax.set_ylabel("Ectopic beats in sent segments (%)", color=INK)
ax.tick_params(labelsize=7, colors=INK, width=0.5, length=2.5)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(loc="upper left", bbox_to_anchor=(0.07, 0.99), fontsize=6.6, frameon=False, handlelength=2)
fig.savefig(f"{FIG}/Fig3_tradeoff.png", bbox_inches="tight", pad_inches=0.02, dpi=600)
plt.close(fig)
print("Fig. 1 and Fig. 3 written")

# ------------------------------------------------------------------ Fig. 4 (S5 payload per record)

x3 = json.load(open(f"{OUT}/extra.json"))["per_record"]
recs = [r for r, v in x3.items() if not v["paced"]]
a = np.array([x3[r]["ectopic_window_pct"] for r in recs])
b = np.array([x3[r]["s5_kB_per_hour"] for r in recs])
pr = stats.pearsonr(a, b)
fit = np.polyfit(a, b, 1)
json.dump({"n": len(recs), "pearson_r": float(pr.statistic), "p": float(pr.pvalue),
           "slope_kB_per_h_per_pct": float(fit[0]), "intercept_kB_per_h": float(fit[1])},
          open(f"{OUT}/fig4_stats.json", "w"), indent=1)
fig, ax = plt.subplots(figsize=(3.3, 2.2), dpi=300)
ax.grid(True, color=GRID, lw=0.4, zorder=0)
ax.scatter(a, b, s=14, color=BLUE, edgecolor="white", linewidth=0.5, zorder=3)
xx = np.array([0, 100])
ax.plot(xx, np.polyval(fit, xx), color=INK2, lw=0.9, ls="--", zorder=2)
ax.text(1, 1620, r"$\mathit{r}$" + f" = {pr.statistic:.2f} (44 non-paced records)", fontsize=6.5, color=INK2)
LABEL = {"122": ((14, -2), "left"), "203": ((-4, 4), "right"), "233": ((4, -3), "left")}
for r, (off, ha) in LABEL.items():
    i = recs.index(r)
    ax.annotate(r, (a[i], b[i]), xytext=off, textcoords="offset points", ha=ha, va="center" if r == "122" else "baseline",
                fontsize=6.5, color=INK2,
                arrowprops=dict(arrowstyle="-", lw=0.4, color=INK2, shrinkA=1, shrinkB=3) if r == "122" else None)
ax.set_xlim(-3, 103)
ax.set_ylim(0, 1900)
ax.set_xlabel("10-s windows containing non-N beats (%)", color=INK)
ax.set_ylabel("S5 payload (kB/h)", color=INK)
ax.tick_params(labelsize=7, colors=INK, width=0.5, length=2.5)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.savefig(f"{FIG}/Fig4_per_record.png", bbox_inches="tight", pad_inches=0.02, dpi=600)
plt.close(fig)
print("fig4", pr, fit)
