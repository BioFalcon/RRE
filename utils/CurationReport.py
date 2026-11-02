#!/usr/bin/env python3
#Report of the alignment cleaning (MSACuration.py) in the AncientMode extension rounds of one repeat:
#for every attempt (accepted rounds and failed attempts), the distribution of the column scores of the raw
#alignment, the peak and valley found in their density, the cap and the threshold used.
#Reads the <...>.Curation.json files in WorkDir/<Repeat>/{Left,Right}/Family_*/{Round_*,Failed_Round*}/
#Then the final model, from the alignments of its copies (consensus columns only; each copy counts only within its
#own span, as most are fragments): depth, agreement with the consensus base, deletions and call margin along the
#model, where copies start and end; and the same metrics, with copy divergence and length, overlaid for the merged
#alignment (WorkDir/<Repeat>/Merge) and every refinement iteration (WorkDir/<Repeat>/FinalModel/Iterations).
#Writes <out>.CurationReport.pdf, <out>.CurationRounds.tsv, <out>.RefinementColumns.tsv and <out>.FinalModelColumns.tsv

import argparse
import filecmp
import glob
import json
import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import gaussian_kde

#Colors: a neutral for the data, one hue for the density and valley, one for the threshold actually used
HistColor = "#c9ccd1"
DensityColor = "#2a6fdb"
ThreshColor = "#e8590c"
CapColor = "#6b7280"
FailColor = "#9aa0a6"
Ink = "#1f2328"
Muted = "#57606a"
#Columns whose two most common bases differ by less than this (fraction of the copies with a base): unstable calls
LowMargin = 0.1

def ReadAttempts(WorkDir):
    Attempts = []
    for Side in ("Left", "Right"):
        for Json in sorted(glob.glob(os.path.join(WorkDir, Side, "Family_*", "*", "*.Curation.json"))):
            RoundDir = os.path.dirname(Json)
            FamDir = os.path.dirname(RoundDir)
            Name = os.path.basename(RoundDir)
            Accepted = re.match(r"^Round_(\d+)$", Name)
            Failed = re.match(r"^Failed_Round(\d+)_Retry(\d+)", Name)
            if Accepted:
                Round = int(Accepted.group(1))
                #The retry level of an accepted round is set by the previous round
                Prev = os.path.join(FamDir, "Round_%02d" % (Round - 1), "NextRetryLevel")
                try:
                    Retry = int(open(Prev).read().strip() or 0)
                except (OSError, ValueError):
                    Retry = 0
                Status = "accepted"
            elif Failed:
                Round, Retry, Status = int(Failed.group(1)), int(Failed.group(2)), "failed"
            else:
                continue
            try:
                Stats = json.load(open(Json))
            except (OSError, ValueError):
                continue
            Attempts.append({"side": Side, "family": os.path.basename(FamDir), "round": Round, "retry": Retry,
                             "status": Status, "dir": Name, "stats": Stats})
    #Failed attempts of a round before the accepted one
    Attempts.sort(key=lambda a: (a["side"], a["family"], a["round"], a["status"] == "accepted", a["retry"], a["dir"]))
    return Attempts

def Cut(Stats):
    return Stats["cuts"][0] if Stats.get("cuts") else None

def WriteTSV(Attempts, Out):
    Cols = ["side", "family", "round", "retry", "status", "dir", "sequences", "columns", "mode", "peaks", "selected_peak",
            "valley", "valley_found", "cap", "threshold", "cap_reached", "columns_after_windows", "columns_kept"]
    with open(Out, "w") as F:
        F.write("\t".join(Cols) + "\n")
        for a in Attempts:
            S, C = a["stats"], Cut(a["stats"])
            Row = [a["side"], a["family"], a["round"], a["retry"], a["status"], a["dir"], S["sequences"], S["columns"],
                   S["mode"], ",".join(str(p) for p in S.get("peaks", [])) or "NA",
                   C["peak"] if C else "NA", C["valley"] if C else "NA", C["valley_found"] if C else "NA",
                   S["max_thresh"], C["threshold"] if C else "NA", C["cap_reached"] if C else "NA",
                   C["columns_after_windows"] if C else "NA", C["columns_kept"] if C else "NA"]
            F.write("\t".join(str(x) for x in Row) + "\n")

def StyleAxis(Ax):
    for Sp in ("top", "right"):
        Ax.spines[Sp].set_visible(False)
    for Sp in ("left", "bottom"):
        Ax.spines[Sp].set_color("#d0d7de")
    Ax.tick_params(colors=Muted, labelsize=7, length=2)
    Ax.grid(axis="y", color="#eaeef2", linewidth=0.6)
    Ax.set_axisbelow(True)

def Label(a):
    L = "%s round %02d" % (a["side"], a["round"])
    if a["retry"] > 0:
        L += ", retry %d" % a["retry"]
    return L

def OverviewPage(Pdf, Attempts, Repeat):
    Sides = [s for s in ("Left", "Right") if any(a["side"] == s for a in Attempts)]
    Fig, Axes = plt.subplots(len(Sides), 1, figsize=(11, 3.4 * len(Sides) + 1.2), squeeze=False)
    Fig.suptitle("%s: threshold of the alignment cleaning in each extension round" % Repeat, x=0.06, ha="left",
                 fontsize=13, color=Ink)
    for Ax, Side in zip(Axes[:, 0], Sides):
        StyleAxis(Ax)
        Sub = [a for a in Attempts if a["side"] == Side]
        Cap = Sub[0]["stats"]["max_thresh"]
        Ax.axhline(Cap, color=CapColor, linestyle="--", linewidth=1)
        for Status in ("accepted", "failed"):
            Pts = [(a["round"], Cut(a["stats"])) for a in Sub if a["status"] == Status and Cut(a["stats"])]
            if not Pts:
                continue
            X = [p[0] for p in Pts]
            Valley = [p[1]["valley"] for p in Pts]
            Thresh = [p[1]["threshold"] for p in Pts]
            if Status == "accepted":
                Ax.scatter(X, Valley, s=26, facecolors="none", edgecolors=DensityColor, linewidths=1, zorder=3)
                Ax.scatter(X, Thresh, s=14, color=ThreshColor, zorder=4)
            else:
                #Failed attempts: the threshold they used
                Ax.scatter(X, Thresh, s=24, marker="x", color=FailColor, linewidths=1, zorder=4)
        NoValley = [a for a in Sub if a["stats"]["mode"] != "valley"]
        for a in NoValley:
            Ax.annotate("no cut", (a["round"], 0), fontsize=6, color=Muted, ha="center", va="bottom")
        Acc = [a for a in Sub if a["status"] == "accepted" and Cut(a["stats"])]
        NCap = sum(Cut(a["stats"])["cap_reached"] for a in Acc)
        NNoValley = sum(not Cut(a["stats"])["valley_found"] for a in Acc)
        NFailed = sum(a["status"] == "failed" for a in Sub)
        Ax.set_title("%s side: %d accepted rounds, %d failed attempt%s; cap reached in %d of the accepted rounds%s" %
                     (Side, len(Acc), NFailed, "" if NFailed == 1 else "s", NCap,
                      ", no valley before 0.5 in %d" % NNoValley if NNoValley else ""),
                     loc="left", fontsize=9, color=Muted)
        Ax.set_xlabel("Round", fontsize=8, color=Muted)
        Ax.set_ylabel("Column score", fontsize=8, color=Muted)
        Top = max([Cap] + [Cut(a["stats"])["valley"] for a in Sub if Cut(a["stats"])])
        Ax.set_ylim(min(0, min([Cut(a["stats"])["valley"] for a in Sub if Cut(a["stats"])] + [0])) - 0.02, Top * 1.15 + 0.02)
        if Sub:
            Ax.set_xlim(min(a["round"] for a in Sub) - 0.7, max(a["round"] for a in Sub) + 0.7)
    Handles = [Line2D([], [], marker="o", linestyle="none", markerfacecolor="none", markeredgecolor=DensityColor, label="Valley after the first peak"),
               Line2D([], [], marker="o", linestyle="none", color=ThreshColor, markersize=4, label="Threshold used (valley, capped)"),
               Line2D([], [], color=CapColor, linestyle="--", label="Cap (--curationMaxThresh)"),
               Line2D([], [], marker="x", linestyle="none", color=FailColor, label="Failed attempt")]
    Fig.legend(handles=Handles, loc="upper left", bbox_to_anchor=(0.06, 0.95), ncol=4, frameon=False, fontsize=8)
    Fig.text(0.06, 0.015, "Columns of the raw alignment are scored (1 - entropy)^a * (1 - gap fraction)^g. The threshold is the valley "
             "after the first peak of their density, capped; windows whose median score is not above it are removed.",
             fontsize=7, color=Muted)
    Fig.tight_layout(rect=(0.03, 0.04, 1, 0.9))
    Pdf.savefig(Fig)
    plt.close(Fig)

def AttemptPanel(Ax, a):
    StyleAxis(Ax)
    S, C = a["stats"], Cut(a["stats"])
    Scores = S["scores"]
    Ax.hist(Scores, bins=[i / 50 for i in range(-5, 51)], density=True, color=HistColor, edgecolor="white", linewidth=0.4)
    if "grid" in S:
        G = [(x, y) for x, y in zip(S["grid"], S["density"]) if -0.1 <= x <= 1.0]
        Ax.plot([g[0] for g in G], [g[1] for g in G], color=DensityColor, linewidth=1.2)
        for p in S.get("peaks", []):
            Y = S["density"][S["grid"].index(p)] if p in S["grid"] else None
            if Y is not None:
                Ax.plot(p, Y, marker="v", markersize=4, color=DensityColor,
                        markerfacecolor=DensityColor if C and p == C["peak"] else "white")
    Ax.axvline(S["max_thresh"], color=CapColor, linestyle="--", linewidth=0.9)
    if C:
        Ax.axvline(C["valley"], color=DensityColor, linestyle=":", linewidth=1.1)
        Ax.axvline(C["threshold"], color=ThreshColor, linewidth=1.3)
    #Scores rarely reach 1: show up to a bit past the highest scores
    High = sorted(Scores)[int(0.99 * (len(Scores) - 1))] if Scores else 1.0
    Ax.set_xlim(-0.05, min(1.0, max(0.5, High + 0.1)))
    Title = Label(a) + ("  (failed)" if a["status"] == "failed" else "")
    Ax.set_title(Title, loc="left", fontsize=8.5, color=FailColor if a["status"] == "failed" else Ink)
    if C:
        Note = "%d seqs, %d cols\nkept %d by score, %d by gaps\nvalley %.2f%s, used %.2f" % (
            S["sequences"], S["columns"], C["columns_after_windows"], C["columns_kept"], C["valley"],
            "" if C["valley_found"] else " (none before 0.5)", C["threshold"])
        if C["cap_reached"]:
            Note += "  CAP"
    else:
        Note = "%d seqs, %d cols\nno peak: gap filter only" % (S["sequences"], S["columns"])
    Ax.text(0.98, 0.95, Note, transform=Ax.transAxes, ha="right", va="top", fontsize=6.5,
            color=ThreshColor if C and C["cap_reached"] else Muted)

def AttemptPages(Pdf, Attempts, Repeat, Rows=4, Cols=3):
    PerPage = Rows * Cols
    for Side in ("Left", "Right"):
        Sub = [a for a in Attempts if a["side"] == Side]
        for Page in range(math.ceil(len(Sub) / PerPage)):
            Chunk = Sub[Page * PerPage:(Page + 1) * PerPage]
            Fig, Axes = plt.subplots(Rows, Cols, figsize=(11, 8.5))
            Fig.suptitle("%s, %s side: column scores of the raw alignment in each round (%d/%d)" %
                         (Repeat, Side, Page + 1, math.ceil(len(Sub) / PerPage)), x=0.06, ha="left", fontsize=12, color=Ink)
            for i, Ax in enumerate(Axes.flat):
                if i < len(Chunk):
                    AttemptPanel(Ax, Chunk[i])
                    if i % Cols == 0:
                        Ax.set_ylabel("Density", fontsize=7, color=Muted)
                    if i >= (Rows - 1) * Cols or i + Cols >= len(Chunk):
                        Ax.set_xlabel("Column score", fontsize=7, color=Muted)
                else:
                    Ax.axis("off")
            Handles = [plt.Rectangle((0, 0), 1, 1, color=HistColor, label="Column scores"),
                       Line2D([], [], color=DensityColor, label="Density (KDE)"),
                       Line2D([], [], marker="v", linestyle="none", color=DensityColor, label="Peaks (filled: used)"),
                       Line2D([], [], color=DensityColor, linestyle=":", label="Valley"),
                       Line2D([], [], color=CapColor, linestyle="--", label="Cap"),
                       Line2D([], [], color=ThreshColor, label="Threshold used")]
            Fig.legend(handles=Handles, loc="upper left", bbox_to_anchor=(0.06, 0.955), ncol=6, frameon=False, fontsize=7.5)
            Fig.tight_layout(rect=(0.02, 0, 1, 0.92))
            Pdf.savefig(Fig)
            plt.close(Fig)

def ReadFasta(Path):
    Seqs, Cur = [], None
    for Line in open(Path):
        Line = Line.strip()
        if Line.startswith(">"):
            if Cur is not None:
                Seqs.append("".join(Cur))
            Cur = []
        elif Cur is not None:
            Cur.append(Line)
    if Cur is not None:
        Seqs.append("".join(Cur))
    return Seqs

def MatchColumns(Seqs, Keep):
    #Sequences restricted to the columns Keep (indexes), without the ones left empty
    Out = []
    for x in Seqs:
        y = np.frombuffer(x.encode(), dtype=np.uint8)[Keep].tobytes().decode()
        if y.strip("-."):
            Out.append(y)
    return Out

def ReadMerge(WorkDir):
    #Merged alignment of the extension, consensus columns only (Merge.Columns.txt: the consensus base, or - for an
    #insert column), as the refinement takes it for its iteration 0
    Aln = os.path.join(WorkDir, "Merge", "Merge.aln.fa")
    Cols = os.path.join(WorkDir, "Merge", "Merge.Columns.txt")
    if not (os.path.isfile(Aln) and os.path.isfile(Cols)):
        return None
    Seqs = ReadFasta(Aln)
    Calls = open(Cols).read().replace("\n", "").strip()
    if not Seqs or len(Calls) != len(Seqs[0]):
        return None
    return MatchColumns(Seqs, np.array([i for i, c in enumerate(Calls) if c != "-"]))

def CopyMetrics(Seqs):
    #Per-column and per-copy metrics of an alignment of copies (consensus columns only). Most copies of an old repeat
    #are fragments, so each copy only counts within its own span (its first to last base): gaps outside the span are
    #missing data, gaps inside it are deletions.
    A = np.array([np.frombuffer(x.encode(), dtype=np.uint8) for x in Seqs])
    A = np.where((A >= 97) & (A <= 122), A - 32, A)
    Base = (A != ord("-")) & (A != ord("."))
    Keep = Base.any(axis=1)
    A, Base = A[Keep], Base[Keep]
    L = A.shape[1]
    First = Base.argmax(axis=1)
    Last = L - 1 - Base[:, ::-1].argmax(axis=1)
    Pos = np.arange(L)
    Span = (Pos[None, :] >= First[:, None]) & (Pos[None, :] <= Last[:, None])
    Depth = Span.sum(axis=0)
    Counts = np.stack([(A == ord(c)).sum(axis=0) for c in "ACGT"])
    Called = Counts.sum(axis=0)
    Sorted = np.sort(Counts, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        Agreement = np.where(Called > 0, Sorted[-1] / Called, np.nan)
        Margin = np.where(Called > 0, (Sorted[-1] - Sorted[-2]) / Called, np.nan)
        Deletions = np.where(Depth > 0, (Span & ~Base).sum(axis=0) / Depth, np.nan)
    #Majority base of each column, and each copy's mismatches to it over the bases it shares with the consensus
    Cons = np.frombuffer(b"ACGT", dtype=np.uint8)[Counts.argmax(axis=0)]
    ACGT = np.isin(A, np.frombuffer(b"ACGT", dtype=np.uint8)) & (Called > 0)[None, :]
    Compared = ACGT.sum(axis=1)
    Mismatch = ((A != Cons[None, :]) & ACGT).sum(axis=1)
    Long = Compared >= 50
    return {"copies": int(A.shape[0]), "columns": int(L), "depth": Depth, "agreement": Agreement, "deletions": Deletions,
            "margin": Margin, "consensus": Cons.tobytes().decode(), "starts": First, "ends": Last,
            "divergence": Mismatch[Long] / Compared[Long], "span": Last - First + 1}

def ReadModels(WorkDir):
    #The merged alignment (iteration 0) and the alignment of every refinement iteration, consensus columns only
    Iters = os.path.join(WorkDir, "FinalModel", "Iterations")
    Files = []
    if os.path.isfile(os.path.join(WorkDir, "Merge", "Merge.aln.fa")):
        Files.append((0, "merge"))
    for F in sorted(glob.glob(os.path.join(Iters, "Iter*.match.aln.fa"))):
        Files.append((int(re.search(r"Iter(\d+)\.match", F).group(1)), F))
    #The iteration kept as the final model (later ones can be discarded)
    Final = None
    for RF in glob.glob(os.path.join(WorkDir, "FinalModel", "*.RFclean.aln.fa")):
        for It, F in Files:
            if F.endswith(".match.aln.fa") and filecmp.cmp(RF, F, shallow=False):
                Final = It
    Out = []
    for It, F in Files:
        Seqs = ReadMerge(WorkDir) if F == "merge" else ReadFasta(F)
        if not Seqs:
            continue
        M = CopyMetrics(Seqs)
        M.update({"iteration": It, "final": It == Final})
        Out.append(M)
    #Not refined: the merged model is the final model
    if Out and Final is None and len(Out) == 1:
        Out[0]["final"] = True
    return Out

def IterName(d):
    if d["iteration"] == 0:
        return "merged alignment" + (" (final model)" if d["final"] else "")
    return "iteration %d%s" % (d["iteration"], " (final model)" if d["final"] else "")

def Summary(d):
    return {"copies": d["copies"], "columns": d["columns"], "median_depth": float(np.median(d["depth"])),
            "median_agreement": float(np.nanmedian(d["agreement"])), "median_deletions": float(np.nanmedian(d["deletions"])),
            "low_margin_columns": int(np.nansum(d["margin"] < LowMargin)),
            "median_copy_divergence": float(np.median(d["divergence"])) if len(d["divergence"]) else float("nan"),
            "median_copy_span": float(np.median(d["span"]))}

def WriteModelTSVs(Models, Out):
    with open(Out + ".RefinementColumns.tsv", "w") as F:
        Keys = list(Summary(Models[0]).keys())
        F.write("iteration\tfinal\t" + "\t".join(Keys) + "\n")
        for d in Models:
            S = Summary(d)
            F.write("%d\t%s\t" % (d["iteration"], d["final"]) + "\t".join(("%.4f" % S[k]) if isinstance(S[k], float) else str(S[k]) for k in Keys) + "\n")
    Fin = [d for d in Models if d["final"]]
    if Fin:
        d = Fin[0]
        Starts = np.bincount(d["starts"], minlength=d["columns"])
        Ends = np.bincount(d["ends"], minlength=d["columns"])
        with open(Out + ".FinalModelColumns.tsv", "w") as F:
            F.write("position\tconsensus\tdepth\tagreement\tdeletions\tmargin\tcopies_starting\tcopies_ending\n")
            for i in range(d["columns"]):
                F.write("%d\t%s\t%d\t%.4f\t%.4f\t%.4f\t%d\t%d\n" % (i + 1, d["consensus"][i], d["depth"][i], d["agreement"][i],
                        d["deletions"][i], d["margin"][i], Starts[i], Ends[i]))

def Rolling(V, W=25):
    #Mean over a centred window of W columns, ignoring missing values
    V = np.asarray(V, dtype=float)
    Ok = ~np.isnan(V)
    K = np.ones(W)
    S = np.convolve(np.where(Ok, V, 0), K, mode="same")
    N = np.convolve(Ok.astype(float), K, mode="same")
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(N > 0, S / N, np.nan)

def FinalModelPage(Pdf, d, Repeat):
    Fig, Axes = plt.subplots(5, 1, figsize=(11, 8.5), sharex=True, gridspec_kw={"height_ratios": (1.1, 1, 1, 1, 1.1)})
    Fig.suptitle("%s: final model along its length (%d consensus columns, %d copies)" % (Repeat, d["columns"], d["copies"]),
                 x=0.06, ha="left", fontsize=13, color=Ink)
    X = np.arange(1, d["columns"] + 1)
    Panels = [("depth", "Copies spanning\nthe column", None),
              ("agreement", "Copies with the\nconsensus base", (0, 1)),
              ("deletions", "Copies with a\ndeletion", (0, None)),
              ("margin", "Call margin\n(top - second base)", (0, 1))]
    for Ax, (Key, Name, Lim) in zip(Axes, Panels):
        StyleAxis(Ax)
        V = d[Key]
        if Key == "depth":
            Ax.fill_between(X, V, step="mid", color=DensityColor, alpha=0.25, linewidth=0)
            Ax.plot(X, V, color=DensityColor, linewidth=0.8, drawstyle="steps-mid")
        else:
            Ax.plot(X, V, color=HistColor, linewidth=0.4)
            Ax.plot(X, Rolling(V), color=DensityColor, linewidth=1.3)
        if Key == "margin":
            Ax.axhline(LowMargin, color=ThreshColor, linestyle="--", linewidth=0.8)
            Ax.text(0.005, 0.97, "%d columns with a margin below %.1f" % (int(np.nansum(V < LowMargin)), LowMargin), transform=Ax.transAxes,
                    fontsize=7, color=ThreshColor, va="top")
        if Lim:
            Ax.set_ylim(Lim[0], Lim[1] if Lim[1] is not None else max(0.05, float(np.nanmax(Rolling(V))) * 1.3))
        Ax.set_ylabel(Name, fontsize=8, color=Muted)
    #Where copies start and end: piles at the edges of the model are copies with an intact end
    Ax = Axes[4]
    StyleAxis(Ax)
    Bin = 10
    Edges = np.arange(0, d["columns"] + Bin, Bin) + 0.5
    Ax.hist(d["starts"] + 1, bins=Edges, color=DensityColor, alpha=0.8, label="copies starting (5' end of the copy)")
    Ax.hist(d["ends"] + 1, bins=Edges, color=ThreshColor, alpha=0.6, label="copies ending (3' end of the copy)")
    Ax.set_ylabel("Copies per\n%d columns" % Bin, fontsize=8, color=Muted)
    Ax.set_xlabel("Position in the final model", fontsize=8, color=Muted)
    Ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0.1, 1.0), ncol=2)
    Ax.set_xlim(0, d["columns"] + 1)
    Fig.text(0.06, 0.012, "Each copy counts only within its own span (first to last base): gaps outside it are missing data, gaps inside "
             "it deletions. Blue lines: mean over 25 columns.\nCopies that run past the model are cut at its edges, so the first and last "
             "bins of starts and ends also count copies that continue beyond the model; piles inside the model are copies ending there.",
             fontsize=7, color=Muted)
    Fig.tight_layout(rect=(0.02, 0.025, 0.97, 0.95))
    Pdf.savefig(Fig)
    plt.close(Fig)

def IterationsPage(Pdf, Models, Repeat):
    Fig = plt.figure(figsize=(11, 8.5))
    Fig.suptitle("%s: merged alignment and refinement iterations" % Repeat, x=0.06, ha="left", fontsize=13, color=Ink)
    Grid = Fig.add_gridspec(3, 3, height_ratios=(1, 1, 0.75), hspace=0.45, wspace=0.25, left=0.06, right=0.98, top=0.86, bottom=0.05)
    #One hue, light to dark with the iteration; the merged alignment (iteration 0) in grey, dashed
    Refined = [d for d in Models if d["iteration"] > 0]
    Blues = plt.get_cmap("Blues")
    Colors = {d["iteration"]: Blues(0.4 + 0.55 * i / max(1, len(Refined) - 1)) for i, d in enumerate(Refined)}
    Colors[0] = CapColor
    def Style(d):
        return {"color": Colors[d["iteration"]], "linewidth": 2.0 if d["final"] else 1.1, "linestyle": "--" if d["iteration"] == 0 else "-"}
    Panels = [("depth", "Copies spanning a column", "column"), ("agreement", "Copies with the consensus base", "column"),
              ("deletions", "Copies with a deletion", "column"), ("margin", "Call margin (top - second base)", "column"),
              ("divergence", "Copy divergence from the consensus", "copy"), ("span", "Copy length (consensus columns)", "copy")]
    for k, (Key, Name, Unit) in enumerate(Panels):
        Ax = Fig.add_subplot(Grid[k // 3, k % 3])
        StyleAxis(Ax)
        Vals = [np.asarray(d[Key], dtype=float) for d in Models]
        Vals = [v[~np.isnan(v)] for v in Vals]
        Lo = 0.0
        Hi = max(float(np.quantile(v, 0.995)) for v in Vals if len(v)) if any(len(v) for v in Vals) else 1.0
        if Key in ("agreement", "margin"):
            Hi = 1.0
        X = np.linspace(Lo, Hi * 1.05 if Hi > 0 else 1, 200)
        for d, V in zip(Models, Vals):
            if len(np.unique(V)) < 2:
                continue
            Ax.plot(X, gaussian_kde(V, bw_method=0.15)(X), **Style(d))
        if Key == "margin":
            Ax.axvline(LowMargin, color=ThreshColor, linestyle="--", linewidth=0.8)
        Ax.set_xlim(X[0], X[-1])
        Ax.set_xlabel(Name, fontsize=8, color=Muted)
        Ax.set_ylabel("Density (per %s)" % Unit, fontsize=7, color=Muted)
    Handles = [Line2D([], [], label=IterName(d), **Style(d)) for d in Models]
    Fig.legend(handles=Handles, loc="upper left", bbox_to_anchor=(0.06, 0.94), ncol=min(4, len(Handles)), frameon=False, fontsize=8)
    Ax = Fig.add_subplot(Grid[2, :])
    Ax.axis("off")
    Rows = []
    for d in Models:
        S = Summary(d)
        Rows.append(["merged" if d["iteration"] == 0 else "%d" % d["iteration"], "yes" if d["final"] else "", S["copies"], S["columns"],
                     "%.0f" % S["median_depth"], "%.3f" % S["median_agreement"], "%.3f" % S["median_deletions"],
                     S["low_margin_columns"], "%.3f" % S["median_copy_divergence"], "%.0f" % S["median_copy_span"]])
    H = 0.17 * (len(Rows) + 1)
    Tab = Ax.table(cellText=Rows, colLabels=["Iteration", "Final", "Copies", "Columns", "Median depth", "Median agreement",
                                              "Median deletions", "Margin < %.1f" % LowMargin, "Median divergence", "Median copy length"],
                   bbox=(0, 1 - H, 1, H), cellLoc="right")
    Tab.auto_set_font_size(False)
    Tab.set_fontsize(7.5)
    for (r, c), Cell in Tab.get_celld().items():
        Cell.set_edgecolor("#d0d7de")
        Cell.set_linewidth(0.5)
        Cell.get_text().set_color(Ink if r > 0 else Muted)
    Fig.text(0.06, 0.012, "Merged: the stitched alignment of the extension rounds (its rows are pieces of copies from each round). "
             "Each iteration realigns the genome-wide hits of the previous model with hmmalign.", fontsize=7, color=Muted)
    Pdf.savefig(Fig)
    plt.close(Fig)

def main():
    Parser = argparse.ArgumentParser(description=__doc__)
    Parser.add_argument("--workdir", required=True, help="WorkDir/<Repeat> of a run")
    Parser.add_argument("--out", required=True, help="Output prefix")
    Parser.add_argument("--name", default=None, help="Repeat name in the titles (default: the workdir name)")
    Args = Parser.parse_args()
    Repeat = Args.name or os.path.basename(os.path.normpath(Args.workdir))
    Attempts = ReadAttempts(Args.workdir)
    Models = ReadModels(Args.workdir)
    os.makedirs(os.path.dirname(os.path.abspath(Args.out)), exist_ok=True)
    if not Attempts:
        print("No Curation.json files under %s: no pages for the extension rounds (rounds run before the scores were saved)" % Args.workdir)
    if not Attempts and not Models:
        return
    if Attempts:
        WriteTSV(Attempts, Args.out + ".CurationRounds.tsv")
    if Models:
        WriteModelTSVs(Models, Args.out)
    with PdfPages(Args.out + ".CurationReport.pdf") as Pdf:
        if Attempts:
            OverviewPage(Pdf, Attempts, Repeat)
        Final = [d for d in Models if d["final"]]
        if Final:
            FinalModelPage(Pdf, Final[0], Repeat)
        if len(Models) > 1:
            IterationsPage(Pdf, Models, Repeat)
        if Attempts:
            AttemptPages(Pdf, Attempts, Repeat)
    print("Curation report of %d extension attempts and %d models (merged and refinement iterations): %s.CurationReport.pdf" %
          (len(Attempts), len(Models), Args.out))

if __name__ == "__main__":
    main()
