from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

run = Path(sys.argv[1])
outputs = [run / "mesh-quality.png", run / "force-flow.png"]
if any(p.exists() for p in outputs):
    raise SystemExit("refusing to overwrite fitted-inlet plots")
fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
for ax, nt in zip(axes, (128, 256, 512)):
    d = np.genfromtxt(run / f"patch-{nt}.tsv", delimiter="\t", names=True)
    quads = np.stack([np.column_stack((d[f"x{i}_mm"], d[f"z{i}_mm"])) for i in range(1, 5)], axis=1)
    pc = PolyCollection(quads, array=d["corner"], cmap="viridis", edgecolors="k", linewidths=0.25, clim=(0.65, 1))
    ax.add_collection(pc)
    a = np.linspace(0, 2*np.pi, 400)
    ax.plot(2.5*np.cos(a), 2.5*np.sin(a), "r-", linewidth=0.9)
    worst = np.argmin(d["corner"])
    xy = quads[worst].mean(axis=0)
    ax.plot(*xy, "rx", markersize=9)
    ax.set(xlim=(-11, 11), ylim=(-11, 11), aspect="equal", xlabel="Drill-plane X (mm)", ylabel="z − midpoint (mm)",
           title=f"Level {nt//128}: min corner = {d['corner'].min():.3f}")
fig.colorbar(pc, ax=axes, label="Signed scaled corner score (NOT 3D Hex8 OQ)", shrink=0.8)
fig.suptitle("Same 5 mm fitted mouth; transition quality worsens on refinement\nRed circle: exact drill rim; cross: poorest displayed cell")
fig.savefig(outputs[0], dpi=160)
plt.close(fig)
r = np.genfromtxt(run / "comparison.tsv", delimiter="\t", names=True, dtype=None, encoding="utf-8")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
for kind, style in (("regular", "o-"), ("fitted", "s-")):
    d = r[r["mesh"] == kind]
    axes[0].plot(d["cells"], d["Fz_N"], style, label=kind)
    axes[1].plot(d["cells"], d["net_feed_m3_s"] * 60000, style, label=kind)
    radial = np.hypot(d["Fx_N"], d["Fy_N"] + 1000)
    axes[2].plot(d["cells"], radial, style, label=kind)
axes[0].axhline(200, color="C3", linestyle="--", label="Axial target")
axes[0].set(ylabel="Axial film force (N)", ylim=(0, None))
axes[1].set(ylabel="Net feed (L/min)", ylim=(0, None))
axes[2].set(ylabel="Transverse force residual at fixed pose (N)", yscale="log")
for ax in axes:
    ax.set(xlabel="Film cells", xscale="log")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
fig.suptitle("Same bush position and physical inputs; no equilibrium re-solve\nFitted candidate is diagnostic, not an accepted production mesh")
fig.savefig(outputs[1], dpi=160)
plt.close(fig)
print(*outputs, sep="\n")
