"""Generate the loss-surface + gradient-descent figure for the
"Gradient descent" section of the ML chapter (P1C3).

Output: ``book/Figures/introduction/gradient-descent.png``.
"""
from __future__ import annotations

import os
import sys

# Avoid the stdlib-``code`` shadow trap (cf. notes in code.py).
_SELF_DIR = os.path.dirname(os.path.abspath(__file__))
if _SELF_DIR in sys.path:
    sys.path.remove(_SELF_DIR)

from pathlib import Path  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


EXP_DIR = Path(__file__).resolve().parent
FIG_DIR = EXP_DIR.parent / "Figures" / "introduction"
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ---- Load data ------------------------------------------------------------

df = pd.read_csv(EXP_DIR / "teengamb.csv")
x = df["income"].values.astype(float)
y = df["gamble"].values.astype(float)
n = len(y)

# Closed-form OLS minimum.
xm, ym = x.mean(), y.mean()
slope_hat = ((x - xm) * (y - ym)).sum() / ((x - xm) ** 2).sum()
intercept_hat = ym - slope_hat * xm


# ---- MSE on a grid of parameter values ------------------------------------

def mse(theta0, theta1):
    pred = theta0 + theta1 * x
    return ((y - pred) ** 2).mean()


t0_grid = np.linspace(-10, 10, 250)
t1_grid = np.linspace(0, 8, 250)
T0, T1 = np.meshgrid(t0_grid, t1_grid)
L = np.zeros_like(T0)
for i in range(T0.shape[0]):
    for j in range(T0.shape[1]):
        L[i, j] = mse(T0[i, j], T1[i, j])


# ---- Gradient descent -----------------------------------------------------

def grad_mse(theta0, theta1):
    pred = theta0 + theta1 * x
    err = pred - y
    g0 = 2 * err.mean()
    g1 = 2 * (err * x).mean()
    return np.array([g0, g1])


theta = np.array([8.0, 1.0])  # poor initial guess in the bottom-right region
lr = 0.01
n_steps = 800
path = [theta.copy()]
for _ in range(n_steps):
    theta = theta - lr * grad_mse(*theta)
    path.append(theta.copy())
path = np.array(path)


# ---- Plot -----------------------------------------------------------------

fig, ax = plt.subplots(figsize=(6.4, 5.0))

levels = np.linspace(L.min(), L.max(), 24)
cf = ax.contourf(T0, T1, L, levels=levels, cmap="YlOrRd", alpha=0.95)
ax.contour(T0, T1, L, levels=levels, colors="white",
           linewidths=0.4, alpha=0.7)
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label("MSE", fontsize=16)
cbar.ax.tick_params(labelsize=13)

# Start.
ax.plot(path[0, 0], path[0, 1], marker="s", markersize=11,
        color="#2E7D32", markeredgecolor="black", markeredgewidth=0.7,
        linestyle="None", label="Start")
# OLS optimum.
ax.plot(intercept_hat, slope_hat, marker="X", markersize=15,
        color="black", linestyle="None", label="OLS")
# GD path.
ax.plot(path[:, 0], path[:, 1], "-", color="#1F77B4",
        linewidth=1.6, label="Path")
# End.
ax.plot(path[-1, 0], path[-1, 1], marker="o", markersize=9,
        color="#D32F2F", markeredgecolor="black", markeredgewidth=0.5,
        linestyle="None", label="End")

ax.set_xlabel(r"$\theta_0$ (intercept)", fontsize=16)
ax.set_ylabel(r"$\theta_1$ (slope)", fontsize=16)
ax.tick_params(axis="both", labelsize=13)
ax.set_xlim(t0_grid.min(), t0_grid.max())
ax.set_ylim(t1_grid.min(), t1_grid.max())
ax.legend(loc="upper right", frameon=True, fontsize=14)

fig.tight_layout()
out = FIG_DIR / "gradient-descent.png"
fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved {out}")
print(f"OLS optimum: theta0={intercept_hat:.3f}, theta1={slope_hat:.3f}")
print(f"GD final:    theta0={path[-1,0]:.3f}, theta1={path[-1,1]:.3f}")
