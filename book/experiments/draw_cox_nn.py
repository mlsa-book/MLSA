"""Fit two Cox-NN variants (proportional hazards) on the tumor data with
the @fig-ffnn-arch (p -> 4 -> 3 -> 1, tanh) architecture and save the
comparison figure.

    C1: Cox-NN on `complications` only (1 input)
    C2: Cox-NN on all 7 features      (p inputs)

Output: book/Figures/neuralnetworks/cox-nn-tumor.png
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE in sys.path:
    sys.path.remove(_HERE)

from pathlib import Path  # noqa: E402

import torch  # noqa: E402
torch.set_num_threads(1)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


EXP = Path(__file__).resolve().parent
FIG = EXP.parent / "Figures" / "neuralnetworks" / "cox-nn-tumor.png"

tumor = pd.read_csv(EXP / "tumor.csv")
tumor["complications_bin"] = (tumor["complications"] == "yes").astype(int)
FACTOR_COLS = ["sex", "transfusion", "complications", "metastases", "resection"]
NUMERIC_COLS = ["charlson_score", "age"]
encoded = pd.get_dummies(tumor[FACTOR_COLS], drop_first=True, dtype=float)
X_full = pd.concat([tumor[NUMERIC_COLS].astype(float), encoded], axis=1)
for c in NUMERIC_COLS:
    X_full[c] = (X_full[c] - X_full[c].mean()) / X_full[c].std()

t = torch.as_tensor(tumor["days"].values, dtype=torch.float32)
delta = torch.as_tensor(tumor["status"].values, dtype=torch.float32)
x_comp = torch.as_tensor(tumor["complications_bin"].values[:, None],
                         dtype=torch.float32)
x_full = torch.as_tensor(X_full.values, dtype=torch.float32)


def cox_nll(eta, t, delta):
    """Cox partial-likelihood NLL (Breslow tie handling)."""
    order = torch.argsort(t, descending=True)
    eta, _, d_s = eta[order], t[order], delta[order]
    log_cum = torch.logcumsumexp(eta, dim=0)
    return -((eta - log_cum) * d_s).sum() / d_s.sum().clamp_min(1.0)


class FFNN43Cox(torch.nn.Module):
    """Reference FFNN of @fig-ffnn-arch: p -> 4 -> 3 -> 1, tanh activations."""
    def __init__(self, p):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(p, 4), torch.nn.Tanh(),
            torch.nn.Linear(4, 3), torch.nn.Tanh(),
            torch.nn.Linear(3, 1, bias=False),
        )
    def forward(self, x):
        return self.net(x).squeeze(-1)


def fit(model, x, epochs=1500, lr=0.02, wd=0.0):
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    for _ in range(epochs):
        opt.zero_grad()
        loss = cox_nll(model(x), t, delta)
        loss.backward()
        opt.step()
    return model.eval()


torch.manual_seed(0); m1 = fit(FFNN43Cox(p=1), x_comp)
torch.manual_seed(0); m3 = fit(FFNN43Cox(p=x_full.shape[1]), x_full, wd=1e-3)


# ---- Breslow baseline + survival functions -------------------------------

def breslow_S(eta_train, eta_eval, time_grid):
    with torch.no_grad():
        order = torch.argsort(t)
        t_s, d_s, e_s = t[order], delta[order], eta_train[order]
    exp_e = torch.exp(e_s).numpy()
    t_arr = t_s.numpy()
    d_arr = d_s.numpy()
    cum_back = np.cumsum(exp_e[::-1])[::-1]
    haz_jumps = np.where(d_arr == 1, 1.0 / cum_back, 0.0)
    H0 = np.zeros_like(time_grid)
    for k, tau in enumerate(time_grid):
        H0[k] = haz_jumps[t_arr <= tau].sum()
    e_eval = np.exp(eta_eval.numpy())
    return np.exp(-np.outer(e_eval, H0))


grid = np.linspace(1.0, float(tumor["days"].max()), 600)

with torch.no_grad():
    eta_train_m1 = m1(x_comp)
    eta_train_m3 = m3(x_full)
    x_eval_m1 = torch.tensor([[0.0], [1.0]], dtype=torch.float32)
    eta_eval_m1 = m1(x_eval_m1)
    eta_eval_m3 = m3(x_full)

S1 = breslow_S(eta_train_m1, eta_eval_m1, grid)
S3 = breslow_S(eta_train_m3, eta_eval_m3, grid)


def km(times, events):
    order = np.argsort(times)
    tt, ee = times[order], events[order]
    uniq = np.unique(tt)
    s, ar = 1.0, len(tt)
    out = []
    for u in uniq:
        m = tt == u; d = ee[m].sum()
        s *= 1.0 - d / ar if ar > 0 else 1.0
        out.append(s); ar -= m.sum()
    return uniq, np.asarray(out)


km_t_no, km_no = km(tumor.loc[tumor["complications"] == "no", "days"].values,
                    tumor.loc[tumor["complications"] == "no", "status"].values)
km_t_yes, km_yes = km(tumor.loc[tumor["complications"] == "yes", "days"].values,
                      tumor.loc[tumor["complications"] == "yes", "status"].values)


# ---- plot ----------------------------------------------------------------

COL_NO, COL_YES = "#2C7FB8", "#C2185B"
yes = (tumor["complications"].values == "yes")
no = ~yes

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)

ax = axes[0]
ax.step(km_t_no, km_no, where="post", color=COL_NO, ls=":", lw=1.6, alpha=0.85)
ax.step(km_t_yes, km_yes, where="post", color=COL_YES, ls=":", lw=1.6, alpha=0.85)
ax.plot(grid, S1[0], color=COL_NO, lw=2.0)
ax.plot(grid, S1[1], color=COL_YES, lw=2.0)
ax.set_title("C1: Cox-NN on complications", fontsize=11)
ax.set_xlabel("days"); ax.set_ylabel(r"$\hat S(t \mid \mathbf{x})$")
ax.set_ylim(0, 1.02); ax.set_xlim(0, float(tumor["days"].max()))
ax.grid(alpha=0.25)

ax = axes[1]
for i in np.where(no)[0]:
    ax.plot(grid, S3[i], color=COL_NO, lw=0.5, alpha=0.18)
for i in np.where(yes)[0]:
    ax.plot(grid, S3[i], color=COL_YES, lw=0.5, alpha=0.18)
ax.step(km_t_no, km_no, where="post", color=COL_NO, ls=":", lw=1.6, alpha=0.85)
ax.step(km_t_yes, km_yes, where="post", color=COL_YES, ls=":", lw=1.6, alpha=0.85)
ax.set_title("C2: Cox-NN on all features", fontsize=11)
ax.set_xlabel("days")
ax.set_ylim(0, 1.02); ax.set_xlim(0, float(tumor["days"].max()))
ax.grid(alpha=0.25)

handles = [
    Line2D([0], [0], color=COL_NO, lw=2, label="no complications"),
    Line2D([0], [0], color=COL_YES, lw=2, label="complications"),
    Line2D([0], [0], color="#555", lw=2, label="Cox-NN"),
    Line2D([0], [0], color="#555", lw=1.6, ls=":", label="Kaplan-Meier"),
]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, -0.04), fontsize=10)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(FIG, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved {FIG}")
