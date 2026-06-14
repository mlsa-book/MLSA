"""Fit two CoxTime variants on the tumor data with the @fig-p3c17-ffnn-arch
(p+1 -> 3 -> 4 -> 1, tanh) architecture and save the comparison figure.

Variants follow the C-labeling parallel to the Weibull M-variants:
    C3: CoxTime on `complications` only (1 + 1 = 2 inputs)
    C4: CoxTime on all 7 features      (7 + 1 = 8 inputs)

Output: book/Figures/neuralnetworks/fig-p3c17-cox-time-tumor.png
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
matplotlib.rcParams.update({"font.size": 12})   # +2 over matplotlib default
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


EXP = Path(__file__).resolve().parent
FIG = EXP.parent / "Figures" / "neuralnetworks" / "fig-p3c17-cox-time-tumor.png"

tumor = pd.read_csv(EXP / "tumor.csv")

# ---- features ------------------------------------------------------------
tumor["complications_bin"] = (tumor["complications"] == "yes").astype(int)
FACTOR_COLS = ["sex", "transfusion", "complications", "metastases", "resection"]
NUMERIC_COLS = ["charlson_score", "age"]
encoded = pd.get_dummies(tumor[FACTOR_COLS], drop_first=True, dtype=float)
X_full = pd.concat([tumor[NUMERIC_COLS].astype(float), encoded], axis=1)
for c in NUMERIC_COLS:
    X_full[c] = (X_full[c] - X_full[c].mean()) / X_full[c].std()

t_np = tumor["days"].values.astype(np.float32)
delta_np = tumor["status"].values.astype(np.float32)
t_max = float(t_np.max())
t_mean, t_std = t_np.mean(), t_np.std()
t_norm = (t_np - t_mean) / t_std

t = torch.tensor(t_np, dtype=torch.float32)
delta = torch.tensor(delta_np, dtype=torch.float32)
tn = torch.tensor(t_norm, dtype=torch.float32)
x_comp = torch.tensor(tumor["complications_bin"].values[:, None],
                      dtype=torch.float32)
x_full = torch.tensor(X_full.values, dtype=torch.float32)


# ---- CoxTime loss (in-loss broadcasting; full O(n^2) on tumor) -----------

def cox_time_nll(g_fn, x: torch.Tensor) -> torch.Tensor:
    ev = torch.nonzero(delta == 1, as_tuple=False).squeeze(-1)
    if ev.numel() == 0:
        return torch.tensor(0.0, requires_grad=True)
    loss = 0.0
    for i in ev:
        ti_n = tn[i:i+1]
        mask = (t >= t[i])
        x_R = x[mask]
        tii = ti_n.expand(x_R.shape[0])
        scores = g_fn(x_R, tii)
        log_sumexp = torch.logsumexp(scores, dim=0)
        score_i = g_fn(x[i:i+1], ti_n)
        loss = loss + (log_sumexp - score_i.squeeze())
    return loss / ev.numel()


class FFNN34CoxTime(torch.nn.Module):
    """Reference FFNN of @fig-ffnn-arch: p+1 -> 3 -> 4 -> 1, tanh activations."""
    def __init__(self, p):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(p + 1, 3), torch.nn.Tanh(),
            torch.nn.Linear(3, 4),     torch.nn.Tanh(),
            torch.nn.Linear(4, 1, bias=False),
        )

    def forward(self, x, t_in):
        z = torch.cat([x, t_in.unsqueeze(-1)], dim=-1)
        return self.net(z).squeeze(-1)


def fit(model, x, epochs=1500, lr=0.02, wd=0.0):
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    for _ in range(epochs):
        opt.zero_grad()
        loss = cox_time_nll(lambda xx, tt: model(xx, tt), x)
        loss.backward()
        opt.step()
    return model.eval()


torch.manual_seed(0)
m_c3 = fit(FFNN34CoxTime(p=1), x_comp)
torch.manual_seed(0)
m_c4 = fit(FFNN34CoxTime(p=x_full.shape[1]), x_full, wd=1e-3)


# ---- CoxTime survival: discrete Breslow over event times -----------------

def predict_survival(model, x_eval, train_x, grid):
    with torch.no_grad():
        ev_idx = np.where(delta_np == 1)[0]
        ev_t = t_np[ev_idx]
        order = np.argsort(ev_t)
        ev_t = ev_t[order]; ev_idx = ev_idx[order]
        ev_tn = (ev_t - t_mean) / t_std

        n_eval = x_eval.shape[0]
        cum_H = np.zeros((n_eval, len(ev_t)))
        for k, (te, ten) in enumerate(zip(ev_t, ev_tn)):
            mask = (t_np >= te)
            x_R = train_x[mask]
            tii_R = torch.full((x_R.shape[0],), ten, dtype=torch.float32)
            denom = torch.logsumexp(model(x_R, tii_R), dim=0)
            tii_E = torch.full((n_eval,), ten, dtype=torch.float32)
            num = model(x_eval, tii_E)
            dH = torch.exp(num - denom)
            if k == 0:
                cum_H[:, k] = dH.numpy()
            else:
                cum_H[:, k] = cum_H[:, k-1] + dH.numpy()

        S = np.ones((n_eval, len(grid)))
        idx = np.searchsorted(ev_t, grid, side="right") - 1
        for j, ix in enumerate(idx):
            if ix >= 0:
                S[:, j] = np.exp(-cum_H[:, ix])
        return S


grid = np.linspace(1.0, t_max, 600)
x_eval_c3 = torch.tensor([[0.0], [1.0]], dtype=torch.float32)
S_c3 = predict_survival(m_c3, x_eval_c3, x_comp, grid)
S_c4 = predict_survival(m_c4, x_full, x_full, grid)


# ---- KM overlays ---------------------------------------------------------

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

COL_NO, COL_YES = "#0072B2", "#D55E00"   # agreed complications palette (no=blue, yes=vermillion)
yes = (tumor["complications"].values == "yes")
no = ~yes

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)

ax = axes[0]
ax.step(km_t_no, km_no, where="post", color=COL_NO, ls=":", lw=1.6, alpha=0.85)
ax.step(km_t_yes, km_yes, where="post", color=COL_YES, ls=":", lw=1.6, alpha=0.85)
ax.plot(grid, S_c3[0], color=COL_NO, lw=2.0)
ax.plot(grid, S_c3[1], color=COL_YES, lw=2.0)
ax.set_title("C3: CoxTime on complications", fontsize=13)
ax.set_xlabel("days")
ax.set_ylabel(r"$\hat S(t \mid \mathbf{x})$")
ax.set_ylim(0, 1.02); ax.set_xlim(0, t_max)
ax.set_box_aspect(1)            # square panel
ax.grid(alpha=0.25)

ax = axes[1]
for i in np.where(no)[0]:
    ax.plot(grid, S_c4[i], color=COL_NO, lw=0.5, alpha=0.18)
for i in np.where(yes)[0]:
    ax.plot(grid, S_c4[i], color=COL_YES, lw=0.5, alpha=0.18)
ax.step(km_t_no, km_no, where="post", color=COL_NO, ls=":", lw=1.6, alpha=0.85)
ax.step(km_t_yes, km_yes, where="post", color=COL_YES, ls=":", lw=1.6, alpha=0.85)
ax.set_title("C4: CoxTime on all features", fontsize=13)
ax.set_xlabel("days")
ax.set_ylim(0, 1.02); ax.set_xlim(0, t_max)
ax.set_box_aspect(1)            # square panel
ax.grid(alpha=0.25)

compl_handles = [
    Line2D([0], [0], color=COL_NO, lw=2, label="no"),
    Line2D([0], [0], color=COL_YES, lw=2, label="yes"),
]
model_handles = [
    Line2D([0], [0], color="#555", lw=2, label="CoxTime"),
    Line2D([0], [0], color="#555", lw=1.6, ls=":", label="Kaplan-Meier"),
]
leg_compl = fig.legend(handles=compl_handles, title="Complications",
                       loc="center left", frameon=False,
                       bbox_to_anchor=(0.86, 0.60), fontsize=12, title_fontsize=12)
fig.add_artist(leg_compl)
fig.legend(handles=model_handles, title="Model",
           loc="center left", frameon=False,
           bbox_to_anchor=(0.86, 0.40), fontsize=12, title_fontsize=12)
fig.tight_layout(rect=(0, 0, 0.84, 1))
fig.savefig(FIG, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved {FIG}")
