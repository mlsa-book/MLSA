"""Preview-only Deep Cox Mixtures (DCM) on the tumor data — full feature
set, K=3 components. Saves both into the book Figures folder and to /tmp.

Shared FFNN g(x; theta) -> R^{2K}:
    K mixture-weight logits  -> softmax -> alpha_k(x)
    K component risk scores  -> eta_k(x)
Each component has its own piecewise-constant baseline hazard
h_{0,k}(t) on the unique event-time grid, parametrised as a free vector
log_h0[k, b]. Trained end-to-end with Adam on the mixture NLL.

Output: book/Figures/neuralnetworks/dcm-tumor.png
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
FIG = EXP.parent / "Figures" / "neuralnetworks" / "dcm-tumor.png"

tumor = pd.read_csv(EXP / "tumor.csv")
tumor["complications_bin"] = (tumor["complications"] == "yes").astype(int)
FACTOR_COLS = ["sex", "transfusion", "complications", "metastases", "resection"]
NUMERIC_COLS = ["charlson_score", "age"]
encoded = pd.get_dummies(tumor[FACTOR_COLS], drop_first=True, dtype=float)
X_full = pd.concat([tumor[NUMERIC_COLS].astype(float), encoded], axis=1)
for c in NUMERIC_COLS:
    X_full[c] = (X_full[c] - X_full[c].mean()) / X_full[c].std()

t_np = tumor["days"].values.astype(np.float32)
delta_np = tumor["status"].values.astype(np.float32)
x_full = torch.tensor(X_full.values, dtype=torch.float32)

event_times = np.unique(t_np[delta_np == 1])
B = len(event_times)
bin_idx = np.clip(np.searchsorted(event_times, t_np, side="right") - 1,
                  0, B - 1)
bin_idx_t = torch.tensor(bin_idx, dtype=torch.long)
delta_t = torch.tensor(delta_np, dtype=torch.float32)


K = 3


class DCM(torch.nn.Module):
    """Shared FFNN(p -> 4 -> 3) -> two linear heads of width K
    (mixture-weight logits and component risk scores).
    """
    def __init__(self, p: int, K: int):
        super().__init__()
        self.K = K
        self.shared = torch.nn.Sequential(
            torch.nn.Linear(p, 4), torch.nn.Tanh(),
            torch.nn.Linear(4, 3), torch.nn.Tanh(),
        )
        self.head_alpha = torch.nn.Linear(3, K, bias=False)
        self.head_eta   = torch.nn.Linear(3, K, bias=False)

    def forward(self, x):
        h = self.shared(x)
        log_alpha = torch.log_softmax(self.head_alpha(h), dim=-1)
        eta = self.head_eta(h)
        return log_alpha, eta


log_h0 = torch.nn.Parameter(torch.full((K, B), -6.0))


def dcm_nll(model, x, bin_idx, delta):
    log_alpha, eta = model(x)
    h0 = torch.exp(log_h0)
    H0 = torch.cumsum(h0, dim=1)
    H0_i = H0[:, bin_idx].T
    exp_eta = torch.exp(eta)
    log_S_k = -H0_i * exp_eta
    log_h0_i = log_h0[:, bin_idx].T
    log_f_k = log_h0_i + eta + log_S_k
    log_L_event = torch.logsumexp(log_alpha + log_f_k, dim=-1)
    log_L_cens  = torch.logsumexp(log_alpha + log_S_k, dim=-1)
    log_L = delta * log_L_event + (1.0 - delta) * log_L_cens
    return -log_L.mean()


torch.manual_seed(0)
model = DCM(p=x_full.shape[1], K=K)
opt = torch.optim.Adam(list(model.parameters()) + [log_h0],
                       lr=0.02, weight_decay=1e-3)
for epoch in range(2500):
    opt.zero_grad()
    loss = dcm_nll(model, x_full, bin_idx_t, delta_t)
    loss.backward()
    opt.step()
    if epoch % 500 == 0:
        print(f"  epoch {epoch:4d}  NLL = {loss.item():.4f}")
print(f"  final NLL = {loss.item():.4f}")

with torch.no_grad():
    log_alpha, eta = model(x_full)
    alpha = torch.exp(log_alpha)
    H0 = torch.cumsum(torch.exp(log_h0), dim=1)
    grid = np.linspace(1.0, float(tumor["days"].max()), 500)
    grid_bin = np.clip(np.searchsorted(event_times, grid, side="right") - 1,
                       0, B - 1)
    H0_grid = H0.numpy()[:, grid_bin]
    exp_eta = torch.exp(eta).numpy()
    S_k = np.exp(-H0_grid[None, :, :] * exp_eta[:, :, None])
    S = (alpha.numpy()[:, :, None] * S_k).sum(axis=1)


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


COL_NO, COL_YES = "#2C7FB8", "#C2185B"
yes = (tumor["complications"].values == "yes")
no = ~yes

fig, ax = plt.subplots(figsize=(7, 4.6))
for i in np.where(no)[0]:
    ax.plot(grid, S[i], color=COL_NO, lw=0.5, alpha=0.18)
for i in np.where(yes)[0]:
    ax.plot(grid, S[i], color=COL_YES, lw=0.5, alpha=0.18)
ax.step(km_t_no, km_no, where="post", color=COL_NO, ls=":", lw=1.6, alpha=0.85)
ax.step(km_t_yes, km_yes, where="post", color=COL_YES, ls=":", lw=1.6, alpha=0.85)
ax.set_title(f"Deep Cox Mixtures (K={K}) — full feature set", fontsize=11)
ax.set_xlabel("days"); ax.set_ylabel(r"$\hat S(t \mid \mathbf{x})$")
ax.set_ylim(0, 1.02); ax.set_xlim(0, float(tumor["days"].max()))
ax.grid(alpha=0.25)
handles = [
    Line2D([0], [0], color=COL_NO, lw=2, label="no complications"),
    Line2D([0], [0], color=COL_YES, lw=2, label="complications"),
    Line2D([0], [0], color="#555", lw=2, label="DCM"),
    Line2D([0], [0], color="#555", lw=1.6, ls=":", label="Kaplan-Meier"),
]
ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=10)
fig.tight_layout()
fig.savefig(FIG, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved {FIG}")
