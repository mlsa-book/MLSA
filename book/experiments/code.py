"""Python counterpart to ``code.R``.

Each block below produces a figure used by the book.  The structure mirrors
``code.R``: a leading comment block names the chapter / figure, followed by
self-contained code that loads the data, fits the model(s) and writes the
PNG into the appropriate ``book/Figures/...`` sub-folder.
"""

# ---------------------------------------------------------------------------
# Neural Networks chapter – Weibull AFT with a neural network.
#
# Tumor data (pammtools), three increasingly flexible parametrisations:
#   M1: scale depends on complications, shape is global
#   M2: scale and shape depend on complications
#   M3: scale and shape depend on all features (individualised prediction)
#
# All panels overlay the Kaplan-Meier estimate stratified by complications
# so the effect of letting the network parametrise more pieces of the
# Weibull distribution is directly visible.
#
# Output: book/Figures/neuralnetworks/weibull-aft-nn.png
# ---------------------------------------------------------------------------

from __future__ import annotations

import os
import sys

# ``code.py`` shadows Python's stdlib ``code`` module, which is used
# transitively by ``pdb`` during ``import torch``.  Drop our own directory
# from ``sys.path`` before importing torch so that ``import code`` resolves
# to the stdlib version (otherwise: segfault during torch's lazy imports).
_SELF_DIR = os.path.dirname(os.path.abspath(__file__))
if _SELF_DIR in sys.path:
    sys.path.remove(_SELF_DIR)

import subprocess  # noqa: E402
from pathlib import Path  # noqa: E402

import torch  # noqa: E402
torch.set_num_threads(1)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


# ---- I/O paths -----------------------------------------------------------

EXP_DIR = Path(__file__).resolve().parent
CSV_PATH = EXP_DIR / "tumor.csv"
FIG_OUT = (
    EXP_DIR.parent / "Figures" / "neuralnetworks" / "weibull-aft-nn.png"
)
FIG_OUT.parent.mkdir(parents=True, exist_ok=True)


def ensure_tumor_csv() -> pd.DataFrame:
    """Load tumor.csv; export from the pammtools R package on first run."""
    if not CSV_PATH.exists():
        subprocess.run(
            [
                "Rscript",
                "-e",
                'data("tumor", package="pammtools");'
                f' write.csv(tumor, "{CSV_PATH}", row.names=FALSE)',
            ],
            check=True,
        )
    return pd.read_csv(CSV_PATH)


# ---- Data preparation -----------------------------------------------------

tumor = ensure_tumor_csv()

# binary / one-hot encode all factor features.  ``complications`` is kept
# both as a stand-alone feature (for M1/M2) and as part of the full feature
# matrix (for M3).
tumor["complications_bin"] = (tumor["complications"] == "yes").astype(int)
FACTOR_COLS = ["sex", "transfusion", "complications", "metastases", "resection"]
NUMERIC_COLS = ["charlson_score", "age"]
encoded = pd.get_dummies(
    tumor[FACTOR_COLS], drop_first=True, dtype=float
)
X_full = pd.concat([tumor[NUMERIC_COLS].astype(float), encoded], axis=1)

# numeric features are standardised; binary dummies are left as-is.
for col in NUMERIC_COLS:
    X_full[col] = (X_full[col] - X_full[col].mean()) / X_full[col].std()

t = torch.as_tensor(tumor["days"].values, dtype=torch.float32)
delta = torch.as_tensor(tumor["status"].values, dtype=torch.float32)
log_t = torch.log(t.clamp(min=1.0))

x_comp = torch.as_tensor(
    tumor["complications_bin"].values[:, None], dtype=torch.float32
)
x_full = torch.as_tensor(X_full.values, dtype=torch.float32)


# ---- Weibull AFT log-likelihood ------------------------------------------

def weibull_nll(
    log_t: torch.Tensor,
    delta: torch.Tensor,
    log_lambda: torch.Tensor,
    log_k: torch.Tensor,
) -> torch.Tensor:
    """Negative Weibull right-censored log-likelihood.

    The density is parametrised as ``f(t) = (k/λ)(t/λ)^(k-1) exp(-(t/λ)^k)``,
    with ``λ = exp(log_lambda)`` (scale) and ``k = exp(log_k)`` (shape).
    Computed in log-space for numerical stability.
    """
    k = torch.exp(log_k)
    z = k * (log_t - log_lambda)              # log of (t / λ)^k
    log_S = -torch.exp(z)
    log_f = log_k - log_lambda + (k - 1.0) * (log_t - log_lambda) - torch.exp(z)
    return -(delta * log_f + (1.0 - delta) * log_S).mean()


# ---- Neural-network heads -------------------------------------------------

class WeibullNet(torch.nn.Module):
    """Two-headed network producing log-scale and log-shape.

    Each head emits a *raw* output that is mapped through a bounded
    activation: log-scale is constrained to ``log_t_init ± scale_range`` and
    log-shape to ``±shape_range``.  Without these bounds an unregularised
    MLP head will produce pathological Weibull parameters on small datasets
    (λ → ∞ or k → ∞) and therefore visibly bad individualised curves.

    A ``None`` head means the corresponding parameter is a single global
    scalar (independent of x).  A ``"linear"`` head is a 1-layer linear map;
    a ``"mlp"`` head adds one hidden layer with ``hidden`` ReLU units.
    """

    def __init__(
        self,
        p: int,
        scale_head: str = "linear",
        shape_head: str | None = None,
        hidden: int = 16,
        log_t_init: float = 7.0,
        scale_range: float = 2.5,
        shape_range: float = 1.5,
    ) -> None:
        super().__init__()
        self.scale_head = self._make_head(p, scale_head, hidden)
        self.shape_head = self._make_head(p, shape_head, hidden)
        self.global_log_k = torch.nn.Parameter(torch.zeros(1))
        self.log_t_init = log_t_init
        self.scale_range = scale_range
        self.shape_range = shape_range

    @staticmethod
    def _make_head(p: int, kind: str | None, hidden: int):
        if kind is None:
            return None
        if kind == "linear":
            return torch.nn.Linear(p, 1)
        if kind == "mlp":
            return torch.nn.Sequential(
                torch.nn.Linear(p, hidden),
                torch.nn.ReLU(),
                torch.nn.Linear(hidden, 1),
            )
        raise ValueError(f"unknown head: {kind}")

    def forward(self, x: torch.Tensor):
        raw_scale = self.scale_head(x).squeeze(-1)
        log_lambda = self.log_t_init + self.scale_range * torch.tanh(raw_scale)
        if self.shape_head is None:
            log_k = self.global_log_k.expand(x.shape[0])
        else:
            raw_shape = self.shape_head(x).squeeze(-1)
            log_k = self.shape_range * torch.tanh(raw_shape)
        return log_lambda, log_k


def fit(model: torch.nn.Module, x: torch.Tensor, *, epochs: int = 3000,
        lr: float = 0.05, weight_decay: float = 0.0) -> torch.nn.Module:
    opt = torch.optim.Adam(model.parameters(), lr=lr,
                           weight_decay=weight_decay)
    for _ in range(epochs):
        opt.zero_grad()
        log_lambda, log_k = model(x)
        loss = weibull_nll(log_t, delta, log_lambda, log_k)
        loss.backward()
        opt.step()
    return model.eval()


log_t_init = float(torch.log(t[delta == 1].mean()).item())

torch.manual_seed(0)
m1 = WeibullNet(p=1, scale_head="linear", shape_head=None,
                log_t_init=log_t_init)

torch.manual_seed(0)
m2 = WeibullNet(p=1, scale_head="linear", shape_head="linear",
                log_t_init=log_t_init)

torch.manual_seed(0)
m3 = WeibullNet(p=x_full.shape[1], scale_head="mlp", shape_head="mlp",
                hidden=16, log_t_init=log_t_init)

fit(m1, x_comp)
fit(m2, x_comp)
# M3 has many more parameters relative to the data; a light weight decay
# keeps both heads from drifting to the bounds of the tanh range.
fit(m3, x_full, weight_decay=1e-3)


# ---- Kaplan-Meier estimator (right-censored) -----------------------------

def kaplan_meier(
    times: np.ndarray, events: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(times)
    times, events = times[order], events[order]
    unique = np.unique(times)
    surv = np.ones_like(unique, dtype=float)
    n = len(times)
    s = 1.0
    at_risk = n
    out = []
    for u in unique:
        mask_at_u = times == u
        d = events[mask_at_u].sum()
        if at_risk > 0:
            s *= 1.0 - d / at_risk
        out.append(s)
        at_risk -= mask_at_u.sum()
    return unique, np.asarray(out)


km_times_no, km_no = kaplan_meier(
    tumor.loc[tumor["complications"] == "no", "days"].values,
    tumor.loc[tumor["complications"] == "no", "status"].values,
)
km_times_yes, km_yes = kaplan_meier(
    tumor.loc[tumor["complications"] == "yes", "days"].values,
    tumor.loc[tumor["complications"] == "yes", "status"].values,
)


# ---- Predicted survival curves -------------------------------------------

grid = torch.linspace(1.0, float(tumor["days"].max()), 600)
log_grid = torch.log(grid)


def weibull_S(log_lambda: torch.Tensor, log_k: torch.Tensor,
              log_grid: torch.Tensor) -> torch.Tensor:
    """Survival function on the time grid for one or many subjects."""
    log_lambda = log_lambda.unsqueeze(-1)
    log_k = log_k.unsqueeze(-1)
    k = torch.exp(log_k)
    z = k * (log_grid.unsqueeze(0) - log_lambda)
    return torch.exp(-torch.exp(z))


def two_group_predictions(model: torch.nn.Module) -> tuple:
    x_eval = torch.tensor([[0.0], [1.0]], dtype=torch.float32)
    with torch.no_grad():
        log_lambda, log_k = model(x_eval)
        S = weibull_S(log_lambda, log_k, log_grid).numpy()
    return S[0], S[1]


S1_no, S1_yes = two_group_predictions(m1)
S2_no, S2_yes = two_group_predictions(m2)

with torch.no_grad():
    log_lambda3, log_k3 = m3(x_full)
    S3_ind = weibull_S(log_lambda3, log_k3, log_grid).numpy()


# ---- Plotting ------------------------------------------------------------

COL_NO = "#2C7FB8"
COL_YES = "#C2185B"
GRID_T = grid.numpy()

fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6), sharey=True)
titles = [
    "M1: $\\log\\lambda(x)$, $k$ constant",
    "M2: $\\log\\lambda(x)$, $\\log k(x)$",
    "M3: full $\\log\\lambda(\\mathbf{x})$, $\\log k(\\mathbf{x})$",
]

for ax, title in zip(axes, titles):
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("days")
    ax.set_xlim(0, float(tumor["days"].max()))
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.25)

axes[0].set_ylabel("$\\hat S(t \\mid x)$")


def overlay_km(ax: plt.Axes) -> None:
    ax.step(km_times_no, km_no, where="post", color=COL_NO,
            linestyle=":", linewidth=1.6, alpha=0.85, label="KM (no compl.)")
    ax.step(km_times_yes, km_yes, where="post", color=COL_YES,
            linestyle=":", linewidth=1.6, alpha=0.85, label="KM (compl.)")


# Panel 1
overlay_km(axes[0])
axes[0].plot(GRID_T, S1_no, color=COL_NO, linewidth=2.0,
             label="AFT (no compl.)")
axes[0].plot(GRID_T, S1_yes, color=COL_YES, linewidth=2.0,
             label="AFT (compl.)")

# Panel 2
overlay_km(axes[1])
axes[1].plot(GRID_T, S2_no, color=COL_NO, linewidth=2.0)
axes[1].plot(GRID_T, S2_yes, color=COL_YES, linewidth=2.0)

# Panel 3: individualised curves coloured by complications group.
yes_mask = tumor["complications"].values == "yes"
no_mask = ~yes_mask
for i in np.where(no_mask)[0]:
    axes[2].plot(GRID_T, S3_ind[i], color=COL_NO, linewidth=0.5, alpha=0.18)
for i in np.where(yes_mask)[0]:
    axes[2].plot(GRID_T, S3_ind[i], color=COL_YES, linewidth=0.5, alpha=0.18)
overlay_km(axes[2])

# Single legend at the bottom.
legend_handles = [
    Line2D([0], [0], color=COL_NO, linewidth=2.0,
           label="Weibull AFT NN, no complications"),
    Line2D([0], [0], color=COL_YES, linewidth=2.0,
           label="Weibull AFT NN, complications"),
    Line2D([0], [0], color=COL_NO, linewidth=1.6, linestyle=":",
           label="Kaplan-Meier, no complications"),
    Line2D([0], [0], color=COL_YES, linewidth=1.6, linestyle=":",
           label="Kaplan-Meier, complications"),
]
fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=2,
    frameon=False,
    fontsize=10,
    bbox_to_anchor=(0.5, -0.04),
)

fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(FIG_OUT, dpi=200, bbox_inches="tight")
print(f"Saved {FIG_OUT}")
