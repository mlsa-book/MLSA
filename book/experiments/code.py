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
FIG_DIR = FIG_OUT.parent
FIG_DIR.mkdir(parents=True, exist_ok=True)


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

# Legend: colour = complications status, linestyle = model type.
legend_handles = [
    Line2D([0], [0], color=COL_NO, linewidth=2.0,
           label="no complications"),
    Line2D([0], [0], color=COL_YES, linewidth=2.0,
           label="complications"),
    Line2D([0], [0], color="#555", linewidth=2.0, linestyle="-",
           label="Weibull AFT NN"),
    Line2D([0], [0], color="#555", linewidth=1.6, linestyle=":",
           label="Kaplan-Meier"),
]
fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=4,
    frameon=False,
    fontsize=10,
    bbox_to_anchor=(0.5, -0.04),
)

fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(FIG_OUT, dpi=200, bbox_inches="tight")
print(f"Saved {FIG_OUT}")


# ---------------------------------------------------------------------------
# Neural Networks chapter – Activation function shape plots.
#
# Small individual plots for the activation table (@tbl-activations).
#
# Output: book/Figures/neuralnetworks/activation-{relu,sigmoid,tanh,softplus}.png
# ---------------------------------------------------------------------------

v = np.linspace(-4, 4, 400)

_activations = [
    ("relu",    np.maximum(0, v)),
    ("sigmoid", 1 / (1 + np.exp(-v))),
    ("tanh",    np.tanh(v)),
    ("softplus", np.log1p(np.exp(v))),
]

for name, vals in _activations:
    fig, ax = plt.subplots(figsize=(3.2, 2.2))
    ax.plot(v, vals, color="#C2185B", lw=2.5)
    ax.axhline(0, color="#888", lw=0.5, ls="--")
    ax.axvline(0, color="#888", lw=0.5, ls="--")
    ax.set_xlim(-4, 4)
    ax.set_xlabel("v", fontsize=13)
    ax.set_ylabel("a(v)", fontsize=13)
    ax.tick_params(labelsize=11)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out = FIG_DIR / f"activation-{name}.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Neural Networks chapter – mcycle: tanh vs ReLU comparison.
#
# Two-panel figure showing the same FFNN architecture (5, 3) fitted to
# mcycle with tanh and ReLU activations respectively.
#
# Output: book/Figures/neuralnetworks/mcycle-tanh-vs-relu.png
# ---------------------------------------------------------------------------

mcycle = pd.read_csv(EXP_DIR / "mcycle.csv")

_x = mcycle["times"].values.astype(np.float32)
_y = mcycle["accel"].values.astype(np.float32)
_xm, _xs, _ym, _ys = _x.mean(), _x.std(), _y.mean(), _y.std()
_xt = torch.as_tensor((_x - _xm) / _xs).unsqueeze(-1)
_yt = torch.as_tensor((_y - _ym) / _ys)


def _fit_ffnn(xt, yt, hidden=(5, 3), activation=torch.nn.Tanh,
              epochs=5000, lr=5e-3, wd=1e-3):
    layers = []
    d = xt.shape[-1]
    for h in hidden:
        layers += [torch.nn.Linear(d, h), activation()]
        d = h
    layers.append(torch.nn.Linear(d, 1))
    model = torch.nn.Sequential(*layers)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    for _ in range(epochs):
        opt.zero_grad()
        loss = ((model(xt).squeeze(-1) - yt) ** 2).mean()
        loss.backward()
        opt.step()
    return model.eval()


_grid = np.linspace(_x.min(), _x.max(), 400)
_gt = torch.as_tensor((_grid - _xm) / _xs, dtype=torch.float32).unsqueeze(-1)

torch.manual_seed(0)
_m_tanh = _fit_ffnn(_xt, _yt, activation=torch.nn.Tanh)
torch.manual_seed(0)
_m_relu = _fit_ffnn(_xt, _yt, activation=torch.nn.ReLU)

with torch.no_grad():
    _p_tanh = _m_tanh(_gt).squeeze(-1).numpy() * _ys + _ym
    _p_relu = _m_relu(_gt).squeeze(-1).numpy() * _ys + _ym

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, pred, title in zip(axes, [_p_tanh, _p_relu],
                           ["tanh activation", "ReLU activation"]):
    ax.scatter(_x, _y, s=14, alpha=0.55, color="#5B7CA8", edgecolor="none")
    ax.plot(_grid, pred, color="#C2185B", lw=2.2)
    ax.set_title(title)
    ax.set_xlabel("time after impact (ms)")
    ax.grid(alpha=0.25)
axes[0].set_ylabel("head acceleration (g)")
fig.tight_layout()
fig.savefig(FIG_DIR / "mcycle-tanh-vs-relu.png", dpi=220, bbox_inches="tight")
plt.close(fig)
print(f"Saved {FIG_DIR / 'mcycle-tanh-vs-relu.png'}")


# ---------------------------------------------------------------------------
# Neural Networks chapter – mcycle: distributional regression.
#
# Two-headed NN (shared trunk, variant (a)) outputting mu(x) and sigma(x),
# trained with heteroscedastic Gaussian NLL.
#
# Output: book/Figures/neuralnetworks/distributional-regression.png
# ---------------------------------------------------------------------------


class _DistributionalFFNN(torch.nn.Module):
    def __init__(self, hidden=32):
        super().__init__()
        self.trunk = torch.nn.Sequential(
            torch.nn.Linear(1, hidden), torch.nn.Tanh(),
            torch.nn.Linear(hidden, hidden), torch.nn.Tanh(),
        )
        self.head_mu = torch.nn.Linear(hidden, 1)
        self.head_log_sigma = torch.nn.Linear(hidden, 1)

    def forward(self, x):
        h = self.trunk(x)
        return self.head_mu(h).squeeze(-1), self.head_log_sigma(h).squeeze(-1)


torch.manual_seed(2)
_dist_model = _DistributionalFFNN(hidden=32)
_opt = torch.optim.Adam(_dist_model.parameters(), lr=5e-3, weight_decay=1e-2)
for _ in range(3000):
    _opt.zero_grad()
    _mu, _ls = _dist_model(_xt)
    _loss = (0.5 * ((_yt - _mu) / torch.exp(_ls)) ** 2 + _ls).mean()
    _loss.backward()
    _opt.step()
_dist_model.eval()

with torch.no_grad():
    _mu_mc, _ls_mc = _dist_model(_gt)
_mu_mc = _mu_mc.numpy() * _ys + _ym
_sigma_mc = np.exp(_ls_mc.numpy()) * _ys

fig, ax = plt.subplots(figsize=(5.6, 3.7))
ax.scatter(_x, _y, s=14, alpha=0.55, color="#5B7CA8", edgecolor="none")
ax.fill_between(_grid, _mu_mc - 1.96 * _sigma_mc, _mu_mc + 1.96 * _sigma_mc,
                color="#C2185B", alpha=0.18,
                label="$\\hat μ(x) \\pm 1.96\\,\\hat σ(x)$")
ax.plot(_grid, _mu_mc, color="#C2185B", lw=2.0, label="$\\hat μ(x)$")
ax.set_xlabel("time after impact (ms)")
ax.set_ylabel("head acceleration (g)")
ax.grid(alpha=0.25)
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(FIG_DIR / "distributional-regression.png",
            dpi=220, bbox_inches="tight")
plt.close(fig)
print(f"Saved {FIG_DIR / 'distributional-regression.png'}")
