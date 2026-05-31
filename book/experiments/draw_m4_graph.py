"""Render the actual forward computation graph of the fitted M4
Weibull AFT network (from the example in @sec-nnet-param), in the
same style as ``draw_m1_graph.py``.

M4 architecture (identical to ``WeibullNet(p=p, scale_head='mlp',
shape_head='mlp', hidden=16)`` in ``code.py``):

    x  (p-dimensional feature vector)
    -> scale head: Linear(p, 16) -> ReLU -> Linear(16, 1)  -> raw_scale
       log_lambda = log_t_init + s * tanh(raw_scale)
    -> shape head: Linear(p, 16) -> ReLU -> Linear(16, 1)  -> raw_shape
       log_k       = r * tanh(raw_shape)

Forward outputs (log_lambda, log_k) feed the Weibull NLL.

Output: book/Figures/neuralnetworks/m4-weibull-graph.png
"""
from __future__ import annotations

import os
import sys

_SELF_DIR = os.path.dirname(os.path.abspath(__file__))
if _SELF_DIR in sys.path:
    sys.path.remove(_SELF_DIR)

from pathlib import Path  # noqa: E402

import graphviz  # noqa: E402
import torch  # noqa: E402
import torch.fx  # noqa: E402


FIG_DIR = Path(__file__).resolve().parent.parent / "Figures" / "neuralnetworks"
FIG_DIR.mkdir(parents=True, exist_ok=True)

P_INPUT = 7  # matches the tumor design matrix after one-hot encoding
HIDDEN = 16


# ---- Minimal M4 (matches WeibullNet in code.py) --------------------------

class _MLPHead(torch.nn.Module):
    def __init__(self, p: int, hidden: int):
        super().__init__()
        self.fc1 = torch.nn.Linear(p, hidden)
        self.act = torch.nn.ReLU()
        self.fc2 = torch.nn.Linear(hidden, 1)

    def forward(self, x):
        return self.fc2(self.act(self.fc1(x)))


class M4(torch.nn.Module):
    """Weibull AFT M4: separate MLP heads for log lambda and log k."""

    def __init__(self, p: int = P_INPUT, hidden: int = HIDDEN,
                 log_t_init: float = 7.0, scale_range: float = 2.5,
                 shape_range: float = 1.5):
        super().__init__()
        self.scale_head = _MLPHead(p, hidden)
        self.shape_head = _MLPHead(p, hidden)
        self.log_t_init = log_t_init
        self.scale_range = scale_range
        self.shape_range = shape_range

    def forward(self, x):
        raw_scale = self.scale_head(x).squeeze(-1)
        log_lambda = self.log_t_init + self.scale_range * torch.tanh(raw_scale)
        raw_shape = self.shape_head(x).squeeze(-1)
        log_k = self.shape_range * torch.tanh(raw_shape)
        return log_lambda, log_k


torch.manual_seed(0)
m4 = M4()
traced = torch.fx.symbolic_trace(m4)


# ---- Manual graphviz render -----------------------------------------------
#
# Same palette as draw_m1_graph.py.  Both heads share the same shape;
# we draw them as two horizontal chains stacked vertically.

dot = graphviz.Digraph("m4", format="png")
dot.attr("graph", rankdir="LR", bgcolor="transparent",
         nodesep="0.30", ranksep="0.55")
dot.attr("node", fontname="DejaVu Sans", fontsize="12",
         style="filled,rounded", shape="box", penwidth="1.2")
dot.attr("edge", fontname="DejaVu Sans", fontsize="10",
         color="#555555")

C_INPUT  = ("#E3F2FD", "#42A5F5")
C_PARAM  = ("#FFF9C4", "#FBC02D")
C_LINEAR = ("#FFF3E0", "#FB8C00")
C_OP     = ("#F3E5F5", "#8E24AA")
C_OUT_LAM = ("#E8F5E9", "#43A047")
C_OUT_K  = ("#FCE4EC", "#C2185B")
C_CONST  = ("#ECEFF1", "#607D8B")


def add(name, label, palette):
    fc, ec = palette
    dot.node(name, label, fillcolor=fc, color=ec)


# Shared input.
add("x", f"x  (p = {P_INPUT})", C_INPUT)

# Shared stabilizer constants.
add("c_init", "log t_init = 7.0", C_CONST)
add("c_scale", "scale_range = 2.5", C_CONST)
add("c_shape", "shape_range = 1.5", C_CONST)


def add_mlp_head(prefix: str, w1_label: str, b1_label: str,
                 w2_label: str, b2_label: str):
    add(f"{prefix}_w1",
        f"{prefix}.fc1.weight\\n({HIDDEN}×{P_INPUT})", C_PARAM)
    add(f"{prefix}_b1",
        f"{prefix}.fc1.bias\\n({HIDDEN})", C_PARAM)
    add(f"{prefix}_lin1",
        f"Linear({P_INPUT} → {HIDDEN})\\nW₁·x + b₁", C_LINEAR)
    add(f"{prefix}_relu", "ReLU", C_OP)
    add(f"{prefix}_w2",
        f"{prefix}.fc2.weight\\n(1×{HIDDEN})", C_PARAM)
    add(f"{prefix}_b2",
        f"{prefix}.fc2.bias\\n(1)", C_PARAM)
    add(f"{prefix}_lin2",
        f"Linear({HIDDEN} → 1)\\nW₂·h + b₂", C_LINEAR)

    dot.edge("x", f"{prefix}_lin1")
    dot.edge(f"{prefix}_w1", f"{prefix}_lin1", style="dashed")
    dot.edge(f"{prefix}_b1", f"{prefix}_lin1", style="dashed")
    dot.edge(f"{prefix}_lin1", f"{prefix}_relu")
    dot.edge(f"{prefix}_relu", f"{prefix}_lin2")
    dot.edge(f"{prefix}_w2", f"{prefix}_lin2", style="dashed")
    dot.edge(f"{prefix}_b2", f"{prefix}_lin2", style="dashed")


# ---------- scale head + lambda tail ----------
add_mlp_head("scale", "W₁", "b₁", "W₂", "b₂")
add("raw_lam", "raw_scale = z_λ", C_OP)
add("tanh_lam", "tanh(z_λ)", C_OP)
add("mul_lam", "× scale_range", C_OP)
add("add_lam", "+ log t_init", C_OP)
add("loglam", "log λ", C_OUT_LAM)
add("lam", "λ = exp(log λ)", C_OUT_LAM)

dot.edge("scale_lin2", "raw_lam")
dot.edge("raw_lam", "tanh_lam")
dot.edge("tanh_lam", "mul_lam")
dot.edge("c_scale", "mul_lam", style="dashed")
dot.edge("mul_lam", "add_lam")
dot.edge("c_init", "add_lam", style="dashed")
dot.edge("add_lam", "loglam")
dot.edge("loglam", "lam")


# ---------- shape head + k tail ----------
add_mlp_head("shape", "U₁", "c₁", "U₂", "c₂")
add("raw_k", "raw_shape = z_k", C_OP)
add("tanh_k", "tanh(z_k)", C_OP)
add("mul_k", "× shape_range", C_OP)
add("logk", "log k", C_OUT_K)
add("k", "k = exp(log k)", C_OUT_K)

dot.edge("shape_lin2", "raw_k")
dot.edge("raw_k", "tanh_k")
dot.edge("tanh_k", "mul_k")
dot.edge("c_shape", "mul_k", style="dashed")
dot.edge("mul_k", "logk")
dot.edge("logk", "k")


# Same-rank groupings.
with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("x")
    s.node("scale_w1"); s.node("scale_b1")
    s.node("shape_w1"); s.node("shape_b1")

with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("scale_w2"); s.node("scale_b2")
    s.node("shape_w2"); s.node("shape_b2")

with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("lam"); s.node("k")


out_base = FIG_DIR / "m4-weibull-graph"
rendered = dot.render(filename=str(out_base), cleanup=True)
print(f"Saved {rendered}")

# Print the textual fx trace for reproducibility.
print("\nForward trace (torch.fx):")
print(traced.graph)
