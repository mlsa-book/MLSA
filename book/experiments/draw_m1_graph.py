"""Render the actual forward computation graph of the fitted M1
Weibull AFT network (from the example in @sec-nnet-param) using
``torch.fx``.

M1 architecture (identical to ``WeibullNet(p=1, scale_head='linear',
shape_head=None)`` in ``code.py``):

    x  (= complications, scalar 0/1)
    -> Linear(1 -> 1)                       (weight, bias)
    -> raw_scale
    -> log_lambda = log_t_init + s * tanh(raw_scale)

    global_log_k                            (parameter, no input)
    -> log_k = global_log_k.expand(...)

Forward outputs (log_lambda, log_k) feed the Weibull NLL.

Output: book/Figures/neuralnetworks/m1-weibull-graph.png
"""
from __future__ import annotations

import os
import sys

_SELF_DIR = os.path.dirname(os.path.abspath(__file__))
if _SELF_DIR in sys.path:
    sys.path.remove(_SELF_DIR)

import subprocess  # noqa: E402
from pathlib import Path  # noqa: E402

import graphviz  # noqa: E402
import torch  # noqa: E402
import torch.fx  # noqa: E402


FIG_DIR = Path(__file__).resolve().parent.parent / "Figures" / "neuralnetworks"
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ---- Minimal M1 (matches WeibullNet in code.py) --------------------------

class M1(torch.nn.Module):
    """Weibull AFT M1: linear scale head on 1 input, global shape bias."""

    def __init__(self, log_t_init: float = 7.0, scale_range: float = 2.5):
        super().__init__()
        self.scale_head = torch.nn.Linear(1, 1)
        self.global_log_k = torch.nn.Parameter(torch.zeros(1))
        self.log_t_init = log_t_init
        self.scale_range = scale_range

    def forward(self, x):
        raw_scale = self.scale_head(x).squeeze(-1)
        log_lambda = self.log_t_init + self.scale_range * torch.tanh(raw_scale)
        log_k = self.global_log_k.expand(x.shape[0])
        return log_lambda, log_k


torch.manual_seed(0)
m1 = M1()
traced = torch.fx.symbolic_trace(m1)


# ---- Manual graphviz render of the forward graph --------------------------
#
# torch.fx's bundled drawer dumps every internal field; for a textbook
# figure we want a small, hand-curated set of nodes.  We translate each
# fx node into a graphviz node with a short, math-style label and skip
# the noisy "get x.shape -> getitem 0 -> expand" sub-chain that exists
# only because ``global_log_k.expand`` needs the batch dimension.

dot = graphviz.Digraph("m1", format="png")
dot.attr("graph", rankdir="LR", bgcolor="transparent",
         nodesep="0.30", ranksep="0.55")
dot.attr("node", fontname="DejaVu Sans", fontsize="12",
         style="filled,rounded", shape="box", penwidth="1.2")
dot.attr("edge", fontname="DejaVu Sans", fontsize="10",
         color="#555555")

# Colour palette consistent with the architecture figure.
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


# Input.
add("x", "x  (complications)", C_INPUT)

# Constants in the bounded scale activation.
add("c_init", "log t_init = 7.0", C_CONST)
add("c_scale", "scale_range = 2.5", C_CONST)

# Linear scale head with explicit parameters.
add("w", "scale_head.weight\\n(1×1)", C_PARAM)
add("b", "scale_head.bias\\n(1)", C_PARAM)
add("lin", "Linear(1 → 1)\\nW·x + b", C_LINEAR)

dot.edge("x", "lin")
dot.edge("w", "lin", style="dashed")
dot.edge("b", "lin", style="dashed")

add("raw", "raw_scale = z", C_OP)
dot.edge("lin", "raw")

add("tanh", "tanh(z)", C_OP)
dot.edge("raw", "tanh")

add("mul", "× scale_range", C_OP)
dot.edge("tanh", "mul")
dot.edge("c_scale", "mul", style="dashed")

add("addlam", "+ log t_init", C_OP)
dot.edge("mul", "addlam")
dot.edge("c_init", "addlam", style="dashed")

add("loglam", "log λ", C_OUT_LAM)
dot.edge("addlam", "loglam")

# Shape: a single scalar parameter, broadcast over the batch.
add("glk", "global_log_k\\n(scalar parameter)", C_PARAM)
add("logk", "log k", C_OUT_K)
dot.edge("glk", "logk", label="expand(batch)")

# Final positivity response.
add("lam", "λ = exp(log λ)", C_OUT_LAM)
add("k", "k = exp(log k)", C_OUT_K)
dot.edge("loglam", "lam")
dot.edge("logk", "k")

# Same-rank groupings (keeps inputs/parameters on the left, outputs right).
with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("x")
    s.node("w")
    s.node("b")
    s.node("glk")

with dot.subgraph() as s:
    s.attr(rank="same")
    s.node("lam")
    s.node("k")

out_base = FIG_DIR / "m1-weibull-graph"
rendered = dot.render(filename=str(out_base), cleanup=True)
print(f"Saved {rendered}")

# Also print a textual fx trace for the script's stdout, so the reader
# can see what was actually traced.
print("\nForward trace (torch.fx):")
print(traced.graph)
