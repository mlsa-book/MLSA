"""
Transformer-based Multi-State Model via PEM Reduction
=====================================================
Combines the PEM reduction (Ch. 22) with the multi-state → single-event
reduction (Ch. 23) and uses a causal transformer to learn non-Markov
transition hazards from the history of the process.

Illness-death model with recovery:
  States: 0 (healthy), 1 (ill), 2 (dead, absorbing)
  Transitions: 0→1, 0→2, 1→0, 1→2

Non-Markov effects (ground truth):
  h_{01} = λ_{01} · (1 + α · n_prev_episodes)
  h_{02} = λ_{02}                               (Markov baseline)
  h_{10} = λ_{10} · exp(−γ · n_episodes)
  h_{12} = λ_{12} · (1 + δ · n_episodes)

The causal attention mask ensures that the hazard prediction at interval j
depends only on information from intervals <= j (past and present),
so that information flows from past to present, respecting the arrow of time.
"""

import os
os.environ["PYTHONUNBUFFERED"] = "1"

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import math
import sys
from collections import defaultdict

def log(msg):
    print(msg, flush=True)

TRANSITIONS = ["0→1", "0→2", "1→0", "1→2"]
TRANS_IDX = {t: i for i, t in enumerate(TRANSITIONS)}
N_TRANS = len(TRANSITIONS)
STATE_TRANSITIONS = {0: [0, 1], 1: [2, 3]}  # indices into TRANSITIONS


# ============================================================
# 1. SIMULATION
# ============================================================

def simulate_illness_death(n_subjects, params, max_time=10.0, seed=42):
    rng = np.random.default_rng(seed)
    trajectories = []

    for _ in range(n_subjects):
        t = 0.0
        state = 0
        n_episodes = 0  # number of times entered state 1
        sojourns = []

        while t < max_time and state != 2:
            if state == 0:
                h01 = params["base_h01"] * (1 + params["alpha"] * n_episodes)
                h02 = params["base_h02"]
                h_total = h01 + h02
                dt = rng.exponential(1.0 / h_total)
                if t + dt >= max_time:
                    sojourns.append((0, t, max_time, None))
                    break
                if rng.random() < h01 / h_total:
                    next_state = 1
                    n_episodes += 1
                else:
                    next_state = 2
                sojourns.append((0, t, t + dt, next_state))
                t += dt
                state = next_state
            else:  # state == 1
                h10 = params["base_h10"] * np.exp(-params["gamma"] * n_episodes)
                h12 = params["base_h12"] * (1 + params["delta"] * n_episodes)
                h_total = h10 + h12
                dt = rng.exponential(1.0 / h_total)
                if t + dt >= max_time:
                    sojourns.append((1, t, max_time, None))
                    break
                if rng.random() < h10 / h_total:
                    next_state = 0
                else:
                    next_state = 2
                sojourns.append((1, t, t + dt, next_state))
                t += dt
                state = next_state

        trajectories.append(sojourns)
    return trajectories


# ============================================================
# 2. PEM EXPANSION
# ============================================================

def create_pem_sequences(trajectories, n_intervals=100, max_time=10.0,
                         params=None):
    """Expand each trajectory into a sequence of PEM interval tokens.

    Each token represents one time interval for one subject. Features:
      - t_j:   interval midpoint (time representation)
      - state: current state (0 or 1)
    Targets per token:
      - delta[k] for k in {0..3}: transition event indicators
    Offset:
      - log(t_ij): log time-at-risk
    Metadata:
      - n_episodes: number of illness episodes so far (for evaluation)
      - true_hazards[k]: ground-truth hazard for each transition
    """
    cut_points = np.linspace(0, max_time, n_intervals + 1)
    interval_width = cut_points[1] - cut_points[0]

    all_sequences = []

    for subj_idx, sojourns in enumerate(trajectories):
        tokens = []
        n_episodes = 0

        for sojourn in sojourns:
            soj_state, t_enter, t_exit, next_state = sojourn
            if soj_state == 1 and (len(tokens) == 0 or
                                    tokens[-1]["state"] != 1 or
                                    tokens[-1].get("_soj_start") != t_enter):
                pass  # n_episodes already incremented during simulation

            # count episodes for ground-truth computation
            episode_count = n_episodes
            if soj_state == 1:
                # n_episodes for state-1 sojourn = episodes up to now
                # we need to recount from sojourns
                episode_count = sum(1 for s in sojourns[:sojourns.index(sojourn) + 1]
                                    if s[0] == 1 and s[1] <= t_enter + 1e-10)

            j_start = max(0, int(np.searchsorted(cut_points, t_enter, side="right")) - 1)
            j_end = min(n_intervals, int(np.searchsorted(cut_points, t_exit, side="right")))

            for j in range(j_start, j_end):
                a_left = cut_points[j]
                a_right = cut_points[j + 1]
                t_j = (a_left + a_right) / 2.0

                eff_start = max(a_left, t_enter)
                eff_end = min(a_right, t_exit)
                t_ij = eff_end - eff_start
                if t_ij < 1e-8:
                    continue

                deltas = np.zeros(N_TRANS, dtype=np.float32)
                is_last = (t_exit <= a_right + 1e-10) and (t_exit >= a_left - 1e-10)
                if is_last and next_state is not None:
                    if soj_state == 0 and next_state == 1:
                        deltas[0] = 1.0
                    elif soj_state == 0 and next_state == 2:
                        deltas[1] = 1.0
                    elif soj_state == 1 and next_state == 0:
                        deltas[2] = 1.0
                    elif soj_state == 1 and next_state == 2:
                        deltas[3] = 1.0

                true_h = np.zeros(N_TRANS, dtype=np.float32)
                if params is not None:
                    if soj_state == 0:
                        true_h[0] = params["base_h01"] * (1 + params["alpha"] * episode_count)
                        true_h[1] = params["base_h02"]
                    else:
                        true_h[2] = params["base_h10"] * np.exp(-params["gamma"] * episode_count)
                        true_h[3] = params["base_h12"] * (1 + params["delta"] * episode_count)

                trans_mask = np.zeros(N_TRANS, dtype=np.float32)
                for idx in STATE_TRANSITIONS[soj_state]:
                    trans_mask[idx] = 1.0

                tokens.append({
                    "t_j": np.float32(t_j),
                    "state": soj_state,
                    "log_offset": np.float32(np.log(max(t_ij, 1e-8))),
                    "deltas": deltas,
                    "trans_mask": trans_mask,
                    "n_episodes": episode_count,
                    "true_hazards": true_h,
                    "_soj_start": t_enter,
                })

        # update n_episodes by counting state-1 sojourns
        if tokens:
            all_sequences.append(tokens)

    return all_sequences


# recount episodes properly during expansion
def create_pem_sequences_v2(trajectories, n_intervals=100, max_time=10.0,
                            params=None):
    cut_points = np.linspace(0, max_time, n_intervals + 1)

    all_sequences = []
    for sojourns in trajectories:
        tokens = []
        n_episodes = 0

        for s_idx, sojourn in enumerate(sojourns):
            soj_state, t_enter, t_exit, next_state = sojourn

            if soj_state == 1:
                n_episodes_current = n_episodes + 1
            else:
                n_episodes_current = n_episodes

            j_start = max(0, int(np.searchsorted(cut_points, t_enter, side="right")) - 1)
            j_end = min(n_intervals, int(np.searchsorted(cut_points, t_exit, side="right")))

            for j in range(j_start, j_end):
                a_left, a_right = cut_points[j], cut_points[j + 1]
                t_j = (a_left + a_right) / 2.0
                eff_start = max(a_left, t_enter)
                eff_end = min(a_right, t_exit)
                t_ij = eff_end - eff_start
                if t_ij < 1e-8:
                    continue

                deltas = np.zeros(N_TRANS, dtype=np.float32)
                is_last = t_exit <= a_right + 1e-10
                if is_last and next_state is not None:
                    trans_key = f"{soj_state}→{next_state}"
                    if trans_key in TRANS_IDX:
                        deltas[TRANS_IDX[trans_key]] = 1.0

                true_h = np.zeros(N_TRANS, dtype=np.float32)
                if params is not None:
                    if soj_state == 0:
                        true_h[0] = params["base_h01"] * (
                            1 + params["alpha"] * n_episodes)
                        true_h[1] = params["base_h02"]
                    else:
                        true_h[2] = params["base_h10"] * np.exp(
                            -params["gamma"] * n_episodes_current)
                        true_h[3] = params["base_h12"] * (
                            1 + params["delta"] * n_episodes_current)

                trans_mask = np.zeros(N_TRANS, dtype=np.float32)
                for idx in STATE_TRANSITIONS[soj_state]:
                    trans_mask[idx] = 1.0

                tokens.append({
                    "t_j": np.float32(t_j),
                    "state": soj_state,
                    "log_offset": np.float32(np.log(max(t_ij, 1e-8))),
                    "deltas": deltas,
                    "trans_mask": trans_mask,
                    "n_episodes": n_episodes_current,
                    "true_hazards": true_h,
                })

            if soj_state == 1 and next_state == 0:
                n_episodes = n_episodes_current

        if tokens:
            all_sequences.append(tokens)

    return all_sequences


# ============================================================
# 3. DATASET
# ============================================================

class MultiStatePEMDataset(Dataset):
    def __init__(self, sequences):
        self.sequences = sequences

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        seq = self.sequences[idx]
        L = len(seq)
        t_j = np.array([tok["t_j"] for tok in seq], dtype=np.float32)
        states = np.array([tok["state"] for tok in seq], dtype=np.int64)
        log_offsets = np.array([tok["log_offset"] for tok in seq], dtype=np.float32)
        deltas = np.stack([tok["deltas"] for tok in seq])
        trans_masks = np.stack([tok["trans_mask"] for tok in seq])
        n_episodes = np.array([tok["n_episodes"] for tok in seq], dtype=np.int64)
        true_h = np.stack([tok["true_hazards"] for tok in seq])

        return {
            "t_j": torch.from_numpy(t_j),
            "states": torch.from_numpy(states),
            "log_offsets": torch.from_numpy(log_offsets),
            "deltas": torch.from_numpy(deltas),
            "trans_masks": torch.from_numpy(trans_masks),
            "n_episodes": torch.from_numpy(n_episodes),
            "true_hazards": torch.from_numpy(true_h),
            "lengths": L,
        }


def collate_fn(batch):
    max_len = max(b["lengths"] for b in batch)
    B = len(batch)

    t_j = torch.zeros(B, max_len)
    states = torch.zeros(B, max_len, dtype=torch.long)
    log_offsets = torch.zeros(B, max_len)
    deltas = torch.zeros(B, max_len, N_TRANS)
    trans_masks = torch.zeros(B, max_len, N_TRANS)
    n_episodes = torch.zeros(B, max_len, dtype=torch.long)
    true_hazards = torch.zeros(B, max_len, N_TRANS)
    padding_mask = torch.ones(B, max_len, dtype=torch.bool)  # True = padded
    lengths = torch.zeros(B, dtype=torch.long)

    for i, b in enumerate(batch):
        L = b["lengths"]
        t_j[i, :L] = b["t_j"]
        states[i, :L] = b["states"]
        log_offsets[i, :L] = b["log_offsets"]
        deltas[i, :L] = b["deltas"]
        trans_masks[i, :L] = b["trans_masks"]
        n_episodes[i, :L] = b["n_episodes"]
        true_hazards[i, :L] = b["true_hazards"]
        padding_mask[i, :L] = False
        lengths[i] = L

    return {
        "t_j": t_j,
        "states": states,
        "log_offsets": log_offsets,
        "deltas": deltas,
        "trans_masks": trans_masks,
        "n_episodes": n_episodes,
        "true_hazards": true_hazards,
        "padding_mask": padding_mask,
        "lengths": lengths,
    }


# ============================================================
# 4. MODELS
# ============================================================

class LearnedTimeEncoding(nn.Module):
    def __init__(self, d_out):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(1, d_out),
            nn.GELU(),
            nn.Linear(d_out, d_out),
        )

    def forward(self, t):
        return self.net(t.unsqueeze(-1))


class TransformerMultiState(nn.Module):
    def __init__(self, d_model=64, nhead=4, num_layers=3, dim_ff=128,
                 dropout=0.1, n_states=2, n_transitions=N_TRANS):
        super().__init__()
        self.d_model = d_model
        self.time_enc = LearnedTimeEncoding(d_model // 2)
        self.state_emb = nn.Embedding(n_states, d_model // 2)
        self.input_proj = nn.Linear(d_model, d_model)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead, dim_feedforward=dim_ff,
            dropout=dropout, batch_first=True, norm_first=True,
        )
        self.transformer = nn.TransformerEncoder(
            encoder_layer, num_layers=num_layers
        )
        self.output_head = nn.Sequential(
            nn.Linear(d_model, dim_ff),
            nn.GELU(),
            nn.Linear(dim_ff, n_transitions),
        )

    def forward(self, t_j, states, padding_mask=None, **kwargs):
        B, L = t_j.shape
        time_feat = self.time_enc(t_j)
        state_feat = self.state_emb(states)
        x = self.input_proj(torch.cat([time_feat, state_feat], dim=-1))

        causal_mask = nn.Transformer.generate_square_subsequent_mask(
            L, device=x.device, dtype=x.dtype
        )
        if padding_mask is not None:
            pad_float = padding_mask.float().masked_fill(padding_mask, float("-inf"))
        else:
            pad_float = None

        h = self.transformer(x, mask=causal_mask, src_key_padding_mask=pad_float)
        log_hazards = self.output_head(h)
        return log_hazards


class MarkovBaseline(nn.Module):
    """Feedforward model that processes each interval independently."""
    def __init__(self, d_model=64, dim_ff=128, n_states=2,
                 n_transitions=N_TRANS):
        super().__init__()
        self.time_enc = LearnedTimeEncoding(d_model // 2)
        self.state_emb = nn.Embedding(n_states, d_model // 2)
        self.net = nn.Sequential(
            nn.Linear(d_model, dim_ff),
            nn.GELU(),
            nn.Linear(dim_ff, dim_ff),
            nn.GELU(),
            nn.Linear(dim_ff, n_transitions),
        )

    def forward(self, t_j, states, padding_mask=None, **kwargs):
        time_feat = self.time_enc(t_j)
        state_feat = self.state_emb(states)
        x = torch.cat([time_feat, state_feat], dim=-1)
        return self.net(x)


class OraclePEM(nn.Module):
    """Simple PEM: log(h) = f(state, n_episodes). No time feature, no NN."""
    def __init__(self, n_states=2, n_transitions=N_TRANS, max_episodes=10):
        super().__init__()
        self.state_episode_logh = nn.Parameter(
            torch.zeros(n_states, max_episodes, n_transitions)
        )
        self.max_episodes = max_episodes

    def forward(self, t_j, states, padding_mask=None, n_episodes=None,
                **kwargs):
        ep = n_episodes.clamp(0, self.max_episodes - 1)
        log_h = self.state_episode_logh[states, ep]
        return log_h


class OracleModel(nn.Module):
    """FF model with access to the true time-dependent covariate n_episodes."""
    def __init__(self, d_model=64, dim_ff=128, n_states=2,
                 n_transitions=N_TRANS):
        super().__init__()
        self.time_enc = LearnedTimeEncoding(d_model // 2)
        self.state_emb = nn.Embedding(n_states, d_model // 4)
        self.episode_enc = nn.Sequential(
            nn.Linear(1, d_model // 4),
            nn.GELU(),
        )
        self.net = nn.Sequential(
            nn.Linear(d_model, dim_ff),
            nn.GELU(),
            nn.Linear(dim_ff, dim_ff),
            nn.GELU(),
            nn.Linear(dim_ff, n_transitions),
        )

    def forward(self, t_j, states, padding_mask=None, n_episodes=None,
                **kwargs):
        time_feat = self.time_enc(t_j)
        state_feat = self.state_emb(states)
        ep_feat = self.episode_enc(n_episodes.float().unsqueeze(-1))
        x = torch.cat([time_feat, state_feat, ep_feat], dim=-1)
        return self.net(x)


# ============================================================
# 5. LOSS
# ============================================================

def poisson_nll_loss(log_h, deltas, log_offsets, trans_masks):
    """Poisson NLL: exp(log_μ) − δ·log_μ, where log_μ = log_h + log_offset."""
    log_mu = log_h + log_offsets.unsqueeze(-1)
    loss = torch.exp(log_mu) - deltas * log_mu
    loss = loss * trans_masks
    n_valid = trans_masks.sum()
    if n_valid > 0:
        return loss.sum() / n_valid
    return torch.tensor(0.0, device=log_h.device)


# ============================================================
# 6. TRAINING
# ============================================================

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def to_device(batch, device=DEVICE):
    return {k: v.to(device) if isinstance(v, torch.Tensor) else v
            for k, v in batch.items()}


def train_model(model, train_loader, val_loader, n_epochs=150, lr=1e-3,
                label="Model", patience=20):
    model = model.to(DEVICE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, n_epochs)

    train_losses, val_losses = [], []
    best_val = float("inf")
    best_state = None
    epochs_no_improve = 0

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        n_batches = 0
        for batch in train_loader:
            batch = to_device(batch)
            log_h = model(batch["t_j"], batch["states"], batch["padding_mask"],
                         n_episodes=batch["n_episodes"])
            mask = batch["trans_masks"] * (~batch["padding_mask"]).unsqueeze(-1).float()
            loss = poisson_nll_loss(log_h, batch["deltas"],
                                    batch["log_offsets"], mask)
            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            epoch_loss += loss.item()
            n_batches += 1

        scheduler.step()
        avg_train = epoch_loss / max(n_batches, 1)
        train_losses.append(avg_train)

        model.eval()
        val_loss = 0.0
        n_val = 0
        with torch.no_grad():
            for batch in val_loader:
                batch = to_device(batch)
                log_h = model(batch["t_j"], batch["states"],
                              batch["padding_mask"],
                              n_episodes=batch["n_episodes"])
                mask = (batch["trans_masks"]
                        * (~batch["padding_mask"]).unsqueeze(-1).float())
                loss = poisson_nll_loss(log_h, batch["deltas"],
                                        batch["log_offsets"], mask)
                val_loss += loss.item()
                n_val += 1
        avg_val = val_loss / max(n_val, 1)
        val_losses.append(avg_val)

        if avg_val < best_val - 1e-6:
            best_val = avg_val
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1

        if (epoch + 1) % 10 == 0 or epoch == 0:
            log(f"  [{label}] Epoch {epoch+1:3d}/{n_epochs}: "
                f"train={avg_train:.4f}  val={avg_val:.4f}"
                f"  (best={best_val:.4f}, wait={epochs_no_improve})")

        if epochs_no_improve >= patience:
            log(f"  [{label}] Early stopping at epoch {epoch+1} "
                f"(best val={best_val:.4f} at epoch {epoch+1-patience})")
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model = model.cpu()
    return train_losses, val_losses


# ============================================================
# 7. EVALUATION
# ============================================================

def collect_predictions(model, loader):
    model.eval()
    results = defaultdict(lambda: defaultdict(list))

    with torch.no_grad():
        for batch in loader:
            log_h = model(batch["t_j"], batch["states"],
                          batch["padding_mask"],
                          n_episodes=batch["n_episodes"])
            pred_h = torch.exp(log_h).numpy()
            true_h = batch["true_hazards"].numpy()
            trans_masks = batch["trans_masks"].numpy()
            n_eps = batch["n_episodes"].numpy()
            pad = batch["padding_mask"].numpy()
            states_np = batch["states"].numpy()

            B, L, _ = pred_h.shape
            for i in range(B):
                for j in range(L):
                    if pad[i, j]:
                        continue
                    ep = int(n_eps[i, j])
                    st = int(states_np[i, j])
                    for k in range(N_TRANS):
                        if trans_masks[i, j, k] > 0:
                            results[(TRANSITIONS[k], ep)]["pred"].append(
                                pred_h[i, j, k])
                            results[(TRANSITIONS[k], ep)]["true"].append(
                                true_h[i, j, k])

    summary = {}
    for key, vals in results.items():
        summary[key] = {
            "pred_mean": np.mean(vals["pred"]),
            "pred_std": np.std(vals["pred"]),
            "true_mean": np.mean(vals["true"]),
            "n_intervals": len(vals["pred"]),
        }
    return summary


def plot_results(transformer_summary, markov_summary, params,
                 oracle_summary=None, save_path=None):
    transitions_with_history = ["0→1", "1→0", "1→2"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False)
    fig.suptitle("Predicted vs True Hazards by Illness Episode Count",
                 fontsize=13, y=1.02)

    n_bars = 4 if oracle_summary else 3
    bar_width = 0.8 / n_bars
    for ax_idx, trans in enumerate(transitions_with_history):
        ax = axes[ax_idx]
        episodes = sorted(set(ep for (t, ep) in transformer_summary.keys()
                               if t == trans and ep > 0))
        if not episodes:
            continue

        x = np.arange(len(episodes))

        def get_vals(summary):
            return [summary.get((trans, ep), {}).get("pred_mean", 0)
                    for ep in episodes]

        true_vals = [transformer_summary.get((trans, ep), {}).get("true_mean", 0)
                     for ep in episodes]

        offset = 0
        ax.bar(x + offset * bar_width, true_vals, bar_width, label="True",
               color="#2c7bb6", alpha=0.85)
        offset += 1
        if oracle_summary:
            ax.bar(x + offset * bar_width, get_vals(oracle_summary), bar_width,
                   label="Oracle", color="#1a9641", alpha=0.85)
            offset += 1
        ax.bar(x + offset * bar_width, get_vals(transformer_summary), bar_width,
               label="Transformer", color="#d7191c", alpha=0.85)
        offset += 1
        ax.bar(x + offset * bar_width, get_vals(markov_summary), bar_width,
               label="Markov (FF)", color="#fdae61", alpha=0.85)

        ax.set_xlabel("Illness Episode Count")
        ax.set_ylabel("Hazard Rate")
        ax.set_title(f"Transition {trans}")
        ax.set_xticks(x + bar_width * (n_bars - 1) / 2)
        ax.set_xticklabels([str(e) for e in episodes])
        if ax_idx == 0:
            ax.legend(fontsize=8)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        log(f"Saved figure to {save_path}")
    plt.show()


def plot_losses(trans_losses, markov_losses, save_path=None):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (train_l, val_l), title in zip(
        axes,
        [(trans_losses, markov_losses)[:1][0], (trans_losses, markov_losses)[1]],
        ["", ""]
    ):
        pass

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    ax1.plot(trans_losses[0], label="Train", color="#2c7bb6")
    ax1.plot(trans_losses[1], label="Val", color="#d7191c")
    ax1.set_title("Transformer")
    ax1.set_xlabel("Epoch")
    ax1.set_ylabel("Poisson NLL")
    ax1.legend()

    ax2.plot(markov_losses[0], label="Train", color="#2c7bb6")
    ax2.plot(markov_losses[1], label="Val", color="#d7191c")
    ax2.set_title("Markov Baseline (FF)")
    ax2.set_xlabel("Epoch")
    ax2.legend()

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()


def plot_hazard_over_time(model, test_loader, params, model_name="Model",
                          save_path=None):
    """Plot predicted hazard as function of time, grouped by episode count."""
    model.eval()
    data = defaultdict(lambda: defaultdict(list))

    with torch.no_grad():
        for batch in test_loader:
            log_h = model(batch["t_j"], batch["states"],
                          batch["padding_mask"],
                          n_episodes=batch["n_episodes"])
            pred_h = torch.exp(log_h).numpy()
            true_h = batch["true_hazards"].numpy()
            t_vals = batch["t_j"].numpy()
            trans_masks = batch["trans_masks"].numpy()
            n_eps = batch["n_episodes"].numpy()
            pad = batch["padding_mask"].numpy()

            B, L, _ = pred_h.shape
            for i in range(B):
                for j in range(L):
                    if pad[i, j]:
                        continue
                    ep = int(n_eps[i, j])
                    for k in range(N_TRANS):
                        if trans_masks[i, j, k] > 0:
                            data[(TRANSITIONS[k], ep)]["t"].append(t_vals[i, j])
                            data[(TRANSITIONS[k], ep)]["pred"].append(
                                pred_h[i, j, k])
                            data[(TRANSITIONS[k], ep)]["true"].append(
                                true_h[i, j, k])

    transitions_to_plot = ["0→1", "1→0", "1→2"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(f"{model_name}: Hazard Over Time by Episode Count",
                 fontsize=13, y=1.02)
    colors = {1: "#2c7bb6", 2: "#d7191c", 3: "#1a9641", 4: "#fdae61"}

    for ax_idx, trans in enumerate(transitions_to_plot):
        ax = axes[ax_idx]
        episodes = sorted(set(ep for (t, ep) in data.keys()
                               if t == trans and ep > 0))

        for ep in episodes[:4]:
            key = (trans, ep)
            if key not in data or len(data[key]["t"]) < 10:
                continue
            t_arr = np.array(data[key]["t"])
            pred_arr = np.array(data[key]["pred"])
            true_arr = np.array(data[key]["true"])

            # bin by time for cleaner plot
            n_bins = 20
            t_bins = np.linspace(t_arr.min(), t_arr.max(), n_bins + 1)
            t_mids, pred_mids, true_mids = [], [], []
            for b in range(n_bins):
                mask = (t_arr >= t_bins[b]) & (t_arr < t_bins[b + 1])
                if mask.sum() > 5:
                    t_mids.append((t_bins[b] + t_bins[b + 1]) / 2)
                    pred_mids.append(pred_arr[mask].mean())
                    true_mids.append(true_arr[mask].mean())

            c = colors.get(ep, "#999999")
            ax.plot(t_mids, pred_mids, '-', color=c, alpha=0.9,
                    label=f"Pred (ep={ep})")
            ax.plot(t_mids, true_mids, '--', color=c, alpha=0.6,
                    label=f"True (ep={ep})")

        ax.set_xlabel("Time")
        ax.set_ylabel("Hazard Rate")
        ax.set_title(f"Transition {trans}")
        ax.legend(fontsize=7)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()


def print_summary_table(trans_summary, markov_summary, oracle_summary=None):
    log("\n" + "=" * 88)
    log(f"{'Transition':<10} {'Episodes':<10} {'True h':<10} "
        f"{'Oracle':<10} {'Transf.':<10} {'Markov':<10} {'N intervals':<12}")
    log("-" * 88)

    all_keys = set(list(trans_summary.keys()) + list(markov_summary.keys()))
    if oracle_summary:
        all_keys |= set(oracle_summary.keys())
    keys = sorted(all_keys)
    for trans_name, ep in keys:
        if ep == 0:
            continue
        true_val = trans_summary.get((trans_name, ep), {}).get("true_mean", float("nan"))
        oracle_val = (oracle_summary or {}).get((trans_name, ep), {}).get("pred_mean", float("nan"))
        trans_val = trans_summary.get((trans_name, ep), {}).get("pred_mean", float("nan"))
        mark_val = markov_summary.get((trans_name, ep), {}).get("pred_mean", float("nan"))
        n_int = trans_summary.get((trans_name, ep), {}).get("n_intervals", 0)
        log(f"{trans_name:<10} {ep:<10} {true_val:<10.4f} "
            f"{oracle_val:<10.4f} {trans_val:<10.4f} {mark_val:<10.4f} {n_int:<12}")
    log("=" * 88)


# ============================================================
# 8. MAIN
# ============================================================

def main():
    log("=" * 60)
    log("Transformer Multi-State Model via PEM Reduction")
    log(f"Device: {DEVICE}")
    log("=" * 60)

    params = {
        "base_h01": 0.5,
        "base_h02": 0.05,
        "base_h10": 1.5,
        "base_h12": 0.2,
        "alpha": 0.5,
        "gamma": 0.4,
        "delta": 0.7,
    }
    MAX_TIME = 10.0
    N_INTERVALS = 50
    N_TRAIN = 2000
    N_TEST = 500

    fig_dir = os.path.join(os.path.dirname(__file__), "..", "Figures",
                            "neuralnetworks")
    os.makedirs(fig_dir, exist_ok=True)
    fig_dir = os.path.abspath(fig_dir)
    log(f"Figures will be saved to: {fig_dir}")

    # --- Simulate ---
    log("\n[1] Simulating non-Markov illness-death data...")
    train_traj = simulate_illness_death(N_TRAIN, params, MAX_TIME, seed=42)
    test_traj = simulate_illness_death(N_TEST, params, MAX_TIME, seed=123)

    def traj_stats(trajs):
        n_events = sum(1 for t in trajs
                       if any(s[3] == 2 for s in t))
        n_episodes = [sum(1 for s in t if s[0] == 1) for t in trajs]
        return n_events, n_episodes

    n_deaths, ep_counts = traj_stats(train_traj)
    log(f"  Train: {N_TRAIN} subjects, {n_deaths} deaths")
    log(f"  Episode distribution: "
        f"mean={np.mean(ep_counts):.1f}, max={np.max(ep_counts)}")

    # --- PEM expansion ---
    log("\n[2] Creating PEM-expanded sequences...")
    train_seqs = create_pem_sequences_v2(
        train_traj, N_INTERVALS, MAX_TIME, params)
    test_seqs = create_pem_sequences_v2(
        test_traj, N_INTERVALS, MAX_TIME, params)

    seq_lens = [len(s) for s in train_seqs]
    log(f"  Train sequences: {len(train_seqs)}, "
        f"avg length: {np.mean(seq_lens):.0f}, "
        f"max: {np.max(seq_lens)}")

    train_ds = MultiStatePEMDataset(train_seqs)
    test_ds = MultiStatePEMDataset(test_seqs)
    bs = 128
    nw = 2 if torch.cuda.is_available() else 0
    train_loader = DataLoader(train_ds, batch_size=bs, shuffle=True,
                              collate_fn=collate_fn, num_workers=nw,
                              pin_memory=torch.cuda.is_available())
    test_loader = DataLoader(test_ds, batch_size=bs * 2, shuffle=False,
                             collate_fn=collate_fn, num_workers=nw,
                             pin_memory=torch.cuda.is_available())

    # --- Train Transformer ---
    log("\n[3] Training Transformer model...")
    n_epochs = 60
    transformer = TransformerMultiState(
        d_model=48, nhead=4, num_layers=2, dim_ff=96, dropout=0.1
    )
    n_params_t = sum(p.numel() for p in transformer.parameters())
    log(f"  Parameters: {n_params_t:,}")

    trans_train_loss, trans_val_loss = train_model(
        transformer, train_loader, test_loader, n_epochs=n_epochs, lr=1e-3,
        label="Transformer"
    )

    # --- Train Markov baseline ---
    log("\n[4] Training Markov baseline (feedforward)...")
    markov = MarkovBaseline(d_model=48, dim_ff=96)
    n_params_m = sum(p.numel() for p in markov.parameters())
    log(f"  Parameters: {n_params_m:,}")

    mark_train_loss, mark_val_loss = train_model(
        markov, train_loader, test_loader, n_epochs=n_epochs, lr=1e-3,
        label="Markov"
    )

    # --- Train Oracle NN ---
    log("\n[5] Training Oracle NN model (knows n_episodes)...")
    oracle = OracleModel(d_model=48, dim_ff=96)
    n_params_o = sum(p.numel() for p in oracle.parameters())
    log(f"  Parameters: {n_params_o:,}")

    oracle_train_loss, oracle_val_loss = train_model(
        oracle, train_loader, test_loader, n_epochs=n_epochs, lr=1e-3,
        label="Oracle NN"
    )

    # --- Train Oracle PEM ---
    log("\n[6] Training Oracle PEM (state × episode intercepts, no time)...")
    oracle_pem = OraclePEM()
    n_params_p = sum(p.numel() for p in oracle_pem.parameters())
    log(f"  Parameters: {n_params_p:,}")

    pem_train_loss, pem_val_loss = train_model(
        oracle_pem, train_loader, test_loader, n_epochs=n_epochs, lr=1e-2,
        label="Oracle PEM"
    )

    # --- Evaluate ---
    log("\n[7] Evaluating...")
    trans_summary = collect_predictions(transformer, test_loader)
    markov_summary = collect_predictions(markov, test_loader)
    oracle_summary = collect_predictions(oracle, test_loader)
    pem_summary = collect_predictions(oracle_pem, test_loader)

    print_summary_table(trans_summary, markov_summary, oracle_summary)
    log("\nOracle PEM:")
    print_summary_table(trans_summary, markov_summary, pem_summary)

    # --- Transition count heatmap ---
    log("\n[7] Generating plots...")

    state_names = ["0 (healthy)", "1 (ill)", "2 (dead)", "censored"]
    trans_counts = np.zeros((3, 4), dtype=int)  # from states 0,1,2 → to 0,1,2,censored
    for traj in train_traj:
        for soj_state, t_enter, t_exit, next_state in traj:
            if next_state is not None:
                trans_counts[soj_state, next_state] += 1
            else:
                trans_counts[soj_state, 3] += 1

    fig, ax = plt.subplots(figsize=(6, 3.5))
    im = ax.imshow(trans_counts[:2, :], cmap="Blues", aspect="auto")
    for i in range(2):
        for j in range(4):
            v = trans_counts[i, j]
            ax.text(j, i, str(v), ha="center", va="center",
                    fontsize=12, color="white" if v > trans_counts[:2].max() * 0.6 else "black")
    ax.set_xticks(range(4))
    ax.set_xticklabels(state_names)
    ax.set_yticks(range(2))
    ax.set_yticklabels(["0 (healthy)", "1 (ill)"])
    ax.set_xlabel("To")
    ax.set_ylabel("From")
    ax.set_title(f"Transition Counts (N={N_TRAIN} subjects)")
    fig.colorbar(im, ax=ax, shrink=0.8)
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/transformer-ms-transition-counts.png",
                dpi=150, bbox_inches="tight")
    log(f"Saved figure to {fig_dir}/transformer-ms-transition-counts.png")
    plt.show()

    plot_results(
        trans_summary, markov_summary, params, oracle_summary=oracle_summary,
        save_path=f"{fig_dir}/transformer-ms-hazards-by-episode.png"
    )
    plot_hazard_over_time(
        transformer, test_loader, params, "Transformer",
        save_path=f"{fig_dir}/transformer-ms-hazard-over-time.png"
    )
    plot_hazard_over_time(
        oracle, test_loader, params, "Oracle NN",
        save_path=f"{fig_dir}/oracle-ms-hazard-over-time.png"
    )
    plot_hazard_over_time(
        oracle_pem, test_loader, params, "Oracle PEM",
        save_path=f"{fig_dir}/oracle-pem-ms-hazard-over-time.png"
    )

    log("\nDone!")


if __name__ == "__main__":
    main()
