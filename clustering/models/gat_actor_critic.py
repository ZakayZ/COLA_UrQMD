"""Policy + scalar value head (mean-pooled GAT embeddings) for A2C / PPO."""

from __future__ import annotations

from typing import Any

import torch.nn as nn
from torch_geometric.data import Data

from models.policy import GATAffinityPolicy


class GATAffinityActorCritic(nn.Module):
    """Policy (:class:`~models.policy.GATAffinityPolicy` edge logits) + scalar ``V(s)``.

    ``V(s)`` is produced from **mean-pooled** GAT node embeddings (permutation-invariant graph
    summary) plus a small MLP — one forward pass for both heads.
    """

    def __init__(
        self,
        policy: GATAffinityPolicy,
        *,
        value_mlp_hidden: int = 128,
    ) -> None:
        super().__init__()
        self.policy = policy
        d = int(policy.node_embed_dim)
        self.value_head = nn.Sequential(
            nn.Linear(d, int(value_mlp_hidden)),
            nn.GELU(),
            nn.Linear(int(value_mlp_hidden), 1),
        )
        for lin in self.value_head:
            if isinstance(lin, nn.Linear):
                nn.init.orthogonal_(lin.weight, gain=1.0)
                if lin.bias is not None:
                    nn.init.zeros_(lin.bias)
        with torch.no_grad():
            self.value_head[-1].weight.mul_(0.01)

    def forward(
        self, data: Data, *, capture: dict[str, Any] | None = None
    ) -> tuple[torch.Tensor, torch.Tensor]:
        logits, h = self.policy.forward_logits_and_h(data, capture=capture)
        if h.size(0) == 0:
            v = self.value_head(h.new_zeros((self.policy.node_embed_dim,)))
        else:
            v = self.value_head(h.mean(dim=0))
        return logits, v.view(-1)[0]
