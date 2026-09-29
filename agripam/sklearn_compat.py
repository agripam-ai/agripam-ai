"""Narrow compatibility loader for the Yu et al. scikit-learn 0.24 forest.

scikit-learn 1.3 added a ``missing_go_to_left`` byte to tree nodes.  The
published model predates that field.  This loader adds the field with its
legacy-equivalent value (zero) while leaving every published split, threshold,
leaf value and forest parameter unchanged.
"""
from __future__ import annotations

import pickle
import warnings
from pathlib import Path

import numpy as np
from sklearn.tree import _tree


class _LegacyTree:
    def __init__(self, *args):
        self.args = args

    def __setstate__(self, state):
        self.state = state


class _LegacyTreeUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "sklearn.tree._tree" and name == "Tree":
            return _LegacyTree
        return super().find_class(module, name)


_NODE_DTYPE = np.dtype({
    "names": ["left_child", "right_child", "feature", "threshold", "impurity",
              "n_node_samples", "weighted_n_node_samples", "missing_go_to_left"],
    "formats": ["<i8", "<i8", "<i8", "<f8", "<f8", "<i8", "<f8", "u1"],
    "offsets": [0, 8, 16, 24, 32, 40, 48, 56],
    "itemsize": 64,
})


def load_legacy_sklearn_forest(path: str | Path):
    """Load the released 0.24 forest without retraining or changing its nodes."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with Path(path).open("rb") as handle:
            model = _LegacyTreeUnpickler(handle).load()
    for estimator in model.estimators_:
        proxy = estimator.tree_
        old_nodes = proxy.state["nodes"]
        new_nodes = np.zeros(old_nodes.shape, dtype=_NODE_DTYPE)
        for field in old_nodes.dtype.names:
            new_nodes[field] = old_nodes[field]
        state = proxy.state.copy()
        state["nodes"] = new_nodes
        tree = _tree.Tree(*proxy.args)
        tree.__setstate__(state)
        estimator.tree_ = tree
    return model
