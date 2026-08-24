"""P6 wall-Cp recovery smoothing (gate G6.1): the surface-Cp sawtooth is a
per-triangle gradient-recovery artifact on the sliver wall triangulation, not
the artificial-density flux (root-caused 2026-07-08: nodal averaging on the
same field drops the G6.1 metric ~330x). `smooth_wall_tangential_gradients`
removes it with a normal-gated edge-neighbour average that preserves the sharp
TE. These are fast, mesh-free unit tests of that operator + its no-op default.
"""

import numpy as np

from pyfp3d.post.surface import smooth_wall_tangential_gradients


def _chain(n=6):
    """A chain of n triangles, consecutive ones sharing edge k=1 (self edge 2
    links back). Returns (grad, normal, area, adj) with a +/-1 sawtooth in the
    x-gradient and all-aligned +z normals."""
    adj = np.full((n, 3), -1, dtype=np.int64)
    for t in range(n):
        if t + 1 < n:
            adj[t, 1] = t + 1
        if t - 1 >= 0:
            adj[t, 2] = t - 1
    grad = np.zeros((n, 3))
    grad[:, 0] = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)   # 2-cell sawtooth
    normal = np.tile([0.0, 0.0, 1.0], (n, 1))
    area = np.ones(n)
    return grad, normal, area, adj


def test_zero_passes_is_identity():
    grad, normal, area, adj = _chain()
    out = smooth_wall_tangential_gradients(grad, normal, area, adj, n_passes=0)
    assert np.array_equal(out, grad)


def test_smoothing_reduces_sawtooth():
    grad, normal, area, adj = _chain(8)
    var0 = grad[:, 0].var()
    out = smooth_wall_tangential_gradients(grad, normal, area, adj, n_passes=1)
    var1 = out[:, 0].var()
    # one Jacobi pass over aligned neighbours collapses the +/-1 oscillation
    assert var1 < 0.2 * var0


def test_normal_gate_preserves_crease():
    """A neighbour whose outward normal is anti-parallel (a sharp TE fold) must
    NOT be averaged in -- otherwise the two surfaces' velocities mix."""
    grad, normal, area, adj = _chain(3)
    normal[2] = [0.0, 0.0, -1.0]      # triangle 2 across a crease from 1
    out = smooth_wall_tangential_gradients(grad, normal, area, adj,
                                           n_passes=1, cos_thresh=0.2)
    # triangle 1's only aligned neighbour is 0 (t2 gated out): grad_x -> mean(-1,1)=0
    assert np.isclose(out[1, 0], 0.0)
    # triangle 2 has no aligned neighbour (1 is anti-parallel) -> unchanged
    assert np.isclose(out[2, 0], grad[2, 0])


def test_smoothing_is_linear_in_grad():
    """Smoothing is a fixed linear operator (differentiable): S(a*g) = a*S(g)."""
    grad, normal, area, adj = _chain()
    s1 = smooth_wall_tangential_gradients(grad, normal, area, adj, n_passes=2)
    s2 = smooth_wall_tangential_gradients(3.0 * grad, normal, area, adj, n_passes=2)
    assert np.allclose(s2, 3.0 * s1, rtol=1e-12, atol=1e-14)
